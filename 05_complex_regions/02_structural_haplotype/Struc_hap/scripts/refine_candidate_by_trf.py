#!/usr/bin/env python3
"""Apply population TRF removal and recursive divergence splitting to an MSA.

The candidate's real source MSA is the only coordinate bridge. PAF and CIGAR
are intentionally absent, and duplicate occurrences of a path remain separate
population observations rather than overwriting one another in a dictionary.
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Mapping, Sequence

try:
    from .build_grouped_bundle_consensus import aligned_consensus_profile
    from .split_bundle_by_trf_divergence import (
        AlignedOccurrence, TRFHit, evaluate_repeat_candidates, merge_trf_hits,
        recursively_decompose_alignment, ungapped_bases_before,
    )
except ImportError:
    from build_grouped_bundle_consensus import aligned_consensus_profile
    from split_bundle_by_trf_divergence import (
        AlignedOccurrence, TRFHit, evaluate_repeat_candidates, merge_trf_hits,
        recursively_decompose_alignment, ungapped_bases_before,
    )


REPORT_COLUMNS = (
    "source_bundle_id", "position", "status", "candidate_msa_start", "candidate_msa_end",
    "msa_start", "msa_end", "source_path_id", "source_row", "repeat_bp", "repeat_start",
    "repeat_end", "repeat_ids", "population_union_source_occurrence_count",
    "population_union_source_path_count", "population_union_source_paths",
    "sequence_count", "full_repeat_count", "full_repeat_fraction",
    "trf_full_coverage", "population_presence_fraction", "structural_divergence",
    "min_tr_length_bp", "median_tr_length_bp", "max_tr_length_bp",
    "tr_length_allele_count", "tr_length_divergence", "length_divergent",
    "nucleotide_divergence", "divergent", "selected_iteration", "required_fragment_count",
    "supporting_sequence_count", "nonempty_sequence_count", "median_sequence_bp", "depth",
)


def _validate(candidate, args):
    width, left, right = candidate["msa_width"], candidate["msa_start"], candidate["msa_end"]
    if not all(isinstance(value, int) for value in (width, left, right)) or not 0 <= left < right <= width:
        raise ValueError("TRF candidate has invalid MSA boundaries")
    rows, sources = candidate["aligned_sequences"], candidate["sources"]
    if not rows or len(rows) != len(sources):
        raise ValueError("TRF candidate must have one aligned sequence per source occurrence")
    if len(candidate["profile"]) != width:
        raise ValueError("TRF candidate profile does not match its MSA width")
    for row, source in zip(rows, sources):
        if len(row) != width or row[:left].strip("-") or row[right:].strip("-"):
            raise ValueError("TRF candidate source sequence has invalid MSA padding")
        if len(source) != 5:
            raise ValueError("TRF candidate source must contain path,start,end,strand,row")
        _, start, end, strand, _ = source
        if start < 0 or end <= start or strand not in {"+", "-"}:
            raise ValueError("TRF candidate source has invalid interval or strand")
        if sum(base != "-" for base in row) != end - start:
            raise ValueError("TRF candidate ungapped source does not match its genomic interval")
    fractions = (
        "trf_full_coverage", "trf_bundle_removal_frequency", "trf_min_sequence_coverage",
        "trf_min_population_presence", "trf_min_length_divergence",
        "trf_min_nucleotide_divergence", "trf_min_fragment_frequency",
        "min_column_occupancy", "min_base_frequency",
    )
    for name in fractions:
        value = getattr(args, name)
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{name} must be a finite fraction between zero and one")
    if args.trf_full_coverage <= 0 or args.trf_bundle_removal_frequency <= 0:
        raise ValueError("TRF full-coverage/removal-frequency thresholds must be positive")
    for name in ("trf_split_min_bp", "trf_min_fragment_sequence_bp", "trf_min_fragment_count", "addition_min_bp"):
        if getattr(args, name) < 1:
            raise ValueError(f"{name} must be positive")


def _decompose_with_rescanning(records, unions, evaluations, args):
    """Remove scored population-MSA TRF union components recursively.

    Native TRFs have already been projected and merged across all occurrences,
    so overlapping annotations form one connected MSA region and cannot leave
    an unscored annotation residual after that region is removed.
    """
    minimum_bp = max(args.trf_min_fragment_sequence_bp, args.addition_min_bp)
    return recursively_decompose_alignment(
        records,
        evaluations,
        minimum_bp,
        args.trf_min_fragment_count,
        args.trf_min_fragment_frequency,
        False,
    )


def refine_candidate_by_trf(candidate: dict, records_by_path: Mapping[str, Sequence], args):
    """Return refined candidates and audit rows without modifying the input.

    Selected TRF records have already passed the caller's length/identity filter.
    A configurable full-repeat occurrence gate (default >=50%) is also used by
    the pre-group filter.
    Each occurrence must reach the coverage threshold with one merged TRF
    interval; disconnected repeat islands never add together. Otherwise long
    repeats are removed only when the population splitter finds divergence.
    Fragment support is tested only after a real divergent repeat is selected.
    """
    _validate(candidate, args)
    left, right, width = candidate["msa_start"], candidate["msa_end"], candidate["msa_width"]
    occurrences, sources_by_key, hits = [], {}, []
    for index, (source, sequence) in enumerate(zip(candidate["sources"], candidate["aligned_sequences"])):
        path, start, end, direction, row = source
        key = f"occurrence_{index:08d}"
        occurrence = AlignedOccurrence(
            header=key, bundle_id=str(candidate["source_bundle_id"]), path_id=key,
            start=start, end=end, direction=direction, bed_row=row,
            sequence=sequence[left:right].upper(),
        )
        occurrences.append(occurrence)
        sources_by_key[key] = source
        for repeat in records_by_path.get(path, ()):
            if repeat.start < 0 or repeat.end <= repeat.start or not math.isfinite(repeat.identity):
                raise ValueError(f"invalid selected TRF record for {path}")
            start_overlap, end_overlap = max(start, repeat.start), min(end, repeat.end)
            if end_overlap > start_overlap:
                hits.append(TRFHit(
                    path_id=key, repeat_id=repeat.repeat_id, start=repeat.start, end=repeat.end,
                    period=repeat.period, copy_number=repeat.copy_number, identity=repeat.identity,
                    overlap_start=start_overlap, overlap_end=end_overlap,
                ))
    unions = merge_trf_hits(hits)
    max_overlap_by_key = defaultdict(int)
    for repeat in unions:
        max_overlap_by_key[repeat.path_id] = max(
            max_overlap_by_key[repeat.path_id], repeat.length
        )
    full_count = sum(
        max_overlap_by_key[record.path_id] / record.ungapped_length
        >= args.trf_full_coverage
        for record in occurrences
    )
    base = {
        "source_bundle_id": candidate["source_bundle_id"], "position": candidate["position"],
        "candidate_msa_start": left, "candidate_msa_end": right,
        "msa_start": left, "msa_end": right, "sequence_count": len(occurrences),
        "full_repeat_count": full_count, "full_repeat_fraction": full_count / len(occurrences),
        "trf_full_coverage": args.trf_full_coverage,
    }
    if not hits:
        return [candidate], [{**base, "status": "no_trf_overlap"}]
    if full_count / len(occurrences) >= args.trf_bundle_removal_frequency:
        return [], [{**base, "status": "removed_full_repeat"}]

    # Native annotations are first projected to MSA columns and unioned across
    # the population. Length and divergence are evaluated on those population
    # union components, rather than on one annotation/occurrence at a time.
    evaluations = evaluate_repeat_candidates(
        unions, {record.path_id: record for record in occurrences}, occurrences,
        args.trf_split_min_bp, args.trf_min_sequence_coverage,
        args.trf_min_population_presence, args.trf_min_nucleotide_divergence,
        args.trf_min_length_divergence,
    ) if unions else []
    for evaluation in evaluations:
        original_paths = sorted({
            sources_by_key[path_id][0]
            for path_id in evaluation.population_region.source_path_ids
        })
        evaluation.metrics["population_union_source_path_count"] = len(original_paths)
        evaluation.metrics["population_union_source_paths"] = original_paths
    report = []
    if not any(evaluation.divergent for evaluation in evaluations):
        report.append({**base, "status": "retained_no_divergent_repeat"})
        for evaluation in evaluations:
            source = sources_by_key[evaluation.repeat.path_id]
            report.append({
                **base, "status": "repeat_conserved", "msa_start": left + evaluation.start_column,
                "msa_end": left + evaluation.end_column, "source_path_id": source[0], "source_row": source[4],
                "repeat_bp": evaluation.repeat_bp, "repeat_start": evaluation.repeat.start,
                "repeat_end": evaluation.repeat.end,
                "repeat_ids": ",".join(evaluation.population_region.source_repeat_ids),
                "population_union_source_occurrence_count": len(evaluation.population_region.members),
                "population_union_source_path_count": evaluation.metrics["population_union_source_path_count"],
                "population_union_source_paths": ",".join(evaluation.metrics["population_union_source_paths"]),
                "population_presence_fraction": evaluation.metrics["population_presence_fraction"],
                "structural_divergence": evaluation.metrics["structural_divergence"],
                "min_tr_length_bp": evaluation.metrics["min_sequence_bp"],
                "median_tr_length_bp": evaluation.metrics["median_sequence_bp"],
                "max_tr_length_bp": evaluation.metrics["max_sequence_bp"],
                "tr_length_allele_count": evaluation.metrics["tr_length_allele_count"],
                "tr_length_divergence": evaluation.metrics["mean_pairwise_tr_length_divergence"],
                "length_divergent": evaluation.length_divergent,
                "nucleotide_divergence": evaluation.metrics["non_gap_consensus_divergence"], "divergent": False,
            })
        return [candidate], report

    retained, decisions, selected, required = _decompose_with_rescanning(
        occurrences, unions, evaluations, args
    )
    selected_ids = {evaluation.candidate_id for evaluation in selected}
    for evaluation in evaluations:
        source = sources_by_key[evaluation.repeat.path_id]
        report.append({
            **base, "status": ("repeat_removed" if evaluation.candidate_id in selected_ids else
                               "repeat_divergent_unselected" if evaluation.divergent else "repeat_conserved"),
            "msa_start": left + evaluation.start_column, "msa_end": left + evaluation.end_column,
            "source_path_id": source[0], "source_row": source[4], "repeat_bp": evaluation.repeat_bp,
            "repeat_start": evaluation.repeat.start, "repeat_end": evaluation.repeat.end,
            "repeat_ids": ",".join(evaluation.population_region.source_repeat_ids),
            "population_union_source_occurrence_count": len(evaluation.population_region.members),
            "population_union_source_path_count": evaluation.metrics["population_union_source_path_count"],
            "population_union_source_paths": ",".join(evaluation.metrics["population_union_source_paths"]),
            "population_presence_fraction": evaluation.metrics["population_presence_fraction"],
            "structural_divergence": evaluation.metrics["structural_divergence"],
            "min_tr_length_bp": evaluation.metrics["min_sequence_bp"],
            "median_tr_length_bp": evaluation.metrics["median_sequence_bp"],
            "max_tr_length_bp": evaluation.metrics["max_sequence_bp"],
            "tr_length_allele_count": evaluation.metrics["tr_length_allele_count"],
            "tr_length_divergence": evaluation.metrics["mean_pairwise_tr_length_divergence"],
            "length_divergent": evaluation.length_divergent,
            "nucleotide_divergence": evaluation.metrics["non_gap_consensus_divergence"],
            "divergent": evaluation.divergent, "selected_iteration": evaluation.selected_iteration,
            "required_fragment_count": required,
        })
    for decision in decisions:
        if not decision.retained:
            report.append({
                **base, "status": "fragment_discarded_low_support",
                "msa_start": left + decision.child_start_column, "msa_end": left + decision.child_end_column,
                "required_fragment_count": required, "supporting_sequence_count": decision.supporting_sequence_count,
                "nonempty_sequence_count": decision.nonempty_sequence_count, "depth": decision.depth + 1,
                "selected_iteration": decision.iteration,
            })

    refined = []
    for segment in retained:
        local_start, local_end = segment.start_column, segment.end_column
        segment_rows, segment_sources = [], []
        for occurrence in occurrences:
            row = occurrence.sequence[local_start:local_end]
            if sum(base != "-" for base in row) < args.addition_min_bp:
                continue
            offset_start = ungapped_bases_before(occurrence.sequence, local_start)
            offset_end = ungapped_bases_before(occurrence.sequence, local_end)
            if occurrence.direction == "+":
                start, end = occurrence.start + offset_start, occurrence.start + offset_end
            else:
                start, end = occurrence.end - offset_end, occurrence.end - offset_start
            original_source = sources_by_key[occurrence.path_id]
            segment_rows.append(row)
            segment_sources.append((original_source[0], start, end, occurrence.direction, original_source[4]))
        try:
            profile, sequence, _, _, _ = aligned_consensus_profile(
                segment_rows, args.min_column_occupancy, args.min_base_frequency)
        except ValueError as error:
            if str(error) not in {"no MSA column passed the candidate consensus threshold",
                                  "cannot build a consensus from zero sequences"}:
                raise
            report.append({**base, "status": "fragment_discarded_no_consensus",
                           "msa_start": left + local_start, "msa_end": left + local_end})
            continue
        if len(sequence) < args.addition_min_bp:
            report.append({**base, "status": "fragment_discarded_short_consensus",
                           "msa_start": left + local_start, "msa_end": left + local_end})
            continue
        global_start, global_end = left + local_start, left + local_end
        refined.append({
            **candidate, "msa_start": global_start, "msa_end": global_end,
            "profile": "-" * global_start + profile + "-" * (width - global_end),
            "aligned_sequences": ["-" * global_start + row + "-" * (width - global_end)
                                  for row in segment_rows],
            "sequence": sequence, "sources": segment_sources, "whole": False, "trf_refined": True,
        })
        report.append({
            **base, "status": "fragment_retained", "msa_start": global_start, "msa_end": global_end,
            "required_fragment_count": required, "supporting_sequence_count": segment.supporting_sequence_count,
            "nonempty_sequence_count": segment.nonempty_sequence_count,
            "median_sequence_bp": segment.median_sequence_bp, "depth": segment.depth,
        })
    return refined, report
