#!/usr/bin/env python3
"""Split consensus bundles at recurrent, BED-gap-supported truncation endpoints.

Consumes the initial refine_v2 PAF (query=genome, target=bundle). A preliminary
retained-mapping BED gates endpoint evidence from incomplete alignments in its
uncovered genome regions. Split sequences must be remapped, not projected back
to genomes. All pieces are retained; no genomic or consensus sequence is deleted.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    from .build_bundle_consensus import iter_fasta, wrap_fasta
    from .select_gap_alignments import GapAlignment, maximum_coverage_selection
    from .refine_bundle_mappings import (
        chain_blocks, chain_to_mapping, deduplicate_same_locus, interval_total,
        merge_intervals, natural_name_key, normalize_bundle_name,
        parse_alias, read_existing_bundle_metadata, read_paf, resolve_alias, subtract_intervals,
    )
except ImportError:
    from build_bundle_consensus import iter_fasta, wrap_fasta
    from select_gap_alignments import GapAlignment, maximum_coverage_selection
    from refine_bundle_mappings import (
        chain_blocks, chain_to_mapping, deduplicate_same_locus, interval_total,
        merge_intervals, natural_name_key, normalize_bundle_name,
        parse_alias, read_existing_bundle_metadata, read_paf, resolve_alias, subtract_intervals,
    )


@dataclass(frozen=True)
class Observation:
    raw_name: str
    path: str
    strand: str
    genome_start: int
    genome_end: int
    aligned_intervals: Tuple[Tuple[int, int], ...]
    unmapped_intervals: Tuple[Tuple[int, int], ...]
    aligned_bp: int
    fraction: float
    source_lines: Tuple[int, ...]


def collect_mapping_context(paf, records, args):
    blocks = read_paf(paf)
    eligible = []
    aligned_by_line = {}
    genome_lengths = {}
    for block in blocks:
        previous_length = genome_lengths.setdefault(block.target_name, block.target_length)
        if previous_length != block.target_length:
            raise ValueError(f"PAF population path {block.target_name!r} has inconsistent lengths")
        if block.raw_query_name not in records:
            raise ValueError(f"PAF bundle {block.raw_query_name!r} is absent from consensus FASTA")
        if block.query_length != len(records[block.raw_query_name]):
            raise ValueError(f"PAF/FASTA length mismatch for {block.raw_query_name}")
        if block.mapq < args.min_mapq or block.identity < args.min_identity:
            continue
        # Only outer PAF bundle coordinates are used. Internal CIGAR deletions
        # are deliberately ignored, including for truncation-length statistics.
        aligned_by_line[block.line_number] = ((block.query_start, block.query_end),)
        eligible.append(block)
    chains = chain_blocks(
        eligible, args.max_genome_gap, args.max_bundle_gap,
        max_overlap_fraction=args.max_overlap_fraction,
        fixed_bundle_min_coverage=args.fixed_bundle_min_coverage,
        max_fixed_overlap_fraction=args.max_fixed_overlap_fraction,
        max_fixed_overlap_bp=args.max_fixed_overlap_bp,
        aliases=args.bundle_aliases,
    )
    # Raw names intentionally remain distinct: allele consensus coordinates
    # are not interchangeable. Deduplicate secondary hits within each raw name.
    mappings = [chain_to_mapping(chain, chain.raw_query_name, f"S{i:07d}")
                for i, chain in enumerate(chains, 1)]
    mappings, duplicates = deduplicate_same_locus(mappings, args.same_locus_overlap)
    observations = defaultdict(list)
    for mapping in mappings:
        aligned = merge_intervals(interval for line in mapping.source_line_numbers
                                  for interval in aligned_by_line[line])
        unmapped = subtract_intervals(((0, mapping.query_length),), aligned)
        aligned_bp = interval_total(aligned)
        observations[mapping.raw_query_name].append(Observation(
            mapping.raw_query_name, mapping.target_name, mapping.strand,
            mapping.target_start, mapping.target_end, aligned, unmapped,
            aligned_bp, aligned_bp / mapping.query_length, mapping.source_line_numbers,
        ))
    return {
        "observations": observations, "input_count": len(blocks), "eligible_count": len(eligible),
        "duplicate_count": len(duplicates), "eligible_blocks": eligible,
        "chains": chains, "genome_lengths": genome_lengths,
    }


def collect_observations(paf, records, args):
    """Compatibility wrapper returning observations and the three PAF counts."""
    context = collect_mapping_context(paf, records, args)
    return (context["observations"], context["input_count"], context["eligible_count"],
            context["duplicate_count"])


def preliminary_bed_rows(context, args):
    """Build final-filter-like BED6 rows, or read an explicitly supplied BED.

    BED spans enclose complete retained chains, including inter-block gaps.
    Class normalization, coverage filtering and same-locus deduplication match
    final refinement; raw allele coordinates remain separate in observations.
    """
    if args.preliminary_bed is not None:
        rows = []
        with args.preliminary_bed.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    raise ValueError(f"{args.preliminary_bed}:{line_number}: BED requires at least 3 fields")
                start, end = int(fields[1]), int(fields[2])
                if start < 0 or end <= start:
                    raise ValueError(f"{args.preliminary_bed}:{line_number}: invalid BED interval")
                if fields[0] in context["genome_lengths"] and end > context["genome_lengths"][fields[0]]:
                    raise ValueError(f"{args.preliminary_bed}:{line_number}: BED exceeds PAF population length")
                rows.append((fields[0], start, end, fields[3] if len(fields) > 3 else ".",
                             fields[4] if len(fields) > 4 else "0", fields[5] if len(fields) > 5 else "."))
        return rows
    retained = []
    for index, chain in enumerate(context["chains"], 1):
        name = normalize_bundle_name(chain.raw_query_name, args.bundle_aliases)
        mapping = chain_to_mapping(chain, name, f"P{index:07d}")
        if mapping.query_coverage >= args.preliminary_min_bundle_coverage and mapping.identity >= args.min_identity:
            retained.append(mapping)
    retained, _ = deduplicate_same_locus(retained, args.same_locus_overlap)
    return [(item.target_name, item.target_start, item.target_end, item.bundle_name,
             round(item.query_coverage * 1000), item.strand) for item in retained]


def collect_gap_evidence(context, args):
    """Select gap coverage first, then retain only safe incomplete-hit endpoints.

    Each candidate must cover >gap_min_mapping_bp within one preliminary-BED
    blank. Selection optimizes the union of PAF outer genome intervals with
    strictly limited pairwise overlap. No endpoint is interpolated from a
    clipped interval: only original PAF ends located on/in the blank qualify.
    A selected complete hit or hit spanning both blank boundaries supplies no
    endpoint and does not cause fallback to a lower-coverage alternative set.
    """
    bed_rows = preliminary_bed_rows(context, args)
    covered = defaultdict(list)
    for path, start, end, *_ in bed_rows:
        covered[path].append((start, end))
    blanks = {path: subtract_intervals(((0, length),), merge_intervals(covered.get(path, ())))
              for path, length in context["genome_lengths"].items()}
    by_path = defaultdict(list)
    observation_by_line = {}
    for name, observations in context["observations"].items():
        for index, item in enumerate(observations):
            for line in item.source_lines:
                observation_by_line[line] = (name, index, item)
    # Same-locus secondary rows removed from observations cannot inflate either
    # maximum-coverage selection or population endpoint support.
    for block in context["eligible_blocks"]:
        if block.line_number in observation_by_line:
            by_path[block.target_name].append(block)
    selected_intervals = defaultdict(list)
    endpoint_rows, candidate_rows, group_rows = [], [], []
    for path in sorted(blanks, key=natural_name_key):
        for gap_start, gap_end in blanks[path]:
            candidates = []
            for block in by_path[path]:
                start = max(gap_start, block.target_start)
                end = min(gap_end, block.target_end)
                if end - start > args.gap_min_mapping_bp:
                    candidates.append(GapAlignment(block=block, gap_start=gap_start, gap_end=gap_end,
                                                   start=start, end=end))
            selected = maximum_coverage_selection(candidates, max_overlap_bp=args.gap_max_overlap_bp)
            selected_lines = {item.block.line_number for item in selected}
            qualified = set()
            for candidate in sorted(candidates, key=lambda item: (item.start, item.end, item.block.line_number)):
                block = candidate.block
                name, index, observation = observation_by_line[block.line_number]
                is_selected = block.line_number in selected_lines
                status = "not_selected"
                safe_endpoints = []
                if is_selected:
                    if observation.fraction >= args.truncated_fraction:
                        status = "selected_not_short"
                    else:
                        # On '-' the left genome boundary is the bundle end.
                        pairs = ((block.target_start, block.query_start, "start"),
                                 (block.target_end, block.query_end, "end")) if block.strand == "+" else (
                                 (block.target_start, block.query_end, "end"),
                                 (block.target_end, block.query_start, "start"))
                        safe_endpoints = [(genome_point, bundle_point, side) for genome_point, bundle_point, side in pairs
                                          if gap_start <= genome_point <= gap_end and 0 < bundle_point < block.query_length]
                        status = "endpoint_supported" if safe_endpoints else "selected_no_original_endpoint_in_blank"
                if safe_endpoints:
                    qualified.add((name, index))
                    # Reconstruct only selected raw bundle intervals. They are
                    # used for diagnostics, never clipped by linear projection.
                    selected_intervals[(name, index)].append((block.query_start, block.query_end))
                    for genome_point, bundle_point, side in safe_endpoints:
                        endpoint_rows.append((name, path, block.strand, gap_start, gap_end, block.line_number,
                                              index, genome_point, bundle_point, side, candidate.end - candidate.start,
                                              observation.fraction))
                candidate_rows.append((path, gap_start, gap_end, block.line_number, name, block.strand,
                                       block.target_start, block.target_end, candidate.start, candidate.end,
                                       candidate.end - candidate.start, block.query_start, block.query_end,
                                       observation.fraction, int(is_selected), status))
            selected_coverage = interval_total(merge_intervals((item.start, item.end) for item in selected))
            group_rows.append((path, gap_start, gap_end, gap_end - gap_start, len(candidates), len(selected),
                               selected_coverage, selected_coverage / (gap_end - gap_start),
                               ",".join(map(str, sorted(selected_lines))), len(qualified)))
    endpoints_by_occurrence = defaultdict(set)
    for name, _, _, _, _, _, index, _, point, *_ in endpoint_rows:
        endpoints_by_occurrence[(name, index)].add(point)
    supported = defaultdict(list)
    for (name, index), positions in endpoints_by_occurrence.items():
        original = context["observations"][name][index]
        supported[name].append({"observation": original, "positions": tuple(sorted(positions)),
                                "selected_bundle_intervals": merge_intervals(selected_intervals[(name, index)])})
    return {"bed_rows": bed_rows, "blanks": blanks, "supported": supported,
            "endpoint_rows": endpoint_rows, "candidate_rows": candidate_rows, "group_rows": group_rows}


def protect_terminal_fragments(cuts, length, min_terminal_bp=1000):
    """Suppress terminal slivers without deleting any sequence.

    Candidate cuts are already PAF-evidence-supported. Reject cuts strictly
    within min_terminal_bp of either original consensus end. Every remaining
    internal cut is retained even if adjacent cuts delimit a short middle child.
    """
    if length <= 0 or min_terminal_bp < 0:
        raise ValueError("bundle length must be positive and terminal threshold non-negative")
    positions = sorted(set(cuts))
    if any(point <= 0 or point >= length for point in positions):
        raise ValueError("cuts must be strictly inside the bundle")
    accepted = tuple(point for point in positions
                     if min_terminal_bp <= point <= length - min_terminal_bp)
    accepted_set = set(accepted)
    blocked = tuple(point for point in positions if point not in accepted_set)
    return accepted, blocked


def supported_cut_plan(supported_short, length, args):
    """Make one joint plan from BED-gap-verified incomplete occurrences.

    Count each deduplicated occurrence once per endpoint cluster. Uncovered
    consensus intervals alone are diagnostic and never mandatory endpoints.
    Cluster medians must be actual observed positions, not invented coordinates.
    """
    endpoints = []
    for index, item in enumerate(supported_short):
        endpoints.extend((point, index) for point in item["positions"] if 0 < point < length)
    clusters = []
    for endpoint in sorted(endpoints):
        if not clusters or endpoint[0] - clusters[-1][0][0] > args.cut_merge_bp:
            clusters.append([endpoint])
        else:
            clusters[-1].append(endpoint)
    evidence = []
    cuts = []
    for cluster in clusters:
        per_occurrence = defaultdict(list)
        for point, occurrence in cluster:
            per_occurrence[occurrence].append(point)
        medians = sorted(values[len(values) // 2] for values in per_occurrence.values())
        support = len(per_occurrence)
        positions = [medians[len(medians) // 2]]
        selected = len(supported_short) >= args.min_truncated_count and support >= args.min_cut_support
        if selected:
            cuts.extend(positions)
        evidence.append({"cluster_start": cluster[0][0], "cluster_end": cluster[-1][0],
                         "positions": tuple(positions), "short_occurrence_support": support,
                         "short_path_support": len({supported_short[index]["observation"].path for index in per_occurrence}),
                         "global_long_gap_boundary": False, "selected": selected,
                         "selected_before_terminal_guard": selected})
    cuts, blocked = protect_terminal_fragments(
        cuts, length, getattr(args, "min_terminal_fragment_bp", 1000))
    for item in evidence:
        suppressed = item["selected"] and any(point in blocked for point in item["positions"])
        item["terminal_fragment_suppressed"] = suppressed
        if suppressed:
            item["selected"] = False
    return cuts, evidence


def plan_splits(records, observations, args, gap_evidence):
    plans = {}
    for name, sequence in records.items():
        length = len(sequence)
        seen = observations.get(name, [])
        short = [item for item in seen if item.fraction < args.truncated_fraction]
        recurrent = len(short) >= args.min_truncated_count
        global_aligned = merge_intervals(interval for item in seen for interval in item.aligned_intervals)
        global_unmapped = subtract_intervals(((0, length),), global_aligned)
        long_gaps = [interval for interval in global_unmapped
                     if interval[1] - interval[0] > args.max_unmapped_bp]
        supported_short = gap_evidence["supported"].get(name, [])
        cuts, evidence = supported_cut_plan(supported_short, length, args)
        plans[name] = {
            "occurrence_count": len(seen), "truncated_count": len(short),
            "gap_supported_truncated_count": len(supported_short),
            "gap_supported_recurrent_truncation": len(supported_short) >= args.min_truncated_count,
            "recurrent_truncation": recurrent, "long_unmapped_count": len(long_gaps),
            "max_unmapped_bp": max((end - start for start, end in global_unmapped), default=0),
            "population_path_count": len({item.path for item in seen}),
            "global_aligned_bp": interval_total(global_aligned),
            "global_coverage": interval_total(global_aligned) / length,
            "global_aligned_intervals": global_aligned, "global_unmapped_intervals": global_unmapped,
            "cut_positions": cuts, "cut_evidence": evidence, "children": [],
            "candidate_cut_positions": tuple(sorted({point for item in evidence
                if item["selected_before_terminal_guard"] for point in item["positions"]})),
            "suppressed_terminal_cuts": tuple(sorted({point for item in evidence
                if item["terminal_fragment_suppressed"] for point in item["positions"]})),
        }
    return plans


def interval_string(intervals):
    return ",".join(f"{start}-{end}" for start, end in intervals)


def write_outputs(records, observations, plans, args, gap_evidence):
    aliases = args.bundle_aliases
    existing_ids, metadata = read_existing_bundle_metadata(args.bundle_metadata_bed)
    used_ids = set(existing_ids)
    classes = {name: normalize_bundle_name(name, aliases) for name in records}
    used_ids.update(value for value in classes.values() if value.isdigit())
    next_id = max((int(value) for value in used_ids), default=-1) + 1
    metadata_sizes = {key: value.bundle_size for key, value in metadata.items()}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    preliminary_output = args.output_dir / "preliminary_bundles.bed"
    # A caller may deliberately reuse this artifact as --preliminary-bed.
    # Preserve that input's exact bytes/mtime instead of rewriting it in place.
    if args.preliminary_bed is None or args.preliminary_bed.resolve() != preliminary_output.resolve():
        with preliminary_output.open("w", encoding="utf-8") as handle:
            handle.write("# Preliminary retained mapping BED6; score=rounded bundle coverage * 1000\n")
            for row in sorted(gap_evidence["bed_rows"], key=lambda item: (natural_name_key(item[0]), item[1], item[2], item[3])):
                handle.write("\t".join(map(str, row)) + "\n")
    with (args.output_dir / "uncovered_genome_regions.bed").open("w", encoding="utf-8") as handle:
        handle.write("# Genome blanks from preliminary BED; lengths available only for PAF-observed paths\n")
        for path in sorted(gap_evidence["blanks"], key=natural_name_key):
            for start, end in gap_evidence["blanks"][path]:
                handle.write(f"{path}\t{start}\t{end}\n")
    audits = (
        ("gap_alignment_candidates.tsv", ("population_path", "gap_start", "gap_end", "paf_line",
         "raw_consensus", "strand", "raw_genome_start", "raw_genome_end", "clipped_genome_start",
         "clipped_genome_end", "mapped_bp_in_blank", "raw_bundle_start", "raw_bundle_end",
         "occurrence_bundle_coverage", "selected", "status"), gap_evidence["candidate_rows"]),
        ("gap_alignment_groups.tsv", ("population_path", "gap_start", "gap_end", "gap_length",
         "candidate_count", "selected_count", "selected_union_bp", "selected_blank_coverage",
         "selected_paf_lines", "qualified_occurrence_count"), gap_evidence["group_rows"]),
        ("potential_bundle_endpoints.tsv", ("raw_consensus", "population_path", "strand", "gap_start",
         "gap_end", "paf_line", "occurrence_index", "genome_endpoint", "bundle_endpoint",
         "bundle_endpoint_side", "mapped_bp_in_blank", "occurrence_bundle_coverage"), gap_evidence["endpoint_rows"]),
    )
    for filename, columns, rows in audits:
        with (args.output_dir / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(rows)
    fasta_path = args.output_dir / "bundle.fa"
    fragment_rows = []
    with fasta_path.open("w", encoding="utf-8") as fasta:
        for name in sorted(records, key=natural_name_key):
            sequence = records[name]
            cuts = plans[name]["cut_positions"]
            boundaries = (0, *cuts, len(sequence))
            parent_size = metadata_sizes.get(classes[name])
            for part, (start, end) in enumerate(zip(boundaries, boundaries[1:]), 1):
                if cuts:
                    child_id = str(next_id)
                    next_id += 1
                    child_name = f"bundle_{child_id}"
                    # Native minimizer coordinates cannot be recovered.
                    size = max(1, math.ceil(parent_size * (end - start) / len(sequence))) if parent_size else 100
                    metadata_sizes[child_id] = size
                else:
                    child_name = name
                    child_id = classes[name]
                    size = parent_size or ""
                plans[name]["children"].append(child_name)
                fasta.write(f">{child_name} parent_consensus={name} parent_interval={start}-{end}\n")
                fasta.writelines(wrap_fasta(sequence[start:end], 80))
                fragment_rows.append((name, classes[name], len(sequence), part,
                                      start, end, child_name, child_id, size, bool(cuts)))
    with (args.output_dir / "bundle_metadata.bed").open("w", encoding="utf-8") as handle:
        handle.write("# ID/size metadata only; synthetic coordinates, NOT genomic mappings\n")
        for bundle_id in sorted(metadata_sizes, key=int):
            size = metadata_sizes[bundle_id]
            handle.write(f"__metadata__\t0\t1\t{bundle_id}:{size}:0:0:{size - 1}:U\n")
    with (args.output_dir / "bundle_fragments.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("parent_consensus", "parent_class", "parent_length", "part",
                         "parent_start", "parent_end", "consensus_name", "numeric_bundle_id",
                         "bundle_size", "is_split"))
        writer.writerows((*row[:-1], int(row[-1])) for row in fragment_rows)
    with (args.output_dir / "bundle_split_stats.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("raw_consensus", "bundle_class", "bundle_length", "occurrence_count",
                         "truncated_count", "recurrent_truncation", "long_unmapped_count",
                         "max_unmapped_bp", "cut_positions", "fragment_names", "status",
                         "population_path_count", "global_aligned_bp", "global_coverage",
                         "global_aligned_intervals", "global_unmapped_intervals",
                         "gap_supported_truncated_count", "gap_supported_recurrent_truncation",
                         "candidate_cut_positions", "suppressed_terminal_cuts", "min_terminal_fragment_bp"))
        for name in sorted(records, key=natural_name_key):
            plan = plans[name]
            writer.writerow((name, classes[name], len(records[name]), plan["occurrence_count"],
                             plan["truncated_count"], int(plan["recurrent_truncation"]),
                             plan["long_unmapped_count"], plan["max_unmapped_bp"],
                             ",".join(map(str, plan["cut_positions"])), ",".join(plan["children"]),
                             "split" if plan["cut_positions"] else
                             "no_eligible_hit" if not plan["occurrence_count"] else "retained",
                             plan["population_path_count"], plan["global_aligned_bp"],
                             f"{plan['global_coverage']:.8g}",
                             interval_string(plan["global_aligned_intervals"]),
                             interval_string(plan["global_unmapped_intervals"]),
                             plan["gap_supported_truncated_count"], int(plan["gap_supported_recurrent_truncation"]),
                             ",".join(map(str, plan["candidate_cut_positions"])),
                             ",".join(map(str, plan["suppressed_terminal_cuts"])), args.min_terminal_fragment_bp))
    with (args.output_dir / "bundle_cut_support.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("raw_consensus", "cluster_start", "cluster_end", "candidate_positions",
                         "short_occurrence_support", "short_path_support", "global_long_gap_boundary",
                         "selected", "selected_before_terminal_guard", "terminal_fragment_suppressed"))
        for name in sorted(records, key=natural_name_key):
            for item in plans[name]["cut_evidence"]:
                writer.writerow((name, item["cluster_start"], item["cluster_end"],
                                 ",".join(map(str, item["positions"])), item["short_occurrence_support"],
                                 item["short_path_support"], int(item["global_long_gap_boundary"]),
                                 int(item["selected"]), int(item["selected_before_terminal_guard"]),
                                 int(item["terminal_fragment_suppressed"])))
    with (args.output_dir / "mapping_observations.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("raw_consensus", "population_path", "strand", "genome_start", "genome_end",
                         "bundle_aligned_bp", "bundle_aligned_fraction", "aligned_bundle_intervals",
                         "unmapped_bundle_intervals", "source_line_numbers"))
        for name in sorted(observations, key=natural_name_key):
            for item in observations[name]:
                writer.writerow((name, item.path, item.strand, item.genome_start, item.genome_end,
                                 item.aligned_bp, f"{item.fraction:.8g}",
                                 interval_string(item.aligned_intervals), interval_string(item.unmapped_intervals),
                                 ",".join(map(str, item.source_lines))))


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle_fasta", type=Path)
    parser.add_argument("paf", type=Path)
    parser.add_argument("--bundle-metadata-bed", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-truncated-count", type=int, default=5)
    parser.add_argument("--truncated-fraction", type=float, default=0.60)
    parser.add_argument("--max-unmapped-bp", type=int, default=10_000,
                        help="diagnostic threshold for long consensus-uncovered intervals; no automatic cuts")
    parser.add_argument("--preliminary-bed", type=Path,
                        help="optional existing preliminary retained mapping BED; empty BED is allowed")
    parser.add_argument("--preliminary-min-bundle-coverage", type=float, default=0.60,
                        help="coverage filter when generating the preliminary mapping BED from PAF")
    parser.add_argument("--gap-min-mapping-bp", type=int, default=10_000,
                        help="each candidate's PAF outer genome span intersected with one BED blank must be strictly longer")
    parser.add_argument("--gap-max-overlap-bp", type=int, default=1_000,
                        help="selected clipped spans must overlap by strictly less than this; zero means no overlap")
    parser.add_argument("--cut-merge-bp", type=int, default=100)
    parser.add_argument("--min-terminal-fragment-bp", type=int, default=1000,
                        help="suppress cuts creating a terminal child shorter than this; short internal children remain; 0 disables")
    parser.add_argument("--min-cut-support", type=int, default=1,
                        help="minimum distinct short occurrences per endpoint cluster; bundle count gate still applies")
    parser.add_argument("--min-identity", type=float, default=0.0)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--max-bundle-gap", type=int, default=1_000)
    parser.add_argument("--max-genome-gap", type=int, default=1_000)
    parser.add_argument("--max-overlap-fraction", type=float, default=0.10,
                        help="maximum overlap on either axis divided by the original consensus length")
    parser.add_argument("--fixed-bundle-min-coverage", type=float, default=0.80,
                        help="single-block bundle coverage required to fix it before chaining")
    parser.add_argument("--max-fixed-overlap-fraction", type=float, default=0.20,
                        help="reject candidate span overlap above this fraction of a fixed consensus length")
    parser.add_argument("--max-fixed-overlap-bp", type=int, default=10_000,
                        help="reject candidate span overlap above this many bp of any fixed occurrence")
    parser.add_argument("--bundle-alias", action="append", type=parse_alias, default=[],
                        metavar="OLD=NEW", help="repeatable class alias; coordinates remain raw-consensus-specific")
    parser.add_argument("--same-locus-overlap", type=float, default=0.80)
    return parser


def main(argv: Optional[Sequence[str]] = None):
    parser = build_parser()
    args = parser.parse_args(argv)
    for name in ("truncated_fraction", "min_identity", "same_locus_overlap",
                 "max_overlap_fraction", "fixed_bundle_min_coverage", "max_fixed_overlap_fraction",
                 "preliminary_min_bundle_coverage"):
        if not 0 <= getattr(args, name) <= 1:
            parser.error(f"--{name.replace('_', '-')} must be in [0, 1]")
    for name in ("min_truncated_count", "min_cut_support", "max_unmapped_bp", "cut_merge_bp", "min_mapq",
                 "max_bundle_gap", "max_genome_gap", "max_fixed_overlap_bp",
                 "gap_min_mapping_bp", "gap_max_overlap_bp", "min_terminal_fragment_bp"):
        if getattr(args, name) < (1 if name in {"min_truncated_count", "min_cut_support"} else 0):
            parser.error(f"invalid --{name.replace('_', '-')}")
    for path in (args.bundle_fasta, args.paf, args.bundle_metadata_bed):
        if not path.is_file():
            parser.error(f"input is missing: {path}")
    if args.preliminary_bed is not None and not args.preliminary_bed.is_file():
        parser.error(f"preliminary BED is missing: {args.preliminary_bed}")
    if args.gap_max_overlap_bp > args.gap_min_mapping_bp:
        parser.error("--gap-max-overlap-bp must not exceed --gap-min-mapping-bp")
    records: Dict[str, str] = {}
    for header, sequence in iter_fasta(args.bundle_fasta):
        name = header.split()[0]
        if name in records or not sequence:
            raise ValueError(f"duplicate name or empty sequence: {name}")
        if set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"consensus {name} contains invalid or gapped bases")
        records[name] = sequence
    if not records:
        raise ValueError("empty bundle consensus FASTA")
    args.bundle_aliases = {"dup1_1": "dup1", "dup1_2": "dup1"}
    args.bundle_aliases.update(dict(args.bundle_alias))
    for alias_name in args.bundle_aliases:
        resolve_alias(alias_name, args.bundle_aliases)
    context = collect_mapping_context(args.paf, records, args)
    observations = context["observations"]
    gap_evidence = collect_gap_evidence(context, args)
    plans = plan_splits(records, observations, args, gap_evidence)
    write_outputs(records, observations, plans, args, gap_evidence)
    summary = {
        "inputs": {"bundle_fasta": str(args.bundle_fasta.resolve()), "paf": str(args.paf.resolve()),
                   "bundle_metadata_bed": str(args.bundle_metadata_bed.resolve()),
                   "preliminary_bed": str(args.preliminary_bed.resolve()) if args.preliminary_bed else None},
        "parameters": {key: getattr(args, key) for key in (
            "min_truncated_count", "truncated_fraction", "max_unmapped_bp", "cut_merge_bp", "min_cut_support",
            "min_identity", "min_mapq", "max_bundle_gap", "max_genome_gap", "same_locus_overlap",
            "max_overlap_fraction", "fixed_bundle_min_coverage", "max_fixed_overlap_fraction",
            "max_fixed_overlap_bp", "preliminary_min_bundle_coverage", "gap_min_mapping_bp", "gap_max_overlap_bp",
            "min_terminal_fragment_bp")},
        "counts": {"input_paf_blocks": context["input_count"], "eligible_paf_blocks": context["eligible_count"],
                   "same_locus_duplicates_removed": context["duplicate_count"], "input_consensus_count": len(records),
                   "preliminary_bed_intervals": len(gap_evidence["bed_rows"]),
                   "paf_observed_population_paths": len(context["genome_lengths"]),
                   "genome_blank_regions": sum(len(items) for items in gap_evidence["blanks"].values()),
                   "gap_alignment_candidates": len(gap_evidence["candidate_rows"]),
                   "selected_gap_alignments": sum(row[14] for row in gap_evidence["candidate_rows"]),
                   "potential_endpoint_observations": sum(len(items) for items in gap_evidence["supported"].values()),
                   "potential_bundle_endpoints": len(gap_evidence["endpoint_rows"]),
                   "split_consensus_count": sum(bool(plan["cut_positions"]) for plan in plans.values()),
                   "output_consensus_count": sum(len(plan["children"]) for plan in plans.values()),
                   "suppressed_terminal_cuts": sum(len(plan["suppressed_terminal_cuts"]) for plan in plans.values())},
        "semantics": {"coordinate_system": "original raw consensus, zero-based half-open",
                      "truncation_test": "aligned bundle bp < original length * truncated_fraction",
                      "unmapped_test": "population-wide uncovered consensus intervals are diagnostic only",
                      "coverage_filter_applied_before_detection": False,
                      "preliminary_bed_coverage_filter_applied": args.preliminary_bed is None,
                      "gap_selection_objective": "maximum union of clipped PAF outer genome intervals per BED blank",
                      "gap_mapping_length_test": "clipped span within one blank strictly > gap_min_mapping_bp",
                      "gap_selected_overlap_test": "pairwise overlap strictly < gap_max_overlap_bp; zero means disjoint",
                      "potential_endpoint_test": "selected deduplicated short occurrence; original PAF end on/in BED blank",
                      "raw_truncation_count_is_diagnostic": True,
                      "recurrent_gate_uses_gap_supported_occurrences": True,
                      "no_cuts_without_paf_gap_evidence": True,
                      "cigar_deletions_count_as_unmapped_bundle": False,
                      "short_cut_minimum_support": args.min_cut_support,
                      "short_cut_support_unit": "deduplicated occurrence, once per endpoint cluster",
                      "one_joint_cut_plan_per_consensus": True,
                      "all_fragment_sequences_retained": True,
                      "trim_is_non_destructive_supported_splitting": True,
                      "terminal_fragment_test": "suppress cuts creating terminal children strictly shorter than min_terminal_fragment_bp",
                      "short_terminal_sequence_fate": "retain within adjacent child by suppressing the cut; never delete",
                      "short_internal_fragments_retained": True,
                      "split_endpoints_use_cigar": False,
                      "different_raw_alleles_split_independently": True,
                      "fixed_occurrences_are_single_eligible_paf_blocks": True,
                      "candidate_whole_genome_span_checks_other_classes": True,
                      "fixed_overlap_fraction_denominator": "fixed original consensus length",
                      "bundle_aliases": args.bundle_aliases,
                      "remap_required_after_split": True},
        "outputs": {name: str((args.output_dir / name).resolve()) for name in (
            "bundle.fa", "bundle_metadata.bed", "bundle_fragments.tsv", "bundle_split_stats.tsv",
            "bundle_cut_support.tsv", "mapping_observations.tsv", "preliminary_bundles.bed",
            "uncovered_genome_regions.bed", "gap_alignment_candidates.tsv", "gap_alignment_groups.tsv",
            "potential_bundle_endpoints.tsv", "mapping_split_summary.json")},
    }
    with (args.output_dir / "mapping_split_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
