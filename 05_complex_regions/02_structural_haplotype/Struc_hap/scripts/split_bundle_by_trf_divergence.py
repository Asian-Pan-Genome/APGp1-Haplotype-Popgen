#!/usr/bin/env python3
"""Recursively remove long, population-variable TRFs from a bundle MSA."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from urllib.parse import unquote


@dataclass(frozen=True)
class AlignedOccurrence:
    header: str
    bundle_id: str
    path_id: str
    start: int
    end: int
    direction: str
    bed_row: Optional[int]
    sequence: str

    @property
    def ungapped_length(self) -> int:
        return sum(base != "-" for base in self.sequence)


@dataclass(frozen=True)
class TRFHit:
    path_id: str
    repeat_id: str
    start: int
    end: int
    period: str
    copy_number: str
    identity: float
    overlap_start: int
    overlap_end: int

    @property
    def overlap_bp(self) -> int:
        return self.overlap_end - self.overlap_start


@dataclass(frozen=True)
class RepeatUnion:
    path_id: str
    start: int
    end: int
    source_repeat_ids: Tuple[str, ...]

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class ProjectedRepeatMember:
    repeat: RepeatUnion
    start_offset: int
    end_offset: int
    start_column: int
    end_column: int


@dataclass(frozen=True)
class PopulationRepeatRegion:
    start_column: int
    end_column: int
    members: Tuple[ProjectedRepeatMember, ...]

    @property
    def source_path_ids(self) -> Tuple[str, ...]:
        return tuple(sorted({member.repeat.path_id for member in self.members}))

    @property
    def source_repeat_ids(self) -> Tuple[str, ...]:
        return tuple(sorted({repeat_id for member in self.members
                             for repeat_id in member.repeat.source_repeat_ids}))


@dataclass(frozen=True)
class BedOccurrence:
    line_number: int
    path_id: str
    start: int
    end: int
    bundle_id: str
    direction: str


@dataclass
class CandidateEvaluation:
    candidate_id: int
    repeat: RepeatUnion
    reference_start: int
    reference_end: int
    start_offset: int
    end_offset: int
    start_column: int
    end_column: int
    population_region: PopulationRepeatRegion
    metrics: Dict[str, object]
    structural_divergent: bool
    length_divergent: bool
    nucleotide_divergent: bool
    divergent: bool
    selected_iteration: Optional[int] = None

    @property
    def repeat_bp(self) -> int:
        return int(self.metrics["reference_repeat_bp"])


@dataclass(frozen=True)
class ColumnPrefixStats:
    gap: Tuple[int, ...]
    canonical_calls: Tuple[int, ...]
    disagreements: Tuple[int, ...]
    informative: Tuple[int, ...]
    variable: Tuple[int, ...]


@dataclass(frozen=True)
class RetainedSegment:
    start_column: int
    end_column: int
    supporting_sequence_count: int
    nonempty_sequence_count: int
    median_sequence_bp: float
    max_sequence_bp: int
    depth: int


@dataclass(frozen=True)
class RecursiveDecision:
    iteration: int
    depth: int
    parent_start_column: int
    parent_end_column: int
    candidate_id: int
    repeat_path_id: str
    repeat_start: int
    repeat_end: int
    repeat_bp: int
    repeat_start_column: int
    repeat_end_column: int
    child_side: str
    child_start_column: int
    child_end_column: int
    supporting_sequence_count: int
    nonempty_sequence_count: int
    required_sequence_count: int
    retained: bool


@dataclass(frozen=True)
class SegmentProjection:
    path_id: str
    direction: str
    parent_start: int
    parent_end: int
    segment_index: int
    new_bundle_id: str
    alignment_start_column: int
    alignment_end_column: int
    ungapped_start_offset: int
    ungapped_end_offset: int
    genomic_start: int
    genomic_end: int


def parse_alignment_header(header: str) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for field in header.split("|"):
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = unquote(value)
    required = {"bundle", "path", "start", "end", "strand"}
    missing = sorted(required - values.keys())
    if missing:
        raise ValueError(f"alignment header lacks {', '.join(missing)}: {header}")
    return values


def read_alignment(path: Path) -> List[AlignedOccurrence]:
    records: List[AlignedOccurrence] = []
    header: Optional[str] = None
    chunks: List[str] = []

    def append_record() -> None:
        if header is None:
            return
        values = parse_alignment_header(header)
        sequence = "".join(chunks).upper()
        record = AlignedOccurrence(
            header=header,
            bundle_id=values["bundle"],
            path_id=values["path"],
            start=int(values["start"]),
            end=int(values["end"]),
            direction=values["strand"],
            bed_row=int(values["bed_row"]) if values.get("bed_row", "").isdigit() else None,
            sequence=sequence,
        )
        if record.direction not in {"+", "-"}:
            raise ValueError(f"unsupported alignment strand {record.direction!r}")
        if record.ungapped_length != record.end - record.start:
            raise ValueError(
                f"{record.path_id}: ungapped alignment length {record.ungapped_length} "
                f"does not match header interval {record.start}-{record.end}"
            )
        records.append(record)

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                append_record()
                header = stripped[1:]
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"{path}:{line_number}: sequence before first header")
                chunks.append(stripped)
    append_record()
    if not records:
        raise ValueError(f"{path}: empty alignment")
    alignment_lengths = {len(record.sequence) for record in records}
    if len(alignment_lengths) != 1:
        raise ValueError(f"{path}: sequences do not have one alignment length")
    bundle_ids = {record.bundle_id for record in records}
    if len(bundle_ids) != 1:
        raise ValueError(f"{path}: alignment contains multiple bundle IDs")
    path_counts = Counter(record.path_id for record in records)
    duplicated_paths = sorted(path_id for path_id, count in path_counts.items() if count > 1)
    if duplicated_paths:
        raise ValueError(
            "this splitter currently requires at most one target-bundle occurrence per path; "
            f"duplicates include {', '.join(duplicated_paths[:5])}"
        )
    return records


def read_trf_hits(
    path: Path,
    occurrences: Mapping[str, AlignedOccurrence],
    min_length: int,
    min_identity: float,
) -> List[TRFHit]:
    hits: List[TRFHit] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 7:
                fields = stripped.split()
            if len(fields) < 7:
                raise ValueError(f"{path}:{line_number}: expected at least seven TRF columns")
            occurrence = occurrences.get(fields[0])
            if occurrence is None:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
                identity = float(fields[6])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid TRF coordinates/identity"
                ) from exc
            if end - start < min_length or identity < min_identity:
                continue
            overlap_start = max(start, occurrence.start)
            overlap_end = min(end, occurrence.end)
            if overlap_end <= overlap_start:
                continue
            hits.append(
                TRFHit(
                    path_id=occurrence.path_id,
                    repeat_id=fields[3],
                    start=start,
                    end=end,
                    period=fields[4] if len(fields) > 4 else ".",
                    copy_number=fields[5] if len(fields) > 5 else ".",
                    identity=identity,
                    overlap_start=overlap_start,
                    overlap_end=overlap_end,
                )
            )
    return hits


def merge_trf_hits(hits: Sequence[TRFHit]) -> List[RepeatUnion]:
    by_path: Dict[str, List[TRFHit]] = defaultdict(list)
    for hit in hits:
        by_path[hit.path_id].append(hit)
    unions: List[RepeatUnion] = []
    for path_id, path_hits in sorted(by_path.items()):
        ordered = sorted(path_hits, key=lambda hit: (hit.overlap_start, hit.overlap_end))
        start = ordered[0].overlap_start
        end = ordered[0].overlap_end
        ids = {ordered[0].repeat_id}
        for hit in ordered[1:]:
            if hit.overlap_start <= end:
                end = max(end, hit.overlap_end)
                ids.add(hit.repeat_id)
            else:
                unions.append(RepeatUnion(path_id, start, end, tuple(sorted(ids))))
                start, end, ids = hit.overlap_start, hit.overlap_end, {hit.repeat_id}
        unions.append(RepeatUnion(path_id, start, end, tuple(sorted(ids))))
    return unions


def ungapped_interval_to_alignment_columns(
    sequence: str,
    start_offset: int,
    end_offset: int,
) -> Tuple[int, int]:
    if not 0 <= start_offset < end_offset <= sum(base != "-" for base in sequence):
        raise ValueError(f"invalid ungapped interval {start_offset}-{end_offset}")
    seen = 0
    start_column: Optional[int] = None
    for column, base in enumerate(sequence):
        if base == "-":
            continue
        if seen == start_offset and start_column is None:
            start_column = column
        seen += 1
        if seen == end_offset:
            if start_column is None:
                raise ValueError("failed to find alignment start column")
            return start_column, column + 1
    raise ValueError("failed to map ungapped interval to alignment columns")


def repeat_alignment_columns(
    repeat: RepeatUnion,
    reference: AlignedOccurrence,
) -> Tuple[int, int, int, int]:
    if reference.direction == "+":
        start_offset = repeat.start - reference.start
        end_offset = repeat.end - reference.start
    else:
        start_offset = reference.end - repeat.end
        end_offset = reference.end - repeat.start
    start_column, end_column = ungapped_interval_to_alignment_columns(
        reference.sequence,
        start_offset,
        end_offset,
    )
    return start_offset, end_offset, start_column, end_column


def population_repeat_regions(
    unions: Sequence[RepeatUnion],
    aligned_by_path: Mapping[str, AlignedOccurrence],
) -> List[PopulationRepeatRegion]:
    """Project all population TRFs to MSA columns and union their positions."""

    projected: List[ProjectedRepeatMember] = []
    for repeat in unions:
        reference = aligned_by_path[repeat.path_id]
        start_offset, end_offset, start_column, end_column = repeat_alignment_columns(
            repeat, reference
        )
        projected.append(
            ProjectedRepeatMember(
                repeat, start_offset, end_offset, start_column, end_column
            )
        )
    if not projected:
        return []
    projected.sort(
        key=lambda item: (
            item.start_column, item.end_column, item.repeat.path_id,
            item.repeat.start, item.repeat.end,
        )
    )
    regions: List[PopulationRepeatRegion] = []
    members = [projected[0]]
    start_column = projected[0].start_column
    end_column = projected[0].end_column
    for member in projected[1:]:
        if member.start_column <= end_column:
            end_column = max(end_column, member.end_column)
            members.append(member)
        else:
            regions.append(
                PopulationRepeatRegion(start_column, end_column, tuple(members))
            )
            start_column, end_column, members = (
                member.start_column, member.end_column, [member]
            )
    regions.append(PopulationRepeatRegion(start_column, end_column, tuple(members)))
    return regions


def ungapped_bases_before(sequence: str, column: int) -> int:
    return column - sequence.count("-", 0, column)


def build_column_prefix_stats(
    records: Sequence[AlignedOccurrence],
) -> ColumnPrefixStats:
    """Precompute alignment-column statistics once for fast repeat scans."""
    gaps = [0]
    canonical_calls = [0]
    disagreements = [0]
    informative = [0]
    variable = [0]
    alignment_length = len(records[0].sequence)
    for column in range(alignment_length):
        counts = Counter(record.sequence[column] for record in records)
        canonical = {base: counts[base] for base in "ACGT" if counts[base]}
        calls = sum(canonical.values())
        gaps.append(gaps[-1] + counts["-"])
        canonical_calls.append(canonical_calls[-1] + calls)
        disagreements.append(
            disagreements[-1] + (calls - max(canonical.values()) if canonical else 0)
        )
        informative.append(informative[-1] + bool(canonical))
        variable.append(variable[-1] + (len(canonical) > 1))
    return ColumnPrefixStats(
        gap=tuple(gaps),
        canonical_calls=tuple(canonical_calls),
        disagreements=tuple(disagreements),
        informative=tuple(informative),
        variable=tuple(variable),
    )


def prefix_difference(values: Sequence[int], start: int, end: int) -> int:
    return values[end] - values[start]


def mean_pairwise_length_divergence(lengths: Sequence[int]) -> float:
    """Mean pairwise relative TR-length difference across population rows.

    Each row contributes once.  A pair is scored as ``abs(a-b)/max(a,b)``;
    two zero-length rows score zero.  This makes the metric population-aware:
    a rare length allele contributes in proportion to its population frequency,
    rather than letting one extreme row define the whole region.
    """
    if len(lengths) < 2:
        return 0.0
    difference = 0.0
    pair_count = 0
    for index, left in enumerate(lengths):
        for right in lengths[index + 1:]:
            denominator = max(left, right)
            difference += abs(left - right) / denominator if denominator else 0.0
            pair_count += 1
    return difference / pair_count


def analyze_region(
    records: Sequence[AlignedOccurrence],
    start_column: int,
    end_column: int,
    reference_repeat_bp: int,
    min_sequence_coverage: float,
    column_stats: Optional[ColumnPrefixStats] = None,
) -> Dict[str, object]:
    if column_stats is None:
        column_stats = build_column_prefix_stats(records)
    alignment_span = end_column - start_column
    sequence_bp = [
        alignment_span - record.sequence.count("-", start_column, end_column)
        for record in records
    ]
    coverage = [value / reference_repeat_bp for value in sequence_bp]
    qualifying = sum(value >= min_sequence_coverage for value in coverage)
    nonzero = sum(value > 0 for value in sequence_bp)
    gap_count = prefix_difference(column_stats.gap, start_column, end_column)
    non_gap_count = prefix_difference(
        column_stats.canonical_calls, start_column, end_column
    )
    disagreement_count = prefix_difference(
        column_stats.disagreements, start_column, end_column
    )
    variable_columns = prefix_difference(
        column_stats.variable, start_column, end_column
    )
    informative_columns = prefix_difference(
        column_stats.informative, start_column, end_column
    )
    ordered_bp = sorted(sequence_bp)
    return {
        "alignment_start_column_0based": start_column,
        "alignment_end_column_0based_exclusive": end_column,
        "alignment_start_position_1based": start_column + 1,
        "alignment_end_position_1based_inclusive": end_column,
        "alignment_span_columns": alignment_span,
        "reference_repeat_bp": reference_repeat_bp,
        "sequence_count": len(records),
        "sequence_with_any_aligned_bp_count": nonzero,
        "sequence_with_any_aligned_bp_fraction": nonzero / len(records),
        "min_sequence_coverage": min_sequence_coverage,
        "sequence_meeting_coverage_count": qualifying,
        "population_presence_fraction": qualifying / len(records),
        "structural_divergence": 1 - qualifying / len(records),
        "gap_fraction": gap_count / (alignment_span * len(records)),
        "non_gap_consensus_divergence": (
            disagreement_count / non_gap_count if non_gap_count else 0.0
        ),
        "variable_column_fraction": (
            variable_columns / informative_columns if informative_columns else 0.0
        ),
        "tr_length_allele_count": len(set(sequence_bp)),
        "mean_pairwise_tr_length_divergence": mean_pairwise_length_divergence(
            sequence_bp
        ),
        "min_sequence_bp": ordered_bp[0],
        "median_sequence_bp": float(statistics.median(ordered_bp)),
        "max_sequence_bp": ordered_bp[-1],
    }


def normalize_bed_direction(value: str) -> str:
    if value in {"0", "+"}:
        return "+"
    if value in {"1", "-"}:
        return "-"
    raise ValueError(f"unsupported PGR direction {value!r}")


def read_target_bed(
    path: Path,
    target_bundle_id: str,
) -> Tuple[List[str], Dict[str, BedOccurrence], List[str], int]:
    lines = path.read_text(encoding="utf-8").splitlines()
    occurrences: Dict[str, BedOccurrence] = {}
    all_bundle_ids: List[str] = []
    for line_number, line in enumerate(lines, 1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"{path}:{line_number}: expected native four-column PGR BED")
        annotation = fields[3].split(":")
        if len(annotation) < 3:
            raise ValueError(f"{path}:{line_number}: invalid PGR annotation")
        bundle_id = annotation[0]
        all_bundle_ids.append(bundle_id)
        if bundle_id != target_bundle_id:
            continue
        if fields[0] in occurrences:
            raise ValueError(
                f"{path}: more than one bundle {target_bundle_id} occurrence on {fields[0]}"
            )
        occurrences[fields[0]] = BedOccurrence(
            line_number=line_number,
            path_id=fields[0],
            start=int(fields[1]),
            end=int(fields[2]),
            bundle_id=bundle_id,
            direction=normalize_bed_direction(annotation[2]),
        )
    numeric_ids = [int(bundle_id) for bundle_id in all_bundle_ids if bundle_id.isdigit()]
    return lines, occurrences, all_bundle_ids, max(numeric_ids, default=-1)


def synthetic_annotation(bundle_id: str, direction: str) -> str:
    digest = hashlib.blake2b(bundle_id.encode("utf-8"), digest_size=8).digest()
    bundle_size = 100 + int.from_bytes(digest, "big") % 900
    if direction == "+":
        bundle_start, bundle_end, encoded_direction = 0, bundle_size - 1, "0"
    else:
        bundle_start, bundle_end, encoded_direction = bundle_size - 1, 0, "1"
    return (
        f"{bundle_id}:{bundle_size}:{encoded_direction}:{bundle_start}:{bundle_end}:U"
    )


def wrap_sequence(sequence: str, width: int = 80) -> Iterable[str]:
    for offset in range(0, len(sequence), width):
        yield sequence[offset : offset + width] + "\n"


def write_reports(
    output_dir: Path,
    hits: Sequence[TRFHit],
    unions: Sequence[RepeatUnion],
) -> None:
    with (output_dir / "trf_intersections.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "path_id", "repeat_id", "trf_start", "trf_end", "trf_length",
                "period", "copy_number", "identity", "overlap_start", "overlap_end",
                "overlap_bp",
            )
        )
        for hit in sorted(hits, key=lambda item: (item.path_id, item.start, item.end)):
            writer.writerow(
                (
                    hit.path_id, hit.repeat_id, hit.start, hit.end, hit.end - hit.start,
                    hit.period, hit.copy_number, f"{hit.identity:g}", hit.overlap_start,
                    hit.overlap_end, hit.overlap_bp,
                )
            )
    with (output_dir / "trf_unions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("path_id", "union_start", "union_end", "union_bp", "source_repeat_ids"))
        for repeat in sorted(unions, key=lambda item: (item.path_id, item.start, item.end)):
            writer.writerow(
                (repeat.path_id, repeat.start, repeat.end, repeat.length, ",".join(repeat.source_repeat_ids))
            )


def write_population_union_report(
    output_dir: Path,
    regions: Sequence[PopulationRepeatRegion],
    records: Sequence[AlignedOccurrence],
) -> None:
    with (output_dir / "population_trf_msa_unions.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "population_union_id", "alignment_start_column",
            "alignment_end_column", "alignment_span_columns",
            "min_tr_length_bp", "median_tr_length_bp", "max_tr_length_bp",
            "tr_length_allele_count", "mean_pairwise_tr_length_divergence",
            "source_occurrence_count",
            "source_path_count", "source_paths", "source_repeat_ids",
            "projected_members",
        ))
        for index, region in enumerate(regions, 1):
            sequence_bp = [
                region.end_column - region.start_column
                - record.sequence.count("-", region.start_column, region.end_column)
                for record in records
            ]
            members = ";".join(
                f"{member.repeat.path_id}:{member.repeat.start}-{member.repeat.end}:"
                f"{member.start_column}-{member.end_column}"
                for member in region.members
            )
            writer.writerow((
                index, region.start_column, region.end_column,
                region.end_column - region.start_column,
                min(sequence_bp), float(statistics.median(sequence_bp)),
                max(sequence_bp), len(set(sequence_bp)),
                mean_pairwise_length_divergence(sequence_bp),
                len(region.members), len(region.source_path_ids),
                ",".join(region.source_path_ids),
                ",".join(region.source_repeat_ids), members,
            ))


def evaluate_repeat_candidates(
    unions: Sequence[RepeatUnion],
    aligned_by_path: Mapping[str, AlignedOccurrence],
    records: Sequence[AlignedOccurrence],
    min_repeat_bp: int,
    min_sequence_coverage: float,
    min_population_presence: float,
    min_nucleotide_divergence: float,
    min_length_divergence: float = 0.05,
) -> List[CandidateEvaluation]:
    """Score connected population unions of TRF-covered MSA columns."""
    column_stats = build_column_prefix_stats(records)
    evaluations: List[CandidateEvaluation] = []
    regions = population_repeat_regions(unions, aligned_by_path)
    for region in regions:
        sequence_bp = [
            region.end_column - region.start_column
            - record.sequence.count("-", region.start_column, region.end_column)
            for record in records
        ]
        reference_repeat_bp = max(sequence_bp)
        if reference_repeat_bp < min_repeat_bp:
            continue
        representative = min(
            region.members,
            key=lambda item: (
                -(item.end_column - item.start_column),
                -item.repeat.length,
                item.repeat.path_id,
                item.repeat.start,
                item.repeat.end,
            ),
        )
        reference = aligned_by_path[representative.repeat.path_id]
        start_offset = ungapped_bases_before(reference.sequence, region.start_column)
        end_offset = ungapped_bases_before(reference.sequence, region.end_column)
        if reference.direction == "+":
            repeat_start = reference.start + start_offset
            repeat_end = reference.start + end_offset
        else:
            repeat_start = reference.end - end_offset
            repeat_end = reference.end - start_offset
        repeat = RepeatUnion(
            reference.path_id,
            repeat_start,
            repeat_end,
            region.source_repeat_ids,
        )
        metrics = analyze_region(
            records,
            region.start_column,
            region.end_column,
            reference_repeat_bp,
            min_sequence_coverage,
            column_stats,
        )
        metrics.update({
            "population_union_source_occurrence_count": len(region.members),
            "population_union_source_path_count": len(region.source_path_ids),
            "population_union_source_paths": list(region.source_path_ids),
            "population_union_source_repeat_ids": list(region.source_repeat_ids),
        })
        structural_divergent = (
            float(metrics["population_presence_fraction"]) < min_population_presence
        )
        nucleotide_divergent = (
            float(metrics["non_gap_consensus_divergence"])
            >= min_nucleotide_divergence
        )
        length_divergent = (
            float(metrics["mean_pairwise_tr_length_divergence"])
            >= min_length_divergence
        )
        evaluations.append(
            CandidateEvaluation(
                candidate_id=len(evaluations) + 1,
                repeat=repeat,
                reference_start=reference.start,
                reference_end=reference.end,
                start_offset=start_offset,
                end_offset=end_offset,
                start_column=region.start_column,
                end_column=region.end_column,
                population_region=region,
                metrics=metrics,
                structural_divergent=structural_divergent,
                length_divergent=length_divergent,
                nucleotide_divergent=nucleotide_divergent,
                divergent=(
                    structural_divergent or length_divergent
                    or nucleotide_divergent
                ),
            )
        )
    return evaluations


def fragment_support(
    records: Sequence[AlignedOccurrence],
    start_column: int,
    end_column: int,
    min_sequence_bp: int,
) -> Tuple[int, int, float, int]:
    """Count sequences with enough real bases in an alignment fragment."""
    span = end_column - start_column
    sequence_bp = sorted(
        span - record.sequence.count("-", start_column, end_column)
        for record in records
    )
    supporting = sum(bp >= min_sequence_bp for bp in sequence_bp)
    nonempty = sum(bp > 0 for bp in sequence_bp)
    return supporting, nonempty, float(statistics.median(sequence_bp)), sequence_bp[-1]


def recursively_decompose_alignment(
    records: Sequence[AlignedOccurrence],
    candidates: Sequence[CandidateEvaluation],
    min_fragment_sequence_bp: int,
    min_fragment_count: int,
    min_fragment_frequency: float,
    force_split: bool,
) -> Tuple[List[RetainedSegment], List[RecursiveDecision], List[CandidateEvaluation], int]:
    """Remove the longest eligible repeat in each child until none remains."""
    required_count = max(
        min_fragment_count,
        int(math.ceil(min_fragment_frequency * len(records))),
    )
    decisions: List[RecursiveDecision] = []
    retained: List[RetainedSegment] = []
    selected: List[CandidateEvaluation] = []
    iteration = 0

    def visit(start_column: int, end_column: int, depth: int) -> None:
        nonlocal iteration
        eligible = [
            candidate
            for candidate in candidates
            if (candidate.divergent or force_split)
            and candidate.selected_iteration is None
            and candidate.start_column >= start_column
            and candidate.end_column <= end_column
        ]
        if not eligible:
            supporting, nonempty, median_bp, max_bp = fragment_support(
                records,
                start_column,
                end_column,
                min_fragment_sequence_bp,
            )
            if end_column > start_column and supporting >= required_count:
                retained.append(
                    RetainedSegment(
                        start_column=start_column,
                        end_column=end_column,
                        supporting_sequence_count=supporting,
                        nonempty_sequence_count=nonempty,
                        median_sequence_bp=median_bp,
                        max_sequence_bp=max_bp,
                        depth=depth,
                    )
                )
            return

        chosen = max(
            eligible,
            key=lambda candidate: (
                candidate.repeat_bp,
                candidate.end_column - candidate.start_column,
                -candidate.candidate_id,
            ),
        )
        iteration += 1
        chosen.selected_iteration = iteration
        selected.append(chosen)
        for side, child_start, child_end in (
            ("left", start_column, chosen.start_column),
            ("right", chosen.end_column, end_column),
        ):
            if child_end > child_start:
                supporting, nonempty, median_bp, max_bp = fragment_support(
                    records,
                    child_start,
                    child_end,
                    min_fragment_sequence_bp,
                )
            else:
                supporting, nonempty, median_bp, max_bp = 0, 0, 0.0, 0
            keep = child_end > child_start and supporting >= required_count
            decisions.append(
                RecursiveDecision(
                    iteration=iteration,
                    depth=depth,
                    parent_start_column=start_column,
                    parent_end_column=end_column,
                    candidate_id=chosen.candidate_id,
                    repeat_path_id=chosen.repeat.path_id,
                    repeat_start=chosen.repeat.start,
                    repeat_end=chosen.repeat.end,
                    repeat_bp=chosen.repeat_bp,
                    repeat_start_column=chosen.start_column,
                    repeat_end_column=chosen.end_column,
                    child_side=side,
                    child_start_column=child_start,
                    child_end_column=child_end,
                    supporting_sequence_count=supporting,
                    nonempty_sequence_count=nonempty,
                    required_sequence_count=required_count,
                    retained=keep,
                )
            )
            if keep:
                visit(child_start, child_end, depth + 1)

    visit(0, len(records[0].sequence), 0)
    retained.sort(key=lambda segment: (segment.start_column, segment.end_column))
    return retained, decisions, selected, required_count


def choose_segment_ids(
    explicit: Optional[str],
    maximum_numeric_id: int,
    existing_ids: Sequence[str],
    segment_count: int,
) -> Tuple[str, ...]:
    if segment_count == 0:
        if explicit:
            raise ValueError("--new-ids was supplied, but no non-repeat segment was retained")
        return ()
    if explicit:
        values = tuple(value.strip() for value in explicit.split(","))
        if len(values) != segment_count or any(not value.isdigit() for value in values):
            raise ValueError(
                "--new-ids must contain exactly one comma-separated numeric ID "
                f"for each of the {segment_count} retained segments"
            )
    else:
        values = tuple(
            str(maximum_numeric_id + offset)
            for offset in range(1, segment_count + 1)
        )
    if len(set(values)) != len(values):
        raise ValueError("retained-segment bundle IDs must be distinct")
    collisions = sorted(set(values) & set(existing_ids))
    if collisions:
        raise ValueError(f"new bundle IDs already exist: {', '.join(collisions)}")
    return values


def project_retained_segment(
    occurrence: BedOccurrence,
    aligned: AlignedOccurrence,
    segment: RetainedSegment,
    segment_index: int,
    new_bundle_id: str,
) -> Optional[SegmentProjection]:
    start_offset = ungapped_bases_before(aligned.sequence, segment.start_column)
    end_offset = ungapped_bases_before(aligned.sequence, segment.end_column)
    if aligned.direction == "+":
        genomic_start = aligned.start + start_offset
        genomic_end = aligned.start + end_offset
    else:
        genomic_start = aligned.end - end_offset
        genomic_end = aligned.end - start_offset
    genomic_start = max(genomic_start, occurrence.start)
    genomic_end = min(genomic_end, occurrence.end)
    if genomic_end <= genomic_start:
        return None
    return SegmentProjection(
        path_id=occurrence.path_id,
        direction=occurrence.direction,
        parent_start=occurrence.start,
        parent_end=occurrence.end,
        segment_index=segment_index,
        new_bundle_id=new_bundle_id,
        alignment_start_column=segment.start_column,
        alignment_end_column=segment.end_column,
        ungapped_start_offset=start_offset,
        ungapped_end_offset=end_offset,
        genomic_start=genomic_start,
        genomic_end=genomic_end,
    )


def write_recursive_split_bed(
    input_lines: Sequence[str],
    target_occurrences: Mapping[str, BedOccurrence],
    aligned_by_path: Mapping[str, AlignedOccurrence],
    target_bundle_id: str,
    retained_segments: Sequence[RetainedSegment],
    new_ids: Sequence[str],
    output_path: Path,
) -> List[SegmentProjection]:
    output_lines: List[str] = []
    projections: List[SegmentProjection] = []
    replaced_paths = set()
    for line in input_lines:
        if not line or line.startswith("#"):
            output_lines.append(line)
            continue
        fields = line.split("\t")
        if fields[3].split(":", 1)[0] != target_bundle_id:
            output_lines.append(line)
            continue
        occurrence = target_occurrences[fields[0]]
        aligned = aligned_by_path.get(fields[0])
        if aligned is None:
            raise ValueError(f"alignment lacks target occurrence for {fields[0]}")
        if occurrence.direction != aligned.direction:
            raise ValueError(f"BED/alignment direction differs for {fields[0]}")
        path_projections: List[SegmentProjection] = []
        for segment_index, (segment, bundle_id) in enumerate(
            zip(retained_segments, new_ids), 1
        ):
            projection = project_retained_segment(
                occurrence, aligned, segment, segment_index, bundle_id
            )
            if projection is not None:
                path_projections.append(projection)
        for projection in sorted(
            path_projections, key=lambda item: (item.genomic_start, item.genomic_end)
        ):
            output_lines.append(
                "\t".join(
                    (
                        projection.path_id,
                        str(projection.genomic_start),
                        str(projection.genomic_end),
                        synthetic_annotation(
                            projection.new_bundle_id, projection.direction
                        ),
                    )
                )
            )
            projections.append(projection)
        replaced_paths.add(fields[0])
    missing = sorted(set(aligned_by_path) - replaced_paths)
    if missing:
        raise ValueError(
            f"BED lacks {len(missing)} aligned target occurrences; "
            f"examples: {', '.join(missing[:5])}"
        )
    output_path.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return projections


def write_retained_alignments(
    records: Sequence[AlignedOccurrence],
    retained_segments: Sequence[RetainedSegment],
    new_ids: Sequence[str],
    output_dir: Path,
) -> Dict[str, str]:
    output_paths: Dict[str, str] = {}
    for segment_index, (segment, bundle_id) in enumerate(
        zip(retained_segments, new_ids), 1
    ):
        path = output_dir / f"bundle_{bundle_id}.segment_{segment_index}.aln.fa"
        with path.open("w", encoding="utf-8") as handle:
            for record in records:
                sequence = record.sequence[segment.start_column : segment.end_column]
                if not sequence or not sequence.strip("-"):
                    continue
                handle.write(
                    f">{record.header}|recursive_segment={segment_index}"
                    f"|new_bundle={bundle_id}\n"
                )
                handle.writelines(wrap_sequence(sequence))
        output_paths[bundle_id] = str(path)
    return output_paths


def write_recursive_reports(
    output_dir: Path,
    candidates: Sequence[CandidateEvaluation],
    decisions: Sequence[RecursiveDecision],
    retained_segments: Sequence[RetainedSegment],
    new_ids: Sequence[str],
    projections: Sequence[SegmentProjection],
) -> None:
    with (output_dir / "candidate_evaluations.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "candidate_id", "selected_iteration", "path_id", "repeat_start",
                "repeat_end", "repeat_bp", "source_repeat_ids", "alignment_start_column",
                "alignment_end_column", "alignment_span_columns",
                "population_union_source_occurrence_count",
                "population_union_source_path_count", "population_union_source_paths",
                "presence_count",
                "presence_fraction", "structural_divergence",
                "tr_length_allele_count", "min_tr_length_bp",
                "median_tr_length_bp", "max_tr_length_bp",
                "tr_length_divergence", "nucleotide_divergence",
                "structural_divergent", "length_divergent",
                "nucleotide_divergent", "divergent",
            )
        )
        for candidate in candidates:
            metrics = candidate.metrics
            writer.writerow(
                (
                    candidate.candidate_id,
                    candidate.selected_iteration if candidate.selected_iteration is not None else ".",
                    candidate.repeat.path_id,
                    candidate.repeat.start,
                    candidate.repeat.end,
                    candidate.repeat_bp,
                    ",".join(candidate.population_region.source_repeat_ids),
                    candidate.start_column,
                    candidate.end_column,
                    candidate.end_column - candidate.start_column,
                    len(candidate.population_region.members),
                    len(candidate.population_region.source_path_ids),
                    ",".join(candidate.population_region.source_path_ids),
                    metrics["sequence_meeting_coverage_count"],
                    metrics["population_presence_fraction"],
                    metrics["structural_divergence"],
                    metrics["tr_length_allele_count"],
                    metrics["min_sequence_bp"],
                    metrics["median_sequence_bp"],
                    metrics["max_sequence_bp"],
                    metrics["mean_pairwise_tr_length_divergence"],
                    metrics["non_gap_consensus_divergence"],
                    int(candidate.structural_divergent),
                    int(candidate.length_divergent),
                    int(candidate.nucleotide_divergent),
                    int(candidate.divergent),
                )
            )
    with (output_dir / "recursive_split_decisions.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        fields = tuple(RecursiveDecision.__dataclass_fields__)
        writer.writerow(fields)
        for decision in decisions:
            writer.writerow(tuple(getattr(decision, field) for field in fields))
    with (output_dir / "retained_segments.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "segment_index", "new_bundle_id", "alignment_start_column",
                "alignment_end_column", "alignment_span_columns", "supporting_sequence_count",
                "nonempty_sequence_count", "median_sequence_bp", "max_sequence_bp", "depth",
            )
        )
        for index, segment in enumerate(retained_segments, 1):
            writer.writerow(
                (
                    index,
                    new_ids[index - 1] if index <= len(new_ids) else ".",
                    segment.start_column,
                    segment.end_column,
                    segment.end_column - segment.start_column,
                    segment.supporting_sequence_count,
                    segment.nonempty_sequence_count,
                    segment.median_sequence_bp,
                    segment.max_sequence_bp,
                    segment.depth,
                )
            )
    with (output_dir / "split_map.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        fields = tuple(SegmentProjection.__dataclass_fields__)
        writer.writerow(fields)
        for projection in projections:
            writer.writerow(tuple(getattr(projection, field) for field in fields))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Recursively remove long, population-variable TRF regions from a bundle MSA; "
            "discard child fragments with insufficient population support."
        )
    )
    parser.add_argument("alignment", type=Path)
    parser.add_argument("trf_bed", type=Path)
    parser.add_argument("bundle_bed", type=Path)
    parser.add_argument("-o", "--output-dir", required=True, type=Path)
    parser.add_argument("--bundle-id", help="default: bundle ID stored in alignment headers")
    parser.add_argument(
        "--trf-min-length", type=int, default=0,
        help="minimum input TRF interval length (default: 0; no length filter)",
    )
    parser.add_argument("--trf-min-identity", type=float, default=80.0)
    parser.add_argument(
        "--split-repeat-min-bp", type=int, default=1_000,
        help="minimum merged TRF length eligible to split a current fragment",
    )
    parser.add_argument(
        "--min-sequence-coverage", type=float, default=0.50,
        help="fraction of the longest repeat length required to call the region present in a sequence",
    )
    parser.add_argument(
        "--min-population-presence", type=float, default=0.50,
        help="presence fraction below which the MSA region is structurally divergent",
    )
    parser.add_argument(
        "--min-length-divergence", type=float, default=0.05,
        help="minimum mean pairwise relative TR-length difference (default: 0.05)",
    )
    parser.add_argument("--min-nucleotide-divergence", type=float, default=0.05)
    parser.add_argument(
        "--min-fragment-sequence-bp", type=int, default=1_000,
        help="real bases a sequence must retain to support a child fragment",
    )
    parser.add_argument(
        "--min-fragment-count", type=int, default=2,
        help="absolute minimum supporting sequence count for a retained fragment",
    )
    parser.add_argument(
        "--min-fragment-frequency", type=float, default=0.02,
        help="minimum fraction of all aligned sequences supporting a retained fragment",
    )
    parser.add_argument(
        "--new-ids",
        help="unused numeric IDs, one per final retained non-repeat segment",
    )
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--force-split", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    for name in (
        "trf_min_identity", "min_sequence_coverage", "min_population_presence",
        "min_length_divergence", "min_nucleotide_divergence",
        "min_fragment_frequency",
    ):
        value = getattr(args, name)
        upper = 100 if name == "trf_min_identity" else 1
        if not 0 <= value <= upper:
            parser.error(f"--{name.replace('_', '-')} must be in [0, {upper}]")
    for name in (
        "trf_min_length", "split_repeat_min_bp", "min_fragment_sequence_bp",
        "min_fragment_count",
    ):
        if getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} must be non-negative")
    if args.min_fragment_sequence_bp == 0:
        parser.error("--min-fragment-sequence-bp must be positive")
    try:
        records = read_alignment(args.alignment)
        inferred_bundle_id = records[0].bundle_id
        bundle_id = args.bundle_id or inferred_bundle_id
        if bundle_id != inferred_bundle_id:
            raise ValueError(
                f"requested bundle {bundle_id}, but alignment contains bundle {inferred_bundle_id}"
            )
        aligned_by_path = {record.path_id: record for record in records}
        hits = read_trf_hits(
            args.trf_bed,
            aligned_by_path,
            args.trf_min_length,
            args.trf_min_identity,
        )
        unions = merge_trf_hits(hits)
        population_regions = population_repeat_regions(unions, aligned_by_path)
        candidates = evaluate_repeat_candidates(
            unions,
            aligned_by_path,
            records,
            args.split_repeat_min_bp,
            args.min_sequence_coverage,
            args.min_population_presence,
            args.min_nucleotide_divergence,
            args.min_length_divergence,
        )
        retained_segments, decisions, selected, required_count = (
            recursively_decompose_alignment(
                records,
                candidates,
                args.min_fragment_sequence_bp,
                args.min_fragment_count,
                args.min_fragment_frequency,
                args.force_split,
            )
        )

        input_lines, target_occurrences, existing_ids, maximum_id = read_target_bed(
            args.bundle_bed,
            bundle_id,
        )
        if set(target_occurrences) != set(aligned_by_path):
            missing_alignment = sorted(set(target_occurrences) - set(aligned_by_path))
            missing_bed = sorted(set(aligned_by_path) - set(target_occurrences))
            raise ValueError(
                "BED/alignment target paths differ: "
                f"missing_alignment={missing_alignment[:5]}, missing_bed={missing_bed[:5]}"
            )
        should_split = bool(selected) and not args.report_only
        new_ids = (
            choose_segment_ids(
                args.new_ids,
                maximum_id,
                existing_ids,
                len(retained_segments),
            )
            if should_split
            else ()
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        projections: List[SegmentProjection] = []
        split_bed_path: Optional[Path] = None
        alignment_outputs: Dict[str, str] = {}
        if should_split:
            split_bed_path = args.output_dir / f"bundle_{bundle_id}.split.bed"
            projections = write_recursive_split_bed(
                input_lines,
                target_occurrences,
                aligned_by_path,
                bundle_id,
                retained_segments,
                new_ids,
                split_bed_path,
            )
            alignment_outputs = write_retained_alignments(
                records,
                retained_segments,
                new_ids,
                args.output_dir,
            )
        write_reports(args.output_dir, hits, unions)
        write_population_union_report(args.output_dir, population_regions, records)
        write_recursive_reports(
            args.output_dir,
            candidates,
            decisions,
            retained_segments,
            new_ids,
            projections,
        )

        primary = selected[0] if selected else (
            max(candidates, key=lambda item: item.repeat_bp) if candidates else None
        )
        reasons: List[str] = []
        if primary is not None and primary.structural_divergent:
            reasons.append("low_population_presence")
        if primary is not None and primary.length_divergent:
            reasons.append("high_tr_length_divergence")
        if primary is not None and primary.nucleotide_divergent:
            reasons.append("high_nucleotide_divergence")

        def candidate_dict(candidate: CandidateEvaluation) -> Dict[str, object]:
            return {
                "candidate_id": candidate.candidate_id,
                "selected_iteration": candidate.selected_iteration,
                "path_id": candidate.repeat.path_id,
                "start": candidate.repeat.start,
                "end": candidate.repeat.end,
                "length_bp": candidate.repeat_bp,
                "source_repeat_ids": list(candidate.population_region.source_repeat_ids),
                "population_union_source_occurrence_count": len(
                    candidate.population_region.members
                ),
                "population_union_source_path_count": len(
                    candidate.population_region.source_path_ids
                ),
                "population_union_source_paths": list(
                    candidate.population_region.source_path_ids
                ),
                "bundle_alignment_header_start": candidate.reference_start,
                "bundle_alignment_header_end": candidate.reference_end,
                "ungapped_start_offset": candidate.start_offset,
                "ungapped_end_offset": candidate.end_offset,
                "alignment_start_column_0based": candidate.start_column,
                "alignment_end_column_0based_exclusive": candidate.end_column,
                "structural_divergent": candidate.structural_divergent,
                "length_divergent": candidate.length_divergent,
                "nucleotide_divergent": candidate.nucleotide_divergent,
                "divergent": candidate.divergent,
                "metrics": candidate.metrics,
            }

        remaining_eligible = [
            candidate
            for candidate in candidates
            if (candidate.divergent or args.force_split)
            and candidate.selected_iteration is None
            and any(
                candidate.start_column >= segment.start_column
                and candidate.end_column <= segment.end_column
                for segment in retained_segments
            )
        ]
        analysis = {
            "alignment": str(args.alignment),
            "bundle_bed": str(args.bundle_bed),
            "bundle_id": bundle_id,
            "trf_bed": str(args.trf_bed),
            "trf_min_length": args.trf_min_length,
            "trf_min_identity": args.trf_min_identity,
            "filtered_trf_hit_count": len(hits),
            "trf_union_count": len(unions),
            "population_msa_union_count": len(population_regions),
            "long_repeat_candidate_count": len(candidates),
            "divergent_repeat_candidate_count": sum(
                candidate.divergent for candidate in candidates
            ),
            "selected_repeat_count": len(selected),
            "remaining_eligible_repeat_count": len(remaining_eligible),
            "longest_repeat": candidate_dict(primary) if primary else None,
            "region_metrics": primary.metrics if primary else None,
            "thresholds": {
                "split_repeat_min_bp": args.split_repeat_min_bp,
                "min_sequence_coverage": args.min_sequence_coverage,
                "min_population_presence": args.min_population_presence,
                "min_length_divergence": args.min_length_divergence,
                "min_nucleotide_divergence": args.min_nucleotide_divergence,
                "min_fragment_sequence_bp": args.min_fragment_sequence_bp,
                "min_fragment_count": args.min_fragment_count,
                "min_fragment_frequency": args.min_fragment_frequency,
                "required_fragment_sequence_count": required_count,
            },
            "structural_divergent": any(
                candidate.structural_divergent for candidate in candidates
            ),
            "length_divergent": any(
                candidate.length_divergent for candidate in candidates
            ),
            "nucleotide_divergent": any(
                candidate.nucleotide_divergent for candidate in candidates
            ),
            "divergent": any(candidate.divergent for candidate in candidates),
            "divergence_reasons": reasons,
            "selected_repeats": [candidate_dict(candidate) for candidate in selected],
            "removed_repeat_alignment_regions": [
                {
                    "iteration": candidate.selected_iteration,
                    "start_column_0based": candidate.start_column,
                    "end_column_0based_exclusive": candidate.end_column,
                }
                for candidate in selected
            ],
            "retained_segments": [
                {
                    "segment_index": index,
                    "new_bundle_id": (
                        new_ids[index - 1] if index <= len(new_ids) else None
                    ),
                    "alignment_start_column_0based": segment.start_column,
                    "alignment_end_column_0based_exclusive": segment.end_column,
                    "alignment_span_columns": segment.end_column - segment.start_column,
                    "supporting_sequence_count": segment.supporting_sequence_count,
                    "nonempty_sequence_count": segment.nonempty_sequence_count,
                    "median_sequence_bp": segment.median_sequence_bp,
                    "max_sequence_bp": segment.max_sequence_bp,
                    "depth": segment.depth,
                }
                for index, segment in enumerate(retained_segments, 1)
            ],
            "dropped_child_fragment_count": sum(
                not decision.retained for decision in decisions
            ),
            "split_performed": should_split,
            "new_bundle_ids": list(new_ids) if should_split else None,
            "split_bed": str(split_bed_path) if split_bed_path else None,
            "split_alignments": alignment_outputs,
            "target_occurrence_count": len(target_occurrences),
            "emitted_segment_occurrence_count": len(projections),
            "repeat_segments_emitted": 0,
        }
        with (args.output_dir / "analysis.json").open("w", encoding="utf-8") as handle:
            json.dump(analysis, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print(json.dumps(analysis, indent=2, sort_keys=True))
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
