#!/usr/bin/env python3
"""Population-graph cleanup of PGR-TK principal-bundle decompositions.

The input is the four-column BED emitted by ``pgr-pbundle-decomp``::

    path  start  end  bundle_id:bundle_size:direction:pb_start:pb_end:R|U

The workflow keeps every observed haplotype walk, builds a bidirected bundle
adjacency graph, computes a deterministic topology on the SCC-condensed signed
graph, and examines the observed subpaths between stable anchor bundles. A
subpath is either retained as a structural alternative or trimmed into a more
supported subpath when its extra nodes are short (and, optionally, rare-long)
and the total removed length is small. Retained non-empty subpaths become new
graph bundles that are projected back onto the original haplotype coordinates.
When a TRF BED is supplied, bundle IDs recurrently overlapped by selected TRF
calls are removed before population statistics or graph operations are computed.
Anchor subgraphs with enough shared interior bundle span can use one joint
retained-target pool for cross-boundary trimming.
The automatic anchor length threshold is Q25 of the bundle-type median lengths.
Unmerged short types are retained here. Optional coverage/haplotype-preserving
deletion is performed only on the final population remapping BED.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple


SignedToken = Tuple[str, str]
Port = Tuple[str, str]
PortEdgeKey = Tuple[Port, Port]
BoundaryKey = Tuple[SignedToken, SignedToken]
RepeatPattern = Tuple[bool, ...]
SubpathKey = Tuple[Tuple[SignedToken, ...], RepeatPattern]


@dataclass(frozen=True)
class Occurrence:
    path_id: str
    start: int
    end: int
    bundle_id: str
    bundle_size: Optional[int]
    direction: str
    bundle_start: Optional[int]
    bundle_end: Optional[int]
    repeat_flag: str
    line_number: int

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def token(self) -> SignedToken:
        return self.bundle_id, self.direction


@dataclass
class BundleStats:
    bundle_id: str
    occurrence_count: int = 0
    units: Set[str] = field(default_factory=set)
    paths: Set[str] = field(default_factory=set)
    repeat_units: Set[str] = field(default_factory=set)
    repeat_paths: Set[str] = field(default_factory=set)
    lengths: List[int] = field(default_factory=list)
    repeat_count: int = 0
    median_length: float = 0.0
    frequency: float = 0.0
    repeat_frequency: float = 0.0
    single_copy_path_fraction: float = 0.0
    rare: bool = False
    short: bool = False
    anchor: bool = False
    duplicated: bool = False
    long_dup_breakpoint: bool = False
    action: str = ""
    short_deletion_action: str = ""


@dataclass(frozen=True)
class AdjacencyObservation:
    path_id: str
    unit_id: str
    left_index: int
    right_index: int
    gap_bp: int


@dataclass
class PopulationEdge:
    key: PortEdgeKey
    observations: List[AdjacencyObservation] = field(default_factory=list)

    @property
    def path_ids(self) -> Set[str]:
        return {observation.path_id for observation in self.observations}

    @property
    def unit_ids(self) -> Set[str]:
        return {observation.unit_id for observation in self.observations}

    @property
    def gaps(self) -> List[int]:
        return [observation.gap_bp for observation in self.observations]


@dataclass(frozen=True)
class SegmentObservation:
    path_id: str
    unit_id: str
    left_index: int
    right_index: int
    boundary: BoundaryKey
    interior: Tuple[SignedToken, ...]
    canonical_orientation: str
    duplication_orientation: str


@dataclass
class SubpathClass:
    boundary: BoundaryKey
    interior: Tuple[SignedToken, ...]
    repeat_pattern: RepeatPattern
    dup_tokens: Tuple[SignedToken, ...]
    duplication_signature: Tuple[SignedToken, ...]
    observations: List[SegmentObservation] = field(default_factory=list)
    support_units: Set[str] = field(default_factory=set)
    support_paths: Set[str] = field(default_factory=set)
    common: bool = False
    action: str = ""
    target_interior: Optional[Tuple[SignedToken, ...]] = None
    target_key: Optional[SubpathKey] = None
    target_boundary: Optional[BoundaryKey] = None
    merged_repeat_pattern: RepeatPattern = ()
    merged_dup_tokens: Tuple[SignedToken, ...] = ()
    merged_duplication_signature: Tuple[SignedToken, ...] = ()
    trim_bp: Optional[int] = None
    trim_bundle_count: Optional[int] = None
    new_bundle_id: Optional[str] = None

    @property
    def key(self) -> SubpathKey:
        return self.interior, self.repeat_pattern


@dataclass(frozen=True)
class GraphBundle:
    new_id: str
    boundaries: Tuple[BoundaryKey, ...]
    canonical_interior: Tuple[SignedToken, ...]
    support_units: int
    support_paths: int
    source_action: str
    duplicated: bool


@dataclass(frozen=True)
class SubgraphGroup:
    group_id: int
    boundaries: Tuple[BoundaryKey, ...]
    interior_bundle_ids: Tuple[str, ...]
    average_total_bp: float


@dataclass(frozen=True)
class SharedBoundaryPair:
    boundary_a: BoundaryKey
    boundary_b: BoundaryKey
    shared_bundle_ids: Tuple[str, ...]
    shared_average_total_bp: float
    merged: bool


@dataclass(frozen=True)
class LongDupBreakpoint:
    bundle_id: str
    candidate_graph_bundle_ids: Tuple[str, ...]
    candidate_boundaries: Tuple[BoundaryKey, ...]


@dataclass(frozen=True)
class OutputOccurrence:
    path_id: str
    start: int
    end: int
    new_bundle_id: str
    direction: str
    source: str
    original_bundle_ids: Tuple[str, ...]
    absorbed_bundle_ids: Tuple[str, ...] = ()
    boundary_text: str = "."
    bundle_size: Optional[int] = None
    bundle_start: Optional[int] = None
    bundle_end: Optional[int] = None
    repeat_flag: str = "U"
    trf_action: str = "none"
    trf_repeat_bp: int = 0
    trf_coverage: float = 0.0


@dataclass(frozen=True)
class CoordinateAdjustment:
    path_id: str
    left_bundle_id: str
    right_bundle_id: str
    original_left_end: int
    original_right_start: int
    overlap_bp: int
    split_coordinate: int
    left_trim_bp: int
    right_trim_bp: int


@dataclass(frozen=True)
class TRFInterval:
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class TRFRecord:
    path_id: str
    start: int
    end: int
    repeat_id: str
    period: str
    copy_number: str
    identity: float


def normalize_direction(value: str) -> str:
    value = value.strip().lower()
    if value in {"0", "+", "f", "forward", "fw"}:
        return "+"
    if value in {"1", "-", "r", "reverse", "rev", "rc"}:
        return "-"
    raise ValueError(f"unsupported PGR-TK bundle direction: {value!r}")


def flip_direction(direction: str) -> str:
    return "-" if direction == "+" else "+"


def reverse_complement_tokens(tokens: Sequence[SignedToken]) -> Tuple[SignedToken, ...]:
    return tuple((bundle_id, flip_direction(direction)) for bundle_id, direction in reversed(tokens))


def canonicalize_tokens(tokens: Sequence[SignedToken]) -> Tuple[Tuple[SignedToken, ...], str]:
    forward = tuple(tokens)
    reverse = reverse_complement_tokens(tokens)
    return (forward, "+") if forward <= reverse else (reverse, "-")


def token_text(token: SignedToken) -> str:
    return f"{token[1]}{token[0]}"


def chain_text(tokens: Sequence[SignedToken]) -> str:
    return ",".join(token_text(token) for token in tokens) if tokens else "."


def repeat_pattern_text(pattern: Sequence[bool]) -> str:
    return ",".join("R" if duplicated else "U" for duplicated in pattern) if pattern else "."


def boundary_text(boundary: BoundaryKey) -> str:
    return f"{token_text(boundary[0])}>{token_text(boundary[1])}"


def pgr_output_annotation(item: OutputOccurrence) -> str:
    """Format an output block using the native PGR-TK annotation schema.

    Original retained bundles preserve their PGR minimizer coordinates. A new
    graph bundle has no valid coordinate system in the original principal
    bundle graph, so it receives reproducible synthetic integer placeholders.
    """
    direction = "0" if item.direction == "+" else "1"
    bundle_size = item.bundle_size
    bundle_start = item.bundle_start
    bundle_end = item.bundle_end
    if bundle_size is None or bundle_start is None or bundle_end is None:
        digest = hashlib.blake2b(item.new_bundle_id.encode("utf-8"), digest_size=8).digest()
        bundle_size = 100 + int.from_bytes(digest, "big") % 900
        if item.direction == "+":
            bundle_start, bundle_end = 0, bundle_size - 1
        else:
            bundle_start, bundle_end = bundle_size - 1, 0
    optional = lambda value: "." if value is None else str(value)
    repeat_flag = "R" if item.repeat_flag.upper() == "R" else "U"
    return ":".join(
        (
            item.new_bundle_id,
            optional(bundle_size),
            direction,
            optional(bundle_start),
            optional(bundle_end),
            repeat_flag,
        )
    )


def parse_optional_int(value: str) -> Optional[int]:
    return None if value in {"", ".", "NA", "na", "None", "none"} else int(value)


def parse_pgr_annotation(annotation: str) -> Tuple[str, Optional[int], str, Optional[int], Optional[int], str]:
    fields = annotation.rstrip().split(":")
    if len(fields) < 6:
        raise ValueError(
            "annotation requires bundle_id:bundle_size:direction:bundle_start:bundle_end:R|U"
        )
    return (
        fields[0],
        parse_optional_int(fields[1]),
        normalize_direction(fields[2]),
        parse_optional_int(fields[3]),
        parse_optional_int(fields[4]),
        fields[5],
    )


def read_pgr_bed(path: Path) -> Dict[str, List[Occurrence]]:
    by_path: Dict[str, List[Occurrence]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 4:
                fields = stripped.split()
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected at least four BED columns")
            path_id = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: start/end must be integers") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval [{start}, {end})")
            try:
                bundle_id, bundle_size, direction, bundle_start, bundle_end, repeat_flag = (
                    parse_pgr_annotation(fields[3])
                )
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc
            by_path[path_id].append(
                Occurrence(
                    path_id,
                    start,
                    end,
                    bundle_id,
                    bundle_size,
                    direction,
                    bundle_start,
                    bundle_end,
                    repeat_flag,
                    line_number,
                )
            )
    if not by_path:
        raise ValueError(f"{path}: no bundle records found")
    for records in by_path.values():
        records.sort(key=lambda record: (record.start, record.end, record.line_number))
    return dict(by_path)


def merge_trf_intervals(intervals: Sequence[Tuple[int, int]]) -> Tuple[TRFInterval, ...]:
    merged: List[TRFInterval] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1].end:
            merged[-1] = TRFInterval(merged[-1].start, max(merged[-1].end, end))
        else:
            merged.append(TRFInterval(start, end))
    return tuple(merged)


def read_trf_bed(
    path: Path,
    min_length: int,
    min_identity: float,
) -> Tuple[
    Dict[str, Tuple[TRFInterval, ...]],
    Dict[str, Tuple[TRFRecord, ...]],
    Dict[str, int],
]:
    """Read TRF BED and union intervals passing total-length/identity filters."""
    selected: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    selected_records: Dict[str, List[TRFRecord]] = defaultdict(list)
    input_rows = 0
    selected_rows = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 7:
                fields = stripped.split()
            if len(fields) < 7:
                raise ValueError(f"{path}:{line_number}: expected at least seven TRF BED columns")
            input_rows += 1
            try:
                start, end = int(fields[1]), int(fields[2])
                identity = float(fields[6])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: TRF start/end/identity must be numeric"
                ) from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid TRF interval [{start}, {end})")
            if end - start < min_length or identity < min_identity:
                continue
            path_id = fields[0]
            selected[path_id].append((start, end))
            selected_records[path_id].append(
                TRFRecord(
                    path_id=path_id,
                    start=start,
                    end=end,
                    repeat_id=fields[3] if len(fields) > 3 else ".",
                    period=fields[4] if len(fields) > 4 else ".",
                    copy_number=fields[5] if len(fields) > 5 else ".",
                    identity=identity,
                )
            )
            selected_rows += 1
    merged = {
        path_id: merge_trf_intervals(intervals)
        for path_id, intervals in selected.items()
    }
    records = {
        path_id: tuple(sorted(items, key=lambda item: (item.start, item.end, item.repeat_id)))
        for path_id, items in selected_records.items()
    }
    return merged, records, {
        "trf_input_row_count": input_rows,
        "trf_selected_row_count": selected_rows,
        "trf_union_interval_count": sum(map(len, merged.values())),
    }


def apply_trf_filter_and_report(
    outputs: Mapping[str, Sequence[OutputOccurrence]],
    intervals_by_path: Mapping[str, Sequence[TRFInterval]],
    records_by_path: Mapping[str, Sequence[TRFRecord]],
    full_coverage: float,
    remove_full_repeat: bool = True,
) -> Tuple[
    Dict[str, List[OutputOccurrence]],
    List[Sequence[object]],
    List[Sequence[object]],
    Dict[str, int],
]:
    """Drop repeat-covered occurrences and report all retained/removed overlaps.

    Selected TRF intervals are merged when they overlap or touch.  Full-repeat
    coverage is the largest overlap with any one merged interval, so separate
    repeat islands cannot be added together to trigger whole-bundle removal.
    Total repeat bases remain available in the report. Internal intervals are
    reported as future break candidates but are not split in this stage.
    """
    processed: Dict[str, List[OutputOccurrence]] = {}
    report_rows: List[Sequence[object]] = []
    bundle_metrics: Dict[str, Dict[str, object]] = defaultdict(
        lambda: {
            "occurrence_count": 0,
            "repeat_occurrence_count": 0,
            "internal_repeat_occurrence_count": 0,
            "dropped_occurrence_count": 0,
            "repeat_bp": [],
            "coverage": [],
        }
    )
    matched_paths = 0
    dropped_count = 0
    repeat_overlap_count = 0
    internal_repeat_count = 0
    reported_count = 0
    for path_id, items in outputs.items():
        intervals = intervals_by_path.get(path_id, ())
        records = records_by_path.get(path_id, ())
        matched_paths += bool(intervals)
        path_outputs: List[OutputOccurrence] = []
        for item in items:
            overlaps = tuple(
                (max(item.start, interval.start), min(item.end, interval.end))
                for interval in intervals
                if interval.end > item.start and interval.start < item.end
            )
            repeat_bp = sum(end - start for start, end in overlaps)
            coverage = max((end - start for start, end in overlaps), default=0) / (
                item.end - item.start
            )
            internal = tuple(
                (start, end)
                for start, end in overlaps
                if start > item.start and end < item.end
            )
            raw_records = tuple(
                record
                for record in records
                if record.end > item.start and record.start < item.end
            )
            full_repeat = coverage >= full_coverage
            action = (
                "drop_full_repeat"
                if full_repeat and remove_full_repeat
                else "keep_reported_full_repeat"
                if full_repeat
                else "keep_partial_repeat"
                if repeat_bp
                else "keep_no_repeat"
            )
            metric = bundle_metrics[item.new_bundle_id]
            metric["occurrence_count"] = int(metric["occurrence_count"]) + 1
            cast_repeat_bp = metric["repeat_bp"]
            cast_coverage = metric["coverage"]
            assert isinstance(cast_repeat_bp, list) and isinstance(cast_coverage, list)
            cast_repeat_bp.append(repeat_bp)
            cast_coverage.append(coverage)
            if repeat_bp:
                repeat_overlap_count += 1
                metric["repeat_occurrence_count"] = int(metric["repeat_occurrence_count"]) + 1
            if internal:
                internal_repeat_count += 1
                metric["internal_repeat_occurrence_count"] = (
                    int(metric["internal_repeat_occurrence_count"]) + 1
                )
            if full_repeat:
                metric["full_repeat_occurrence_count"] = (
                    int(metric.get("full_repeat_occurrence_count", 0)) + 1
                )
            if action == "drop_full_repeat":
                dropped_count += 1
                metric["dropped_occurrence_count"] = (
                    int(metric["dropped_occurrence_count"]) + 1
                )
            else:
                path_outputs.append(
                    replace(
                        item,
                        trf_action=action,
                        trf_repeat_bp=repeat_bp,
                        trf_coverage=coverage,
                    )
                )

            internal_bp = sum(end - start for start, end in internal)
            break_candidates = ",".join(
                f"{start}-{end}:left={start - item.start}:right={item.end - end}"
                for start, end in internal
            )
            raw_text = ";".join(
                f"{record.repeat_id}:{record.start}-{record.end}:period={record.period}:"
                f"copies={record.copy_number}:identity={record.identity:g}"
                for record in raw_records
            )
            report_rows.append(
                (
                    path_id, item.start, item.end, item.end - item.start,
                    item.new_bundle_id, item.direction, item.source, action,
                    repeat_bp, f"{coverage:.6f}", len(overlaps),
                    ",".join(f"{start}-{end}" for start, end in overlaps) or ".",
                    len(internal), internal_bp,
                    max((end - start for start, end in internal), default=0),
                    ",".join(f"{start}-{end}" for start, end in internal) or ".",
                    break_candidates or ".", len(raw_records), raw_text or ".",
                )
            )
            reported_count += 1
        path_outputs.sort(key=lambda item: (item.start, item.end, item.new_bundle_id))
        processed[path_id] = path_outputs

    bundle_summary_rows: List[Sequence[object]] = []
    for bundle_id, metric in sorted(bundle_metrics.items()):
        occurrence_count = int(metric["occurrence_count"])
        repeat_occurrence_count = int(metric["repeat_occurrence_count"])
        dropped = int(metric["dropped_occurrence_count"])
        repeat_bps = metric["repeat_bp"]
        coverages = metric["coverage"]
        assert isinstance(repeat_bps, list) and isinstance(coverages, list)
        bundle_summary_rows.append(
            (
                bundle_id, occurrence_count, repeat_occurrence_count,
                f"{repeat_occurrence_count / occurrence_count:.6f}",
                int(metric["internal_repeat_occurrence_count"]),
                int(metric.get("full_repeat_occurrence_count", 0)), dropped,
                occurrence_count - dropped, f"{statistics.median(repeat_bps):.1f}",
                max(repeat_bps), f"{statistics.median(coverages):.6f}",
                f"{max(coverages):.6f}",
                "removed_all_occurrences"
                if dropped == occurrence_count
                else "removed_some_occurrences"
                if dropped
                else "retained",
            )
        )
    return processed, report_rows, bundle_summary_rows, {
        "trf_matched_path_count": matched_paths,
        "trf_reported_occurrence_count": reported_count,
        "trf_repeat_overlapping_occurrence_count": repeat_overlap_count,
        "trf_internal_repeat_occurrence_count": internal_repeat_count,
        "trf_full_repeat_occurrence_count": sum(
            int(metric.get("full_repeat_occurrence_count", 0))
            for metric in bundle_metrics.values()
        ),
        "trf_dropped_occurrence_count": dropped_count,
        "trf_dropped_bundle_type_count": sum(
            int(metric["dropped_occurrence_count"]) == int(metric["occurrence_count"])
            for metric in bundle_metrics.values()
        ),
    }


def filter_input_pgr_by_trf(
    by_path: Mapping[str, Sequence[Occurrence]],
    intervals_by_path: Mapping[str, Sequence[TRFInterval]],
    records_by_path: Mapping[str, Sequence[TRFRecord]],
    full_coverage: float,
    bundle_removal_frequency: float,
) -> Tuple[
    Dict[str, List[Occurrence]],
    List[Sequence[object]],
    List[Sequence[object]],
    Dict[str, int],
]:
    """Remove a whole bundle ID when one repeat island covers enough occurrences.

    A bundle occurrence qualifies only when one merged length/identity-filtered
    TRF interval covers at least ``full_coverage`` of it. Separate merged
    intervals never add their coverage. If the qualifying occurrence fraction
    reaches ``bundle_removal_frequency``, every occurrence of that bundle ID is
    removed; otherwise every occurrence is retained.
    """
    wrapped: Dict[str, List[OutputOccurrence]] = {
        path_id: [
            OutputOccurrence(
                path_id=record.path_id,
                start=record.start,
                end=record.end,
                new_bundle_id=record.bundle_id,
                direction=record.direction,
                source="premerge_pgr_bundle",
                original_bundle_ids=(record.bundle_id,),
                bundle_size=record.bundle_size,
                bundle_start=record.bundle_start,
                bundle_end=record.bundle_end,
                repeat_flag=record.repeat_flag,
            )
            for record in items
        ]
        for path_id, items in by_path.items()
    }
    _, report_rows, bundle_summary_rows, action_summary = (
        apply_trf_filter_and_report(
            wrapped,
            intervals_by_path,
            records_by_path,
            full_coverage,
            remove_full_repeat=False,
        )
    )
    occurrence_count = Counter(str(row[4]) for row in report_rows)
    covered_count = Counter(
        str(row[4]) for row in report_rows if row[7] == "keep_reported_full_repeat"
    )
    removed_bundle_ids = {
        bundle_id
        for bundle_id, count in occurrence_count.items()
        if covered_count[bundle_id] / count >= bundle_removal_frequency
    }

    report_rows = [
        (
            *row[:7],
            "drop_bundle_id_trf_frequency"
            if str(row[4]) in removed_bundle_ids
            else "retain_bundle_id_trf_frequency",
            *row[8:],
        )
        for row in report_rows
    ]
    revised_summary_rows: List[Sequence[object]] = []
    for row in bundle_summary_rows:
        values = list(row)
        bundle_id = str(values[0])
        total = int(values[1])
        removed = bundle_id in removed_bundle_ids
        values[6] = total if removed else 0
        values[7] = 0 if removed else total
        values[12] = "removed_all_occurrences" if removed else "retained_all_occurrences"
        revised_summary_rows.append(tuple(values))

    filtered = {
        path_id: [record for record in records if record.bundle_id not in removed_bundle_ids]
        for path_id, records in by_path.items()
    }
    action_summary["trf_dropped_occurrence_count"] = sum(
        occurrence_count[bundle_id] for bundle_id in removed_bundle_ids
    )
    action_summary["trf_dropped_bundle_type_count"] = len(removed_bundle_ids)
    return filtered, report_rows, revised_summary_rows, action_summary


def read_sample_map(path: Optional[Path], path_ids: Iterable[str]) -> Dict[str, str]:
    observed = set(path_ids)
    if path is None:
        return {path_id: path_id for path_id in observed}
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected path_id and sample_id")
            path_id, sample_id = fields[:2]
            if path_id in mapping and mapping[path_id] != sample_id:
                raise ValueError(f"{path}:{line_number}: conflicting sample for {path_id}")
            mapping[path_id] = sample_id
    missing = sorted(observed - mapping.keys())
    if missing:
        raise ValueError(f"sample map lacks {len(missing)} paths, e.g. {', '.join(missing[:5])}")
    return {path_id: mapping[path_id] for path_id in observed}


def read_anchor_file(path: Optional[Path]) -> Optional[Set[str]]:
    if path is None:
        return None
    anchors: Set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip().split("#", 1)[0].strip()
            if value:
                anchors.add(value.split()[0])
    if not anchors:
        raise ValueError(f"{path}: anchor file is empty")
    return anchors


def compute_bundle_stats(
    by_path: Mapping[str, Sequence[Occurrence]],
    unit_by_path: Mapping[str, str],
    min_count: int,
    min_frequency: float,
    min_dup_count: int,
    min_dup_frequency: float,
    short_bp: int,
) -> Tuple[Dict[str, BundleStats], int, int, int]:
    """Census original bundle types and assign rare/short labels."""
    stats: Dict[str, BundleStats] = {}
    copy_by_bundle_path: Dict[str, Counter[str]] = defaultdict(Counter)
    for path_id, records in by_path.items():
        unit_id = unit_by_path[path_id]
        for record in records:
            stat = stats.setdefault(record.bundle_id, BundleStats(record.bundle_id))
            stat.occurrence_count += 1
            stat.units.add(unit_id)
            stat.paths.add(path_id)
            stat.lengths.append(record.length)
            stat.repeat_count += record.repeat_flag.upper() == "R"
            if record.repeat_flag.upper() == "R":
                stat.repeat_units.add(unit_id)
                stat.repeat_paths.add(path_id)
            copy_by_bundle_path[record.bundle_id][path_id] += 1
    total_units = len(set(unit_by_path.values()))
    required_units = max(min_count, math.ceil(min_frequency * total_units))
    required_dup_units = max(min_dup_count, math.ceil(min_dup_frequency * total_units))
    for bundle_id, stat in stats.items():
        stat.median_length = float(statistics.median(stat.lengths))
        stat.frequency = len(stat.units) / total_units if total_units else 0.0
        stat.repeat_frequency = len(stat.repeat_units) / total_units if total_units else 0.0
        counts = copy_by_bundle_path[bundle_id]
        stat.single_copy_path_fraction = (
            sum(count == 1 for count in counts.values()) / len(counts) if counts else 0.0
        )
        stat.rare = len(stat.units) < required_units
        stat.short = stat.median_length < short_bp
        stat.duplicated = bool(stat.repeat_units) and len(stat.repeat_units) >= required_dup_units
    return stats, total_units, required_units, required_dup_units


def choose_anchors(
    stats: MutableMapping[str, BundleStats],
    explicit: Optional[Set[str]],
    min_frequency: float,
    min_bp: float,
    min_single_copy_fraction: float,
) -> Set[str]:
    if explicit is not None:
        missing = sorted(explicit - stats.keys())
        if missing:
            raise ValueError(f"anchor file contains unknown bundle IDs: {', '.join(missing[:10])}")
        anchors = set(explicit)
    else:
        anchors = {
            bundle_id
            for bundle_id, stat in stats.items()
            if (
                not stat.rare
                and stat.frequency >= min_frequency
                and stat.median_length >= min_bp
                and stat.single_copy_path_fraction >= min_single_copy_fraction
            )
        }
    for bundle_id, stat in stats.items():
        stat.anchor = bundle_id in anchors
    return anchors


def find_long_dup_breakpoints_in_candidates(
    stats: Mapping[str, BundleStats],
    candidate_subgraphs: Mapping[BoundaryKey, Mapping[SubpathKey, SubpathClass]],
    min_bp: int,
    min_repeat_frequency: float,
) -> List[LongDupBreakpoint]:
    """Find long, recurrent duplicated IDs only inside first-pass candidates."""
    candidate_ids: Dict[str, Set[str]] = defaultdict(set)
    candidate_boundaries: Dict[str, Set[BoundaryKey]] = defaultdict(set)
    for boundary, classes in candidate_subgraphs.items():
        for subpath in classes.values():
            if subpath.new_bundle_id is None:
                continue
            for bundle_id in {token[0] for token in subpath.interior}:
                stat = stats[bundle_id]
                if (
                    stat.duplicated
                    and stat.median_length > min_bp
                    and stat.repeat_frequency >= min_repeat_frequency
                ):
                    candidate_ids[bundle_id].add(subpath.new_bundle_id)
                    candidate_boundaries[bundle_id].add(boundary)
    return [
        LongDupBreakpoint(
            bundle_id=bundle_id,
            candidate_graph_bundle_ids=tuple(
                sorted(candidate_ids[bundle_id], key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value))
            ),
            candidate_boundaries=tuple(sorted(candidate_boundaries[bundle_id])),
        )
        for bundle_id in sorted(candidate_ids)
    ]


def promote_long_dup_breakpoints(
    stats: MutableMapping[str, BundleStats],
    anchors: Set[str],
    breakpoints: Sequence[LongDupBreakpoint],
) -> None:
    """Promote candidate-derived breakpoint IDs before the final merge pass."""
    for breakpoint in breakpoints:
        anchors.add(breakpoint.bundle_id)
        stats[breakpoint.bundle_id].anchor = True
        stats[breakpoint.bundle_id].long_dup_breakpoint = True


def left_entry_port(record: Occurrence) -> Port:
    return record.bundle_id, ("L" if record.direction == "+" else "R")


def right_exit_port(record: Occurrence) -> Port:
    return record.bundle_id, ("R" if record.direction == "+" else "L")


def port_edge(left: Occurrence, right: Occurrence) -> PortEdgeKey:
    return tuple(sorted((right_exit_port(left), left_entry_port(right))))  # type: ignore[return-value]


def build_population_graph(
    by_path: Mapping[str, Sequence[Occurrence]],
    unit_by_path: Mapping[str, str],
) -> Dict[PortEdgeKey, PopulationEdge]:
    edges: Dict[PortEdgeKey, PopulationEdge] = {}
    for path_id, records in by_path.items():
        for index, (left, right) in enumerate(zip(records, records[1:])):
            key = port_edge(left, right)
            edge = edges.setdefault(key, PopulationEdge(key))
            edge.observations.append(
                AdjacencyObservation(
                    path_id=path_id,
                    unit_id=unit_by_path[path_id],
                    left_index=index,
                    right_index=index + 1,
                    gap_bp=max(0, right.start - left.end),
                )
            )
    return edges


def normalized_walk(records: Sequence[Occurrence]) -> Tuple[SignedToken, ...]:
    return canonicalize_tokens([record.token for record in records])[0]


def tarjan_scc(
    nodes: Iterable[SignedToken],
    adjacency: Mapping[SignedToken, Set[SignedToken]],
) -> List[Tuple[SignedToken, ...]]:
    index = 0
    stack: List[SignedToken] = []
    on_stack: Set[SignedToken] = set()
    indices: Dict[SignedToken, int] = {}
    lowlink: Dict[SignedToken, int] = {}
    components: List[Tuple[SignedToken, ...]] = []

    def visit(node: SignedToken) -> None:
        nonlocal index
        indices[node] = lowlink[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)
        for neighbor in sorted(adjacency.get(node, set())):
            if neighbor not in indices:
                visit(neighbor)
                lowlink[node] = min(lowlink[node], lowlink[neighbor])
            elif neighbor in on_stack:
                lowlink[node] = min(lowlink[node], indices[neighbor])
        if lowlink[node] == indices[node]:
            component: List[SignedToken] = []
            while True:
                member = stack.pop()
                on_stack.remove(member)
                component.append(member)
                if member == node:
                    break
            components.append(tuple(sorted(component)))

    for node in sorted(set(nodes)):
        if node not in indices:
            visit(node)
    return components


def compute_topology(
    by_path: Mapping[str, Sequence[Occurrence]],
) -> Tuple[Dict[SignedToken, Tuple[int, int, int]], List[Tuple[SignedToken, SignedToken, int]]]:
    """Return signed token -> (SCC id, topological rank, SCC size)."""
    nodes: Set[SignedToken] = set()
    adjacency: Dict[SignedToken, Set[SignedToken]] = defaultdict(set)
    directed_support: Dict[Tuple[SignedToken, SignedToken], Set[str]] = defaultdict(set)
    for path_id, records in by_path.items():
        walk = normalized_walk(records)
        nodes.update(walk)
        for left, right in zip(walk, walk[1:]):
            adjacency[left].add(right)
            directed_support[(left, right)].add(path_id)

    components = tarjan_scc(nodes, adjacency)
    components.sort()
    component_of: Dict[SignedToken, int] = {}
    component_size: Dict[int, int] = {}
    for component_id, component in enumerate(components, 1):
        component_size[component_id] = len(component)
        for token in component:
            component_of[token] = component_id
    dag: Dict[int, Set[int]] = defaultdict(set)
    indegree = {component_id: 0 for component_id in component_size}
    for left, neighbors in adjacency.items():
        for right in neighbors:
            source, target = component_of[left], component_of[right]
            if source != target and target not in dag[source]:
                dag[source].add(target)
                indegree[target] += 1
    queue = sorted(component_id for component_id, degree in indegree.items() if degree == 0)
    ranks: Dict[int, int] = {}
    while queue:
        component_id = queue.pop(0)
        ranks[component_id] = len(ranks)
        for target in sorted(dag.get(component_id, set())):
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(target)
                queue.sort()
    topology = {
        token: (component_of[token], ranks[component_of[token]], component_size[component_of[token]])
        for token in nodes
    }
    directed_rows = [
        (left, right, len(paths))
        for (left, right), paths in sorted(directed_support.items())
    ]
    return topology, directed_rows


def collect_anchor_subpaths(
    by_path: Mapping[str, Sequence[Occurrence]],
    unit_by_path: Mapping[str, str],
    anchors: Set[str],
    stats: Mapping[str, BundleStats],
) -> Dict[BoundaryKey, Dict[SubpathKey, SubpathClass]]:
    """Collect only subpaths actually observed between consecutive anchors."""
    result: Dict[BoundaryKey, Dict[SubpathKey, SubpathClass]] = defaultdict(dict)
    for path_id, records in by_path.items():
        anchor_indices = [index for index, record in enumerate(records) if record.bundle_id in anchors]
        for left_index, right_index in zip(anchor_indices, anchor_indices[1:]):
            segment = records[left_index : right_index + 1]
            chain = tuple(record.token for record in segment)
            canonical_chain, orientation = canonicalize_tokens(chain)
            repeat_pattern = tuple(
                record.repeat_flag.upper() == "R" and stats[record.bundle_id].duplicated
                for record in segment
            )
            if orientation == "-":
                repeat_pattern = tuple(reversed(repeat_pattern))
            boundary = canonical_chain[0], canonical_chain[-1]
            interior = canonical_chain[1:-1]
            interior_repeat_pattern = repeat_pattern[1:-1]
            dup_tokens = tuple(
                token
                for token, duplicated in zip(interior, interior_repeat_pattern)
                if duplicated
            )
            duplication_signature, _ = canonicalize_tokens(dup_tokens)
            raw_dup_tokens = tuple(
                record.token
                for record in segment[1:-1]
                if record.repeat_flag.upper() == "R" and stats[record.bundle_id].duplicated
            )
            _, duplication_orientation = canonicalize_tokens(raw_dup_tokens)
            subpath_key = interior, interior_repeat_pattern
            classes = result[boundary]
            subpath = classes.setdefault(
                subpath_key,
                SubpathClass(
                    boundary,
                    interior,
                    interior_repeat_pattern,
                    dup_tokens,
                    duplication_signature,
                ),
            )
            observation = SegmentObservation(
                path_id=path_id,
                unit_id=unit_by_path[path_id],
                left_index=left_index,
                right_index=right_index,
                boundary=boundary,
                interior=interior,
                canonical_orientation=orientation,
                duplication_orientation=duplication_orientation,
            )
            subpath.observations.append(observation)
            subpath.support_units.add(observation.unit_id)
            subpath.support_paths.add(path_id)
    return dict(result)


def group_subgraphs_by_shared_bundles(
    subgraphs: Mapping[BoundaryKey, Mapping[SubpathKey, SubpathClass]],
    stats: Mapping[str, BundleStats],
    min_shared_bp: float,
    topology: Mapping[SignedToken, Tuple[int, int, int]],
) -> Tuple[Dict[BoundaryKey, int], List[SubgraphGroup], List[SharedBoundaryPair]]:
    """Join boundary subgraphs whose shared interior bundle span is sufficient.

    Each shared bundle ID contributes its mean observed occurrence length once
    to a boundary-pair score. Qualifying pairs are connected, and transitive
    connected components become joint trim groups.
    """
    boundaries = sorted(
        subgraphs,
        key=lambda boundary: (
            topology.get(boundary[0], (0, 10**9, 0))[1],
            topology.get(boundary[1], (0, 10**9, 0))[1],
            boundary,
        ),
    )
    boundary_rank = {boundary: index for index, boundary in enumerate(boundaries)}
    bundle_ids_by_boundary = {
        boundary: {
            token[0]
            for subpath in subgraphs[boundary].values()
            for token in subpath.interior
        }
        for boundary in boundaries
    }
    mean_length_by_bundle = {
        bundle_id: float(statistics.mean(stat.lengths))
        for bundle_id, stat in stats.items()
        if stat.lengths
    }
    parent = {boundary: boundary for boundary in boundaries}

    def find(boundary: BoundaryKey) -> BoundaryKey:
        root = boundary
        while parent[root] != root:
            root = parent[root]
        while parent[boundary] != boundary:
            next_boundary = parent[boundary]
            parent[boundary] = root
            boundary = next_boundary
        return root

    def union(left: BoundaryKey, right: BoundaryKey) -> None:
        left_root, right_root = find(left), find(right)
        if left_root == right_root:
            return
        if boundary_rank[left_root] <= boundary_rank[right_root]:
            parent[right_root] = left_root
        else:
            parent[left_root] = right_root

    pair_rows: List[SharedBoundaryPair] = []
    for left_index, left in enumerate(boundaries):
        for right in boundaries[left_index + 1 :]:
            shared = tuple(sorted(bundle_ids_by_boundary[left] & bundle_ids_by_boundary[right]))
            if not shared:
                continue
            shared_bp = sum(mean_length_by_bundle[bundle_id] for bundle_id in shared)
            merged = shared_bp >= min_shared_bp
            pair_rows.append(SharedBoundaryPair(left, right, shared, shared_bp, merged))
            if merged:
                union(left, right)

    boundaries_by_root: Dict[BoundaryKey, List[BoundaryKey]] = defaultdict(list)
    for boundary in boundaries:
        boundaries_by_root[find(boundary)].append(boundary)
    ordered_groups = sorted(
        (tuple(items) for items in boundaries_by_root.values()),
        key=lambda items: min(boundary_rank[boundary] for boundary in items),
    )
    group_by_boundary: Dict[BoundaryKey, int] = {}
    groups: List[SubgraphGroup] = []
    for group_id, group_boundaries in enumerate(ordered_groups, 1):
        interior_ids = tuple(
            sorted(set().union(*(bundle_ids_by_boundary[boundary] for boundary in group_boundaries)))
        )
        groups.append(
            SubgraphGroup(
                group_id=group_id,
                boundaries=group_boundaries,
                interior_bundle_ids=interior_ids,
                average_total_bp=sum(mean_length_by_bundle[bundle_id] for bundle_id in interior_ids),
            )
        )
        for boundary in group_boundaries:
            group_by_boundary[boundary] = group_id
    return group_by_boundary, groups, pair_rows


def is_trimmable_token(
    token: SignedToken,
    stats: Mapping[str, BundleStats],
    trim_rare_long: bool,
) -> bool:
    stat = stats[token[0]]
    return stat.short or (trim_rare_long and stat.rare)


def deletion_trim_cost(
    source: Sequence[SignedToken],
    target: Sequence[SignedToken],
    stats: Mapping[str, BundleStats],
    trim_rare_long: bool,
) -> Optional[Tuple[int, int]]:
    """Cost of reducing source to target by deleting trimmable tokens only."""
    source_tuple, target_tuple = tuple(source), tuple(target)
    infinity = (10**30, 10**30)
    dp: List[List[Tuple[int, int]]] = [
        [infinity for _ in range(len(target_tuple) + 1)]
        for _ in range(len(source_tuple) + 1)
    ]
    dp[0][0] = (0, 0)
    for i in range(1, len(source_tuple) + 1):
        token = source_tuple[i - 1]
        if is_trimmable_token(token, stats, trim_rare_long) and dp[i - 1][0] != infinity:
            dp[i][0] = (
                dp[i - 1][0][0] + max(1, round(stats[token[0]].median_length)),
                dp[i - 1][0][1] + 1,
            )
    for i in range(1, len(source_tuple) + 1):
        token = source_tuple[i - 1]
        for j in range(1, len(target_tuple) + 1):
            if token == target_tuple[j - 1]:
                dp[i][j] = min(dp[i][j], dp[i - 1][j - 1])
            if is_trimmable_token(token, stats, trim_rare_long) and dp[i - 1][j] != infinity:
                deletion = (
                    dp[i - 1][j][0] + max(1, round(stats[token[0]].median_length)),
                    dp[i - 1][j][1] + 1,
                )
                dp[i][j] = min(dp[i][j], deletion)
    return None if dp[-1][-1] == infinity else dp[-1][-1]


def deleted_tokens(
    source: Sequence[SignedToken], target: Sequence[SignedToken]
) -> Tuple[SignedToken, ...]:
    """Return the multiset of source tokens absent from its target subsequence."""
    remaining = Counter(target)
    deleted: List[SignedToken] = []
    for token in source:
        if remaining[token] > 0:
            remaining[token] -= 1
        else:
            deleted.append(token)
    return tuple(deleted)


def project_repeat_pattern_to_subsequence(
    source: Sequence[SignedToken],
    source_repeat_pattern: Sequence[bool],
    target: Sequence[SignedToken],
) -> RepeatPattern:
    """Project source R/U states onto a retained target subsequence.

    Trimming is based only on token deletion. When repeated token identities
    make several alignments possible, the deterministic leftmost alignment is
    used for the later duplication-aware merge decision.
    """
    if len(source) != len(source_repeat_pattern):
        raise ValueError("token and repeat-pattern lengths differ")
    projected: List[bool] = []
    source_index = 0
    for target_token in target:
        while source_index < len(source) and source[source_index] != target_token:
            source_index += 1
        if source_index == len(source):
            raise ValueError("trim target is not a subsequence of its source")
        projected.append(source_repeat_pattern[source_index])
        source_index += 1
    return tuple(projected)


def classify_subpaths(
    subgraphs: MutableMapping[BoundaryKey, Dict[SubpathKey, SubpathClass]],
    stats: Mapping[str, BundleStats],
    min_path_count: int,
    min_path_frequency: float,
    max_trim_bp: int,
    max_trim_bundles: int,
    trim_rare_long: bool,
    preserve_common_paths: bool,
    topology: Mapping[SignedToken, Tuple[int, int, int]],
    group_by_boundary: Mapping[BoundaryKey, int],
) -> List[GraphBundle]:
    retained_classes: List[SubpathClass] = []
    ordered_boundaries = sorted(
        subgraphs,
        key=lambda boundary: (
            topology.get(boundary[0], (0, 10**9, 0))[1],
            topology.get(boundary[1], (0, 10**9, 0))[1],
            boundary,
        ),
    )
    boundary_order = {boundary: index for index, boundary in enumerate(ordered_boundaries)}
    classes_by_group: Dict[int, List[SubpathClass]] = defaultdict(list)
    for boundary in ordered_boundaries:
        classes = list(subgraphs[boundary].values())
        spanning_units = set().union(*(subpath.support_units for subpath in classes))
        required = max(min_path_count, math.ceil(min_path_frequency * len(spanning_units)))
        for subpath in classes:
            subpath.common = len(subpath.support_units) >= required
        classes_by_group[group_by_boundary[boundary]].extend(classes)

    ordered_group_ids = sorted(
        classes_by_group,
        key=lambda group_id: min(
            boundary_order[subpath.boundary] for subpath in classes_by_group[group_id]
        ),
    )
    for group_id in ordered_group_ids:
        classes = classes_by_group[group_id]
        classes.sort(
            key=lambda subpath: (
                -len(subpath.support_units),
                -len(subpath.support_paths),
                len(subpath.interior),
                boundary_order[subpath.boundary],
                subpath.interior,
                subpath.repeat_pattern,
            )
        )
        local_retained: List[SubpathClass] = []
        for index, subpath in enumerate(classes):
            best_target: Optional[SubpathClass] = None
            best_cost: Optional[Tuple[int, int]] = None
            for target in local_retained:
                cost = deletion_trim_cost(subpath.interior, target.interior, stats, trim_rare_long)
                if cost is None or cost[0] > max_trim_bp or cost[1] > max_trim_bundles:
                    continue
                target_key = (
                    cost, -len(target.support_units), target.interior,
                    boundary_order[target.boundary], target.repeat_pattern,
                )
                if best_target is None:
                    best_target, best_cost = target, cost
                else:
                    best_key = (
                        best_cost, -len(best_target.support_units), best_target.interior,
                        boundary_order[best_target.boundary], best_target.repeat_pattern,
                    )
                    if target_key < best_key:
                        best_target, best_cost = target, cost
            if best_target is not None and (not subpath.common or not preserve_common_paths):
                subpath.action = "trimmed"
                subpath.target_interior = best_target.interior
                subpath.target_key = best_target.key
                subpath.target_boundary = best_target.boundary
                subpath.trim_bp, subpath.trim_bundle_count = best_cost  # type: ignore[misc]
            else:
                if index == 0:
                    subpath.action = "retain_backbone"
                elif subpath.common:
                    subpath.action = "retain_common"
                else:
                    subpath.action = "retain_protected"
                subpath.target_interior = subpath.interior
                subpath.target_key = subpath.key
                subpath.target_boundary = subpath.boundary
                subpath.trim_bp = subpath.trim_bundle_count = 0
                local_retained.append(subpath)
                retained_classes.append(subpath)

    retained_lookup = {
        (subpath.boundary, subpath.key): subpath for subpath in retained_classes
    }
    all_subpaths = [
        subpath
        for boundary in ordered_boundaries
        for subpath in subgraphs[boundary].values()
    ]
    for subpath in all_subpaths:
        target = retained_lookup[(subpath.target_boundary, subpath.target_key)]  # type: ignore[index]
        subpath.merged_repeat_pattern = project_repeat_pattern_to_subsequence(
            subpath.interior,
            subpath.repeat_pattern,
            target.interior,
        )
        subpath.merged_dup_tokens = tuple(
            token
            for token, duplicated in zip(target.interior, subpath.merged_repeat_pattern)
            if duplicated
        )
        subpath.merged_duplication_signature, _ = canonicalize_tokens(
            subpath.merged_dup_tokens
        )

    used_bundle_ids = set(stats)
    numeric_bundle_ids = [int(bundle_id) for bundle_id in used_bundle_ids if bundle_id.isdigit()]
    next_bundle_id = max(numeric_bundle_ids, default=-1) + 1
    identity_to_id: Dict[Tuple[object, ...], str] = {}

    def graph_identity(subpath: SubpathClass) -> Tuple[object, ...]:
        if subpath.merged_duplication_signature:
            return "dup", subpath.merged_duplication_signature
        return "path", group_by_boundary[subpath.boundary], subpath.target_interior

    retained_ids = {id(subpath) for subpath in retained_classes}
    assignment_order = retained_classes + [
        subpath for subpath in all_subpaths if id(subpath) not in retained_ids
    ]
    for subpath in assignment_order:
        if not subpath.target_interior:
            continue
        identity = graph_identity(subpath)
        if identity not in identity_to_id:
            while str(next_bundle_id) in used_bundle_ids:
                next_bundle_id += 1
            identity_to_id[identity] = str(next_bundle_id)
            used_bundle_ids.add(str(next_bundle_id))
            next_bundle_id += 1
        subpath.new_bundle_id = identity_to_id[identity]

    assigned_units: Dict[str, Set[str]] = defaultdict(set)
    assigned_paths: Dict[str, Set[str]] = defaultdict(set)
    classes_by_id: Dict[str, List[SubpathClass]] = defaultdict(list)
    for subpath in all_subpaths:
        if subpath.new_bundle_id is None:
            continue
        assigned_units[subpath.new_bundle_id].update(subpath.support_units)
        assigned_paths[subpath.new_bundle_id].update(subpath.support_paths)
        classes_by_id[subpath.new_bundle_id].append(subpath)

    graph_bundles: List[GraphBundle] = []
    for new_id in sorted(classes_by_id, key=int):
        group = classes_by_id[new_id]
        representative = group[0]
        duplicated = bool(representative.merged_duplication_signature)
        graph_bundles.append(
            GraphBundle(
                new_id,
                tuple(sorted({subpath.boundary for subpath in group})),
                (
                    representative.merged_duplication_signature
                    if duplicated
                    else representative.target_interior or ()
                ),
                len(assigned_units[new_id]),
                len(assigned_paths[new_id]),
                "dup_order_merged" if duplicated else representative.action,
                duplicated,
            )
        )
    return graph_bundles


def target_subpath(
    subgraphs: Mapping[BoundaryKey, Mapping[SubpathKey, SubpathClass]],
    subpath: SubpathClass,
) -> SubpathClass:
    if subpath.target_key is None or subpath.target_boundary is None:
        raise ValueError("subpath target has not been assigned")
    return subgraphs[subpath.target_boundary][subpath.target_key]


def project_graph_bundles(
    by_path: Mapping[str, Sequence[Occurrence]],
    anchors: Set[str],
    subgraphs: Mapping[BoundaryKey, Mapping[SubpathKey, SubpathClass]],
    stats: Mapping[str, BundleStats],
    keep_unmerged_short: bool,
) -> Dict[str, List[OutputOccurrence]]:
    """Project graph bundles, retaining non-rare short records for coverage trials.

    ``keep_unmerged_short`` is retained for API compatibility; short deletion
    is not performed by the graph-merge workflow.
    """
    outputs: Dict[str, List[OutputOccurrence]] = {}
    observation_to_class: Dict[Tuple[str, int, int], SubpathClass] = {}
    observation_lookup: Dict[Tuple[str, int, int], SegmentObservation] = {}
    for classes in subgraphs.values():
        for subpath in classes.values():
            for observation in subpath.observations:
                key = observation.path_id, observation.left_index, observation.right_index
                observation_to_class[key] = subpath
                observation_lookup[key] = observation

    for path_id, records in by_path.items():
        consumed: Set[int] = set()
        path_outputs: List[OutputOccurrence] = []
        anchor_indices = [index for index, record in enumerate(records) if record.bundle_id in anchors]
        for left_index, right_index in zip(anchor_indices, anchor_indices[1:]):
            key = path_id, left_index, right_index
            subpath = observation_to_class[key]
            observation = observation_lookup[key]
            interior_indices = tuple(range(left_index + 1, right_index))
            consumed.update(interior_indices)
            if subpath.new_bundle_id is None:
                continue
            if not interior_indices:
                raise ValueError(f"cannot project {subpath.new_bundle_id} onto empty interior in {path_id}")
            target = target_subpath(subgraphs, subpath)
            target_members = (
                subpath.merged_dup_tokens
                if subpath.merged_duplication_signature
                else target.interior
            )
            target_ids = Counter(token[0] for token in target_members)
            observed_ids = [records[index].bundle_id for index in interior_indices]
            remaining = target_ids.copy()
            absorbed: List[str] = []
            for bundle_id in observed_ids:
                if remaining[bundle_id] > 0:
                    remaining[bundle_id] -= 1
                else:
                    absorbed.append(bundle_id)
            if subpath.merged_duplication_signature:
                raw_dup_tokens = (
                    subpath.merged_dup_tokens
                    if observation.canonical_orientation == "+"
                    else reverse_complement_tokens(subpath.merged_dup_tokens)
                )
                _, output_direction = canonicalize_tokens(raw_dup_tokens)
            else:
                output_direction = observation.canonical_orientation
            path_outputs.append(
                OutputOccurrence(
                    path_id=path_id,
                    start=records[interior_indices[0]].start,
                    end=records[interior_indices[-1]].end,
                    new_bundle_id=subpath.new_bundle_id,
                    direction=output_direction,
                    source="graph_bundle" if subpath.action != "trimmed" else "trimmed_to_graph_bundle",
                    original_bundle_ids=tuple(observed_ids),
                    absorbed_bundle_ids=tuple(absorbed),
                    boundary_text=boundary_text(subpath.boundary),
                    repeat_flag="R" if subpath.merged_duplication_signature else "U",
                )
            )

        for index, record in enumerate(records):
            if index in consumed:
                continue
            stat = stats[record.bundle_id]
            if stat.anchor:
                source = (
                    "long_dup_breakpoint" if stat.long_dup_breakpoint else "anchor"
                )
            elif stat.rare:
                continue
            else:
                source = "retained_short" if stat.short else "retained_long"
            path_outputs.append(
                OutputOccurrence(
                    path_id,
                    record.start,
                    record.end,
                    record.bundle_id,
                    record.direction,
                    source,
                    (record.bundle_id,),
                    bundle_size=record.bundle_size,
                    bundle_start=record.bundle_start,
                    bundle_end=record.bundle_end,
                    repeat_flag=record.repeat_flag,
                )
            )
        path_outputs.sort(key=lambda item: (item.start, item.end, item.new_bundle_id))
        outputs[path_id] = path_outputs
    return outputs


def resolve_output_coordinate_overlaps(
    outputs: Mapping[str, Sequence[OutputOccurrence]],
    seed: int,
) -> Tuple[Dict[str, List[OutputOccurrence]], List[CoordinateAdjustment]]:
    """Assign each adjacent overlap to one side using a reproducible random cut.

    For an overlap ``right.start < left.end``, a pseudo-random boundary is
    selected from the inclusive interval ``[right.start, left.end]``. The left
    end and right start are both set to that boundary, preserving the covered
    union while ensuring that the two final BED intervals do not overlap.
    """
    adjusted_outputs: Dict[str, List[OutputOccurrence]] = {}
    adjustments: List[CoordinateAdjustment] = []
    for path_id, path_items in sorted(outputs.items()):
        items = sorted(path_items, key=lambda item: (item.start, item.end, item.new_bundle_id))
        adjusted = list(items)
        for index in range(len(adjusted) - 1):
            left, right = adjusted[index], adjusted[index + 1]
            if right.start >= left.end:
                continue
            lower = max(right.start, left.start + 1)
            upper = min(left.end, right.end - 1)
            if lower > upper:
                raise ValueError(
                    f"cannot resolve coordinate overlap on {path_id}: "
                    f"{left.new_bundle_id}:{left.start}-{left.end} and "
                    f"{right.new_bundle_id}:{right.start}-{right.end}"
                )
            digest = hashlib.blake2b(
                (
                    f"{seed}\t{path_id}\t{index}\t{left.new_bundle_id}\t"
                    f"{right.new_bundle_id}\t{left.end}\t{right.start}"
                ).encode("utf-8"),
                digest_size=8,
            ).digest()
            split_coordinate = lower + int.from_bytes(digest, "big") % (upper - lower + 1)
            overlap_bp = left.end - right.start
            adjustments.append(
                CoordinateAdjustment(
                    path_id=path_id,
                    left_bundle_id=left.new_bundle_id,
                    right_bundle_id=right.new_bundle_id,
                    original_left_end=left.end,
                    original_right_start=right.start,
                    overlap_bp=overlap_bp,
                    split_coordinate=split_coordinate,
                    left_trim_bp=left.end - split_coordinate,
                    right_trim_bp=split_coordinate - right.start,
                )
            )
            adjusted[index] = replace(left, end=split_coordinate)
            adjusted[index + 1] = replace(right, start=split_coordinate)

        for left, right in zip(adjusted, adjusted[1:]):
            if left.end > right.start:
                raise ValueError(
                    f"coordinate overlap remains on {path_id}: "
                    f"{left.new_bundle_id}:{left.start}-{left.end} and "
                    f"{right.new_bundle_id}:{right.start}-{right.end}"
                )
        adjusted_outputs[path_id] = adjusted
    return adjusted_outputs, adjustments


def quantile(values: Sequence[float], probability: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return float(ordered[lower])
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def union_intervals(intervals: Iterable[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    """Return a half-open interval union, without counting overlapping copies twice."""
    result: List[Tuple[int, int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if result and start <= result[-1][1]:
            result[-1] = result[-1][0], max(result[-1][1], end)
        else:
            result.append((start, end))
    return tuple(result)


def covered_interval_bp(
    records: Sequence[OutputOccurrence],
    allowed_intervals: Sequence[Tuple[int, int]],
) -> int:
    covered = union_intervals((record.start, record.end) for record in records)
    total = 0
    left = right = 0
    while left < len(covered) and right < len(allowed_intervals):
        start = max(covered[left][0], allowed_intervals[right][0])
        end = min(covered[left][1], allowed_intervals[right][1])
        total += max(0, end - start)
        if covered[left][1] <= allowed_intervals[right][1]:
            left += 1
        else:
            right += 1
    return total


def read_sequence_lengths(path: Path) -> Dict[str, int]:
    """Read FASTA lengths without loading the population sequences into memory."""
    lengths: Dict[str, int] = {}
    current: Optional[str] = None
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                fields = stripped[1:].split()
                if not fields:
                    raise ValueError(f"{path}:{line_number}: empty FASTA sequence ID")
                current = fields[0]
                if current in lengths:
                    raise ValueError(f"{path}:{line_number}: duplicate FASTA sequence ID {current}")
                lengths[current] = 0
            else:
                if current is None:
                    raise ValueError(f"{path}:{line_number}: sequence before first FASTA header")
                lengths[current] += len("".join(stripped.split()))
    if not lengths:
        raise ValueError(f"{path}: no FASTA records found")
    return lengths


def build_coverage_intervals(
    input_by_path: Mapping[str, Sequence[Occurrence]],
    fasta: Optional[Path],
    denominator: str = "auto",
) -> Tuple[Dict[str, Tuple[Tuple[int, int], ...]], str]:
    """Freeze the denominator before TRF, rare filtering, merging, or deletions.

    Only paths present in the input BED are the evaluated sequence population;
    sample-map aggregation is deliberately not used for sequence coverage.
    """
    if denominator not in {"auto", "fasta-length", "input-bed-union"}:
        raise ValueError(f"unknown coverage denominator: {denominator}")
    mode = ("fasta-length" if fasta is not None else "input-bed-union") if denominator == "auto" else denominator
    if mode == "input-bed-union":
        return {
            path_id: union_intervals((record.start, record.end) for record in records)
            for path_id, records in input_by_path.items()
        }, mode
    if fasta is None:
        raise ValueError("--coverage-denominator fasta-length requires --fasta")
    lengths = read_sequence_lengths(fasta)
    result: Dict[str, Tuple[Tuple[int, int], ...]] = {}
    for path_id, records in input_by_path.items():
        if path_id not in lengths:
            raise ValueError(f"{fasta}: input BED sequence {path_id} is missing from FASTA")
        length = lengths[path_id]
        if length <= 0:
            raise ValueError(f"{fasta}: input BED sequence {path_id} has zero length")
        if any(record.end > length for record in records):
            raise ValueError(f"input BED coordinates exceed FASTA length for {path_id}: {length}")
        result[path_id] = ((0, length),)
    return result, mode


def prune_unmerged_short_by_coverage(
    outputs: Mapping[str, Sequence[OutputOccurrence]],
    stats: MutableMapping[str, BundleStats],
    coverage_intervals: Mapping[str, Sequence[Tuple[int, int]]],
    min_coverage: float = 0.90,
    max_low_coverage_fraction: float = 0.05,
    keep_unmerged_short: bool = False,
) -> Tuple[Dict[str, List[OutputOccurrence]], List[Sequence[object]], List[Sequence[object]], Dict[str, object]]:
    """Greedily delete unmerged short types with cumulative sequence protection.

    All standalone occurrences of a candidate type are tested together. Anchors,
    rare types, and intervals already incorporated into graph bundles are not
    candidates. A failed candidate leaves the working library unchanged. The
    denominator and evaluated paths stay fixed for every trial.
    """
    if not 0 <= min_coverage <= 1 or not 0 <= max_low_coverage_fraction <= 1:
        raise ValueError("sequence coverage thresholds must be in [0, 1]")
    unknown = set(outputs) - set(coverage_intervals)
    if unknown:
        raise ValueError(f"coverage denominator is missing paths: {', '.join(sorted(unknown)[:10])}")
    allowed = {path_id: union_intervals(intervals) for path_id, intervals in coverage_intervals.items()}
    denominators = {path_id: sum(end - start for start, end in intervals) for path_id, intervals in allowed.items()}
    if any(length <= 0 for length in denominators.values()):
        raise ValueError("sequence coverage denominators must be positive")
    working = {path_id: list(outputs.get(path_id, ())) for path_id in allowed}
    current_bp = {path_id: covered_interval_bp(records, allowed[path_id]) for path_id, records in working.items()}
    baseline_bp = current_bp.copy()
    total_paths = len(working)
    def is_low(path_id: str, covered_bp: int) -> bool:
        return covered_bp / denominators[path_id] < min_coverage
    low_count = sum(is_low(path_id, bp) for path_id, bp in current_bp.items())
    baseline_low_count = low_count
    candidate_paths: Dict[str, Set[str]] = defaultdict(set)
    for path_id, records in working.items():
        for record in records:
            stat = stats.get(record.new_bundle_id)
            if (record.source == "retained_short" and stat is not None
                    and stat.short and not stat.rare and not stat.anchor):
                candidate_paths[record.new_bundle_id].add(path_id)
    candidate_ids = sorted(candidate_paths, key=lambda bundle_id: (
        stats[bundle_id].median_length, not bundle_id.isdigit(),
        int(bundle_id) if bundle_id.isdigit() else bundle_id, bundle_id,
    ))
    decision_rows: List[Sequence[object]] = []
    deleted_count = protected_count = deleted_occurrences = 0
    for attempt_index, bundle_id in enumerate(candidate_ids, 1):
        stat = stats[bundle_id]
        trial_records: Dict[str, List[OutputOccurrence]] = {}
        trial_bp: Dict[str, int] = {}
        removed_count = 0
        trial_low_count = low_count
        for path_id in sorted(candidate_paths[bundle_id]):
            remaining = [record for record in working[path_id]
                         if not (record.new_bundle_id == bundle_id and record.source == "retained_short")]
            removed_count += len(working[path_id]) - len(remaining)
            trial_records[path_id] = remaining
            trial_bp[path_id] = covered_interval_bp(remaining, allowed[path_id])
            trial_low_count += int(is_low(path_id, trial_bp[path_id])) - int(is_low(path_id, current_bp[path_id]))
        trial_fraction = trial_low_count / total_paths if total_paths else 0.0
        lost_bp = sum(current_bp[path_id] - bp for path_id, bp in trial_bp.items())
        before_low_count = low_count
        if keep_unmerged_short:
            action = "retain_keep_unmerged_short"
        elif trial_fraction > max_low_coverage_fraction:
            action = "retain_coverage_protected"
            protected_count += 1
        else:
            action = "delete_coverage_safe"
            deleted_count += 1
            deleted_occurrences += removed_count
            working.update(trial_records)
            current_bp.update(trial_bp)
            low_count = trial_low_count
        stat.short_deletion_action = action
        if stat.action == "retain_short":
            stat.action = {
                "delete_coverage_safe": "drop_unmerged_short_coverage_safe",
                "retain_coverage_protected": "retain_short_coverage_protected",
                "retain_keep_unmerged_short": "retain_short",
            }[action]
        decision_rows.append((
            attempt_index, bundle_id, f"{stat.median_length:.3f}", removed_count,
            len(trial_records), lost_bp, before_low_count, trial_low_count,
            f"{trial_fraction:.9f}", low_count, action,
            ",".join(path_id for path_id in sorted(trial_bp) if is_low(path_id, trial_bp[path_id])) or ".",
        ))
    coverage_rows = [(
        path_id, denominators[path_id], baseline_bp[path_id],
        f"{baseline_bp[path_id] / denominators[path_id]:.9f}", current_bp[path_id],
        f"{current_bp[path_id] / denominators[path_id]:.9f}",
        baseline_bp[path_id] - current_bp[path_id], int(is_low(path_id, current_bp[path_id])),
    ) for path_id in sorted(working)]
    return working, decision_rows, coverage_rows, {
        "short_deletion_candidate_bundle_count": len(candidate_ids),
        "short_deleted_bundle_count": deleted_count,
        "short_coverage_protected_bundle_count": protected_count,
        "short_deleted_occurrence_count": deleted_occurrences,
        "coverage_sequence_count": total_paths,
        "coverage_baseline_low_sequence_count": baseline_low_count,
        "coverage_final_low_sequence_count": low_count,
        "coverage_baseline_low_sequence_fraction": baseline_low_count / total_paths if total_paths else 0.0,
        "coverage_final_low_sequence_fraction": low_count / total_paths if total_paths else 0.0,
        "coverage_baseline_covered_bp": sum(baseline_bp.values()),
        "coverage_final_covered_bp": sum(current_bp.values()),
    }


def median_absolute_deviation(values: Sequence[int]) -> float:
    if not values:
        return 0.0
    median = statistics.median(values)
    return float(statistics.median(abs(value - median) for value in values))


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def assign_bundle_actions(
    stats: MutableMapping[str, BundleStats],
    graph_bundles: Sequence[GraphBundle],
    subgraphs: Mapping[BoundaryKey, Mapping[SubpathKey, SubpathClass]],
    keep_unmerged_short: bool,
) -> None:
    merged_members = {
        token[0] for bundle in graph_bundles for token in bundle.canonical_interior
    }
    absorbed_members = {
        token[0]
        for boundary_classes in subgraphs.values()
        for subpath in boundary_classes.values()
        if subpath.new_bundle_id is not None and subpath.target_key is not None
        for token in deleted_tokens(
            subpath.interior,
            (
                subpath.merged_dup_tokens
                if subpath.merged_duplication_signature
                else target_subpath(subgraphs, subpath).interior
            ),
        )
    }
    for bundle_id, stat in stats.items():
        if stat.long_dup_breakpoint:
            stat.action = "retain_long_dup_breakpoint"
        elif stat.anchor:
            stat.action = "retain_anchor"
        elif bundle_id in merged_members:
            stat.action = "graph_bundle_member"
        elif bundle_id in absorbed_members:
            stat.action = "absorbed_trimmed_branch"
        elif stat.rare:
            stat.action = "drop_rare_unrepresented"
        elif stat.short:
            stat.action = "retain_short"
        else:
            stat.action = "retain_long"


def write_outputs(
    prefix: Path,
    outputs: Mapping[str, Sequence[OutputOccurrence]],
    stats: Mapping[str, BundleStats],
    graph_edges: Mapping[PortEdgeKey, PopulationEdge],
    topology: Mapping[SignedToken, Tuple[int, int, int]],
    directed_edges: Sequence[Tuple[SignedToken, SignedToken, int]],
    subgraphs: Mapping[BoundaryKey, Mapping[SubpathKey, SubpathClass]],
    graph_bundles: Sequence[GraphBundle],
    group_by_boundary: Mapping[BoundaryKey, int],
    subgraph_groups: Sequence[SubgraphGroup],
    shared_boundary_pairs: Sequence[SharedBoundaryPair],
    long_dup_breakpoints: Sequence[LongDupBreakpoint],
    coordinate_adjustments: Sequence[CoordinateAdjustment],
    trf_premerge_report_rows: Sequence[Sequence[object]],
    trf_premerge_summary_rows: Sequence[Sequence[object]],
    trf_report_rows: Sequence[Sequence[object]],
    trf_bundle_summary_rows: Sequence[Sequence[object]],
    summary: Mapping[str, object],
    short_deletion_rows: Sequence[Sequence[object]] = (),
    sequence_coverage_rows: Sequence[Sequence[object]] = (),
) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with Path(f"{prefix}.bed").open("w", encoding="utf-8") as handle:
        handle.write(
            "# format: path_id\tstart\tend\t"
            "bundle_id:bundle_size:direction:bundle_pos_start:bundle_pos_end:R|U\n"
        )
        for path_id in sorted(outputs):
            for item in outputs[path_id]:
                handle.write(
                    f"{path_id}\t{item.start}\t{item.end}\t{pgr_output_annotation(item)}\n"
                )
    with Path(f"{prefix}.bed6").open("w", encoding="utf-8") as handle:
        handle.write("#path_id\tstart\tend\tnew_bundle_id\tscore\tstrand\n")
        for path_id in sorted(outputs):
            for item in outputs[path_id]:
                handle.write(
                    f"{path_id}\t{item.start}\t{item.end}\t{item.new_bundle_id}\t0\t{item.direction}\n"
                )
    write_tsv(
        Path(f"{prefix}.paths.tsv"),
        ("path_id", "block_count", "structural_haplotype"),
        (
            (path_id, len(items), ",".join(f"{item.direction}{item.new_bundle_id}" for item in items))
            for path_id, items in sorted(outputs.items())
        ),
    )
    write_tsv(
        Path(f"{prefix}.occurrences.tsv"),
        (
            "path_id", "start", "end", "length", "new_bundle_id", "direction", "source",
            "original_bundle_ids", "absorbed_bundle_ids", "boundary", "trf_action",
            "trf_repeat_bp", "trf_coverage",
        ),
        (
            (
                path_id, item.start, item.end, item.end - item.start, item.new_bundle_id,
                item.direction, item.source, ",".join(item.original_bundle_ids),
                ",".join(item.absorbed_bundle_ids), item.boundary_text,
                item.trf_action, item.trf_repeat_bp, f"{item.trf_coverage:.6f}",
            )
            for path_id, items in sorted(outputs.items())
            for item in items
        ),
    )
    write_tsv(
        Path(f"{prefix}.coordinate_adjustments.tsv"),
        (
            "path_id", "left_bundle_id", "right_bundle_id", "original_left_end",
            "original_right_start", "overlap_bp", "split_coordinate",
            "left_trim_bp", "right_trim_bp",
        ),
        (
            (
                adjustment.path_id, adjustment.left_bundle_id,
                adjustment.right_bundle_id, adjustment.original_left_end,
                adjustment.original_right_start, adjustment.overlap_bp,
                adjustment.split_coordinate, adjustment.left_trim_bp,
                adjustment.right_trim_bp,
            )
            for adjustment in coordinate_adjustments
        ),
    )
    write_tsv(
        Path(f"{prefix}.short_bundle_deletions.tsv"),
        (
            "attempt_index", "bundle_id", "median_length_bp", "candidate_occurrence_count",
            "affected_sequence_count", "trial_lost_coverage_bp", "low_sequence_count_before",
            "trial_low_sequence_count", "trial_low_sequence_fraction", "low_sequence_count_after",
            "action", "trial_low_sequences",
        ),
        short_deletion_rows,
    )
    write_tsv(
        Path(f"{prefix}.sequence_coverage.tsv"),
        (
            "path_id", "denominator_bp", "baseline_covered_bp", "baseline_coverage",
            "final_covered_bp", "final_coverage", "lost_coverage_bp", "below_min_coverage",
        ),
        sequence_coverage_rows,
    )
    write_tsv(
        Path(f"{prefix}.graph_nodes.tsv"),
        (
            "bundle_id", "occurrence_count", "path_count", "population_unit_count", "frequency",
            "median_length_bp", "single_copy_path_fraction", "repeat_occurrence_count",
            "repeat_path_count", "repeat_population_unit_count", "repeat_frequency", "duplicated",
            "rare", "short", "anchor", "long_dup_breakpoint", "action", "short_deletion_action",
        ),
        (
            (
                bundle_id, stat.occurrence_count, len(stat.paths), len(stat.units),
                f"{stat.frequency:.6f}", f"{stat.median_length:.1f}",
                f"{stat.single_copy_path_fraction:.6f}", stat.repeat_count,
                len(stat.repeat_paths), len(stat.repeat_units), f"{stat.repeat_frequency:.6f}",
                int(stat.duplicated), int(stat.rare), int(stat.short), int(stat.anchor),
                int(stat.long_dup_breakpoint), stat.action, stat.short_deletion_action or ".",
            )
            for bundle_id, stat in sorted(stats.items())
        ),
    )
    write_tsv(
        Path(f"{prefix}.long_dup_breakpoints.tsv"),
        (
            "bundle_id", "median_length_bp", "occurrence_count", "path_count",
            "population_unit_count", "repeat_occurrence_count", "repeat_path_count",
            "repeat_population_unit_count", "repeat_frequency",
            "candidate_graph_bundle_ids", "candidate_boundaries", "action",
        ),
        (
            (
                breakpoint.bundle_id,
                f"{stats[breakpoint.bundle_id].median_length:.1f}",
                stats[breakpoint.bundle_id].occurrence_count,
                len(stats[breakpoint.bundle_id].paths),
                len(stats[breakpoint.bundle_id].units),
                stats[breakpoint.bundle_id].repeat_count,
                len(stats[breakpoint.bundle_id].repeat_paths),
                len(stats[breakpoint.bundle_id].repeat_units),
                f"{stats[breakpoint.bundle_id].repeat_frequency:.6f}",
                ",".join(breakpoint.candidate_graph_bundle_ids),
                ";".join(boundary_text(boundary) for boundary in breakpoint.candidate_boundaries),
                stats[breakpoint.bundle_id].action,
            )
            for breakpoint in long_dup_breakpoints
        ),
    )
    write_tsv(
        Path(f"{prefix}.graph_edges.tsv"),
        (
            "bundle_a", "port_a", "bundle_b", "port_b", "occurrence_count", "path_count",
            "population_unit_count", "mean_gap_bp", "median_gap_bp", "gap_mad_bp",
            "gap_q10_bp", "gap_q90_bp", "supporting_paths",
        ),
        (
            (
                key[0][0], key[0][1], key[1][0], key[1][1], len(edge.observations),
                len(edge.path_ids), len(edge.unit_ids), f"{statistics.mean(edge.gaps):.3f}",
                f"{statistics.median(edge.gaps):.3f}", f"{median_absolute_deviation(edge.gaps):.3f}",
                f"{quantile(edge.gaps, 0.10):.3f}", f"{quantile(edge.gaps, 0.90):.3f}",
                ",".join(sorted(edge.path_ids)),
            )
            for key, edge in sorted(graph_edges.items())
        ),
    )
    self_loop_tokens = {
        left for left, right, _ in directed_edges if left == right
    }
    write_tsv(
        Path(f"{prefix}.topology.tsv"),
        ("bundle_id", "direction", "scc_id", "topology_rank", "scc_size", "in_cycle"),
        (
            (
                token[0], token[1], data[0], data[1], data[2],
                int(data[2] > 1 or token in self_loop_tokens),
            )
            for token, data in sorted(topology.items(), key=lambda item: (item[1][1], item[0]))
        ),
    )
    write_tsv(
        Path(f"{prefix}.topology_edges.tsv"),
        ("source_bundle", "source_direction", "target_bundle", "target_direction", "support_paths"),
        ((left[0], left[1], right[0], right[1], support) for left, right, support in directed_edges),
    )
    write_tsv(
        Path(f"{prefix}.subgraph_groups.tsv"),
        (
            "trim_group_id", "boundary_count", "boundaries", "interior_bundle_count",
            "interior_bundle_ids", "average_total_bp",
        ),
        (
            (
                group.group_id, len(group.boundaries),
                ";".join(boundary_text(boundary) for boundary in group.boundaries),
                len(group.interior_bundle_ids), ",".join(group.interior_bundle_ids) or ".",
                f"{group.average_total_bp:.3f}",
            )
            for group in subgraph_groups
        ),
    )
    write_tsv(
        Path(f"{prefix}.subgraph_shared_pairs.tsv"),
        (
            "boundary_a", "boundary_b", "shared_bundle_count", "shared_bundle_ids",
            "shared_average_total_bp", "qualifies_for_group_union",
        ),
        (
            (
                boundary_text(pair.boundary_a), boundary_text(pair.boundary_b),
                len(pair.shared_bundle_ids), ",".join(pair.shared_bundle_ids),
                f"{pair.shared_average_total_bp:.3f}", int(pair.merged),
            )
            for pair in shared_boundary_pairs
        ),
    )
    spanning_by_boundary = {
        boundary: len(set().union(*(subpath.support_units for subpath in classes.values())))
        for boundary, classes in subgraphs.items()
    }
    write_tsv(
        Path(f"{prefix}.subgraph_paths.tsv"),
        (
            "boundary", "trim_group_id", "canonical_interior", "support_units", "support_paths",
            "spanning_units", "frequency_within_boundary", "common", "action",
            "target_boundary", "target_interior", "trim_bp", "trim_bundle_count",
            "new_bundle_id", "supporting_paths", "effective_repeat_pattern",
            "duplication_signature", "target_class_duplication_signature",
            "post_trim_repeat_pattern", "post_trim_duplication_signature", "merged_as_dup",
        ),
        (
            (
                boundary_text(boundary), group_by_boundary[boundary], chain_text(subpath.interior),
                len(subpath.support_units),
                len(subpath.support_paths), spanning_by_boundary[boundary],
                f"{len(subpath.support_units) / spanning_by_boundary[boundary]:.6f}",
                int(subpath.common), subpath.action,
                boundary_text(subpath.target_boundary) if subpath.target_boundary else ".",
                chain_text(subpath.target_interior or ()),
                subpath.trim_bp, subpath.trim_bundle_count, subpath.new_bundle_id or ".",
                ",".join(sorted(subpath.support_paths)),
                repeat_pattern_text(subpath.repeat_pattern), chain_text(subpath.duplication_signature),
                chain_text(
                    target_subpath(subgraphs, subpath).duplication_signature
                    if subpath.target_key is not None
                    else ()
                ),
                repeat_pattern_text(subpath.merged_repeat_pattern),
                chain_text(subpath.merged_duplication_signature),
                int(bool(subpath.merged_duplication_signature)),
            )
            for boundary, classes in sorted(subgraphs.items())
            for subpath in sorted(classes.values(), key=lambda item: item.interior)
        ),
    )
    catalog_rows: List[Sequence[object]] = []
    for bundle in graph_bundles:
        catalog_rows.append(
            (
                bundle.new_id, "graph_bundle",
                ";".join(boundary_text(boundary) for boundary in bundle.boundaries),
                chain_text(bundle.canonical_interior), len(bundle.canonical_interior),
                bundle.support_units, bundle.support_paths, bundle.source_action,
                int(bundle.duplicated),
            )
        )
    for bundle_id, stat in sorted(stats.items()):
        if stat.action.startswith("retain"):
            catalog_rows.append(
                (
                    bundle_id, stat.action, ".", f"+{bundle_id}", 1, len(stat.units),
                    len(stat.paths), stat.action, int(stat.duplicated),
                )
            )
    write_tsv(
        Path(f"{prefix}.bundle_catalog.tsv"),
        (
            "new_bundle_id", "type", "boundary", "canonical_members", "member_count",
            "support_units", "support_paths", "source_action", "duplicated",
        ),
        catalog_rows,
    )
    bundle_map_rows: List[Sequence[object]] = [
        (
            token[0], bundle.new_id,
            ";".join(boundary_text(boundary) for boundary in bundle.boundaries), token[1],
            "duplicated_graph_bundle_member" if bundle.duplicated else "graph_bundle_member",
        )
        for bundle in graph_bundles
        for token in bundle.canonical_interior
    ]
    for boundary, classes in subgraphs.items():
        for subpath in classes.values():
            if subpath.action != "trimmed" or subpath.new_bundle_id is None:
                continue
            target = target_subpath(subgraphs, subpath)
            for token in deleted_tokens(subpath.interior, target.interior):
                bundle_map_rows.append(
                    (
                        token[0], subpath.new_bundle_id, boundary_text(boundary), token[1],
                        "absorbed_trimmed_branch",
                    )
                )

    for boundary, classes in subgraphs.items():
        for subpath in classes.values():
            if (
                subpath.new_bundle_id is None
                or not subpath.merged_duplication_signature
                or subpath.target_key is None
            ):
                continue
            target = target_subpath(subgraphs, subpath)
            for token in deleted_tokens(target.interior, subpath.merged_dup_tokens):
                bundle_map_rows.append(
                    (
                        token[0], subpath.new_bundle_id, boundary_text(boundary), token[1],
                        "nondup_interrupter_not_merged",
                    )
                )
    write_tsv(
        Path(f"{prefix}.bundle_map.tsv"),
        ("original_bundle_id", "new_bundle_id", "boundary", "canonical_direction", "role"),
        bundle_map_rows,
    )
    write_tsv(
        Path(f"{prefix}.trf_premerge_report.tsv"),
        (
            "path_id", "bundle_start", "bundle_end", "bundle_length", "bundle_id",
            "direction", "source", "action", "repeat_bp", "repeat_coverage",
            "repeat_union_interval_count", "repeat_union_intervals",
            "internal_repeat_count", "internal_repeat_bp", "max_internal_repeat_bp",
            "internal_repeat_intervals", "future_break_candidates", "raw_trf_count",
            "raw_trf_records",
        ),
        trf_premerge_report_rows,
    )
    trf_summary_header = (
        "bundle_id", "occurrence_count", "repeat_occurrence_count",
        "repeat_occurrence_fraction", "internal_repeat_occurrence_count",
        "full_repeat_occurrence_count", "dropped_occurrence_count",
        "remaining_occurrence_count",
        "median_repeat_bp", "max_repeat_bp", "median_repeat_coverage",
        "max_repeat_coverage", "action",
    )
    write_tsv(
        Path(f"{prefix}.trf_premerge_summary.tsv"),
        trf_summary_header,
        trf_premerge_summary_rows,
    )
    write_tsv(
        Path(f"{prefix}.trf_bundle_report.tsv"),
        (
            "path_id", "bundle_start", "bundle_end", "bundle_length", "bundle_id",
            "direction", "source", "action", "repeat_bp", "repeat_coverage",
            "repeat_union_interval_count", "repeat_union_intervals",
            "internal_repeat_count", "internal_repeat_bp", "max_internal_repeat_bp",
            "internal_repeat_intervals", "future_break_candidates", "raw_trf_count",
            "raw_trf_records",
        ),
        trf_report_rows,
    )
    write_tsv(
        Path(f"{prefix}.trf_bundle_summary.tsv"),
        trf_summary_header,
        trf_bundle_summary_rows,
    )
    with Path(f"{prefix}.summary.json").open("w", encoding="utf-8") as handle:
        json.dump(dict(summary), handle, indent=2, sort_keys=True)
        handle.write("\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prune and merge PGR-TK bundles using population-supported anchor-to-anchor subgraphs."
    )
    parser.add_argument("input_bed", type=Path)
    parser.add_argument("-o", "--output-prefix", required=True, type=Path)
    parser.add_argument("--sample-map", type=Path)
    parser.add_argument("--fasta", type=Path, help="population FASTA providing sequence lengths for coverage checks")
    parser.add_argument(
        "--coverage-denominator", choices=("auto", "fasta-length", "input-bed-union"), default="auto",
        help="auto: full FASTA sequence lengths when --fasta is provided, otherwise original input BED interval unions",
    )
    parser.add_argument("--anchor-file", type=Path, help="optional one-bundle-ID-per-line anchor list")
    parser.add_argument("--min-frequency", type=float, default=0.01)
    parser.add_argument("--min-count", type=int, default=2)
    parser.add_argument(
        "--min-dup-frequency", type=float, default=0.01,
        help="minimum population-unit frequency of R observations for a bundle to be treated as duplicated",
    )
    parser.add_argument(
        "--min-dup-count", type=int, default=2,
        help="minimum population-unit count of R observations for a bundle to be treated as duplicated",
    )
    parser.add_argument(
        "--long-dup-break-bp", type=int, default=10_000,
        help="make duplicated bundles with median length greater than this value hard merge boundaries",
    )
    parser.add_argument(
        "--long-dup-break-frequency", type=float, default=0.50,
        help="minimum population-unit R frequency for a long duplicated bundle to break merging",
    )
    parser.add_argument("--short-bp", type=int, default=10_000, help="median-length threshold for short-bundle deletion candidates and graph trimming")
    parser.add_argument("--anchor-min-frequency", type=float, default=0.95)
    parser.add_argument(
        "--anchor-min-bp", type=int, default=None,
        help="default: Q25 of per-bundle population median lengths; explicit value overrides the automatic threshold",
    )
    parser.add_argument("--anchor-min-single-copy", type=float, default=0.95)
    parser.add_argument("--min-path-frequency", type=float, default=0.02)
    parser.add_argument("--min-path-count", type=int, default=2)
    parser.add_argument("--max-trim-bp", type=int, default=10_000)
    parser.add_argument("--max-trim-bundles", type=int, default=5)
    parser.add_argument(
        "--subgraph-merge-min-shared-bp", type=float, default=10_000.0,
        help="jointly trim boundary subgraphs when shared bundle mean lengths sum to at least this many bp",
    )
    parser.add_argument(
        "--trim-rare-long", action="store_true",
        help="allow rare bundles >= --short-bp to be removed from a subpath",
    )
    parser.add_argument(
        "--merge-common-short-paths", action="store_true",
        help="also merge common alternative paths when they differ only by trimmable short nodes",
    )
    parser.add_argument("--keep-unmerged-short", action="store_true", help="compatibility option; graph merging now always retains unmerged non-rare short bundles")
    parser.add_argument(
        "--short-min-sequence-coverage", type=float, default=0.90,
        help="report-only sequence coverage threshold; graph merging no longer deletes short bundles",
    )
    parser.add_argument(
        "--short-max-low-coverage-fraction", type=float, default=0.05,
        help="compatibility option; bundle deletion has moved to final refinement",
    )
    parser.add_argument(
        "--coordinate-overlap-seed", type=int, default=17,
        help="seed for reproducible random split points used to remove final BED overlaps",
    )
    parser.add_argument(
        "--keep-coordinate-overlaps", action="store_true",
        help="leave adjacent final bundle coordinates unchanged even when they overlap",
    )
    parser.add_argument("--trf-bed", type=Path, help="optional 9-column TRF BED annotation")
    parser.add_argument(
        "--trf-min-length", type=int, default=0,
        help="minimum total tandem-repeat interval length in bp (default: 0; no length filter)",
    )
    parser.add_argument(
        "--trf-min-identity", type=float, default=80.0,
        help="minimum TRF percent identity",
    )
    parser.add_argument(
        "--trf-full-coverage", type=float, default=0.90,
        help=("minimum fraction covered by one merged TRF interval for an occurrence "
              "to count toward whole-ID removal"),
    )
    parser.add_argument(
        "--trf-bundle-removal-frequency", type=float, default=0.50,
        help=("remove every occurrence of a bundle ID when this fraction are covered "
              "to --trf-full-coverage by one merged TRF interval"),
    )
    return parser


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    for name in (
        "min_frequency", "min_dup_frequency", "anchor_min_frequency",
        "long_dup_break_frequency", "anchor_min_single_copy", "min_path_frequency",
        "trf_bundle_removal_frequency",
        "short_min_sequence_coverage", "short_max_low_coverage_fraction",
    ):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            parser.error(f"--{name.replace('_', '-')} must be in [0, 1]")
    for name in (
        "min_count", "min_dup_count", "long_dup_break_bp", "short_bp", "min_path_count",
        "max_trim_bp", "max_trim_bundles", "trf_min_length",
    ):
        if getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} must be non-negative")
    if args.anchor_min_bp is not None and args.anchor_min_bp < 0:
        parser.error("--anchor-min-bp must be non-negative")
    if args.subgraph_merge_min_shared_bp < 0:
        parser.error("--subgraph-merge-min-shared-bp must be non-negative")
    if not 0 <= args.trf_min_identity <= 100:
        parser.error("--trf-min-identity must be in [0, 100]")
    if not 0 <= args.trf_full_coverage <= 1:
        parser.error("--trf-full-coverage must be in [0, 1]")


def run(args: argparse.Namespace) -> Mapping[str, object]:
    input_by_path = read_pgr_bed(args.input_bed)
    coverage_intervals, coverage_denominator = build_coverage_intervals(
        input_by_path, getattr(args, "fasta", None), getattr(args, "coverage_denominator", "auto"),
    )
    unit_by_path = read_sample_map(args.sample_map, input_by_path)
    by_path = input_by_path
    trf_intervals: Dict[str, Tuple[TRFInterval, ...]] = {}
    trf_records: Dict[str, Tuple[TRFRecord, ...]] = {}
    trf_premerge_report_rows: List[Sequence[object]] = []
    trf_premerge_summary_rows: List[Sequence[object]] = []
    trf_summary: Dict[str, int] = {
        "trf_input_row_count": 0,
        "trf_selected_row_count": 0,
        "trf_union_interval_count": 0,
        "trf_matched_path_count": 0,
        "trf_reported_occurrence_count": 0,
        "trf_repeat_overlapping_occurrence_count": 0,
        "trf_internal_repeat_occurrence_count": 0,
        "trf_full_repeat_occurrence_count": 0,
        "trf_dropped_occurrence_count": 0,
        "trf_dropped_bundle_type_count": 0,
        "trf_postmerge_matched_path_count": 0,
        "trf_postmerge_reported_occurrence_count": 0,
        "trf_postmerge_repeat_overlapping_occurrence_count": 0,
        "trf_postmerge_internal_repeat_occurrence_count": 0,
        "trf_postmerge_full_repeat_occurrence_count": 0,
    }
    if args.trf_bed is not None:
        trf_intervals, trf_records, trf_read_summary = read_trf_bed(
            args.trf_bed,
            args.trf_min_length,
            args.trf_min_identity,
        )
        if not set(input_by_path) & set(trf_intervals):
            raise ValueError("TRF BED has no selected intervals matching an input path ID")
        (
            by_path,
            trf_premerge_report_rows,
            trf_premerge_summary_rows,
            trf_premerge_action_summary,
        ) = filter_input_pgr_by_trf(
            input_by_path,
            trf_intervals,
            trf_records,
            args.trf_full_coverage,
            args.trf_bundle_removal_frequency,
        )
        trf_summary.update(trf_read_summary)
        trf_summary.update(trf_premerge_action_summary)

    stats, total_units, required_units, required_dup_units = compute_bundle_stats(
        by_path,
        unit_by_path,
        args.min_count,
        args.min_frequency,
        args.min_dup_count,
        args.min_dup_frequency,
        args.short_bp,
    )
    anchor_length_q25 = quantile([stat.median_length for stat in stats.values()], 0.25)
    anchor_min_bp = args.anchor_min_bp if args.anchor_min_bp is not None else anchor_length_q25
    anchors = choose_anchors(
        stats,
        read_anchor_file(args.anchor_file),
        args.anchor_min_frequency,
        anchor_min_bp,
        args.anchor_min_single_copy,
    )
    graph_edges = build_population_graph(by_path, unit_by_path)
    topology, directed_edges = compute_topology(by_path)

    # Pass 1 determines merge candidates without allowing long-duplication
    # breakpoints to alter the initial anchor graph.
    candidate_subgraphs = collect_anchor_subpaths(by_path, unit_by_path, anchors, stats)
    (
        candidate_group_by_boundary,
        candidate_subgraph_groups,
        candidate_shared_boundary_pairs,
    ) = (
        group_subgraphs_by_shared_bundles(
            candidate_subgraphs,
            stats,
            args.subgraph_merge_min_shared_bp,
            topology,
        )
    )
    candidate_graph_bundles = classify_subpaths(
        candidate_subgraphs,
        stats,
        args.min_path_count,
        args.min_path_frequency,
        args.max_trim_bp,
        args.max_trim_bundles,
        args.trim_rare_long,
        not args.merge_common_short_paths,
        topology,
        candidate_group_by_boundary,
    )
    long_dup_breakpoints = find_long_dup_breakpoints_in_candidates(
        stats,
        candidate_subgraphs,
        args.long_dup_break_bp,
        args.long_dup_break_frequency,
    )

    # Pass 2 is needed only when a first-pass candidate contains a qualifying
    # long duplicated ID. Promoting that ID to an anchor splits every occurrence
    # consistently and prevents the final graph bundle from spanning it.
    if long_dup_breakpoints:
        promote_long_dup_breakpoints(stats, anchors, long_dup_breakpoints)
        subgraphs = collect_anchor_subpaths(by_path, unit_by_path, anchors, stats)
        group_by_boundary, subgraph_groups, shared_boundary_pairs = (
            group_subgraphs_by_shared_bundles(
                subgraphs,
                stats,
                args.subgraph_merge_min_shared_bp,
                topology,
            )
        )
        graph_bundles = classify_subpaths(
            subgraphs,
            stats,
            args.min_path_count,
            args.min_path_frequency,
            args.max_trim_bp,
            args.max_trim_bundles,
            args.trim_rare_long,
            not args.merge_common_short_paths,
            topology,
            group_by_boundary,
        )
    else:
        subgraphs = candidate_subgraphs
        group_by_boundary = candidate_group_by_boundary
        subgraph_groups = candidate_subgraph_groups
        shared_boundary_pairs = candidate_shared_boundary_pairs
        graph_bundles = candidate_graph_bundles
    assign_bundle_actions(stats, graph_bundles, subgraphs, args.keep_unmerged_short)
    outputs = project_graph_bundles(by_path, anchors, subgraphs, stats, args.keep_unmerged_short)
    coordinate_adjustments: List[CoordinateAdjustment] = []
    if not args.keep_coordinate_overlaps:
        outputs, coordinate_adjustments = resolve_output_coordinate_overlaps(
            outputs,
            args.coordinate_overlap_seed,
        )
    short_min_coverage = getattr(args, "short_min_sequence_coverage", 0.90)
    short_max_low_fraction = getattr(args, "short_max_low_coverage_fraction", 0.05)
    outputs, short_deletion_rows, sequence_coverage_rows, short_coverage_summary = prune_unmerged_short_by_coverage(
        outputs, stats, coverage_intervals, short_min_coverage, short_max_low_fraction,
        keep_unmerged_short=True,
    )
    trf_report_rows: List[Sequence[object]] = []
    trf_bundle_summary_rows: List[Sequence[object]] = []
    if args.trf_bed is not None:
        outputs, trf_report_rows, trf_bundle_summary_rows, postmerge_action_summary = (
            apply_trf_filter_and_report(
                outputs,
                trf_intervals,
                trf_records,
                args.trf_full_coverage,
                remove_full_repeat=False,
            )
        )
        trf_summary.update(
            {
                key.replace("trf_", "trf_postmerge_", 1): value
                for key, value in postmerge_action_summary.items()
                if key not in {"trf_dropped_occurrence_count", "trf_dropped_bundle_type_count"}
            }
        )
    summary: Dict[str, object] = {
        "input": str(args.input_bed),
        "path_count": len(input_by_path),
        "population_unit_count": total_units,
        "rare_required_unit_count": required_units,
        "dup_required_unit_count": required_dup_units,
        "duplicated_bundle_type_count": sum(stat.duplicated for stat in stats.values()),
        "long_dup_breakpoint_bundle_count": len(long_dup_breakpoints),
        "prebreak_candidate_anchor_boundary_count": len(candidate_subgraphs),
        "prebreak_candidate_subpath_class_count": sum(
            len(classes) for classes in candidate_subgraphs.values()
        ),
        "prebreak_candidate_graph_bundle_count": len(candidate_graph_bundles),
        "input_bundle_type_count": len(
            {record.bundle_id for records in input_by_path.values() for record in records}
        ),
        "input_occurrence_count": sum(len(records) for records in input_by_path.values()),
        "premerge_bundle_type_count": len(stats),
        "premerge_occurrence_count": sum(len(records) for records in by_path.values()),
        "graph_edge_count": len(graph_edges),
        "topology_scc_count": len({value[0] for value in topology.values()}),
        "cyclic_signed_node_count": sum(
            value[2] > 1
            or any(left == token and right == token for left, right, _ in directed_edges)
            for token, value in topology.items()
        ),
        "anchor_bundle_count": len(anchors),
        "anchor_boundary_count": len(subgraphs),
        "subgraph_trim_group_count": len(subgraph_groups),
        "joint_subgraph_trim_group_count": sum(
            len(group.boundaries) > 1 for group in subgraph_groups
        ),
        "qualifying_shared_boundary_pair_count": sum(
            pair.merged for pair in shared_boundary_pairs
        ),
        "cross_subgraph_trimmed_class_count": sum(
            subpath.action == "trimmed"
            and subpath.target_boundary is not None
            and subpath.target_boundary != subpath.boundary
            for classes in subgraphs.values()
            for subpath in classes.values()
        ),
        "observed_subpath_class_count": sum(len(classes) for classes in subgraphs.values()),
        "trimmed_subpath_class_count": sum(
            subpath.action == "trimmed" for classes in subgraphs.values() for subpath in classes.values()
        ),
        "retained_subpath_class_count": sum(
            subpath.action.startswith("retain") for classes in subgraphs.values() for subpath in classes.values()
        ),
        "new_graph_bundle_count": len(graph_bundles),
        "duplicated_graph_bundle_count": sum(bundle.duplicated for bundle in graph_bundles),
        "coordinate_overlap_adjustment_count": len(coordinate_adjustments),
        "coordinate_overlap_adjusted_bp": sum(
            adjustment.overlap_bp for adjustment in coordinate_adjustments
        ),
        "coordinate_overlap_max_bp": max(
            (adjustment.overlap_bp for adjustment in coordinate_adjustments),
            default=0,
        ),
        "output_occurrence_count": sum(len(items) for items in outputs.values()),
        "output_bundle_type_count": len(
            {item.new_bundle_id for items in outputs.values() for item in items}
        ),
        "anchor_length_q25_bp": anchor_length_q25,
        **short_coverage_summary,
        **trf_summary,
        "parameters": {
            "min_frequency": args.min_frequency,
            "min_count": args.min_count,
            "min_dup_frequency": args.min_dup_frequency,
            "min_dup_count": args.min_dup_count,
            "long_dup_break_bp": args.long_dup_break_bp,
            "long_dup_break_frequency": args.long_dup_break_frequency,
            "short_bp": args.short_bp,
            "anchor_min_frequency": args.anchor_min_frequency,
            "anchor_min_bp": anchor_min_bp,
            "anchor_min_bp_mode": "explicit" if args.anchor_min_bp is not None else "bundle_median_q25",
            "anchor_min_single_copy": args.anchor_min_single_copy,
            "min_path_frequency": args.min_path_frequency,
            "min_path_count": args.min_path_count,
            "max_trim_bp": args.max_trim_bp,
            "max_trim_bundles": args.max_trim_bundles,
            "subgraph_merge_min_shared_bp": args.subgraph_merge_min_shared_bp,
            "trim_rare_long": args.trim_rare_long,
            "preserve_common_paths": not args.merge_common_short_paths,
            "keep_unmerged_short": True,
            "keep_unmerged_short_requested": args.keep_unmerged_short,
            "short_bundle_pruning_enabled": False,
            "coverage_denominator_requested": getattr(args, "coverage_denominator", "auto"),
            "coverage_denominator": coverage_denominator,
            "coverage_fasta": str(args.fasta) if getattr(args, "fasta", None) is not None else None,
            "short_min_sequence_coverage": short_min_coverage,
            "short_max_low_coverage_fraction": short_max_low_fraction,
            "resolve_coordinate_overlaps": not args.keep_coordinate_overlaps,
            "coordinate_overlap_seed": args.coordinate_overlap_seed,
            "trf_bed": str(args.trf_bed) if args.trf_bed is not None else None,
            "trf_min_length": args.trf_min_length,
            "trf_min_identity": args.trf_min_identity,
            "trf_full_coverage": args.trf_full_coverage,
            "trf_full_coverage_semantics": "maximum overlap with one merged TRF interval",
            "trf_bundle_removal_frequency": args.trf_bundle_removal_frequency,
        },
    }
    write_outputs(
        args.output_prefix,
        outputs,
        stats,
        graph_edges,
        topology,
        directed_edges,
        subgraphs,
        graph_bundles,
        group_by_boundary,
        subgraph_groups,
        shared_boundary_pairs,
        long_dup_breakpoints,
        coordinate_adjustments,
        trf_premerge_report_rows,
        trf_premerge_summary_rows,
        trf_report_rows,
        trf_bundle_summary_rows,
        summary,
        short_deletion_rows,
        sequence_coverage_rows,
    )
    return summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)
    try:
        summary = run(args)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
