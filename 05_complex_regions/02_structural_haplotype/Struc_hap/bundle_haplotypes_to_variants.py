#!/usr/bin/env python3
"""Decompose signed PGR-TK bundle haplotypes into anchored variant regions.

The input is a native four-column PGR-TK-like BED.  Each BED path is treated as
one haplotype sequence; when a FAI is supplied, only paths found in both files
are analyzed.  The workflow is:

1. Check that the signed left and right terminal bundle signatures are
   conserved.  Paths differing from either dominant terminal signature are
   flank outliers.  If their union is greater than the configured fraction
   (1% by default), return the original PGR-TK BED unchanged for the whole
   region as both ``original_pgrtk.bed`` and
   ``subvariation_regions/R0001.pgr.bed``, then stop decomposition.
2. On the remaining core paths, retain bundle IDs occurring exactly once on
   every path and having one common orientation.  When fewer than two
   order-consistent anchors remain, return the original BED unchanged and
   retain the anchor-selection audit.
3. Find a maximum-cardinality subset whose relative order is identical on all
   core paths.  Because every candidate is single-copy, this is a longest chain
   in the all-path precedence partial order.
4. Collapse consecutive anchors and their intervening bundles when the full
   signed order is identical on every core path.  Variable intervals between
   these composite anchors become the regions.  Empty interiors are retained.
   If no variable interval remains, write the complete original BED as
   ``subvariation_regions/R0001.pgr.bed`` without counting a variant region.
5. Combine alleles from variable regions into population structural-haplotype
   types and report a coarse event class relative to each region's most common
   allele.

Coordinates are zero-based, half-open.  Directions are reported as ``+/-``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import re
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:
    from .merge_pgr_bundles import (
        Occurrence,
        flip_direction,
        read_pgr_bed,
        write_tsv,
    )
except ImportError:
    from merge_pgr_bundles import (
        Occurrence,
        flip_direction,
        read_pgr_bed,
        write_tsv,
    )


SignedToken = Tuple[str, str]
Signature = Tuple[SignedToken, ...]


OUTPUT_FILES = (
    "bundle_variant_summary.json",
    "flank_check.tsv",
    "flank_summary.json",
    "skipped_paths.tsv",
    "anchor_candidates.tsv",
    "anchors.tsv",
    "merged_anchors.tsv",
    "merged_anchors.bed",
    "variant_regions.tsv",
    "subvariation_regions.ref.bed",
    "variant_alleles.tsv",
    "path_variant_regions.tsv",
    "structural_variant_types.tsv",
    "path_variant_types.tsv",
    "original_pgrtk.bed",
)


@dataclass(frozen=True)
class MergedAnchor:
    anchor_id: str
    first_bundle: str
    last_bundle: str
    signature: Signature
    reference_start: int
    reference_end: int


def natural_key(value: str):
    return (0, int(value), value) if value.isdigit() else (1, value)


def token_key(token: SignedToken):
    return natural_key(token[0]), token[1]


def signature_key(signature: Signature):
    return tuple(token_key(token) for token in signature)


def signed_text(signature: Sequence[SignedToken]) -> str:
    return ",".join(f"{direction}{bundle_id}" for bundle_id, direction in signature) or "."


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fai_paths(path: Path) -> List[str]:
    paths = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2 or not fields[0]:
                raise ValueError(f"{path}:{line_number}: expected FASTA-index name and length")
            if fields[0] in seen:
                raise ValueError(f"{path}:{line_number}: duplicate sequence {fields[0]!r}")
            try:
                length = int(fields[1])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid sequence length") from error
            if length <= 0:
                raise ValueError(f"{path}:{line_number}: sequence length must be positive")
            seen.add(fields[0])
            paths.append(fields[0])
    return paths


def read_source_rows(path: Path) -> Dict[int, str]:
    """Keep the input's first four fields verbatim for regional PGR BEDs."""
    rows = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 4:
                fields = stripped.split()
            rows[line_number] = "\t".join(fields[:4])
    return rows


def occurrence_signature(records: Sequence[Occurrence]) -> Signature:
    return tuple((record.bundle_id, record.direction) for record in records)


def dominant_signature(values: Mapping[str, Signature]) -> Tuple[Signature, int]:
    counts = Counter(values.values())
    if not counts:
        return (), 0
    signature, count = min(
        counts.items(), key=lambda item: (-item[1], signature_key(item[0]))
    )
    return signature, count


def terminal_signature(records: Sequence[Occurrence], count: int, side: str) -> Signature:
    if len(records) < count:
        return occurrence_signature(records)
    selected = records[:count] if side == "left" else records[-count:]
    return occurrence_signature(selected)


def write_flank_audit(
    output_dir: Path,
    by_path: Mapping[str, Sequence[Occurrence]],
    path_ids: Sequence[str],
    flank_bundles: int,
    max_outlier_fraction: float,
):
    left = {path: terminal_signature(by_path[path], flank_bundles, "left") for path in path_ids}
    right = {path: terminal_signature(by_path[path], flank_bundles, "right") for path in path_ids}
    dominant_left, left_count = dominant_signature(left)
    dominant_right, right_count = dominant_signature(right)
    outliers = {
        path for path in path_ids
        if left[path] != dominant_left or right[path] != dominant_right
    }
    outlier_fraction = len(outliers) / len(path_ids) if path_ids else 0.0
    rows = []
    for path in path_ids:
        left_match = left[path] == dominant_left
        right_match = right[path] == dominant_right
        rows.append((
            path,
            len(by_path[path]),
            signed_text(left[path]),
            signed_text(right[path]),
            int(left_match),
            int(right_match),
            int(path in outliers),
            "core" if path not in outliers else "flank_outlier",
        ))
    write_tsv(
        output_dir / "flank_check.tsv",
        (
            "path_id", "bundle_count", "left_signature", "right_signature",
            "left_matches_dominant", "right_matches_dominant", "is_outlier", "status",
        ),
        rows,
    )
    summary = {
        "flank_bundle_count": flank_bundles,
        "dominant_left_signature": signed_text(dominant_left),
        "dominant_right_signature": signed_text(dominant_right),
        "left_dominant_path_count": left_count,
        "right_dominant_path_count": right_count,
        "total_path_count": len(path_ids),
        "outlier_path_count": len(outliers),
        "outlier_fraction": outlier_fraction,
        "max_outlier_fraction": max_outlier_fraction,
        "passes": bool(path_ids) and outlier_fraction <= max_outlier_fraction,
        "outlier_paths": sorted(outliers),
    }
    (output_dir / "flank_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return [path for path in path_ids if path not in outliers], outliers, summary


def is_precedent(
    left: str,
    right: str,
    positions: Mapping[str, Mapping[str, int]],
    paths: Sequence[str],
) -> bool:
    return all(positions[path][left] < positions[path][right] for path in paths)


def maximum_consistent_anchor_chain(
    candidates: Sequence[str],
    positions: Mapping[str, Mapping[str, int]],
    paths: Sequence[str],
) -> Tuple[str, ...]:
    """Return a deterministic maximum-cardinality all-path order chain."""
    if not candidates:
        return ()
    best_ending: List[Tuple[int, ...]] = []
    for current in range(len(candidates)):
        best = (current,)
        for previous in range(current):
            if not is_precedent(candidates[previous], candidates[current], positions, paths):
                continue
            proposal = best_ending[previous] + (current,)
            if len(proposal) > len(best) or (len(proposal) == len(best) and proposal < best):
                best = proposal
        best_ending.append(best)
    selected = min(best_ending, key=lambda chain: (-len(chain), chain))
    return tuple(candidates[index] for index in selected)


def select_anchors(
    output_dir: Path,
    by_path: Mapping[str, Sequence[Occurrence]],
    core_paths: Sequence[str],
    reference_path: Optional[str],
):
    if not core_paths:
        raise ValueError("no core paths remain after terminal-flank filtering")
    if reference_path is None:
        reference_path = random.SystemRandom().choice(sorted(core_paths))
    if reference_path not in core_paths:
        raise ValueError("--ref-id/--reference-path must be one of the conserved core paths")

    counts = {
        path: Counter(record.bundle_id for record in by_path[path])
        for path in core_paths
    }
    all_ids = sorted(
        {bundle_id for path in core_paths for bundle_id in counts[path]},
        key=natural_key,
    )
    directions = {
        path: defaultdict(set)
        for path in core_paths
    }
    for path in core_paths:
        for record in by_path[path]:
            directions[path][record.bundle_id].add(record.direction)

    positions: Dict[str, Dict[str, int]] = {path: {} for path in core_paths}
    orientations: Dict[str, str] = {}
    ordered_candidates = []
    for index, record in enumerate(by_path[reference_path]):
        bundle_id = record.bundle_id
        if bundle_id in ordered_candidates:
            continue
        if not all(counts[path][bundle_id] == 1 for path in core_paths):
            continue
        observed = {next(iter(directions[path][bundle_id])) for path in core_paths}
        if len(observed) != 1:
            continue
        ordered_candidates.append(bundle_id)
        orientations[bundle_id] = next(iter(observed))
    for path in core_paths:
        for index, record in enumerate(by_path[path]):
            if record.bundle_id in orientations:
                positions[path][record.bundle_id] = index

    selected = maximum_consistent_anchor_chain(ordered_candidates, positions, core_paths)
    selected_set = set(selected)
    candidate_rows = []
    for bundle_id in all_ids:
        presence_count = sum(counts[path][bundle_id] > 0 for path in core_paths)
        single_count = sum(counts[path][bundle_id] == 1 for path in core_paths)
        orientation_values = sorted(
            {direction for path in core_paths for direction in directions[path].get(bundle_id, ())}
        )
        if presence_count < len(core_paths):
            reason = "missing_from_core_path"
        elif single_count < len(core_paths):
            reason = "not_single_copy_on_every_core_path"
        elif len(orientation_values) != 1:
            reason = "orientation_inconsistent"
        elif bundle_id not in selected_set:
            reason = "excluded_by_maximum_order_consistent_chain"
        else:
            reason = "selected_anchor"
        candidate_rows.append((
            bundle_id,
            presence_count,
            single_count,
            ",".join(orientation_values) or ".",
            by_path[reference_path].index(next(
                record for record in by_path[reference_path] if record.bundle_id == bundle_id
            )) + 1 if counts[reference_path][bundle_id] else ".",
            int(bundle_id in selected_set),
            reason,
        ))
    write_tsv(
        output_dir / "anchor_candidates.tsv",
        (
            "bundle_id", "core_path_presence_count", "core_path_single_copy_count",
            "observed_directions", "reference_position_1based", "selected", "reason",
        ),
        candidate_rows,
    )
    write_tsv(
        output_dir / "anchors.tsv",
        (
            "anchor_rank", "bundle_id", "direction", "reference_path",
            "reference_position_1based", "core_path_count",
        ),
        (
            (
                rank,
                bundle_id,
                orientations[bundle_id],
                reference_path,
                positions[reference_path][bundle_id] + 1,
                len(core_paths),
            )
            for rank, bundle_id in enumerate(selected, 1)
        ),
    )
    return selected, orientations, positions, reference_path


def merge_invariant_anchors(
    output_dir: Path,
    by_path: Mapping[str, Sequence[Occurrence]],
    core_paths: Sequence[str],
    anchors: Sequence[str],
    positions: Mapping[str, Mapping[str, int]],
    reference_path: str,
) -> Tuple[MergedAnchor, ...]:
    """Collapse maximal runs of anchors separated only by invariant walks."""
    groups = []
    first = 0
    for index, (left, right) in enumerate(zip(anchors, anchors[1:])):
        observed = {
            occurrence_signature(by_path[path][positions[path][left] + 1:positions[path][right]])
            for path in core_paths
        }
        if len(observed) > 1:
            groups.append((first, index))
            first = index + 1
    groups.append((first, len(anchors) - 1))

    merged = []
    rows = []
    bed_rows = []
    for rank, (first, last) in enumerate(groups, 1):
        left, right = anchors[first], anchors[last]
        ref_records = by_path[reference_path]
        start_index, end_index = positions[reference_path][left], positions[reference_path][right]
        signature = occurrence_signature(ref_records[start_index:end_index + 1])
        # Every internal gap in a group was checked for exact signed equality.
        for path in core_paths:
            records = by_path[path]
            observed = occurrence_signature(
                records[positions[path][left]:positions[path][right] + 1]
            )
            if observed != signature:
                raise ValueError(f"merged anchor order differs on {path}: {left}..{right}")
        anchor_id = f"ANCHOR{rank:04d}"
        reference_span = ref_records[start_index:end_index + 1]
        ref_start = min(record.start for record in reference_span)
        ref_end = max(record.end for record in reference_span)
        merged.append(MergedAnchor(anchor_id, left, right, signature, ref_start, ref_end))
        rows.append((
            rank, anchor_id, left, right, len(signature), signed_text(signature),
            reference_path, ref_start, ref_end, len(core_paths),
        ))
        bed_rows.append((reference_path, ref_start, ref_end, anchor_id))
    write_tsv(
        output_dir / "merged_anchors.tsv",
        (
            "anchor_rank", "anchor_id", "first_boundary_bundle", "last_boundary_bundle",
            "bundle_count", "signed_bundle_order", "reference_path", "reference_start",
            "reference_end", "core_path_count",
        ),
        rows,
    )
    with (output_dir / "merged_anchors.bed").open("w", encoding="utf-8") as handle:
        for row in bed_rows:
            handle.write("\t".join(map(str, row)) + "\n")
    return tuple(merged)


def is_subsequence(short: Signature, long: Signature) -> bool:
    iterator = iter(long)
    return all(any(item == candidate for item in iterator) for candidate in short)


def reverse_complement_signature(signature: Signature) -> Signature:
    return tuple((bundle_id, flip_direction(direction)) for bundle_id, direction in reversed(signature))


def has_single_inversion(reference: Signature, alternate: Signature) -> bool:
    if len(reference) != len(alternate) or reference == alternate:
        return False
    for left in range(len(reference)):
        for right in range(left + 1, len(reference) + 1):
            proposal = (
                reference[:left]
                + reverse_complement_signature(reference[left:right])
                + reference[right:]
            )
            if proposal == alternate:
                return True
    return False


def classify_allele(reference: Signature, alternate: Signature) -> str:
    """Return a conservative coarse class; exact signed alleles remain authoritative."""
    if alternate == reference:
        return "REF"
    if has_single_inversion(reference, alternate):
        return "INV"
    if len(alternate) < len(reference) and is_subsequence(alternate, reference):
        return "DEL"
    if len(alternate) > len(reference) and is_subsequence(reference, alternate):
        reference_ids = Counter(bundle_id for bundle_id, _ in reference)
        alternate_ids = Counter(bundle_id for bundle_id, _ in alternate)
        added_ids = {
            bundle_id for bundle_id, count in alternate_ids.items()
            if count > reference_ids.get(bundle_id, 0)
        }
        return "DUP" if added_ids and added_ids <= set(reference_ids) else "INS"
    reference_ids = Counter(bundle_id for bundle_id, _ in reference)
    alternate_ids = Counter(bundle_id for bundle_id, _ in alternate)
    if reference_ids == alternate_ids:
        return "REARRANGEMENT"
    if set(reference_ids) == set(alternate_ids):
        return "CNV"
    return "COMPLEX"


@dataclass(frozen=True)
class RegionalAllele:
    allele_id: str
    signature: Signature
    count: int
    event_class: str


def decompose_regions(
    output_dir: Path,
    by_path: Mapping[str, Sequence[Occurrence]],
    source_rows: Mapping[int, str],
    core_paths: Sequence[str],
    outlier_paths: Iterable[str],
    anchors: Sequence[MergedAnchor],
    reference_path: str,
):
    outlier_paths = set(outlier_paths)
    position_by_path = {
        path: {record.bundle_id: index for index, record in enumerate(by_path[path])}
        for path in core_paths
    }
    region_summaries = []
    allele_rows = []
    path_region_rows = []
    path_vectors = {path: [] for path in core_paths}
    variable_region_ids = []
    region_bed_files = {}
    reference_bed_rows = []
    for region_index, (left_group, right_group) in enumerate(zip(anchors, anchors[1:]), 1):
        region_id = f"R{region_index:04d}"
        left_anchor = left_group.last_bundle
        right_anchor = right_group.first_bundle
        reference_records = by_path[reference_path]
        reference_start = reference_records[
            position_by_path[reference_path][left_anchor]
        ].end
        reference_end = reference_records[
            position_by_path[reference_path][right_anchor]
        ].start
        reference_interval_status = (
            "gap" if reference_start < reference_end else
            "touching" if reference_start == reference_end else "overlap"
        )
        if reference_interval_status == "gap":
            reference_bed_rows.append((reference_path, reference_start, reference_end, region_id))
        signatures = {}
        coordinates = {}
        for path in core_paths:
            left_index = position_by_path[path][left_anchor]
            right_index = position_by_path[path][right_anchor]
            if left_index >= right_index:
                raise ValueError(
                    f"selected anchor order is inconsistent on {path}: {left_anchor}, {right_anchor}"
                )
            internal = tuple(
                (record.bundle_id, record.direction)
                for record in by_path[path][left_index + 1:right_index]
            )
            left_record = by_path[path][left_index]
            right_record = by_path[path][right_index]
            signatures[path] = internal
            coordinates[path] = (
                left_record.end,
                right_record.start,
                max(0, right_record.start - left_record.end),
                max(0, left_record.end - right_record.start),
            )
        counts = Counter(signatures.values())
        ordered = sorted(counts, key=lambda signature: (-counts[signature], signature_key(signature)))
        reference = ordered[0]
        alleles = []
        for allele_index, signature in enumerate(ordered, 1):
            allele = RegionalAllele(
                f"{region_id}.A{allele_index:03d}",
                signature,
                counts[signature],
                classify_allele(reference, signature),
            )
            alleles.append(allele)
            allele_rows.append((
                region_id,
                allele.allele_id,
                allele.count,
                f"{allele.count / len(core_paths):.9f}",
                int(signature == reference),
                allele.event_class,
                len(signature),
                signed_text(signature),
            ))
        allele_by_signature = {allele.signature: allele for allele in alleles}
        is_variable = len(alleles) > 1
        if not is_variable:
            raise ValueError(f"{region_id}: interval between merged anchors is unexpectedly invariant")
        variable_region_ids.append(region_id)
        region_summaries.append((
            region_id,
            left_group.anchor_id,
            right_group.anchor_id,
            left_anchor,
            right_anchor,
            reference_path,
            reference_start,
            reference_end,
            reference_interval_status,
            len(core_paths),
            len(alleles),
            int(is_variable),
            alleles[0].allele_id,
            signed_text(reference),
        ))
        for path in core_paths:
            allele = allele_by_signature[signatures[path]]
            start, end, gap_bp, overlap_bp = coordinates[path]
            path_region_rows.append((
                path,
                region_id,
                left_group.anchor_id,
                right_group.anchor_id,
                left_anchor,
                right_anchor,
                start,
                end,
                gap_bp,
                overlap_bp,
                allele.allele_id,
                allele.event_class,
                signed_text(signatures[path]),
            ))
            path_vectors[path].append(allele.allele_id)

        # Keep native four-column PGR-TK records, including both composite
        # flanks.  An empty interior allele remains visible through its flanks.
        region_bed = output_dir / "subvariation_regions" / f"{region_id}.pgr.bed"
        with region_bed.open("w", encoding="utf-8") as handle:
            for path in sorted(core_paths):
                records = by_path[path]
                start_index = position_by_path[path][left_group.first_bundle]
                end_index = position_by_path[path][right_group.last_bundle]
                for record in records[start_index:end_index + 1]:
                    handle.write(source_rows[record.line_number] + "\n")
        region_bed_files[region_id] = str(region_bed.resolve())

    write_tsv(
        output_dir / "variant_regions.tsv",
        (
            "region_id", "left_merged_anchor", "right_merged_anchor",
            "left_anchor", "right_anchor", "reference_path", "reference_start",
            "reference_end", "reference_interval_status", "core_path_count",
            "allele_count", "is_variable", "reference_allele", "reference_internal_order",
        ),
        region_summaries,
    )
    with (output_dir / "subvariation_regions.ref.bed").open("w", encoding="utf-8") as handle:
        for row in reference_bed_rows:
            handle.write("\t".join(map(str, row)) + "\n")
    write_tsv(
        output_dir / "variant_alleles.tsv",
        (
            "region_id", "allele_id", "path_count", "path_frequency", "is_reference",
            "event_class", "internal_bundle_count", "signed_internal_bundle_order",
        ),
        allele_rows,
    )
    write_tsv(
        output_dir / "path_variant_regions.tsv",
        (
            "path_id", "region_id", "left_merged_anchor", "right_merged_anchor",
            "left_anchor", "right_anchor",
            "inter_anchor_start", "inter_anchor_end", "inter_anchor_gap_bp",
            "anchor_overlap_bp", "allele_id", "event_class", "signed_internal_bundle_order",
        ),
        path_region_rows,
    )

    vectors = {path: tuple(path_vectors[path]) for path in core_paths}
    vector_counts = Counter(vectors.values())
    ordered_vectors = sorted(vector_counts, key=lambda vector: (-vector_counts[vector], vector))
    type_by_vector = {vector: f"SVH{index:04d}" for index, vector in enumerate(ordered_vectors, 1)}
    type_paths = defaultdict(list)
    for path, vector in vectors.items():
        type_paths[type_by_vector[vector]].append(path)
    write_tsv(
        output_dir / "structural_variant_types.tsv",
        (
            "variant_type", "path_count", "path_frequency", "variable_region_count",
            "regional_alleles", "path_ids",
        ),
        (
            (
                type_id,
                len(paths),
                f"{len(paths) / len(core_paths):.9f}",
                len(vector),
                ",".join(vector) or ".",
                ",".join(sorted(paths)),
            )
            for vector, type_id in sorted(type_by_vector.items(), key=lambda item: item[1])
            for paths in [type_paths[type_id]]
        ),
    )
    path_rows = []
    for path in sorted(by_path):
        if path in outlier_paths:
            path_rows.append((
                path, "flank_outlier", ".", ".", signed_text(occurrence_signature(by_path[path]))
            ))
        elif path in vectors:
            vector = vectors[path]
            path_rows.append((
                path,
                "core",
                type_by_vector[vector],
                ",".join(vector) or ".",
                signed_text(occurrence_signature(by_path[path])),
            ))
        else:
            path_rows.append((
                path, "not_analyzed", ".", ".", signed_text(occurrence_signature(by_path[path]))
            ))
    write_tsv(
        output_dir / "path_variant_types.tsv",
        (
            "path_id", "status", "variant_type", "regional_alleles",
            "full_signed_bundle_haplotype",
        ),
        path_rows,
    )
    return {
        "region_count": len(region_summaries),
        "variable_region_count": len(variable_region_ids),
        "variable_region_ids": variable_region_ids,
        "variant_type_count": len(type_by_vector),
        "region_bed_files": region_bed_files,
        "reference_bed_region_count": len(reference_bed_rows),
    }


def convert_bundle_haplotypes(
    pgr_bed: Path,
    output_dir: Path,
    fai: Optional[Path] = None,
    max_nonconserved_flank_fraction: float = 0.01,
    flank_bundles: int = 1,
    reference_path: Optional[str] = None,
):
    pgr_bed = Path(pgr_bed)
    output_dir = Path(output_dir)
    fai = Path(fai) if fai is not None else None
    if not pgr_bed.is_file():
        raise ValueError(f"PGR BED does not exist: {pgr_bed}")
    if fai is not None and not fai.is_file():
        raise ValueError(f"FASTA index does not exist: {fai}")
    if not math.isfinite(max_nonconserved_flank_fraction) or not (
        0 <= max_nonconserved_flank_fraction <= 1
    ):
        raise ValueError("maximum non-conserved flank fraction must be in [0, 1]")
    if not isinstance(flank_bundles, int) or flank_bundles < 1:
        raise ValueError("flank bundle count must be a positive integer")

    by_path = read_pgr_bed(pgr_bed)
    source_rows = read_source_rows(pgr_bed)
    bed_path_ids = set(by_path)
    if fai is None:
        path_ids = sorted(bed_path_ids)
        fai_paths = None
        missing_from_bed = []
        missing_from_fai = []
    else:
        fai_paths = read_fai_paths(fai)
        fai_path_ids = set(fai_paths)
        path_ids = sorted(bed_path_ids & fai_path_ids)
        missing_from_bed = sorted(fai_path_ids - bed_path_ids)
        missing_from_fai = sorted(bed_path_ids - fai_path_ids)
    protected_inputs = {pgr_bed.resolve()}
    if fai is not None:
        protected_inputs.add(fai.resolve())
    region_dir = output_dir / "subvariation_regions"
    if (pgr_bed.resolve().parent == region_dir.resolve()
            and re.fullmatch(r"R\d{4,}\.pgr\.bed", pgr_bed.name)):
        raise ValueError("input BED would be removed while refreshing regional BED outputs")
    for name in OUTPUT_FILES:
        if (output_dir / name).resolve() in protected_inputs:
            raise ValueError(f"input may not be a workflow output: {output_dir / name}")
    output_dir.mkdir(parents=True, exist_ok=True)
    # This command is safe to rerun in the same directory.  Remove only files
    # owned by this workflow so a failed rerun cannot leave apparently valid
    # downstream tables from an earlier successful run.
    for name in OUTPUT_FILES:
        owned_output = output_dir / name
        if owned_output.is_file() or owned_output.is_symlink():
            owned_output.unlink()
    if region_dir.is_symlink():
        raise ValueError(f"regional BED output directory may not be a symlink: {region_dir}")
    region_dir.mkdir(exist_ok=True)
    for old_bed in region_dir.glob("R*.pgr.bed"):
        if re.fullmatch(r"R\d{4,}\.pgr\.bed", old_bed.name) and (
            old_bed.is_file() or old_bed.is_symlink()
        ):
            old_bed.unlink()
    write_tsv(
        output_dir / "skipped_paths.tsv",
        ("path_id", "reason"),
        [(path, "missing_from_bed") for path in missing_from_bed]
        + [(path, "missing_from_fai") for path in missing_from_fai],
    )
    summary_path = output_dir / "bundle_variant_summary.json"
    base_summary = {
        "version": 7,
        "completed": False,
        "inputs": {
            "pgr_bed": str(pgr_bed.resolve()),
            "pgr_bed_sha256": file_sha256(pgr_bed),
            "fai": str(fai.resolve()) if fai is not None else None,
            "fai_sha256": file_sha256(fai) if fai is not None else None,
        },
        "parameters": {
            "max_nonconserved_flank_fraction": max_nonconserved_flank_fraction,
            "flank_bundles": flank_bundles,
            "reference_path": reference_path,
            "anchor_rule": "present exactly once with one orientation on every core path",
            "anchor_objective": "maximum-cardinality all-path order-consistent chain",
        },
        "path_roster": {
            "bed_path_count": len(bed_path_ids),
            "fai_path_count": len(fai_paths) if fai_paths is not None else None,
            "shared_path_count": len(path_ids),
            "missing_from_bed_count": len(missing_from_bed),
            "missing_from_fai_count": len(missing_from_fai),
            "missing_from_bed": missing_from_bed,
            "missing_from_fai": missing_from_fai,
        },
    }
    summary_path.write_text(json.dumps(base_summary, indent=2, sort_keys=True) + "\n")

    def finish_original_fallback(
        reason: str,
        flank_summary: Mapping[str, object],
        core_paths: Sequence[str],
        outlier_paths: Iterable[str],
        selected_anchors: Sequence[str] = (),
        selected_reference: Optional[str] = None,
    ):
        fallback_bed = output_dir / "original_pgrtk.bed"
        shutil.copyfile(pgr_bed, fallback_bed)
        fallback_region_bed = region_dir / "R0001.pgr.bed"
        shutil.copyfile(pgr_bed, fallback_region_bed)
        outputs = {
            name: str((output_dir / name).resolve())
            for name in (
                "flank_check.tsv", "flank_summary.json", "skipped_paths.tsv",
                "anchor_candidates.tsv", "anchors.tsv",
            )
            if (output_dir / name).is_file()
        }
        outputs["original_pgrtk.bed"] = str(fallback_bed.resolve())
        outputs["subvariation_regions"] = {
            "R0001": str(fallback_region_bed.resolve()),
        }
        summary = dict(base_summary)
        summary.update({
            "completed": True,
            "mode": "original_pgrtk_fallback",
            "fallback_reason": reason,
            "flanks": flank_summary,
            "parameters": {
                **base_summary["parameters"],
                "reference_path": selected_reference or reference_path,
            },
            "counts": {
                "input_path_count": len(path_ids),
                "core_path_count": len(core_paths),
                "flank_outlier_path_count": len(outlier_paths),
                "skipped_path_count": len(missing_from_bed) + len(missing_from_fai),
                "selected_anchor_count": len(selected_anchors),
                "merged_anchor_count": 0,
                "region_count": 1,
                "variable_region_count": 0,
                "variant_type_count": 0,
            },
            "outputs": outputs,
        })
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        return summary

    try:
        core_paths, outlier_paths, flank_summary = write_flank_audit(
            output_dir,
            by_path,
            path_ids,
            flank_bundles,
            max_nonconserved_flank_fraction,
        )
        fallback_reason = None
        if not path_ids:
            fallback_reason = "no_shared_paths"
        elif not flank_summary["passes"]:
            fallback_reason = "nonconserved_flank_fraction_exceeded"
        elif not core_paths:
            fallback_reason = "no_conserved_core_paths"
        elif reference_path is not None and reference_path not in core_paths:
            fallback_reason = (
                "requested_reference_missing_from_shared_paths"
                if reference_path not in path_ids else "requested_reference_flank_outlier"
            )
        if fallback_reason is not None:
            return finish_original_fallback(
                fallback_reason, flank_summary, core_paths, outlier_paths,
            )
        anchors, orientations, positions, selected_reference = select_anchors(
            output_dir, by_path, core_paths, reference_path,
        )
        if len(anchors) < 2:
            return finish_original_fallback(
                "insufficient_anchors", flank_summary, core_paths, outlier_paths,
                anchors, selected_reference,
            )
        merged_anchors = merge_invariant_anchors(
            output_dir, by_path, core_paths, anchors, positions, selected_reference,
        )
        region_summary = decompose_regions(
            output_dir, by_path, source_rows, core_paths, outlier_paths,
            merged_anchors, selected_reference,
        )
        if region_summary["region_count"] == 0:
            # An invariant population has no inferred variant region, but
            # downstream consumers still need one native BED for the whole
            # original region, as in the original-BED fallback modes.
            whole_region_bed = region_dir / "R0001.pgr.bed"
            shutil.copyfile(pgr_bed, whole_region_bed)
            region_summary["region_bed_files"]["R0001"] = str(whole_region_bed.resolve())
    except (OSError, ValueError) as error:
        failed = dict(base_summary)
        failed["error"] = str(error)
        summary_path.write_text(json.dumps(failed, indent=2, sort_keys=True) + "\n")
        raise

    outputs = {
        name: str((output_dir / name).resolve())
        for name in (
            "flank_check.tsv",
            "flank_summary.json",
            "skipped_paths.tsv",
            "anchor_candidates.tsv",
            "anchors.tsv",
            "merged_anchors.tsv",
            "merged_anchors.bed",
            "variant_regions.tsv",
            "subvariation_regions.ref.bed",
            "variant_alleles.tsv",
            "path_variant_regions.tsv",
            "structural_variant_types.tsv",
            "path_variant_types.tsv",
        )
    }
    outputs["subvariation_regions"] = region_summary["region_bed_files"]
    summary = dict(base_summary)
    summary.update({
        "completed": True,
        "mode": "decomposed",
        "parameters": {
            **base_summary["parameters"],
            "reference_path": selected_reference,
        },
        "flanks": flank_summary,
        "counts": {
            "input_path_count": len(path_ids),
            "core_path_count": len(core_paths),
            "flank_outlier_path_count": len(outlier_paths),
            "skipped_path_count": len(missing_from_bed) + len(missing_from_fai),
            "selected_anchor_count": len(anchors),
            "merged_anchor_count": len(merged_anchors),
            "region_count": region_summary["region_count"],
            "reference_bed_region_count": region_summary["reference_bed_region_count"],
            "variable_region_count": region_summary["variable_region_count"],
            "variant_type_count": region_summary["variant_type_count"],
        },
        "variable_region_ids": region_summary["variable_region_ids"],
        "selected_anchors": [
            {"bundle_id": bundle_id, "direction": orientations[bundle_id]}
            for bundle_id in anchors
        ],
        "merged_anchors": [
            {
                "anchor_id": anchor.anchor_id,
                "signed_bundle_order": signed_text(anchor.signature),
                "reference_start": anchor.reference_start,
                "reference_end": anchor.reference_end,
            }
            for anchor in merged_anchors
        ],
        "outputs": outputs,
    })
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pgr_bed", type=Path, help="four-column PGR-TK-like BED")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--fai", type=Path,
        help="optional FASTA index; analyze only chromosome IDs shared with the BED and report missing IDs",
    )
    parser.add_argument(
        "--max-nonconserved-flank-fraction", type=float, default=0.01,
        help="return the original PGR-TK BED when the union of left/right terminal outliers exceeds this fraction (default: 0.01)",
    )
    parser.add_argument(
        "--flank-bundles", type=int, default=1,
        help="number of signed terminal bundles compared on each side (default: 1)",
    )
    parser.add_argument(
        "--ref-id", "--reference-path", dest="reference_path",
        help="core path used for anchor coordinates and tie-breaking; random core path if omitted",
    )
    args = parser.parse_args(argv)
    try:
        summary = convert_bundle_haplotypes(
            args.pgr_bed,
            args.output_dir,
            args.fai,
            args.max_nonconserved_flank_fraction,
            args.flank_bundles,
            args.reference_path,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
