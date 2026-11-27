#!/usr/bin/env python3
"""Refine minimap2 bundle-to-genome mappings and find newly covered sequence.

The expected input follows the historical ``scripts/refine_v2.sh`` command, where the
bundle consensuses are the minimap2 target and population paths are queries::

    minimap2 -x asm20 -c bundle.fa population.fa > bundle_to_genome.paf

PAF axes are converted on input to the logical names used by this module:
``query`` means the bundle consensus and ``target`` means the population path.
This keeps the existing TSV/BED6 output semantics stable while matching the
actual minimap2 invocation.

Filtering by MAPQ and identity is applied to individual PAF records.  Query
coverage is deliberately applied *after* adjacent, collinear records have been
chained, so two partial records can together pass the coverage threshold.
Final-remapping callers can enable cross-bundle locus deduplication. Its stricter
coverage threshold is applied to the union coverage of each completed chain,
not to native PAF records, before repeated local-neighborhood maximum selection
is applied to original genome intervals.

Coordinates in all BED outputs are zero-based, half-open. An optional merged
PGR BED supplies old coverage for novelty checks. Without it, refinement does
not assess or report novelty; a separate metadata BED supplies IDs and sizes.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import sys
from bisect import bisect_left
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, Tuple

try:
    from .final_bundle_filter import filter_final_bundle_classes, read_sequence_lengths, write_final_filter_audits
except ImportError:
    from final_bundle_filter import filter_final_bundle_classes, read_sequence_lengths, write_final_filter_audits


Interval = Tuple[int, int]


@dataclass(frozen=True)
class PafBlock:
    raw_query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    matching_bases: int
    alignment_block_length: int
    mapq: int
    line_number: int

    @property
    def identity(self) -> float:
        if self.alignment_block_length == 0:
            return 0.0
        return self.matching_bases / self.alignment_block_length


@dataclass
class BlockChain:
    raw_query_name: str
    target_name: str
    strand: str
    blocks: List[PafBlock]


@dataclass(frozen=True)
class RefinedMapping:
    mapping_id: str
    bundle_name: str
    raw_query_name: str
    query_length: int
    query_start: int
    query_end: int
    query_aligned_bp: int
    query_coverage: float
    strand: str
    target_name: str
    target_start: int
    target_end: int
    target_footprint: Tuple[Interval, ...]
    target_footprint_bp: int
    matching_bases: int
    alignment_block_length: int
    identity: float
    min_mapq: int
    max_mapq: int
    block_count: int
    source_line_numbers: Tuple[int, ...]

    @property
    def target_span_bp(self) -> int:
        return self.target_end - self.target_start

    @property
    def target_gap_bp(self) -> int:
        return self.target_span_bp - self.target_footprint_bp


@dataclass(frozen=True)
class DuplicateMapping:
    dropped: RefinedMapping
    retained: RefinedMapping
    overlap_bp: int
    overlap_fraction_dropped: float
    overlap_fraction_retained: float


@dataclass(frozen=True)
class CrossBundleLocusDecision:
    pass_number: int
    reference: RefinedMapping
    selected: RefinedMapping
    reference_retained: bool
    overlap_bp: int
    overlap_fraction_reference: float
    overlap_fraction_selected: float
    score_field: str


@dataclass(frozen=True)
class NovelAssessment:
    mapping: RefinedMapping
    old_covered_bp: int
    novel_intervals: Tuple[Interval, ...]
    novel_bp: int
    required_novel_bp: int
    novel_fraction_of_bundle: float
    is_novel: bool


@dataclass(frozen=True)
class ExistingBundleMetadata:
    bundle_id: str
    bundle_size: int


@dataclass(frozen=True)
class BundleClassAssignment:
    bundle_name: str
    numeric_bundle_id: str
    id_source: str
    bundle_size: int
    bundle_size_source: str
    raw_query_names: Tuple[str, ...]
    mapping_count: int


def open_text(path: Path):
    """Open plain text or gzip-compressed text based on the file suffix."""

    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_nonnegative_int(value: str, field: str, path: Path, line_number: int) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_number}: {field} must be an integer, got {value!r}"
        ) from exc
    if parsed < 0:
        raise ValueError(
            f"{path}:{line_number}: {field} must be non-negative, got {parsed}"
        )
    return parsed


def read_paf(path: Path) -> List[PafBlock]:
    """Read historical-refine PAF and normalize it to bundle/query semantics.

    Raw minimap2 PAF columns describe ``population path`` as the query
    (columns 1--4) and ``bundle consensus`` as the target (columns 6--9).
    ``PafBlock`` deliberately exposes the reverse logical terminology so the
    rest of the refinement code and its established outputs consistently call
    the bundle the query and the population path the target.
    """

    records: List[PafBlock] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError(
                    f"{path}:{line_number}: PAF requires at least 12 tab-separated "
                    f"fields, found {len(fields)}"
                )

            # Raw PAF query: population/genome path.
            genome_length = parse_nonnegative_int(fields[1], "PAF query length", path, line_number)
            genome_start = parse_nonnegative_int(fields[2], "PAF query start", path, line_number)
            genome_end = parse_nonnegative_int(fields[3], "PAF query end", path, line_number)
            # Raw PAF target: bundle consensus.
            bundle_length = parse_nonnegative_int(fields[6], "PAF target length", path, line_number)
            bundle_start = parse_nonnegative_int(fields[7], "PAF target start", path, line_number)
            bundle_end = parse_nonnegative_int(fields[8], "PAF target end", path, line_number)
            matching_bases = parse_nonnegative_int(fields[9], "matching bases", path, line_number)
            alignment_length = parse_nonnegative_int(fields[10], "alignment block length", path, line_number)
            mapq = parse_nonnegative_int(fields[11], "MAPQ", path, line_number)

            if fields[4] not in {"+", "-"}:
                raise ValueError(
                    f"{path}:{line_number}: strand must be '+' or '-', got {fields[4]!r}"
                )
            if genome_length == 0 or bundle_length == 0:
                raise ValueError(
                    f"{path}:{line_number}: PAF query and target lengths must be positive"
                )
            if not (0 <= genome_start < genome_end <= genome_length):
                raise ValueError(
                    f"{path}:{line_number}: invalid population-path interval "
                    f"[{genome_start}, {genome_end}) for length {genome_length}"
                )
            if not (0 <= bundle_start < bundle_end <= bundle_length):
                raise ValueError(
                    f"{path}:{line_number}: invalid bundle interval "
                    f"[{bundle_start}, {bundle_end}) for length {bundle_length}"
                )
            if alignment_length == 0:
                raise ValueError(f"{path}:{line_number}: alignment block length must be positive")
            if matching_bases > alignment_length:
                raise ValueError(
                    f"{path}:{line_number}: matching bases ({matching_bases}) exceed "
                    f"alignment block length ({alignment_length})"
                )

            records.append(
                PafBlock(
                    raw_query_name=fields[5],
                    query_length=bundle_length,
                    query_start=bundle_start,
                    query_end=bundle_end,
                    strand=fields[4],
                    target_name=fields[0],
                    target_length=genome_length,
                    target_start=genome_start,
                    target_end=genome_end,
                    matching_bases=matching_bases,
                    alignment_block_length=alignment_length,
                    mapq=mapq,
                    line_number=line_number,
                )
            )
    return records


def merge_intervals(intervals: Iterable[Interval]) -> Tuple[Interval, ...]:
    ordered = sorted((start, end) for start, end in intervals if end > start)
    if not ordered:
        return ()
    merged: List[List[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        previous = merged[-1]
        if start <= previous[1]:
            previous[1] = max(previous[1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def interval_total(intervals: Iterable[Interval]) -> int:
    return sum(end - start for start, end in intervals)


def intersect_intervals(left: Sequence[Interval], right: Sequence[Interval]) -> Tuple[Interval, ...]:
    answer: List[Interval] = []
    i = 0
    j = 0
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        end = min(left[i][1], right[j][1])
        if end > start:
            answer.append((start, end))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(answer)


def subtract_intervals(source: Sequence[Interval], mask: Sequence[Interval]) -> Tuple[Interval, ...]:
    """Return source - mask for sorted, non-overlapping half-open intervals."""

    answer: List[Interval] = []
    mask_index = 0
    for source_start, source_end in source:
        cursor = source_start
        while mask_index < len(mask) and mask[mask_index][1] <= cursor:
            mask_index += 1
        index = mask_index
        while index < len(mask) and mask[index][0] < source_end:
            mask_start, mask_end = mask[index]
            if mask_start > cursor:
                answer.append((cursor, min(mask_start, source_end)))
            cursor = max(cursor, mask_end)
            if cursor >= source_end:
                break
            index += 1
        if cursor < source_end:
            answer.append((cursor, source_end))
    return tuple(answer)


def read_merged_bed(path: Path) -> Dict[str, Tuple[Interval, ...]]:
    by_target: MutableMapping[str, List[Interval]] = defaultdict(list)
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{line_number}: BED requires at least three tab-separated fields"
                )
            start = parse_nonnegative_int(fields[1], "BED start", path, line_number)
            end = parse_nonnegative_int(fields[2], "BED end", path, line_number)
            if end <= start:
                raise ValueError(
                    f"{path}:{line_number}: BED end must be greater than start, got [{start}, {end})"
                )
            by_target[fields[0]].append((start, end))
    return {target: merge_intervals(intervals) for target, intervals in by_target.items()}


def read_existing_bundle_metadata(
    path: Path,
) -> Tuple[Tuple[str, ...], Dict[str, ExistingBundleMetadata]]:
    """Read reusable numeric IDs and sizes from a PGR-like fourth column.

    The existing BED3 coverage behavior remains valid: files without a native
    PGR fourth column simply provide no reusable IDs.  Malformed non-PGR names
    are ignored here because they may be ordinary BED6 labels; malformed BED
    coordinates have already been rejected by :func:`read_merged_bed`.
    """

    numeric_ids = set()
    sizes: Dict[str, ExistingBundleMetadata] = {}
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            annotation = fields[3].split(":")
            if len(annotation) < 2 or not annotation[0].isdigit():
                continue
            bundle_id = str(int(annotation[0]))
            numeric_ids.add(bundle_id)
            if not annotation[1].isdigit() or int(annotation[1]) <= 0:
                continue
            bundle_size = int(annotation[1])
            previous = sizes.get(bundle_id)
            if previous is not None and previous.bundle_size != bundle_size:
                raise ValueError(
                    f"{path}:{line_number}: existing bundle {bundle_id} has "
                    f"inconsistent sizes {previous.bundle_size} and {bundle_size}"
                )
            sizes[bundle_id] = ExistingBundleMetadata(bundle_id, bundle_size)
    return tuple(sorted(numeric_ids, key=int)), sizes


def parse_alias(value: str) -> Tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            f"bundle alias must have OLD=NEW form, got {value!r}"
        )
    old, new = value.split("=", 1)
    if not old or not new:
        raise argparse.ArgumentTypeError(
            f"bundle alias must have non-empty OLD and NEW names, got {value!r}"
        )
    return old, new


def resolve_alias(name: str, aliases: Mapping[str, str]) -> str:
    """Resolve exact aliases transitively while detecting accidental cycles."""

    current = name
    seen = set()
    while current in aliases and aliases[current] != current:
        if current in seen:
            chain = " -> ".join([*sorted(seen), current])
            raise ValueError(f"cycle in bundle aliases while resolving {name!r}: {chain}")
        seen.add(current)
        current = aliases[current]
    return current


def normalize_bundle_name(name: str, aliases: Mapping[str, str]) -> str:
    """Apply exact aliases and recover numeric IDs from pipeline FASTA names."""

    normalized = resolve_alias(name, aliases)
    if normalized.startswith("bundle_"):
        normalized = normalized[len("bundle_") :]
        allele_marker = normalized.rfind(".allele_")
        if allele_marker >= 0 and normalized[allele_marker + 8 :].isdigit():
            normalized = normalized[:allele_marker]
        # Also support aliases written for the class name rather than for the
        # full pipeline FASTA record (for example bundle_dup1_1 -> dup1).
        normalized = resolve_alias(normalized, aliases)
    if normalized.isdigit():
        normalized = str(int(normalized))
    return normalized


def natural_name_key(value: str) -> Tuple[Tuple[int, object], ...]:
    """Return a deterministic numeric-aware key (``dup2`` before ``dup10``)."""

    return tuple(
        (0, int(part)) if part.isdigit() else (1, part)
        for part in re.split(r"(\d+)", value.casefold())
        if part
    )


def deterministic_synthetic_bundle_size(bundle_name: str) -> int:
    """Match the 100--999 deterministic placeholder convention of merge output."""

    digest = hashlib.blake2b(bundle_name.encode("utf-8"), digest_size=8).digest()
    return 100 + int.from_bytes(digest, "big") % 900


def assign_bundle_class_ids(
    mappings: Sequence[RefinedMapping],
    existing_numeric_ids: Sequence[str],
    existing_metadata: Mapping[str, ExistingBundleMetadata],
) -> Dict[str, BundleClassAssignment]:
    """Assign stable numeric IDs while preserving recognizable existing IDs."""

    classes = sorted({mapping.bundle_name for mapping in mappings}, key=natural_name_key)
    existing_ids = set(existing_numeric_ids)
    used_ids = set(existing_ids)

    chosen_ids: Dict[str, Tuple[str, str]] = {}
    for bundle_name in classes:
        if bundle_name.isdigit():
            numeric_id = str(int(bundle_name))
            source = (
                "reused_existing_numeric_id"
                if numeric_id in existing_ids
                else "preserved_numeric_class_name"
            )
            chosen_ids[bundle_name] = numeric_id, source
            used_ids.add(numeric_id)

    next_numeric_id = max((int(value) for value in used_ids), default=0) + 1
    for bundle_name in classes:
        if bundle_name in chosen_ids:
            continue
        while str(next_numeric_id) in used_ids:
            next_numeric_id += 1
        numeric_id = str(next_numeric_id)
        chosen_ids[bundle_name] = numeric_id, "allocated_after_existing_max"
        used_ids.add(numeric_id)
        next_numeric_id += 1

    raw_names: MutableMapping[str, set] = defaultdict(set)
    mapping_counts: Counter[str] = Counter()
    for mapping in mappings:
        raw_names[mapping.bundle_name].add(mapping.raw_query_name)
        mapping_counts[mapping.bundle_name] += 1

    assignments: Dict[str, BundleClassAssignment] = {}
    for bundle_name in classes:
        numeric_id, id_source = chosen_ids[bundle_name]
        metadata = existing_metadata.get(numeric_id)
        if metadata is not None and bundle_name == numeric_id:
            bundle_size = metadata.bundle_size
            size_source = "reused_existing_pgr_bundle_size"
        else:
            bundle_size = deterministic_synthetic_bundle_size(bundle_name)
            size_source = "deterministic_synthetic_placeholder"
        assignments[bundle_name] = BundleClassAssignment(
            bundle_name=bundle_name,
            numeric_bundle_id=numeric_id,
            id_source=id_source,
            bundle_size=bundle_size,
            bundle_size_source=size_source,
            raw_query_names=tuple(sorted(raw_names[bundle_name], key=natural_name_key)),
            mapping_count=mapping_counts[bundle_name],
        )
    return assignments


def blocks_are_adjacent(
    previous: PafBlock,
    following: PafBlock,
    max_genome_gap: int = 1_000,
    max_bundle_gap: int = 1_000,
    max_overlap_fraction: float = 0.10,
) -> Optional[Tuple[int, int]]:
    """Return two-axis gaps for a collinear, locally chainable pair.

    Negative gaps are overlaps. Their magnitude must satisfy both the absolute
    axis-specific gap limit and ``max_overlap_fraction`` of the original bundle
    consensus length. The fraction is checked independently on both axes. This
    local predicate does not check overlap with fixed other-class mappings;
    ``chain_blocks`` applies that additional whole-candidate-span check.
    """

    if (
        previous.raw_query_name != following.raw_query_name
        or previous.target_name != following.target_name
        or previous.strand != following.strand
        or previous.query_length != following.query_length
        or previous.target_length != following.target_length
    ):
        return None

    # Both axes must make genuine progress.  This prevents two mappings of the
    # complete query to tandem target copies from being mistaken for one chain.
    if following.target_start < previous.target_start or following.target_end <= previous.target_end:
        return None
    genome_gap = following.target_start - previous.target_end

    if previous.strand == "+":
        if (
            following.query_start < previous.query_start
            or following.query_end <= previous.query_end
        ):
            return None
        bundle_gap = following.query_start - previous.query_end
    else:
        if (
            following.query_start >= previous.query_start
            or following.query_end > previous.query_end
        ):
            return None
        bundle_gap = previous.query_start - following.query_end

    # Small overlaps are allowed as well as gaps.  A very large overlap is not
    # "adjacent" and commonly represents competing secondary alignments.
    if not (-max_genome_gap <= genome_gap <= max_genome_gap):
        return None
    if not (-max_bundle_gap <= bundle_gap <= max_bundle_gap):
        return None
    max_overlap_bp = max_overlap_fraction * previous.query_length
    if max(0, -genome_gap) > max_overlap_bp:
        return None
    if max(0, -bundle_gap) > max_overlap_bp:
        return None
    return genome_gap, bundle_gap


@dataclass(frozen=True)
class _PathFixedBundleIndex:
    """Static start-sorted interval index with subtree maximum endpoints.

    A small segment tree prunes footprints ending before the candidate and limits
    the search to footprints starting before its other boundary. It avoids
    scanning every bundle on the population path for each candidate extension.
    Each entry is one fixed PAF occurrence with its class and original bundle
    length. Occurrences are tested individually, never cumulatively summed.
    """

    intervals: Tuple[Tuple[int, int, str, int], ...]
    starts: Tuple[int, ...]
    leaf_count: int
    max_ends: Tuple[int, ...]

    @classmethod
    def build(cls, intervals: Iterable[Tuple[int, int, str, int]]) -> "_PathFixedBundleIndex":
        ordered = tuple(sorted(set(intervals)))
        leaf_count = 1
        while leaf_count < len(ordered):
            leaf_count *= 2
        max_ends = [0] * (2 * leaf_count)
        for index, (_, end, _, _) in enumerate(ordered):
            max_ends[leaf_count + index] = end
        for index in range(leaf_count - 1, 0, -1):
            max_ends[index] = max(max_ends[2 * index], max_ends[2 * index + 1])
        return cls(
            ordered,
            tuple(start for start, _, _, _ in ordered),
            leaf_count,
            tuple(max_ends),
        )

    def exceeds_overlap(
        self,
        start: int,
        end: int,
        bundle_class: str,
        max_overlap_fraction: float,
        max_overlap_bp: int,
    ) -> bool:
        """Reject when one different-class fixed occurrence exceeds either limit.

        ``[start, end)`` is the entire proposed chain's genome outer span, not
        just its joining gap. Fractional overlap uses the fixed occurrence's
        original consensus length as denominator. Both thresholds are strict:
        exact equality passes, and merely touching an interval does not count.
        """

        # Exclude fixed footprints beginning exactly at the right boundary.
        limit = bisect_left(self.starts, end)
        if limit == 0:
            return False
        pending = [(1, 0, self.leaf_count)]
        while pending:
            node, left, right = pending.pop()
            if left >= limit or self.max_ends[node] <= start:
                continue
            if right - left == 1:
                fixed_start, fixed_end, fixed_class, fixed_length = self.intervals[left]
                if fixed_class == bundle_class:
                    continue
                overlap_bp = min(end, fixed_end) - max(start, fixed_start)
                if (
                    overlap_bp > max_overlap_bp
                    or overlap_bp > max_overlap_fraction * fixed_length
                ):
                    return True
                continue
            middle = (left + right) // 2
            pending.append((2 * node + 1, middle, right))
            pending.append((2 * node, left, middle))
        return False


def _greedy_chain_blocks(
    blocks: Sequence[PafBlock],
    max_genome_gap: int,
    max_bundle_gap: int,
    max_overlap_fraction: float,
    aliases: Mapping[str, str],
    fixed_indexes: Mapping[str, _PathFixedBundleIndex],
    max_fixed_overlap_fraction: float,
    max_fixed_overlap_bp: int,
) -> List[BlockChain]:
    """Chain partial raw-consensus blocks while honoring an immutable fixed set."""

    grouped: MutableMapping[Tuple[str, str, str], List[PafBlock]] = defaultdict(list)
    for block in blocks:
        # Important: chain on raw query name.  Aliases describe a bundle class,
        # but do not make coordinates from distinct consensus sequences shared.
        grouped[(block.raw_query_name, block.target_name, block.strand)].append(block)

    output: List[BlockChain] = []
    for key in sorted(grouped):
        raw_query_name, target_name, strand = key
        bundle_class = normalize_bundle_name(raw_query_name, aliases)
        fixed_index = fixed_indexes.get(target_name)
        ordered = sorted(
            grouped[key],
            key=lambda block: (
                block.target_start,
                block.target_end,
                block.query_start if strand == "+" else -block.query_end,
                block.line_number,
            ),
        )
        chains: List[BlockChain] = []
        for block in ordered:
            compatible: List[Tuple[int, int, int, int]] = []
            for chain_index, chain in enumerate(chains):
                gaps = blocks_are_adjacent(
                    chain.blocks[-1], block, max_genome_gap, max_bundle_gap,
                    max_overlap_fraction,
                )
                if gaps is None:
                    continue
                if fixed_index is not None and fixed_index.exceeds_overlap(
                    chain.blocks[0].target_start, block.target_end, bundle_class,
                    max_fixed_overlap_fraction, max_fixed_overlap_bp,
                ):
                    continue
                genome_gap, bundle_gap = gaps
                # Prefer the closest extension; then the chain whose last block
                # is nearest on the target, for deterministic handling of ties.
                compatible.append(
                    (
                        abs(genome_gap) + abs(bundle_gap),
                        abs(genome_gap),
                        -chain.blocks[-1].target_end,
                        chain_index,
                    )
                )
            if compatible:
                chosen_index = min(compatible)[-1]
                chains[chosen_index].blocks.append(block)
            else:
                chains.append(BlockChain(raw_query_name, target_name, strand, [block]))
        output.extend(chains)
    return output


def chain_blocks(
    blocks: Sequence[PafBlock],
    max_genome_gap: int = 1_000,
    max_bundle_gap: int = 1_000,
    *,
    max_overlap_fraction: float = 0.10,
    fixed_bundle_min_coverage: float = 0.80,
    max_fixed_overlap_fraction: float = 0.20,
    max_fixed_overlap_bp: int = 10_000,
    aliases: Optional[Mapping[str, str]] = None,
) -> List[BlockChain]:
    """Fix sufficiently covered PAF blocks, then greedily chain only partials.

    The caller supplies PAF blocks already eligible by MAPQ and identity. Each
    individual block covering at least ``fixed_bundle_min_coverage`` of its
    original consensus becomes an immutable, singleton output chain. Fixed
    blocks never join another block, and newly merged partial chains never
    become fixed mappings during this call.

    Remaining partial blocks are extended only if locally collinear and within
    the absolute gap and pairwise-overlap limits. Before each extension, test
    the entire proposed chain's genome outer span against every intersecting
    different-class fixed occurrence on the same path, irrespective of strand.
    Reject if its overlap with any one fixed occurrence is strictly greater
    than ``max_fixed_overlap_fraction`` of that fixed consensus's original
    length OR strictly greater than ``max_fixed_overlap_bp``. Independent fixed
    occurrences do not have their overlaps added together.

    Chaining itself always uses raw consensus names: aliases and allele names
    identify classes for fixed checks but never combine distinct coordinate axes.
    Built-in dup1 aliases are retained, and explicit aliases can override them.
    """

    if max_genome_gap < 0 or max_bundle_gap < 0:
        raise ValueError("chaining gap limits must be non-negative")
    if not 0.0 <= max_overlap_fraction <= 1.0:
        raise ValueError("max_overlap_fraction must be in [0, 1]")
    if not 0.0 <= fixed_bundle_min_coverage <= 1.0:
        raise ValueError("fixed_bundle_min_coverage must be in [0, 1]")
    if not 0.0 <= max_fixed_overlap_fraction <= 1.0:
        raise ValueError("max_fixed_overlap_fraction must be in [0, 1]")
    if max_fixed_overlap_bp < 0:
        raise ValueError("max_fixed_overlap_bp must be non-negative")
    class_aliases: Dict[str, str] = {"dup1_1": "dup1", "dup1_2": "dup1"}
    class_aliases.update(aliases or {})
    fixed_chains: List[BlockChain] = []
    partial_blocks: List[PafBlock] = []
    fixed_intervals: MutableMapping[str, List[Tuple[int, int, str, int]]] = defaultdict(list)
    for block in blocks:
        coverage = (block.query_end - block.query_start) / block.query_length
        if coverage < fixed_bundle_min_coverage:
            partial_blocks.append(block)
            continue
        fixed_chains.append(
            BlockChain(block.raw_query_name, block.target_name, block.strand, [block])
        )
        bundle_class = normalize_bundle_name(block.raw_query_name, class_aliases)
        fixed_intervals[block.target_name].append(
            (block.target_start, block.target_end, bundle_class, block.query_length)
        )
    fixed_indexes = {
        path: _PathFixedBundleIndex.build(intervals)
        for path, intervals in fixed_intervals.items()
    }
    partial_chains = _greedy_chain_blocks(
        partial_blocks, max_genome_gap, max_bundle_gap, max_overlap_fraction,
        class_aliases, fixed_indexes, max_fixed_overlap_fraction, max_fixed_overlap_bp,
    )
    return sorted(
        fixed_chains + partial_chains,
        key=lambda chain: (
            chain.raw_query_name, chain.target_name, chain.strand,
            chain.blocks[0].target_start, chain.blocks[-1].target_end,
            chain.blocks[0].line_number,
        ),
    )


def chain_to_mapping(
    chain: BlockChain,
    bundle_name: str,
    temporary_id: str,
) -> RefinedMapping:
    blocks = chain.blocks
    query_lengths = {block.query_length for block in blocks}
    if len(query_lengths) != 1:
        raise ValueError(
            f"bundle consensus {chain.raw_query_name!r} has inconsistent lengths "
            f"in one chain: "
            f"{sorted(query_lengths)}"
        )
    query_length = next(iter(query_lengths))
    query_intervals = merge_intervals(
        (block.query_start, block.query_end) for block in blocks
    )
    target_footprint = merge_intervals(
        (block.target_start, block.target_end) for block in blocks
    )
    query_aligned_bp = interval_total(query_intervals)
    target_footprint_bp = interval_total(target_footprint)
    alignment_block_length = sum(block.alignment_block_length for block in blocks)
    matching_bases = sum(block.matching_bases for block in blocks)
    return RefinedMapping(
        mapping_id=temporary_id,
        bundle_name=bundle_name,
        raw_query_name=chain.raw_query_name,
        query_length=query_length,
        query_start=min(start for start, _ in query_intervals),
        query_end=max(end for _, end in query_intervals),
        query_aligned_bp=query_aligned_bp,
        query_coverage=query_aligned_bp / query_length,
        strand=chain.strand,
        target_name=chain.target_name,
        target_start=target_footprint[0][0],
        target_end=target_footprint[-1][1],
        target_footprint=target_footprint,
        target_footprint_bp=target_footprint_bp,
        matching_bases=matching_bases,
        alignment_block_length=alignment_block_length,
        identity=matching_bases / alignment_block_length,
        min_mapq=min(block.mapq for block in blocks),
        max_mapq=max(block.mapq for block in blocks),
        block_count=len(blocks),
        source_line_numbers=tuple(block.line_number for block in blocks),
    )


def mapping_quality_key(mapping: RefinedMapping) -> Tuple[float, float, int, int, int, str, int, int, str]:
    """Sort best same-locus representative first, deterministically."""

    return (
        -mapping.query_coverage,
        -mapping.identity,
        -mapping.min_mapq,
        -mapping.query_aligned_bp,
        -mapping.target_footprint_bp,
        mapping.raw_query_name,
        mapping.target_start,
        mapping.target_end,
        mapping.strand,
    )


def deduplicate_same_locus(
    mappings: Sequence[RefinedMapping],
    minimum_overlap_fraction: float,
) -> Tuple[List[RefinedMapping], List[DuplicateMapping]]:
    """Deduplicate aliases mapping to substantially the same genome footprint.

    Two mappings are the same locus when their normalized bundle name and target
    agree, they overlap, and overlap covers at least ``minimum_overlap_fraction``
    of either mapping's aligned genome footprint.  Strand is intentionally not a
    grouping key: opposite-strand alignments of one bundle at one locus compete.
    """

    grouped: MutableMapping[Tuple[str, str], List[RefinedMapping]] = defaultdict(list)
    for mapping in mappings:
        grouped[(mapping.bundle_name, mapping.target_name)].append(mapping)

    retained_all: List[RefinedMapping] = []
    duplicate_all: List[DuplicateMapping] = []
    for group_key in sorted(grouped):
        retained: List[RefinedMapping] = []
        for candidate in sorted(grouped[group_key], key=mapping_quality_key):
            duplicate_of: Optional[Tuple[RefinedMapping, int, float, float]] = None
            for selected in retained:
                overlap_bp = interval_total(
                    intersect_intervals(candidate.target_footprint, selected.target_footprint)
                )
                if overlap_bp == 0:
                    continue
                candidate_fraction = overlap_bp / candidate.target_footprint_bp
                selected_fraction = overlap_bp / selected.target_footprint_bp
                if max(candidate_fraction, selected_fraction) >= minimum_overlap_fraction:
                    duplicate_of = (
                        selected,
                        overlap_bp,
                        candidate_fraction,
                        selected_fraction,
                    )
                    break
            if duplicate_of is None:
                retained.append(candidate)
            else:
                selected, overlap_bp, candidate_fraction, selected_fraction = duplicate_of
                duplicate_all.append(
                    DuplicateMapping(
                        dropped=candidate,
                        retained=selected,
                        overlap_bp=overlap_bp,
                        overlap_fraction_dropped=candidate_fraction,
                        overlap_fraction_retained=selected_fraction,
                    )
                )
        retained_all.extend(retained)
    return retained_all, duplicate_all


def locus_maximum_key(mapping: RefinedMapping, score_field: str = "matching_bases") -> Tuple:
    """The awk-selected fifth column is the score; lexical names resolve ties."""

    return (
        -getattr(mapping, score_field),
        mapping.raw_query_name,
        mapping.target_start,
        mapping.target_end,
        mapping.bundle_name,
        mapping.strand,
        mapping.source_line_numbers,
        mapping.mapping_id,
    )


def deduplicate_cross_bundle_loci(
    mappings: Sequence[RefinedMapping],
    minimum_overlap_fraction: float = 0.80,
    passes: int = 2,
) -> Tuple[List[RefinedMapping], List[CrossBundleLocusDecision]]:
    """Repeated ``bedmap --max-element --fraction-either`` neighborhood passes.

    Every input genome interval independently selects the highest-scoring input
    overlapping at least the specified fraction of *either* outer BED span.
    The distinct selected original mappings become the next pass' input.  This
    is not greedy suppression or transitive connected-component collapse: an
    interval can survive because a neighbor selects it even if it selects a
    better neighbor itself.  Only the population path partitions competition;
    raw bundle name, normalized class and strand do not.

    The user's first awk puts matching bp in the score column; their second awk
    removes it and moves identity into that column.  Pass 1 therefore maximizes
    matching bp, and subsequent passes maximize identity.  Ties prefer lexical
    raw bundle names followed by stable coordinates, names and source lines.
    """

    if not 0.0 < minimum_overlap_fraction <= 1.0:
        raise ValueError("minimum_overlap_fraction must be in (0, 1]")
    if passes < 1:
        raise ValueError("passes must be at least 1")
    current = list(mappings)
    decisions: List[CrossBundleLocusDecision] = []
    for pass_number in range(1, passes + 1):
        score_field = "matching_bases" if pass_number == 1 else "identity"
        grouped: MutableMapping[str, List[RefinedMapping]] = defaultdict(list)
        for mapping in current:
            grouped[mapping.target_name].append(mapping)
        next_pass: List[RefinedMapping] = []
        for path in sorted(grouped):
            intervals = sorted(
                grouped[path],
                key=lambda item: (item.target_start, item.target_end, locus_maximum_key(item, score_field)),
            )
            starts = [item.target_start for item in intervals]
            prefix_max_ends: List[int] = []
            for item in intervals:
                prefix_max_ends.append(max(item.target_end, prefix_max_ends[-1]
                                          if prefix_max_ends else item.target_end))
            winners: List[Tuple[RefinedMapping, RefinedMapping, int]] = []
            selected_ids = set()
            for reference in intervals:
                left = bisect_left(prefix_max_ends, reference.target_start + 1)
                right = bisect_left(starts, reference.target_end)
                neighbors: List[Tuple[RefinedMapping, int]] = []
                for candidate in intervals[left:right]:
                    overlap = (min(reference.target_end, candidate.target_end)
                               - max(reference.target_start, candidate.target_start))
                    if overlap > 0 and max(overlap / reference.target_span_bp,
                                           overlap / candidate.target_span_bp) >= minimum_overlap_fraction:
                        neighbors.append((candidate, overlap))
                # The reference always overlaps itself at 100%.
                selected, overlap = min(neighbors, key=lambda pair: locus_maximum_key(pair[0], score_field))
                selected_ids.add(selected.mapping_id)
                winners.append((reference, selected, overlap))
            next_pass.extend(item for item in intervals if item.mapping_id in selected_ids)
            for reference, selected, overlap in winners:
                decisions.append(CrossBundleLocusDecision(
                    pass_number=pass_number,
                    reference=reference,
                    selected=selected,
                    reference_retained=reference.mapping_id in selected_ids,
                    overlap_bp=overlap,
                    overlap_fraction_reference=overlap / reference.target_span_bp,
                    overlap_fraction_selected=overlap / selected.target_span_bp,
                    score_field=score_field,
                ))
        current = next_pass
    return current, decisions


def assign_mapping_ids(mappings: Sequence[RefinedMapping]) -> List[RefinedMapping]:
    ordered = sorted(
        mappings,
        key=lambda mapping: (
            mapping.target_name,
            mapping.target_start,
            mapping.target_end,
            mapping.bundle_name,
            mapping.raw_query_name,
            mapping.strand,
        ),
    )
    answer: List[RefinedMapping] = []
    for index, mapping in enumerate(ordered, 1):
        values = dict(mapping.__dict__)
        values["mapping_id"] = f"M{index:07d}"
        answer.append(RefinedMapping(**values))
    return answer


def assess_novel_coverage(
    mappings: Sequence[RefinedMapping],
    old_coverage_by_target: Mapping[str, Tuple[Interval, ...]],
    novel_fraction: float,
) -> List[NovelAssessment]:
    answer: List[NovelAssessment] = []
    for mapping in mappings:
        old_intervals = old_coverage_by_target.get(mapping.target_name, ())
        novel_intervals = subtract_intervals(mapping.target_footprint, old_intervals)
        novel_bp = interval_total(novel_intervals)
        old_covered_bp = mapping.target_footprint_bp - novel_bp
        # Integer required length makes the inclusive >= rule explicit.  Thus a
        # bundle with exactly 25% newly covered sequence passes at 0.25.
        required_novel_bp = math.ceil(mapping.query_length * novel_fraction)
        answer.append(
            NovelAssessment(
                mapping=mapping,
                old_covered_bp=old_covered_bp,
                novel_intervals=novel_intervals,
                novel_bp=novel_bp,
                required_novel_bp=required_novel_bp,
                novel_fraction_of_bundle=novel_bp / mapping.query_length,
                is_novel=novel_bp >= required_novel_bp,
            )
        )
    return answer


def score_from_identity(identity: float) -> int:
    return max(0, min(1000, round(identity * 1000)))


def format_float(value: float) -> str:
    return f"{value:.10g}"


def write_refined_tsv(path: Path, mappings: Sequence[RefinedMapping]) -> None:
    columns = (
        "mapping_id",
        "bundle_name",
        "raw_consensus_name",
        "bundle_length",
        "bundle_start",
        "bundle_end",
        "bundle_aligned_bp",
        "bundle_coverage",
        "strand",
        "population_path",
        "genome_start",
        "genome_end",
        "genome_footprint_bp",
        "genome_span_bp",
        "genome_gap_bp",
        "matching_bases",
        "alignment_block_length",
        "identity",
        "min_mapq",
        "max_mapq",
        "block_count",
        "source_line_numbers",
        "genome_footprint_intervals",
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for mapping in mappings:
            handle.write(
                "\t".join(
                    (
                        mapping.mapping_id,
                        mapping.bundle_name,
                        mapping.raw_query_name,
                        str(mapping.query_length),
                        str(mapping.query_start),
                        str(mapping.query_end),
                        str(mapping.query_aligned_bp),
                        format_float(mapping.query_coverage),
                        mapping.strand,
                        mapping.target_name,
                        str(mapping.target_start),
                        str(mapping.target_end),
                        str(mapping.target_footprint_bp),
                        str(mapping.target_span_bp),
                        str(mapping.target_gap_bp),
                        str(mapping.matching_bases),
                        str(mapping.alignment_block_length),
                        format_float(mapping.identity),
                        str(mapping.min_mapq),
                        str(mapping.max_mapq),
                        str(mapping.block_count),
                        ",".join(map(str, mapping.source_line_numbers)),
                        ",".join(f"{start}-{end}" for start, end in mapping.target_footprint),
                    )
                )
                + "\n"
            )


def write_refined_bed(path: Path, mappings: Sequence[RefinedMapping]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# population_path\tstart\tend\tbundle|mapping_id|raw_consensus\t"
            "score(identity*1000)\tstrand\n"
        )
        for mapping in mappings:
            name = f"{mapping.bundle_name}|{mapping.mapping_id}|raw={mapping.raw_query_name}"
            handle.write(
                f"{mapping.target_name}\t{mapping.target_start}\t{mapping.target_end}\t"
                f"{name}\t{score_from_identity(mapping.identity)}\t{mapping.strand}\n"
            )


def project_pgr_bundle_positions(
    mapping: RefinedMapping, bundle_size: int
) -> Tuple[int, int]:
    """Project bundle bp coordinates onto synthetic/retained PGR coordinates.

    Native PGR bundle positions are inclusive minimizer indices and therefore
    cannot be reconstructed from a consensus sequence.  Scaling the observed
    bundle interval preserves partial-vs-full mapping information while making
    no claim that the projected positions are original minimizer coordinates.
    """

    start_index = (mapping.query_start * bundle_size) // mapping.query_length
    end_exclusive = (
        mapping.query_end * bundle_size + mapping.query_length - 1
    ) // mapping.query_length
    start_index = min(bundle_size - 1, max(0, start_index))
    end_index = min(bundle_size - 1, max(start_index, end_exclusive - 1))
    if mapping.strand == "+":
        return start_index, end_index
    return end_index, start_index


def write_refined_pgr_bed(
    path: Path,
    mappings: Sequence[RefinedMapping],
    assignments: Mapping[str, BundleClassAssignment],
) -> None:
    """Write a minimal native-schema four-column PGR-like BED."""

    occurrence_counts = Counter(
        (mapping.target_name, mapping.bundle_name) for mapping in mappings
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# format: path_id\tstart\tend\t"
            "bundle_id:bundle_size:direction:bundle_pos_start:bundle_pos_end:R|U\n"
        )
        for mapping in mappings:
            assignment = assignments[mapping.bundle_name]
            bundle_start, bundle_end = project_pgr_bundle_positions(
                mapping, assignment.bundle_size
            )
            direction = "0" if mapping.strand == "+" else "1"
            repeat_flag = (
                "R"
                if occurrence_counts[(mapping.target_name, mapping.bundle_name)] > 1
                else "U"
            )
            annotation = ":".join(
                (
                    assignment.numeric_bundle_id,
                    str(assignment.bundle_size),
                    direction,
                    str(bundle_start),
                    str(bundle_end),
                    repeat_flag,
                )
            )
            handle.write(
                f"{mapping.target_name}\t{mapping.target_start}\t"
                f"{mapping.target_end}\t{annotation}\n"
            )


def write_bundle_class_id_map(
    path: Path, assignments: Mapping[str, BundleClassAssignment]
) -> None:
    columns = (
        "bundle_class",
        "numeric_bundle_id",
        "id_source",
        "bundle_size",
        "bundle_size_source",
        "raw_consensus_names",
        "refined_mapping_count",
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for assignment in sorted(
            assignments.values(),
            key=lambda item: (int(item.numeric_bundle_id), natural_name_key(item.bundle_name)),
        ):
            handle.write(
                "\t".join(
                    (
                        assignment.bundle_name,
                        assignment.numeric_bundle_id,
                        assignment.id_source,
                        str(assignment.bundle_size),
                        assignment.bundle_size_source,
                        ",".join(assignment.raw_query_names),
                        str(assignment.mapping_count),
                    )
                )
                + "\n"
            )


def write_novel_tsv(path: Path, assessments: Sequence[NovelAssessment]) -> None:
    columns = (
        "mapping_id",
        "bundle_name",
        "raw_consensus_name",
        "bundle_length",
        "population_path",
        "genome_start",
        "genome_end",
        "strand",
        "genome_footprint_bp",
        "old_covered_bp",
        "novel_bp",
        "novel_fraction_of_bundle",
        "required_novel_bp",
        "is_novel",
        "novel_intervals",
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for item in assessments:
            mapping = item.mapping
            handle.write(
                "\t".join(
                    (
                        mapping.mapping_id,
                        mapping.bundle_name,
                        mapping.raw_query_name,
                        str(mapping.query_length),
                        mapping.target_name,
                        str(mapping.target_start),
                        str(mapping.target_end),
                        mapping.strand,
                        str(mapping.target_footprint_bp),
                        str(item.old_covered_bp),
                        str(item.novel_bp),
                        format_float(item.novel_fraction_of_bundle),
                        str(item.required_novel_bp),
                        "1" if item.is_novel else "0",
                        ",".join(f"{start}-{end}" for start, end in item.novel_intervals),
                    )
                )
                + "\n"
            )


def write_novel_bed(path: Path, assessments: Sequence[NovelAssessment]) -> None:
    """Write actual uncovered footprint pieces for mappings passing the threshold."""

    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            "# population_path\tstart\tend\tbundle|mapping_id|novel_part\t"
            "score(identity*1000)\tstrand\n"
        )
        for item in assessments:
            if not item.is_novel:
                continue
            mapping = item.mapping
            for part_index, (start, end) in enumerate(item.novel_intervals, 1):
                name = f"{mapping.bundle_name}|{mapping.mapping_id}|novel_part={part_index}"
                handle.write(
                    f"{mapping.target_name}\t{start}\t{end}\t{name}\t"
                    f"{score_from_identity(mapping.identity)}\t{mapping.strand}\n"
                )


def write_duplicate_tsv(path: Path, duplicates: Sequence[DuplicateMapping]) -> None:
    columns = (
        "dropped_temporary_mapping_id",
        "retained_temporary_mapping_id",
        "bundle_name",
        "population_path",
        "dropped_raw_consensus_name",
        "retained_raw_consensus_name",
        "overlap_bp",
        "overlap_fraction_dropped",
        "overlap_fraction_retained",
        "reason",
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for item in duplicates:
            handle.write(
                "\t".join(
                    (
                        item.dropped.mapping_id,
                        item.retained.mapping_id,
                        item.dropped.bundle_name,
                        item.dropped.target_name,
                        item.dropped.raw_query_name,
                        item.retained.raw_query_name,
                        str(item.overlap_bp),
                        format_float(item.overlap_fraction_dropped),
                        format_float(item.overlap_fraction_retained),
                        "same_normalized_bundle_and_target_locus",
                    )
                )
                + "\n"
            )


def write_cross_bundle_locus_tsv(
    path: Path, decisions: Sequence[CrossBundleLocusDecision],
) -> None:
    columns = (
        "pass", "decision", "reference_temporary_mapping_id", "selected_temporary_mapping_id",
        "population_path", "reference_start", "reference_end", "selected_start", "selected_end",
        "reference_bundle_name", "selected_bundle_name", "reference_raw_consensus_name",
        "selected_raw_consensus_name", "reference_strand", "selected_strand",
        "reference_source_paf_lines", "selected_source_paf_lines", "reference_matching_bases",
        "selected_matching_bases", "reference_identity", "selected_identity",
        "overlap_bp", "overlap_fraction_reference", "overlap_fraction_selected", "score_field",
        "reference_score", "selected_score",
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for item in decisions:
            reference, selected = item.reference, item.selected
            values = (
                item.pass_number, "retained" if item.reference_retained else "dropped",
                reference.mapping_id, selected.mapping_id, reference.target_name,
                reference.target_start, reference.target_end, selected.target_start, selected.target_end,
                reference.bundle_name, selected.bundle_name, reference.raw_query_name,
                selected.raw_query_name, reference.strand, selected.strand,
                ",".join(map(str, reference.source_line_numbers)),
                ",".join(map(str, selected.source_line_numbers)), reference.matching_bases,
                selected.matching_bases, format_float(reference.identity), format_float(selected.identity),
                item.overlap_bp, format_float(item.overlap_fraction_reference),
                format_float(item.overlap_fraction_selected), item.score_field,
                format_float(getattr(reference, item.score_field)),
                format_float(getattr(selected, item.score_field)),
            )
            handle.write("\t".join(map(str, values)) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Chain adjacent minimap2 bundle-to-genome PAF records, filter mappings, "
            "write a PGR-like remapping, and optionally report sequence absent "
            "from a previous merged PGR BED when that baseline is supplied."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "paf",
        type=Path,
        help="historical refine PAF: raw query=population path, raw target=bundle",
    )
    parser.add_argument(
        "merged_bed",
        type=Path,
        nargs="?",
        help="optional previous merged BED; omit to refine without coverage comparison",
    )
    parser.add_argument("-o", "--output-dir", type=Path, required=True)
    parser.add_argument(
        "--bundle-metadata-bed",
        type=Path,
        help=(
            "PGR BED supplying bundle IDs/sizes, independently of old coverage; "
            "defaults to merged_bed; required when merged_bed is omitted. "
            "Use the original PGR BED when remapping raw bundles"
        ),
    )
    parser.add_argument(
        "--min-bundle-coverage",
        "--min-query-coverage",
        dest="min_bundle_coverage",
        type=float,
        default=0.60,
        help=(
            "minimum union of chained bundle intervals divided by consensus length; "
            "--min-query-coverage is a compatibility alias"
        ),
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.0,
        help="minimum PAF matching_bases/alignment_block_length per input block",
    )
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument(
        "--max-genome-gap",
        "--max-query-gap",
        dest="max_genome_gap",
        type=int,
        default=1_000,
        help=(
            "largest population-genome gap (and overlap magnitude) while chaining; "
            "overlaps must also satisfy --max-overlap-fraction; "
            "raw-PAF --max-query-gap is a compatibility alias"
        ),
    )
    parser.add_argument(
        "--max-bundle-gap",
        "--max-target-gap",
        dest="max_bundle_gap",
        type=int,
        default=1_000,
        help=(
            "largest strand-oriented bundle gap (and overlap magnitude) while chaining; "
            "overlaps must also satisfy --max-overlap-fraction; "
            "raw-PAF --max-target-gap is a compatibility alias"
        ),
    )
    parser.add_argument(
        "--max-overlap-fraction",
        type=float,
        default=0.10,
        help=(
            "maximum pairwise overlap magnitude divided by original consensus length, "
            "checked separately on both genome and bundle axes"
        ),
    )
    parser.add_argument(
        "--fixed-bundle-min-coverage",
        type=float,
        default=0.80,
        help=(
            "minimum individual PAF bundle coverage to fix an occurrence as a "
            "singleton mapping excluded from chaining"
        ),
    )
    parser.add_argument(
        "--max-fixed-overlap-fraction",
        type=float,
        default=0.20,
        help=(
            "reject a partial-chain extension if its whole genome outer span overlaps "
            "any different-class fixed occurrence by more than this fraction of that "
            "fixed bundle's original consensus length"
        ),
    )
    parser.add_argument(
        "--max-fixed-overlap-bp",
        type=int,
        default=10_000,
        help=(
            "also reject when the overlap with any one different-class fixed "
            "occurrence strictly exceeds this many bp"
        ),
    )
    parser.add_argument(
        "--novel-fraction",
        type=float,
        default=0.25,
        help="minimum uncovered genome-footprint bp divided by bundle length; inclusive",
    )
    parser.add_argument(
        "--same-locus-overlap",
        type=float,
        default=0.80,
        help=(
            "deduplicate normalized bundle mappings when target-footprint overlap covers "
            "this fraction of either mapping"
        ),
    )
    parser.add_argument(
        "--cross-bundle-locus-dedup",
        action="store_true",
        help=(
            "require post-chaining bundle coverage then select local genome-locus maxima "
            "across all bundle classes; intended for final remapping, not ranked addition masks"
        ),
    )
    parser.add_argument(
        "--filter-final-bundles", action="store_true",
        help=(
            "prune final classes only when deletion merges presence/absence haplotypes, "
            "with full-FASTA sequence-coverage and signed-haplotype protection"
        ),
    )
    parser.add_argument("--genome-fasta", type=Path, help="population FASTA for final bundle filtering; includes unmapped sequences")
    parser.add_argument("--final-min-sequence-coverage", type=float, default=0.95)
    parser.add_argument("--final-max-low-coverage-fraction", type=float, default=0.05)
    parser.add_argument("--final-max-bundle-bp", type=int, default=500, help="only classes with median consensus length at most this many bp can be deleted")
    parser.add_argument(
        "--locus-dedup-min-bundle-coverage", type=float, default=0.70,
        help="minimum union coverage of a completed bundle chain before locus deduplication",
    )
    parser.add_argument(
        "--locus-dedup-overlap", type=float, default=0.80,
        help="inclusive fraction of either original genome BED span required for local competition",
    )
    parser.add_argument(
        "--locus-dedup-passes", type=int, default=2,
        help="number of local maximum passes: first matching bp, subsequent passes identity",
    )
    parser.add_argument(
        "--bundle-alias",
        action="append",
        type=parse_alias,
        default=[],
        metavar="OLD=NEW",
        help=(
            "exact bundle-name replacement used for class barriers and after chaining; repeatable. "
            "Built-ins are dup1_1=dup1 and dup1_2=dup1"
        ),
    )
    return parser


def validate_arguments(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    for option in (
        "min_bundle_coverage", "min_identity", "novel_fraction", "same_locus_overlap",
        "max_overlap_fraction", "fixed_bundle_min_coverage", "max_fixed_overlap_fraction",
        "locus_dedup_min_bundle_coverage",
        "final_min_sequence_coverage", "final_max_low_coverage_fraction",
    ):
        value = getattr(args, option)
        if not 0.0 <= value <= 1.0:
            parser.error(f"--{option.replace('_', '-')} must be in [0, 1], got {value}")
    if args.min_mapq < 0:
        parser.error("--min-mapq must be non-negative")
    if args.max_genome_gap < 0:
        parser.error("--max-genome-gap must be non-negative")
    if args.max_bundle_gap < 0:
        parser.error("--max-bundle-gap must be non-negative")
    if args.max_fixed_overlap_bp < 0:
        parser.error("--max-fixed-overlap-bp must be non-negative")
    if args.final_max_bundle_bp < 0:
        parser.error("--final-max-bundle-bp must be non-negative")
    if not 0.0 < args.locus_dedup_overlap <= 1.0:
        parser.error("--locus-dedup-overlap must be in (0, 1]")
    if args.locus_dedup_passes < 1:
        parser.error("--locus-dedup-passes must be at least 1")
    if not args.paf.is_file():
        parser.error(f"PAF does not exist or is not a file: {args.paf}")
    if args.merged_bed is not None and not args.merged_bed.is_file():
        parser.error(f"merged BED does not exist or is not a file: {args.merged_bed}")
    if args.merged_bed is None and args.bundle_metadata_bed is None:
        parser.error("--bundle-metadata-bed is required when merged_bed is omitted")
    if args.bundle_metadata_bed is not None and not args.bundle_metadata_bed.is_file():
        parser.error(f"bundle metadata BED does not exist or is not a file: {args.bundle_metadata_bed}")
    if args.filter_final_bundles and (args.genome_fasta is None or not args.genome_fasta.is_file()):
        parser.error("--filter-final-bundles requires an existing --genome-fasta")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_arguments(parser, args)

    aliases: Dict[str, str] = {"dup1_1": "dup1", "dup1_2": "dup1"}
    aliases.update(dict(args.bundle_alias))
    # Resolve all aliases now, so cycles are detected even if a name is absent
    # from the current PAF.
    for alias_name in aliases:
        resolve_alias(alias_name, aliases)

    paf_blocks = read_paf(args.paf)
    compare_coverage = args.merged_bed is not None
    old_coverage = read_merged_bed(args.merged_bed) if compare_coverage else {}
    metadata_bed = args.bundle_metadata_bed or args.merged_bed
    existing_numeric_ids, existing_bundle_metadata = read_existing_bundle_metadata(
        metadata_bed
    )

    paf_target_names = {block.target_name for block in paf_blocks}
    paf_query_names = {block.raw_query_name for block in paf_blocks}
    old_target_names = set(old_coverage)
    target_matches = len(paf_target_names & old_target_names)
    query_matches = len(paf_query_names & old_target_names)
    if paf_blocks and old_target_names and target_matches == 0:
        if query_matches > 0:
            raise ValueError(
                "PAF orientation appears reversed: raw PAF target names, rather than raw "
                "query names, match the merged BED paths. Run minimap2 with bundle.fa first "
                "and population.fa second."
            )
        raise ValueError(
            "no population-path name in the raw PAF query column matches the merged BED; "
            "check that both files use the same path identifiers"
        )

    filter_counts: Counter[str] = Counter()
    eligible_blocks: List[PafBlock] = []
    for block in paf_blocks:
        filter_counts["input_paf_blocks"] += 1
        if block.mapq < args.min_mapq:
            filter_counts["blocks_below_min_mapq"] += 1
            continue
        if block.identity < args.min_identity:
            filter_counts["blocks_below_min_identity"] += 1
            continue
        eligible_blocks.append(block)
    filter_counts["blocks_entering_chaining"] = len(eligible_blocks)
    filter_counts["fixed_paf_blocks"] = sum(
        (block.query_end - block.query_start) / block.query_length
        >= args.fixed_bundle_min_coverage
        for block in eligible_blocks
    )
    filter_counts["partial_paf_blocks_entering_chaining"] = (
        len(eligible_blocks) - filter_counts["fixed_paf_blocks"]
    )

    chains = chain_blocks(
        eligible_blocks, args.max_genome_gap, args.max_bundle_gap,
        max_overlap_fraction=args.max_overlap_fraction,
        fixed_bundle_min_coverage=args.fixed_bundle_min_coverage,
        max_fixed_overlap_fraction=args.max_fixed_overlap_fraction,
        max_fixed_overlap_bp=args.max_fixed_overlap_bp,
        aliases=aliases,
    )
    filter_counts["chains_before_bundle_coverage_filter"] = len(chains)
    chained_mappings: List[RefinedMapping] = []
    for index, chain in enumerate(chains, 1):
        normalized_name = normalize_bundle_name(chain.raw_query_name, aliases)
        mapping = chain_to_mapping(chain, normalized_name, f"T{index:07d}")
        if (args.cross_bundle_locus_dedup
                and mapping.query_coverage < args.locus_dedup_min_bundle_coverage):
            filter_counts["chains_below_locus_dedup_min_bundle_coverage"] += 1
            continue
        if mapping.query_coverage < args.min_bundle_coverage:
            filter_counts["chains_below_min_bundle_coverage"] += 1
            continue
        # Retest the weighted chain identity for clarity if future filtering
        # strategies are changed.  With per-block filtering this normally cannot fail.
        if mapping.identity < args.min_identity:
            filter_counts["chains_below_min_identity"] += 1
            continue
        chained_mappings.append(mapping)
    filter_counts["chains_passing_filters"] = len(chained_mappings)
    filter_counts["coverage_only_mappings_before_locus_dedup"] = len(chained_mappings)

    # A class-local identity-first preselection would change the native awk /
    # bedmap matching-bp winner, especially for unequal-length alleles.  Final
    # locus mode therefore sends every eligible mapping through its two passes;
    # legacy ranked masking keeps the original class-local policy.
    if args.cross_bundle_locus_dedup:
        deduplicated, duplicates = list(chained_mappings), []
    else:
        deduplicated, duplicates = deduplicate_same_locus(
            chained_mappings, args.same_locus_overlap
        )
    filter_counts["same_locus_duplicates_removed"] = len(duplicates)
    locus_decisions: List[CrossBundleLocusDecision] = []
    if args.cross_bundle_locus_dedup:
        filter_counts["cross_bundle_locus_mappings_before_dedup"] = len(deduplicated)
        deduplicated, locus_decisions = deduplicate_cross_bundle_loci(
            deduplicated, args.locus_dedup_overlap, args.locus_dedup_passes,
        )
        filter_counts["cross_bundle_locus_duplicates_removed"] = (
            filter_counts["cross_bundle_locus_mappings_before_dedup"] - len(deduplicated)
        )
        for pass_number in range(1, args.locus_dedup_passes + 1):
            rows = [row for row in locus_decisions if row.pass_number == pass_number]
            filter_counts[f"locus_dedup_pass_{pass_number}_input_mappings"] = len(rows)
            filter_counts[f"locus_dedup_pass_{pass_number}_retained_mappings"] = sum(
                row.reference_retained for row in rows
            )
            filter_counts[f"locus_dedup_pass_{pass_number}_removed_mappings"] = sum(
                not row.reference_retained for row in rows
            )
    mappings = assign_mapping_ids(deduplicated)
    class_assignments = assign_bundle_class_ids(
        mappings, existing_numeric_ids, existing_bundle_metadata
    )
    before_final_mappings = list(mappings)
    before_final_assignments = dict(class_assignments)
    final_filter_result = None
    if args.filter_final_bundles:
        mappings, final_filter_result = filter_final_bundle_classes(
            mappings, class_assignments, read_sequence_lengths(args.genome_fasta),
            args.final_min_sequence_coverage, args.final_max_low_coverage_fraction,
            args.final_max_bundle_bp,
        )
        final_name_counts = Counter(mapping.bundle_name for mapping in mappings)
        final_raw_names: Dict[str, set] = defaultdict(set)
        for mapping in mappings:
            final_raw_names[mapping.bundle_name].add(mapping.raw_query_name)
        class_assignments = {
            name: replace(assignment, mapping_count=final_name_counts[name],
                          raw_query_names=tuple(sorted(final_raw_names[name])))
            for name, assignment in class_assignments.items() if final_name_counts[name]
        }
    filter_counts["refined_mappings"] = len(mappings)
    filter_counts["refined_bundle_classes"] = len(class_assignments)
    filter_counts["classes_reusing_existing_numeric_id"] = sum(
        assignment.id_source == "reused_existing_numeric_id"
        for assignment in class_assignments.values()
    )

    assessments = []
    if compare_coverage:
        assessments = assess_novel_coverage(mappings, old_coverage, args.novel_fraction)
        filter_counts["novel_mappings"] = sum(item.is_novel for item in assessments)
        filter_counts["novel_uncovered_intervals"] = sum(
            len(item.novel_intervals) for item in assessments if item.is_novel
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    refined_tsv = args.output_dir / "refined_bundle_mappings.tsv"
    coverage_only_tsv = args.output_dir / "coverage_filtered_bundle_mappings.tsv"
    refined_bed = args.output_dir / "refined_bundle_mappings.bed"
    refined_pgr_bed = args.output_dir / "refined_bundle_mappings.pgr.bed"
    class_id_map = args.output_dir / "refined_bundle_class_ids.tsv"
    novel_tsv = args.output_dir / "novel_bundle_coverage.tsv"
    novel_bed = args.output_dir / "novel_bundle_mappings.bed"
    duplicate_tsv = args.output_dir / "duplicate_same_locus_mappings.tsv"
    locus_dedup_tsv = args.output_dir / "cross_bundle_locus_dedup.tsv"
    summary_json = args.output_dir / "refinement_summary.json"

    final_filter_outputs = {}
    if final_filter_result is not None:
        before_final_bed = args.output_dir / "refined_bundle_mappings.before_filter.pgr.bed"
        write_refined_pgr_bed(before_final_bed, before_final_mappings, before_final_assignments)
        final_filter_outputs = write_final_filter_audits(args.output_dir, final_filter_result)
        final_filter_outputs["before_filter_pgr_bed"] = str(before_final_bed.resolve())

    if args.cross_bundle_locus_dedup:
        write_refined_tsv(coverage_only_tsv, chained_mappings)
    write_refined_tsv(refined_tsv, mappings)
    write_refined_bed(refined_bed, mappings)
    write_refined_pgr_bed(refined_pgr_bed, mappings, class_assignments)
    write_bundle_class_id_map(class_id_map, class_assignments)
    if compare_coverage:
        write_novel_tsv(novel_tsv, assessments)
        write_novel_bed(novel_bed, assessments)
    write_duplicate_tsv(duplicate_tsv, duplicates)
    write_cross_bundle_locus_tsv(locus_dedup_tsv, locus_decisions)

    summary = {
        "coverage_comparison_enabled": compare_coverage,
        "final_bundle_filter": {
            "enabled": args.filter_final_bundles,
            "min_sequence_coverage": args.final_min_sequence_coverage,
            "max_low_coverage_fraction": args.final_max_low_coverage_fraction,
            "max_bundle_bp": args.final_max_bundle_bp,
            "denominator": "full_population_fasta_sequence_length",
            "coverage": "union of final BED outer intervals per sequence",
            "deleted_bundle_position_protection": False,
            "haplotype_rule": (
                "retain a class when deletion leaves the haplotype partition unchanged; "
                "for changed partitions, any newly collapsed original pair must have "
                "compatible signed order of shared classes; copy-number differences are "
                "protected only for classes occurring in an adjacent tandem run in at "
                "least one original haplotype; placement of a deleted class among repeated "
                "units is ignored while its orientation remains protected"
            ),
            "counts": final_filter_result["counts"] if final_filter_result is not None else {},
            "outputs": final_filter_outputs,
        },
        "inputs": {
            "paf": str(args.paf.resolve()),
            "merged_bed": str(args.merged_bed.resolve()) if compare_coverage else None,
            "bundle_metadata_bed": str(metadata_bed.resolve()),
            "genome_fasta": str(args.genome_fasta.resolve()) if args.genome_fasta is not None else None,
        },
        "parameters": {
            "min_bundle_coverage": args.min_bundle_coverage,
            "min_identity": args.min_identity,
            "min_mapq": args.min_mapq,
            "max_genome_gap": args.max_genome_gap,
            "max_bundle_gap": args.max_bundle_gap,
            "max_overlap_fraction": args.max_overlap_fraction,
            "fixed_bundle_min_coverage": args.fixed_bundle_min_coverage,
            "max_fixed_overlap_fraction": args.max_fixed_overlap_fraction,
            "max_fixed_overlap_bp": args.max_fixed_overlap_bp,
            "novel_fraction": args.novel_fraction if compare_coverage else None,
            "same_locus_overlap": args.same_locus_overlap,
            "cross_bundle_locus_dedup": args.cross_bundle_locus_dedup,
            "locus_dedup_min_bundle_coverage": args.locus_dedup_min_bundle_coverage,
            "locus_dedup_overlap": args.locus_dedup_overlap,
            "locus_dedup_passes": args.locus_dedup_passes,
            "locus_dedup_score_by_pass": [
                "matching_bases" if number == 1 else "identity"
                for number in range(1, args.locus_dedup_passes + 1)
            ],
            "bundle_aliases": aliases,
        },
        "counts": dict(sorted(filter_counts.items())),
        "coordinate_semantics": {
            "bed": "zero-based half-open",
            "raw_paf_query": "population/genome path",
            "raw_paf_target": "bundle consensus",
            "bundle_coverage": "union(chained raw-PAF target intervals) / bundle length",
            "cross_bundle_locus_dedup": (
                "when enabled, filter completed chains by union bundle coverage; "
                "native partial PAF blocks remain eligible to chain; "
                "bypass class-local identity-first deduplication so all eligible mappings compete; "
                "repeat independent local maximum selection (pass 1 matching bp, later passes "
                "identity) on the same population "
                "path without class/name/strand barriers; overlap is inclusive fraction of "
                "either outer genome BED span; retain distinct selected original mappings, "
                "not merged coordinates or transitive connected components; ties prefer "
                "lexical raw names followed by stable coordinates and source lines"
            ),
            "coverage_only_mapping_table": (
                "native PAF MAPQ/identity, chaining, and configured post-chain union bundle "
                "coverage only; written before same-locus/cross-bundle locus deduplication "
                "and before optional final class pruning"
            ),
            "chain_overlap_fraction": (
                "pairwise overlap bp / original bundle consensus length, checked "
                "separately on both axes in addition to absolute gap limits"
            ),
            "chain_fixed_bundle_selection": (
                "eligible individual PAF blocks with bundle coverage at least "
                "fixed_bundle_min_coverage; immutable singleton chains, excluded "
                "from partial chaining; newly merged chains do not become fixed"
            ),
            "chain_fixed_bundle_overlap": (
                "entire proposed partial chain genome outer span intersected with "
                "each different normalized-class fixed occurrence's genome interval "
                "on the same path, without strand restriction; reject if overlap "
                "strictly exceeds max_fixed_overlap_fraction of the fixed original "
                "consensus length OR max_fixed_overlap_bp; independent occurrences "
                "are not summed"
            ),
            "pgr_bundle_positions": (
                "inclusive deterministic projection from bundle bp to retained or "
                "synthetic PGR bundle_size; not original minimizer coordinates"
            ),
            "pgr_direction": "0=PAF '+', 1=PAF '-'",
            "pgr_repeat_flag": (
                "R when a normalized bundle class occurs more than once on one path; "
                "otherwise U"
            ),
        },
        "outputs": {
            "refined_tsv": str(refined_tsv.resolve()),
            "coverage_only_tsv": (
                str(coverage_only_tsv.resolve()) if args.cross_bundle_locus_dedup else None
            ),
            "refined_bed": str(refined_bed.resolve()),
            "refined_pgr_bed": str(refined_pgr_bed.resolve()),
            "bundle_class_id_map": str(class_id_map.resolve()),
            "same_locus_duplicates_tsv": str(duplicate_tsv.resolve()),
            "cross_bundle_locus_dedup_tsv": str(locus_dedup_tsv.resolve()),
        },
    }
    if compare_coverage:
        summary["coordinate_semantics"].update({
            "novel_bp": "union(member genome footprints) minus union(old BED intervals)",
            "chain_genome_gaps_count_as_novel": False,
            "novel_threshold_is_inclusive": True,
        })
        summary["outputs"].update({
            "novel_tsv": str(novel_tsv.resolve()),
            "novel_bed": str(novel_bed.resolve()),
        })
    with summary_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    json.dump(summary, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
