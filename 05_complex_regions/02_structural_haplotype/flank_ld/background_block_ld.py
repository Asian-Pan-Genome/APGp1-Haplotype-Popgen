#!/usr/bin/env python3
"""Distance-conditioned block D' and matched-region background tests.

Pair backgrounds require both windows, but NOT their intervening gap, to be
inside the allowed BED union. By default their center distance may differ by
10% from the tested pair and their SNP counts are not matched.
Region backgrounds require the entire central
gap and both flanks to fit inside that union. Only cross-flank pairs enter
the regional maximum and significant-pair-fraction statistics. Region controls
match available block counts and per-side block haplotype entropy distributions
with optional left/right swap, not the target's rank-pair mask;
central lengths may differ by up to --region-length-tolerance-bp.
VCF/BCF genotypes are retrieved on demand through bcftools and require an index
by default; --vcf-access preload explicitly enables the old full-file scan.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
import random
import sys
import threading
import time
from bisect import bisect_left, bisect_right
from collections import Counter, OrderedDict, defaultdict
from contextlib import nullcontext
from dataclasses import dataclass, replace
from datetime import datetime
from concurrent.futures import Future, ThreadPoolExecutor
from fractions import Fraction
from functools import lru_cache
from itertools import islice
from pathlib import Path
from typing import Callable, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple, TypeVar

try:
    from .fixed_window_block_ld import (
        OUTPUT_FIELDS, FixedBlock, cap_variants_by_entropy, entropy_score, pair_blocks,
        summarize_block_pair,
    )
    from .flanking_block_ld import (
        Region, VCFData, common_haplotype_keys, haplotypes_for_keys,
        read_regions, read_vcf,
    )
    from .indexed_vcf import BCFtoolsVCFSource
except ImportError:  # Direct script execution.
    from fixed_window_block_ld import (  # type: ignore
        OUTPUT_FIELDS, FixedBlock, cap_variants_by_entropy, entropy_score, pair_blocks,
        summarize_block_pair,
    )
    from flanking_block_ld import (  # type: ignore
        Region, VCFData, common_haplotype_keys, haplotypes_for_keys,
        read_regions, read_vcf,
    )
    from indexed_vcf import BCFtoolsVCFSource  # type: ignore


VERSION = "1.17"
ANALYSIS_SEMANTICS_VERSION = "pair_specific_fixed_complete_cohort_v1"
PAIR_SNP_SELECTION_MODE = "pair_specific_fixed_complete_cohort_deferred_to_dprime"
STATISTIC_TIE_TOLERANCE = 1e-12
Mask = Sequence[Tuple[int, int]]
PairBackgroundKey = Tuple[int, int, int]  # Reference distance, counts; (0, 0) means no count matching.
BlockSNPCounts = Mapping[Tuple[str, int], int]
VCF_QC_FIELDS = ("vcf_maf_missing_qc_skipped", "window_missing_qc_status")
T = TypeVar("T")
U = TypeVar("U")


def stable_subseed(seed: int, namespace: str, values: Sequence[object]) -> int:
    """Stable task-specific RNG seed, independent of scheduling and PYTHONHASHSEED."""
    payload = json.dumps([seed, namespace, *values], separators=(",", ":"), ensure_ascii=True)
    return int.from_bytes(hashlib.blake2b(payload.encode("utf-8"), digest_size=8).digest(), "big")


def ordered_thread_map(
    function: Callable[[T], U], items: Iterable[T], workers: int,
) -> Iterator[U]:
    """Bounded parallel map whose yielded order is identical for every worker count."""
    if workers <= 1:
        for item in items:
            yield function(item)
        return
    iterator = iter(items)
    executor = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="block-ld-bg")
    pending: List[Future] = []
    try:
        for _ in range(workers * 2):
            try:
                pending.append(executor.submit(function, next(iterator)))
            except StopIteration:
                break
        while pending:
            future = pending.pop(0)
            yield future.result()
            try:
                pending.append(executor.submit(function, next(iterator)))
            except StopIteration:
                pass
    finally:
        for future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)


def resolve_snp_limits(args: argparse.Namespace) -> Tuple[int, int]:
    """Keep legacy fixed-k commands while allowing an explicit count range."""
    fixed = getattr(args, "snps_per_block", None)
    lower = getattr(args, "min_snps_per_block", None)
    upper = getattr(args, "max_snps_per_block", None)
    if fixed is not None:
        if lower is not None or upper is not None:
            raise ValueError("--snps-per-block cannot be combined with --min-snps-per-block or --max-snps-per-block")
        lower = upper = fixed
    else:
        lower = 10 if lower is None else lower
        upper = 10 if upper is None else upper
    if lower < 1:
        raise ValueError("--min-snps-per-block must be >= 1")
    if upper < lower:
        raise ValueError("--max-snps-per-block must be >= --min-snps-per-block (unspecified bounds default to 10)")
    return lower, upper


@dataclass(frozen=True)
class StartDomain:
    """Inclusive legal integer starts; unlike BED, end IS included."""

    chrom: str
    display_chrom: str
    start: int
    end: int

    @property
    def size(self) -> int:
        return self.end - self.start + 1


def precise_float(value: Optional[float]) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return format(value, ".17g")


@dataclass(frozen=True)
class RegionEntropyProfile:
    """One marginal entropy per available block, in genomic side/rank order."""

    left: Tuple[float, ...]
    right: Tuple[float, ...]
    left_complete_haplotypes: Tuple[int, ...] = ()
    right_complete_haplotypes: Tuple[int, ...] = ()
    left_ranks: Tuple[int, ...] = ()
    right_ranks: Tuple[int, ...] = ()


@dataclass(frozen=True)
class RegionEntropyMatch:
    orientation: str
    left_wasserstein: Optional[float]
    right_wasserstein: Optional[float]
    max_wasserstein: Optional[float]
    direct_max_wasserstein: Optional[float]
    swapped_max_wasserstein: Optional[float]


@dataclass(frozen=True)
class RegionCandidateMatch:
    """Symmetric direct/swapped covariate match for strict or nearest selection."""

    orientation: str
    score: float
    strict: bool
    left_count_error: int
    right_count_error: int
    entropy: RegionEntropyMatch


@dataclass(frozen=True)
class RegionCandidateScreen:
    """D'-free candidate screen retained until final background selection."""

    left_blocks: Tuple[FixedBlock, ...]
    right_blocks: Tuple[FixedBlock, ...]
    eligible_pairs: Tuple[Tuple[int, int], ...]
    entropy: Optional[RegionEntropyProfile]

    @property
    def available_counts(self) -> Tuple[int, int]:
        return (
            len({left for left, _ in self.eligible_pairs}),
            len({right for _, right in self.eligible_pairs}),
        )


ENTROPY_PROFILE_FIELDS = (
    "left_block_haplotype_entropies", "right_block_haplotype_entropies",
    "left_block_entropy_complete_haplotype_counts", "right_block_entropy_complete_haplotype_counts",
    "left_entropy_block_ranks", "right_entropy_block_ranks",
)
ENTROPY_CONTROL_FIELDS = ("region_entropy_matching_enabled", "region_entropy_max_wasserstein")
ENTROPY_MATCH_FIELDS = (
    "matched_left_entropy_wasserstein", "matched_right_entropy_wasserstein",
    "region_entropy_wasserstein", "direct_entropy_max_wasserstein", "swapped_entropy_max_wasserstein",
)
REGION_FALLBACK_CONTROL_FIELDS = (
    "region_match_fallback", "region_fallback_max_score", "background_region_attempt_limit",
)
REGION_FALLBACK_SUMMARY_FIELDS = (
    "region_strict_match_count", "region_fallback_candidate_count",
    "region_fallback_selected_count", "selected_region_match_score_max",
)
REGION_FALLBACK_MATCH_FIELDS = (
    "region_match_selection", "region_match_score", "region_match_strict",
)


def entropy_profile_fields(profile: Optional[RegionEntropyProfile]) -> Dict[str, str]:
    if profile is None:
        return {field: "NA" for field in ENTROPY_PROFILE_FIELDS}
    values = (profile.left, profile.right, profile.left_complete_haplotypes,
              profile.right_complete_haplotypes, profile.left_ranks, profile.right_ranks)
    return {field: json.dumps(value, separators=(",", ":"), allow_nan=False)
            for field, value in zip(ENTROPY_PROFILE_FIELDS, values)}


def wasserstein_distance(first: Sequence[float], second: Sequence[float]) -> float:
    """Exact 1D W1 for equally weighted empirical samples; no SciPy required.

    Integrate the absolute difference between empirical CDFs. Each sample
    has total mass one, including when their block counts differ.
    """
    if not first or not second or any(not math.isfinite(value) for value in (*first, *second)):
        raise ValueError("Wasserstein distance requires nonempty finite samples")
    left, right = sorted(first), sorted(second)
    support = sorted(set(left + right))
    return math.fsum((end - start) * abs(
        bisect_right(left, start) / len(left) - bisect_right(right, start) / len(right)
    ) for start, end in zip(support, support[1:]))


def merge_background_intervals(regions: Sequence[Region]) -> List[Region]:
    merged: List[Region] = []
    for region in sorted(regions, key=lambda item: (item.chrom, item.start, item.end)):
        if merged and merged[-1].chrom == region.chrom and region.start <= merged[-1].end:
            merged[-1] = replace(merged[-1], end=max(merged[-1].end, region.end))
        else:
            merged.append(replace(region, region_id=f"background_{len(merged) + 1}"))
    return merged


def subtract_exclusions(
    background: Sequence[Region], exclusions: Sequence[Region],
) -> List[Region]:
    excluded: Dict[str, List[Region]] = defaultdict(list)
    for region in merge_background_intervals(exclusions):
        excluded[region.chrom].append(region)
    pieces: List[Region] = []
    for source in merge_background_intervals(background):
        cursor = source.start
        for exclusion in excluded[source.chrom]:
            if exclusion.end <= cursor:
                continue
            if exclusion.start >= source.end:
                break
            if exclusion.start > cursor:
                pieces.append(replace(source, start=cursor, end=exclusion.start))
            cursor = max(cursor, exclusion.end)
            if cursor >= source.end:
                break
        if cursor < source.end:
            pieces.append(replace(source, start=cursor))
    return merge_background_intervals(pieces)


def pair_sampling_domains(
    background: Sequence[Region], window_bp: int, distance_bp: int,
) -> List[StartDomain]:
    """Legal first-window starts; second start = first start + distance.

    Equal-width windows have the same start and center distance. Intersect
    the window-start universe with its distance-shifted copy, chromosome by
    chromosome. Separate BED components can contain the two windows.
    """
    if window_bp < 1 or distance_bp < window_bp:
        raise ValueError("pair windows must have positive width and cannot overlap")
    by_chrom: Dict[str, List[Region]] = defaultdict(list)
    for region in merge_background_intervals(background):
        if region.end - region.start >= window_bp:
            by_chrom[region.chrom].append(region)
    domains: List[StartDomain] = []
    for chrom, intervals in sorted(by_chrom.items()):
        starts = [(item.start, item.end - window_bp) for item in intervals]
        shifted = [(start - distance_bp, end - distance_bp) for start, end in starts]
        i = j = 0
        while i < len(starts) and j < len(shifted):
            low = max(starts[i][0], shifted[j][0])
            high = min(starts[i][1], shifted[j][1])
            if low <= high:
                domains.append(StartDomain(chrom, intervals[0].display_chrom, low, high))
            if starts[i][1] < shifted[j][1]:
                i += 1
            else:
                j += 1
    return domains


def pair_distance_bounds(distance_bp: int, tolerance: float, window_bp: int = 1) -> Tuple[int, int]:
    """Inclusive integer distance band; windows may touch but not overlap."""
    if distance_bp < window_bp or window_bp < 1:
        raise ValueError("pair center distance must be >= the positive window width")
    if not math.isfinite(tolerance) or not 0 <= tolerance < 1:
        raise ValueError("--pair-distance-tolerance must be finite and in [0, 1)")
    delta = distance_bp * Fraction(str(tolerance))
    lower, upper = distance_bp - delta, distance_bp + delta
    return max(window_bp, -(-lower.numerator // lower.denominator)), upper.numerator // upper.denominator


def pair_configuration_count(
    first_start: int, first_end: int, second_start: int, second_end: int,
    min_distance: int, max_distance: int,
) -> int:
    """Count integer pairs in a rectangle clipped by a distance band, O(1)."""
    if first_start > first_end or second_start > second_end:
        return 0
    first_n, second_n = first_end - first_start + 1, second_end - second_start + 1
    base = second_start - first_end

    def triangle(value: int) -> int:
        return (value + 1) * (value + 2) // 2 if value >= 0 else 0

    def prefix(distance: int) -> int:
        value = distance - base
        return (triangle(value) - triangle(value - first_n) - triangle(value - second_n)
                + triangle(value - first_n - second_n))

    return prefix(max_distance) - prefix(min_distance - 1)


@dataclass(frozen=True)
class PairDistanceDomain:
    """Legal (first start, second start) configurations; all bounds inclusive."""

    chrom: str
    display_chrom: str
    start: int
    end: int
    second_start: int
    second_end: int
    min_distance: int
    max_distance: int

    @property
    def size(self) -> int:
        return self.count_through(self.end)

    def count_through(self, end: int) -> int:
        return pair_configuration_count(self.start, end, self.second_start, self.second_end,
                                        self.min_distance, self.max_distance)

    def decode(self, index: int) -> Tuple[int, int]:
        """Invert cumulative first-start counts without enumerating bp distances."""
        if not 0 <= index < self.size:
            raise ValueError("pair configuration index is out of range")
        low, high = self.start, self.end
        while low < high:
            middle = (low + high) // 2
            if self.count_through(middle) > index:
                high = middle
            else:
                low = middle + 1
        second_low = max(self.second_start, low + self.min_distance)
        return low, second_low + index - self.count_through(low - 1)


def pair_distance_sampling_domains(
    background: Sequence[Region], window_bp: int, distance_bp: int, tolerance: float,
) -> List[PairDistanceDomain]:
    lower, upper = pair_distance_bounds(distance_bp, tolerance, window_bp)
    by_chrom: Dict[str, List[Region]] = defaultdict(list)
    for region in merge_background_intervals(background):
        if region.end - region.start >= window_bp:
            by_chrom[region.chrom].append(region)
    domains: List[PairDistanceDomain] = []
    for chrom, intervals in sorted(by_chrom.items()):
        starts = [item.start for item in intervals]
        ends = [item.end - window_bp for item in intervals]
        for first, first_end in zip(intervals, ends):
            begin = bisect_left(ends, first.start + lower)
            stop = bisect_right(starts, first_end + upper)
            for second_index in range(begin, stop):
                second_start, second_end = starts[second_index], ends[second_index]
                low = max(first.start, second_start - upper)
                high = min(first_end, second_end - lower)
                if low <= high:
                    domain = PairDistanceDomain(chrom, first.display_chrom, low, high,
                                                second_start, second_end, lower, upper)
                    if domain.size:
                        domains.append(domain)
    return domains


def random_pair_windows(
    domains: Sequence[PairDistanceDomain], limit: int, rng: random.Random,
) -> Iterator[Tuple[PairDistanceDomain, int, int]]:
    """Uniform legal window-pair configurations, without replacement."""
    cumulative: List[int] = []
    total = 0
    for domain in domains:
        total += domain.size
        cumulative.append(total)
    for index in random_indices(total, limit, rng):
        domain_index = bisect_right(cumulative, index)
        previous = cumulative[domain_index - 1] if domain_index else 0
        domain = domains[domain_index]
        first_start, second_start = domain.decode(index - previous)
        yield domain, first_start, second_start


def region_sampling_domains(
    background: Sequence[Region], central_bp: int, flank_span: int,
) -> List[StartDomain]:
    if central_bp < 1 or flank_span < 1:
        raise ValueError("central length and flank span must be positive")
    domains: List[StartDomain] = []
    for region in merge_background_intervals(background):
        low = region.start + flank_span
        high = region.end - central_bp - flank_span
        if low <= high:
            domains.append(StartDomain(region.chrom, region.display_chrom, low, high))
    return domains


def region_length_bounds(central_bp: int, tolerance_bp: int) -> Tuple[int, int]:
    if central_bp < 1 or tolerance_bp < 0:
        raise ValueError("central length must be positive; --region-length-tolerance-bp must be >= 0")
    return max(1, central_bp - tolerance_bp), central_bp + tolerance_bp


@dataclass(frozen=True)
class RegionLengthDomain:
    """Legal (central start, length) configurations inside a BED component."""

    chrom: str
    display_chrom: str
    start: int
    end: int  # Inclusive legal start for the shortest length.
    min_length: int
    max_length: int

    @property
    def start_count(self) -> int:
        return self.end - self.start + 1

    @property
    def size(self) -> int:
        width = self.max_length - self.min_length + 1
        ramp = min(width, self.start_count)
        return (self.start_count - ramp) * width + ramp * (ramp + 1) // 2

    def decode(self, index: int) -> Tuple[int, int]:
        """Invert a plateau + triangular count in O(log length-range) time."""
        if not 0 <= index < self.size:
            raise ValueError("region configuration index is out of range")
        width = self.max_length - self.min_length + 1
        ramp = min(width, self.start_count)
        plateau = self.start_count - ramp
        if index < plateau * width:
            offset, extra_length = divmod(index, width)
            return self.start + offset, self.min_length + extra_length
        index -= plateau * width
        low, high = 0, ramp
        while low + 1 < high:
            middle = (low + high) // 2
            prefix = middle * (2 * ramp - middle + 1) // 2
            if prefix <= index:
                low = middle
            else:
                high = middle
        prefix = low * (2 * ramp - low + 1) // 2
        return self.start + plateau + low, self.min_length + index - prefix


def region_length_sampling_domains(
    background: Sequence[Region], central_bp: int, flank_span: int, tolerance_bp: int,
) -> List[RegionLengthDomain]:
    lower, upper = region_length_bounds(central_bp, tolerance_bp)
    return [RegionLengthDomain(item.chrom, item.display_chrom, item.start, item.end, lower, upper)
            for item in region_sampling_domains(background, lower, flank_span)]


def random_indices(total: int, limit: int, rng: random.Random) -> Iterator[int]:
    """Uniform integer indices without replacement in a bounded-memory shuffle."""
    swaps: Dict[int, int] = {}
    for consumed in range(min(limit, total)):
        remaining = total - consumed
        draw = rng.randrange(remaining)
        index = swaps.get(draw, draw)
        last = remaining - 1
        replacement = swaps.get(last, last)
        if draw != last:
            swaps[draw] = replacement
        swaps.pop(last, None)
        yield index


def random_starts(
    domains: Sequence[StartDomain], limit: int, rng: random.Random,
) -> Iterator[Tuple[StartDomain, int]]:
    """Uniform legal positions without replacement, without a genome-sized list."""
    cumulative: List[int] = []
    total = 0
    for domain in domains:
        total += domain.size
        cumulative.append(total)
    for index in random_indices(total, limit, rng):
        domain_index = bisect_right(cumulative, index)
        previous = cumulative[domain_index - 1] if domain_index else 0
        domain = domains[domain_index]
        yield domain, domain.start + index - previous


def random_region_intervals(
    domains: Sequence[RegionLengthDomain], limit: int, rng: random.Random,
) -> Iterator[Tuple[RegionLengthDomain, int, int]]:
    cumulative: List[int] = []
    total = 0
    for domain in domains:
        total += domain.size
        cumulative.append(total)
    for index in random_indices(total, limit, rng):
        domain_index = bisect_right(cumulative, index)
        previous = cumulative[domain_index - 1] if domain_index else 0
        domain = domains[domain_index]
        start, length = domain.decode(index - previous)
        yield domain, start, length


def empirical_tail_p(value: float, sorted_background: Sequence[float]) -> float:
    if not sorted_background:
        raise ValueError("cannot calculate an empirical p value from an empty background")
    if not math.isfinite(value):
        raise ValueError("empirical-test statistic must be finite")
    # Multiallelic sums can return theoretical ties such as 1 +/- 2e-16.
    # Count numerical ties conservatively rather than turn roundoff into LD
    # enrichment. This tolerance is far below meaningful D' differences.
    greater_equal = len(sorted_background) - bisect_left(
        sorted_background, value - STATISTIC_TIE_TOLERANCE,
    )
    return (greater_equal + 1) / (len(sorted_background) + 1)


class BlockStore:
    """Indexed slicing and bounded reuse of maximum-entropy selected windows."""

    def __init__(self, data: VCFData, args: argparse.Namespace,
                 logger: Optional[RunLogger] = None, source: Optional[BCFtoolsVCFSource] = None):
        self.data = data
        self.args = args
        self.logger = logger
        self.source = source
        self.window_cache_size = getattr(args, "vcf_window_cache_size", 8192)
        if self.window_cache_size < 2 * args.block_count:
            raise ValueError("--vcf-window-cache-size must be >= 2 * --block-count")
        self.window_cache: OrderedDict = OrderedDict()
        self.window_cache_hits = 0
        self._cache_lock = threading.RLock()
        # A region batch can contain more windows than the persistent LRU.
        # Keep its raw candidate blocks in a short-lived mapping so the whole batch
        # is served by one indexed VCF query without evict/refetch churn.
        self._active_window_batch: Optional[Mapping[Tuple[str, int, str], FixedBlock]] = None
        self.min_snps, self.max_snps = resolve_snp_limits(args)
        self.positions = {
            chrom: tuple(variant.bed_pos for variant in variants)
            for chrom, variants in data.variants_by_chrom.items()
        }
        # Instance-owned cache: a module-global method cache would keep whole
        # VCFData objects alive across separate calls to run().
        self.window_entropy = lru_cache(maxsize=self.window_cache_size)(self._window_entropy)

    def _cache_window(self, key: Tuple[str, int, str], block: FixedBlock) -> None:
        with self._cache_lock:
            self.window_cache[key] = block
            self.window_cache.move_to_end(key)
            while len(self.window_cache) > self.window_cache_size:
                self.window_cache.popitem(last=False)

    def prefetch_windows(self, requested: Sequence[Tuple[str, int, str]]) -> None:
        if self.source is None:
            return
        missing = []
        with self._cache_lock:
            for key in dict.fromkeys(requested):
                if key in self.window_cache:
                    self.window_cache.move_to_end(key)
                else:
                    missing.append(key)
        if not missing:
            return
        intervals = [Region(chrom, chrom, start, start + self.args.window_bp, "window", 0)
                     for chrom, start, _ in missing]
        data = self.source.fetch(intervals)
        positions = {chrom: tuple(variant.bed_pos for variant in variants)
                     for chrom, variants in data.variants_by_chrom.items()}
        for key in missing:
            self._cache_window(key, self._window(*key, local_data=data, local_positions=positions))
        # Only per-window candidate SNPs survive in the block cache; the larger
        # merged fetched data object is released rather than accumulating batches.

    def load_window_batch(
        self, requested: Sequence[Tuple[str, int, str]],
    ) -> Dict[Tuple[str, int, str], FixedBlock]:
        """Build many random windows with at most one indexed query.

        Newly built blocks are returned in a batch-local mapping instead of
        flooding the persistent LRU. Existing cached blocks are reused.
        """
        unique = list(dict.fromkeys(requested))
        blocks: Dict[Tuple[str, int, str], FixedBlock] = {}
        missing: List[Tuple[str, int, str]] = []
        with self._cache_lock:
            for key in unique:
                block = self.window_cache.get(key)
                if block is None:
                    missing.append(key)
                else:
                    self.window_cache_hits += 1
                    self.window_cache.move_to_end(key)
                    blocks[key] = block
        if not missing:
            return blocks
        if self.source is None:
            for key in missing:
                blocks[key] = self._window(*key)
            return blocks
        intervals = [Region(chrom, chrom, start, start + self.args.window_bp, "window_batch", 0)
                     for chrom, start, _ in missing]
        data = self.source.fetch(intervals)
        positions = {chrom: tuple(variant.bed_pos for variant in variants)
                     for chrom, variants in data.variants_by_chrom.items()}
        for key in missing:
            blocks[key] = self._window(*key, local_data=data, local_positions=positions)
        return blocks

    def window(self, chrom: str, start: int, side: str) -> FixedBlock:
        key = (chrom, start, side)
        with self._cache_lock:
            if self._active_window_batch is not None:
                block = self._active_window_batch.get(key)
                if block is not None:
                    return block
            block = self.window_cache.get(key)
            if block is not None:
                self.window_cache_hits += 1
                self.window_cache.move_to_end(key)
                return block
        if block is None:
            if self.source is None:
                self._cache_window(key, self._window(chrom, start, side))
            else:
                self.prefetch_windows([key])
        # Another thread can evict a freshly prefetched entry only when the
        # configured cache is smaller than concurrent working sets. Refetching
        # is safe and bounded; the final lookup/copy remains atomic.
        with self._cache_lock:
            block = self.window_cache.get(key)
        if block is None:
            if self.source is None:
                block = self._window(chrom, start, side)
            else:
                intervals = [Region(chrom, chrom, start, start + self.args.window_bp, "window", 0)]
                local_data = self.source.fetch(intervals)
                local_positions = {name: tuple(item.bed_pos for item in variants)
                                   for name, variants in local_data.variants_by_chrom.items()}
                block = self._window(chrom, start, side, local_data, local_positions)
            self._cache_window(key, block)
        return block

    def _window_entropy(self, chrom: str, start: int, side: str) -> Tuple[float, int]:
        # Unlike pair-row entropy, this never restricts to a partner's calls.
        # The cached block keeps all QC-retained candidates; marginal region
        # matching performs its own block-only entropy cap lazily.
        block = self.window(chrom, start, side)
        selected, _method = cap_variants_by_entropy(
            block.variants, self.max_snps, self.args.max_exact_entropy_combinations,
        )
        return entropy_score(selected)

    def region_entropy_profile(
        self, region: Region, rows: Sequence[Mapping[str, str]],
    ) -> RegionEntropyProfile:
        by_side: Dict[str, Dict[int, int]] = {"left": {}, "right": {}}
        for row in rows:
            if valid_row(row):
                by_side["left"][int(row["block1_rank"])] = int(row["left_window_start"])
                by_side["right"][int(row["block2_rank"])] = int(row["right_window_start"])
        return self._region_entropy_profile_from_starts(region, by_side)

    def _region_entropy_profile_from_starts(
        self, region: Region, by_side: Mapping[str, Mapping[int, int]],
    ) -> RegionEntropyProfile:
        values: Dict[str, Tuple[float, ...]] = {}
        counts: Dict[str, Tuple[int, ...]] = {}
        ranks: Dict[str, Tuple[int, ...]] = {}
        for side in ("left", "right"):
            ranks[side] = tuple(sorted(by_side[side]))
            scores = [self.window_entropy(region.chrom, by_side[side][rank], side)
                      for rank in ranks[side]]
            if any(not math.isfinite(entropy) or entropy < 0 for entropy, _ in scores):
                raise ValueError("background region has invalid marginal block entropy")
            values[side] = tuple(entropy for entropy, _ in scores)
            counts[side] = tuple(count for _, count in scores)
        return RegionEntropyProfile(values["left"], values["right"], counts["left"], counts["right"],
                                    ranks["left"], ranks["right"])

    def _window(self, chrom: str, start: int, side: str,
                local_data: Optional[VCFData] = None,
                local_positions: Optional[Mapping[str, Sequence[int]]] = None) -> FixedBlock:
        end = start + self.args.window_bp
        data = self.data if local_data is None else local_data
        positions = (self.positions if local_positions is None else local_positions).get(chrom, ())
        variants = data.variants_by_chrom.get(chrom, ())
        candidates = variants[bisect_left(positions, start):bisect_left(positions, end)]
        if data.maf_missing_qc_skipped:
            # Counts before external filtering are unknown. These placeholders
            # are never used for filtering and are rendered as NA in output.
            pre_count = missed_count = 0
            fraction, passes = 0.0, True
        else:
            pre = data.pre_missing_snp_positions_by_chrom.get(chrom, ())
            missed = data.missing_filtered_snp_positions_by_chrom.get(chrom, ())
            pre_count = bisect_left(pre, end) - bisect_left(pre, start)
            missed_count = bisect_left(missed, end) - bisect_left(missed, start)
            fraction = missed_count / pre_count if pre_count else 0.0
            passes = fraction <= self.args.max_window_missing_filtered_fraction
        if passes:
            ordered = tuple(sorted(candidates, key=lambda item: (
                -item.pos if side == "left" else item.pos, item.variant_id,
            )))
            if self.logger is not None:
                self.logger.periodic("block_work", f"WORK WINDOW chrom={chrom} start={start} end={end} "
                                     f"side={side} eligible_snps={len(candidates)} "
                                     f"entropy_cap={self.max_snps} "
                                     "selection=deferred_to_dprime")
            selected, method = ordered, "deferred_pair_specific_fixed_cohort"
        else:
            selected, method = (), "not_selected_failed_missing_qc"
        return FixedBlock(
            side, 1, start, end, pre_count, missed_count, fraction, passes,
            len(candidates), method, selected,
        )

    def region_blocks(
        self, region: Region,
    ) -> Tuple[Tuple[FixedBlock, ...], Tuple[FixedBlock, ...]]:
        """Build and cache every window around a region, without computing D'."""
        requested = self.region_window_keys(region)
        if not (self._active_window_batch is not None and
                all(key in self._active_window_batch for key in requested)):
            self.prefetch_windows(requested)
        left_keys = requested[:self.args.block_count]
        right_keys = requested[self.args.block_count:]
        left = tuple(replace(self.window(*key), rank=rank)
                     for rank, key in enumerate(left_keys, 1))
        right = tuple(replace(self.window(*key), rank=rank)
                      for rank, key in enumerate(right_keys, 1))
        return left, right

    def final_snp_count(self, block: FixedBlock) -> int:
        """Count retained SNPs without forcing deferred pair-specific selection."""
        return min(len(block.variants), self.max_snps)

    def region_window_keys(self, region: Region) -> Tuple[Tuple[str, int, str], ...]:
        """Deterministic left-then-right block keys for one region."""
        step = self.args.window_step_bp or self.args.window_bp
        return tuple([
            (region.chrom, region.start - (rank - 1) * step - self.args.window_bp, "left")
            for rank in range(1, self.args.block_count + 1)
        ] + [
            (region.chrom, region.end + (rank - 1) * step, "right")
            for rank in range(1, self.args.block_count + 1)
        ])

    def pair_is_eligible_without_dprime(
        self, region: Region, left: FixedBlock, right: FixedBlock,
    ) -> bool:
        """Reproduce pair validity without calculating D' or pair entropy."""
        if (not left.passes_missing_qc or not right.passes_missing_qc or
                len(left.variants) < self.min_snps or len(right.variants) < self.min_snps):
            return False
        keys = common_haplotype_keys(left.variants + right.variants)
        fraction = getattr(self.args, "min_pair_complete_haplotype_fraction", 0.20)
        minimum = max(
            2,
            getattr(self.args, "min_complete_haplotypes", 2),
            math.ceil(2 * self.data.sample_count * fraction),
        )
        if len(keys) < minimum:
            return False
        left_haplotypes = haplotypes_for_keys(left.variants, keys)
        right_haplotypes = haplotypes_for_keys(right.variants, keys)
        return len(set(left_haplotypes)) >= 2 and len(set(right_haplotypes)) >= 2

    def screen_region(self, region: Region) -> RegionCandidateScreen:
        """Screen controls without D'; entropy is conditional on its match flag."""
        left, right = self.region_blocks(region)
        eligible_pairs = tuple(
            (first.rank, second.rank)
            for first, second in pair_blocks(left, right, "all")
            if self.pair_is_eligible_without_dprime(region, first, second)
        )
        entropy = None
        if not getattr(self.args, "skip_region_entropy_matching", False):
            left_by_rank = {block.rank: block for block in left}
            right_by_rank = {block.rank: block for block in right}
            # Lightweight rows preserve the public entropy-profile hook used by
            # callers/tests while carrying no calculated D' or pair entropy.
            eligibility_rows = [
                {
                    "status": "ok", "multiallelic_dprime": "0",
                    "block1_rank": str(left_rank), "block2_rank": str(right_rank),
                    "left_window_start": str(left_by_rank[left_rank].start),
                    "right_window_start": str(right_by_rank[right_rank].start),
                }
                for left_rank, right_rank in eligible_pairs
            ]
            entropy = self.region_entropy_profile(region, eligibility_rows)
        return RegionCandidateScreen(left, right, eligible_pairs, entropy)

    def region_pairs_from_blocks(
        self, region: Region, left: Sequence[FixedBlock], right: Sequence[FixedBlock],
    ) -> List[Dict[str, str]]:
        """Calculate full D' rows only after a candidate has been selected."""
        rows = []
        for first, second in pair_blocks(left, right, "all"):
            rows.append(self.summarize(region, first, second))
            if self.logger is not None:
                self.logger.periodic("block_work", f"WORK BLOCK_D region={region.region_id} "
                                     f"computed_pairs={len(rows)}/{len(left) * len(right)}")
        return rows

    def region_pairs(
        self, region: Region, required_snp_counts: Optional[BlockSNPCounts] = None,
    ) -> List[Dict[str, str]]:
        left, right = self.region_blocks(region)
        if required_snp_counts is not None:
            by_key = {(block.side, block.rank): block for block in left + right}
            for key, expected_count in required_snp_counts.items():
                block = by_key.get(key)
                if block is None or not block.passes_missing_qc or len(block.variants) < self.min_snps:
                    raise ValueError("background region fails target eligible window QC")
                if self.final_snp_count(block) != expected_count:
                    raise ValueError("background region fails target selected SNP count matching")
        return self.region_pairs_from_blocks(region, left, right)

    def summarize(self, region: Region, left: FixedBlock, right: FixedBlock) -> Dict[str, str]:
        row = summarize_block_pair(
            region, self.data.sample_count, left, right,
            min_snps_per_block=self.min_snps,
            permutations=0, rng=random.Random(0), block_count=self.args.block_count,
            float_formatter=precise_float,
            pair_specific_entropy_selection=True,
            max_snps_per_block=self.max_snps,
            max_exact_entropy_combinations=self.args.max_exact_entropy_combinations,
            min_complete_haplotypes=getattr(self.args, "min_complete_haplotypes", 2),
            min_complete_haplotype_fraction=getattr(
                self.args, "min_pair_complete_haplotype_fraction", 0.20,
            ),
        )
        row["vcf_maf_missing_qc_skipped"] = str(self.data.maf_missing_qc_skipped).lower()
        row["window_missing_qc_status"] = ("not_assessable_prefiltered_vcf"
                                           if self.data.maf_missing_qc_skipped else "enabled")
        row["pair_snp_selection_mode"] = PAIR_SNP_SELECTION_MODE
        if self.data.maf_missing_qc_skipped:
            for field in ("left_pre_missing_snp_count", "right_pre_missing_snp_count",
                          "left_missing_filtered_snp_count", "right_missing_filtered_snp_count",
                          "left_missing_filtered_fraction", "right_missing_filtered_fraction",
                          "block1_missed_snp_number", "block2_missed_snp_number", "total_missed_snp_number"):
                row[field] = "NA"
        return row


def row_key(row: Mapping[str, str]) -> Tuple[int, int]:
    return int(row["block1_rank"]), int(row["block2_rank"])


def valid_row(row: Mapping[str, str]) -> bool:
    return row["status"] == "ok" and math.isfinite(float(row["multiallelic_dprime"]))


def pair_background_key(row: Mapping[str, str], match_snp_counts: bool = False) -> PairBackgroundKey:
    # Background output rows belong to their reference-distance bank, not a bank
    # centered on their own (potentially different) sampled distance.
    distance = row.get("distance_stratum_bp", row["block_center_distance_bp"])
    return (
        int(float(distance)),
        int(row["block1_snp_number_pass_threshold"]) if match_snp_counts else 0,
        int(row["block2_snp_number_pass_threshold"]) if match_snp_counts else 0,
    )


def region_snp_counts(pair_rows: Sequence[Mapping[str, str]], mask: Mask) -> Dict[Tuple[str, int], int]:
    # Retained for callers of the earlier strict-rank helper API. The workflow
    # no longer uses this vector to accept/reject background regions.
    by_key = {row_key(row): row for row in pair_rows}
    counts: Dict[Tuple[str, int], int] = {}
    for key in mask:
        row = by_key[key]
        for side, rank, field in (("left", key[0], "block1_snp_number_pass_threshold"),
                                  ("right", key[1], "block2_snp_number_pass_threshold")):
            count = int(row[field])
            block_key = (side, rank)
            if block_key in counts and counts[block_key] != count:
                raise ValueError("inconsistent selected SNP counts for the same ranked block")
            counts[block_key] = count
    return counts


def format_block_snp_counts(counts: BlockSNPCounts) -> str:
    return ",".join(f"{'L' if side == 'left' else 'R'}{rank}={count}"
                    for (side, rank), count in sorted(counts.items()))


def available_block_counts(rows: Sequence[Mapping[str, str]]) -> Tuple[int, int]:
    """A usable block participates in at least one QC-valid cross-flank D' pair."""
    valid = [row for row in rows if valid_row(row)]
    return (len({int(row["block1_rank"]) for row in valid}),
            len({int(row["block2_rank"]) for row in valid}))


def region_valid_pairs(rows: Sequence[Mapping[str, str]]) -> Tuple[Tuple[int, int], ...]:
    return tuple(row_key(row) for row in rows if valid_row(row))


def available_block_count_alignments(
    rows: Sequence[Mapping[str, str]], required_counts: Tuple[int, int],
    tolerance: int = 1, allow_swap: bool = True,
) -> List[Tuple[int, int, bool, str]]:
    """Return every count-valid orientation, sorted by count discrepancy."""
    if tolerance < 0:
        raise ValueError("--region-block-count-tolerance must be >= 0")
    actual = available_block_counts(rows)
    alignments = [("direct", actual)]
    if allow_swap:
        alignments.append(("swapped", (actual[1], actual[0])))
    matches = []
    # Even a large tolerance must not admit a region without any usable pair.
    if min(actual) > 0:
        for orientation, aligned in alignments:
            errors = tuple(abs(count - expected) for count, expected in zip(aligned, required_counts))
            if max(errors) <= tolerance:
                matches.append((max(errors), sum(errors), orientation == "swapped", orientation))
    if not matches:
        raise ValueError("background region fails target available block count matching")
    return sorted(matches)


def check_available_block_counts(
    rows: Sequence[Mapping[str, str]], required_counts: Tuple[int, int],
    tolerance: int = 1, allow_swap: bool = True,
) -> str:
    """Return the best count alignment; never flip genomic coordinates or pairs."""
    return available_block_count_alignments(rows, required_counts, tolerance, allow_swap)[0][-1]


def match_region_entropy(
    rows: Sequence[Mapping[str, str]], required_counts: Tuple[int, int],
    profile: RegionEntropyProfile, target_profile: RegionEntropyProfile,
    tolerance: int = 1, allow_swap: bool = True,
    max_wasserstein: float = 0.05, enabled: bool = True,
) -> RegionEntropyMatch:
    """Counts and entropy must pass in the SAME orientation; no rank matching."""
    if not math.isfinite(max_wasserstein) or max_wasserstein < 0:
        raise ValueError("--region-entropy-max-wasserstein must be finite and >= 0")
    alignments = available_block_count_alignments(rows, required_counts, tolerance, allow_swap)
    if ((len(profile.left), len(profile.right)) != available_block_counts(rows) or
            (len(target_profile.left), len(target_profile.right)) != required_counts):
        raise ValueError("entropy profile does not represent available blocks")
    if any(not math.isfinite(value) or value < 0
           for values in (profile.left, profile.right, target_profile.left, target_profile.right)
           for value in values):
        raise ValueError("region entropy matching requires finite nonnegative marginal entropies")
    distances = {}
    for _, _, _, orientation in alignments:
        left, right = ((profile.left, profile.right) if orientation == "direct"
                       else (profile.right, profile.left))
        dl = wasserstein_distance(target_profile.left, left)
        dr = wasserstein_distance(target_profile.right, right)
        distances[orientation] = (dl, dr, max(dl, dr))
    if enabled:
        best = min(alignments, key=lambda alignment: (distances[alignment[-1]][2], *alignment))
    else:
        # Opt-out restores the previous count-only choice, but reports W1.
        best = alignments[0]
    orientation = best[-1]
    dl, dr, score = distances[orientation]
    if enabled and score > max_wasserstein + STATISTIC_TIE_TOLERANCE:
        raise ValueError("background region fails block haplotype entropy distribution matching")
    return RegionEntropyMatch(
        orientation, dl, dr, score,
        distances["direct"][2] if "direct" in distances else None,
        distances["swapped"][2] if "swapped" in distances else None,
    )


def score_region_candidate(
    rows: Sequence[Mapping[str, str]], required_counts: Tuple[int, int],
    profile: Optional[RegionEntropyProfile], target_profile: Optional[RegionEntropyProfile],
    central_length_delta: int, length_tolerance: int,
    count_tolerance: int = 1, allow_swap: bool = True,
    max_wasserstein: float = 0.05, entropy_enabled: bool = True,
) -> RegionCandidateMatch:
    """Compatibility wrapper for callers that already have full D' rows."""
    return score_region_candidate_counts(
        available_block_counts(rows), required_counts, profile, target_profile,
        central_length_delta, length_tolerance, count_tolerance, allow_swap,
        max_wasserstein, entropy_enabled,
    )


def score_region_candidate_counts(
    actual: Tuple[int, int], required_counts: Tuple[int, int],
    profile: Optional[RegionEntropyProfile], target_profile: Optional[RegionEntropyProfile],
    central_length_delta: int, length_tolerance: int,
    count_tolerance: int = 1, allow_swap: bool = True,
    max_wasserstein: float = 0.05, entropy_enabled: bool = True,
) -> RegionCandidateMatch:
    """Score both orientations and choose one orientation for every covariate.

    The score is the maximum normalized discrepancy. Direct and left/right
    swapped mappings are evaluated symmetrically; count and entropy can never
    silently choose different orientations.
    """
    if min(actual) < 1:
        # A zero-block flank is a block-count mismatch, rather than an
        # unrelated QC failure, and cannot form a cross-flank comparison.
        raise ValueError("background region fails target available block count matching")
    if entropy_enabled:
        if profile is None or target_profile is None:
            raise ValueError("region entropy matching requires entropy profiles")
        if ((len(profile.left), len(profile.right)) != actual or
                (len(target_profile.left), len(target_profile.right)) != required_counts):
            count_alignments = [actual]
            if allow_swap:
                count_alignments.append(actual[::-1])
            if not any(max(abs(value - target) for value, target in zip(aligned, required_counts))
                       <= count_tolerance for aligned in count_alignments):
                raise ValueError("background region fails target available block count matching")
            raise ValueError("entropy profile does not represent available blocks")
        if any(not math.isfinite(value) or value < 0
               for values in (profile.left, profile.right, target_profile.left, target_profile.right)
               for value in values):
            raise ValueError("region entropy matching requires finite nonnegative marginal entropies")

    length_score = abs(central_length_delta) / max(1, length_tolerance)
    candidates = []
    distances: Dict[str, Tuple[float, float, float]] = {}
    for orientation, aligned_counts in (("direct", actual), ("swapped", actual[::-1])):
        if orientation == "swapped" and not allow_swap:
            continue
        count_errors = tuple(abs(value - target) for value, target in zip(aligned_counts, required_counts))
        entropy_score: Optional[float] = None
        if entropy_enabled:
            assert profile is not None and target_profile is not None
            aligned_entropy = ((profile.left, profile.right) if orientation == "direct"
                               else (profile.right, profile.left))
            dl = wasserstein_distance(target_profile.left, aligned_entropy[0])
            dr = wasserstein_distance(target_profile.right, aligned_entropy[1])
            entropy_score = max(dl, dr)
            distances[orientation] = (dl, dr, entropy_score)
        strict = (max(count_errors) <= count_tolerance and
                  (not entropy_enabled or
                   entropy_score <= max_wasserstein + STATISTIC_TIE_TOLERANCE))
        components = [length_score, max(count_errors) / max(1, count_tolerance)]
        if entropy_enabled:
            assert entropy_score is not None
            components.append(entropy_score / max_wasserstein if max_wasserstein > 0
                              else (0.0 if entropy_score <= STATISTIC_TIE_TOLERANCE else math.inf))
        score = max(components)
        candidates.append({
            "orientation": orientation,
            "swapped": orientation == "swapped",
            "count_errors": count_errors,
            "max_count_error": max(count_errors),
            "sum_count_error": sum(count_errors),
            "entropy_score": entropy_score,
            "score": score,
            "strict": strict,
        })

    strict_candidates = [candidate for candidate in candidates if candidate["strict"]]
    if strict_candidates:
        # Preserve the legacy strict-match orientation choice.  Entropy-enabled
        # matching first minimizes W1 and then the count discrepancy; otherwise
        # it minimizes the count discrepancy directly.  The combined normalized
        # score is used to rank only relaxed fallback candidates.
        if entropy_enabled:
            best = min(
                strict_candidates,
                key=lambda candidate: (
                    candidate["entropy_score"],
                    candidate["max_count_error"],
                    candidate["sum_count_error"],
                    candidate["swapped"],
                ),
            )
        else:
            best = min(
                strict_candidates,
                key=lambda candidate: (
                    candidate["max_count_error"],
                    candidate["sum_count_error"],
                    candidate["swapped"],
                ),
            )
    else:
        best = min(
            candidates,
            key=lambda candidate: (
                candidate["score"],
                candidate["max_count_error"],
                candidate["sum_count_error"],
                candidate["swapped"],
            ),
        )
    score = best["score"]
    orientation = best["orientation"]
    count_errors = best["count_errors"]
    if entropy_enabled:
        dl, dr, entropy_score = distances[orientation]
        entropy = RegionEntropyMatch(
            orientation, dl, dr, entropy_score,
            distances["direct"][2] if "direct" in distances else None,
            distances["swapped"][2] if "swapped" in distances else None,
        )
    else:
        entropy = RegionEntropyMatch(orientation, None, None, None, None, None)
    return RegionCandidateMatch(
        orientation, score, best["strict"], count_errors[0], count_errors[1], entropy,
    )


def attach_pair_p(
    rows: Sequence[Dict[str, str]], backgrounds: Mapping[PairBackgroundKey, Sequence[float]],
    alpha: float, match_snp_counts: bool = False, distance_tolerance: float = 0.10, window_bp: int = 1,
) -> None:
    for row in rows:
        row.update(pair_nominal_p="NA", pair_significant="NA", pair_distance_background_n="0",
                   pair_snp_count_stratum="NA", pair_reference_distance_bp="NA",
                   pair_distance_min_bp="NA", pair_distance_max_bp="NA",
                   pair_distance_tolerance=precise_float(distance_tolerance),
                   pair_snp_count_matching="exact_ordered_selected_counts" if match_snp_counts else "none")
        if valid_row(row):
            key = pair_background_key(row, match_snp_counts)
            lower, upper = pair_distance_bounds(key[0], distance_tolerance, window_bp)
            row.update(pair_reference_distance_bp=str(key[0]), pair_distance_min_bp=str(lower),
                       pair_distance_max_bp=str(upper))
            if match_snp_counts:
                row["pair_snp_count_stratum"] = f"{key[1]},{key[2]}"
            reference = backgrounds.get(key)
            if reference:
                p = empirical_tail_p(float(row["multiallelic_dprime"]), reference)
                row.update(
                    pair_nominal_p=precise_float(p), pair_significant=str(p < alpha).lower(),
                    pair_distance_background_n=str(len(reference)),
                )


def union_length(intervals: Sequence[Tuple[int, int]]) -> int:
    length = 0
    previous_end = -1
    for start, end in sorted(intervals):
        length += max(0, end - max(start, previous_end))
        previous_end = max(previous_end, end)
    return length


def analyze_region(
    pair_rows: Sequence[Mapping[str, str]], mask: Mask, alpha: float,
    required_snp_counts: Optional[BlockSNPCounts] = None,
    calculate_breadth: bool = True,
) -> Dict[str, object]:
    """Summarize specified valid pairs; workflow uses each region's own set."""
    if not mask:
        raise ValueError("regional test requires at least one eligible cross-flank pair")
    by_key = {row_key(row): row for row in pair_rows}
    if len(by_key) != len(pair_rows) or len(set(mask)) != len(mask):
        raise ValueError("duplicate rank pairs in regional rows or mask")
    compared = []
    for key in mask:
        row = by_key.get(key)
        if row is None or not valid_row(row):
            raise ValueError("background region fails the target's eligible pair mask")
        if calculate_breadth and row.get("pair_nominal_p", "NA") == "NA":
            raise ValueError("regional breadth requires a calibrated pair background")
        if required_snp_counts is not None:
            for side, rank, field in (("left", key[0], "block1_snp_number_pass_threshold"),
                                      ("right", key[1], "block2_snp_number_pass_threshold")):
                block_key = (side, rank)
                if block_key not in required_snp_counts or int(row[field]) != required_snp_counts[block_key]:
                    raise ValueError("background region fails target selected SNP count matching")
        compared.append(row)
    significant = (
        [row for row in compared if float(row["pair_nominal_p"]) < alpha]
        if calculate_breadth else []
    )
    peak = max(compared, key=lambda row: float(row["multiallelic_dprime"]))
    result: Dict[str, object] = {
        "max_dprime": float(peak["multiallelic_dprime"]),
        "n_valid_pairs": len(compared),
        "n_significant_pairs": len(significant) if calculate_breadth else None,
        "significant_pair_fraction": len(significant) / len(compared) if calculate_breadth else None,
        "peak_block1_label": peak["block1_label"],
        "peak_block2_label": peak["block2_label"],
        "peak_pair_nominal_p": float(peak["pair_nominal_p"]) if calculate_breadth else None,
        "peak_pair_distance_bp": float(peak["block_center_distance_bp"]),
        "peak_left_window_start": int(peak["left_window_start"]),
        "peak_left_window_end": int(peak["left_window_end"]),
        "peak_right_window_start": int(peak["right_window_start"]),
        "peak_right_window_end": int(peak["right_window_end"]),
    }
    for side, rank_field in (("left", "block1_rank"), ("right", "block2_rank")):
        eligible = {row[rank_field] for row in compared}
        supported = {row[rank_field] for row in significant}
        all_intervals = [(int(row[f"{side}_window_start"]), int(row[f"{side}_window_end"])) for row in compared]
        supported_intervals = [(int(row[f"{side}_window_start"]), int(row[f"{side}_window_end"])) for row in significant]
        result[f"n_eligible_{side}_blocks"] = len(eligible)
        result[f"n_supported_{side}_blocks"] = len(supported) if calculate_breadth else None
        result[f"supported_{side}_block_fraction"] = len(supported) / len(eligible) if calculate_breadth else None
        result[f"eligible_{side}_window_coverage_bp"] = union_length(all_intervals)
        result[f"supported_{side}_window_coverage_bp"] = union_length(supported_intervals) if calculate_breadth else None
        result[f"supported_{side}_window_coverage_fraction"] = (
            union_length(supported_intervals) / union_length(all_intervals) if calculate_breadth else None
        )
    return result


PAIR_CONTROL_FIELDS = ("pair_distance_tolerance", "pair_snp_count_matching")
PAIR_SELECTION_FIELDS = ("pair_snp_selection_mode",)
TEST_CONTROL_FIELDS = ("pair_test_enabled", "region_test_enabled")
PAIR_DISTANCE_FIELDS = ("pair_distance_min_bp", "pair_distance_max_bp")
PAIR_EXTRA_FIELDS = ("pair_nominal_p", "pair_significant", "pair_distance_background_n", "pair_in_region_test", "pair_snp_count_stratum", "pair_reference_distance_bp") + PAIR_DISTANCE_FIELDS + PAIR_CONTROL_FIELDS
TARGET_PAIR_FIELDS = OUTPUT_FIELDS[:14] + PAIR_EXTRA_FIELDS + OUTPUT_FIELDS[14:] + PAIR_SELECTION_FIELDS + VCF_QC_FIELDS
BACKGROUND_PAIR_FIELDS = ("background_pair_id", "distance_stratum_bp", "snp_number_stratum_block1", "snp_number_stratum_block2") + OUTPUT_FIELDS + PAIR_DISTANCE_FIELDS + PAIR_CONTROL_FIELDS + PAIR_SELECTION_FIELDS + VCF_QC_FIELDS
STAT_FIELDS = (
    "max_dprime", "significant_pair_fraction", "n_significant_pairs", "n_valid_pairs",
    "n_supported_left_blocks", "n_supported_right_blocks",
    "supported_left_block_fraction", "supported_right_block_fraction",
    "supported_left_window_coverage_fraction", "supported_right_window_coverage_fraction",
    "peak_pair_nominal_p", "peak_block1_label", "peak_block2_label", "peak_pair_distance_bp",
    "peak_left_window_start", "peak_left_window_end", "peak_right_window_start", "peak_right_window_end",
    "n_eligible_left_blocks", "n_eligible_right_blocks",
    "eligible_left_window_coverage_bp", "eligible_right_window_coverage_bp",
    "supported_left_window_coverage_bp", "supported_right_window_coverage_bp",
)
SUMMARY_FIELDS = (
    "region_id", "max_dprime", "region_max_empirical_p",
    "significant_pair_fraction", "region_breadth_empirical_p",
) + tuple(field for field in STAT_FIELDS if field not in ("max_dprime", "significant_pair_fraction")) + (
    "chrom", "region_start", "region_end", "central_length_bp", "flank_span_bp",
    "n_total_cross_pairs", "n_target_qc_failed_pairs", "n_background_regions",
    "background_region_attempts", "background_region_qc_rejections",
    "background_region_snp_count_rejections", "background_region_rejections",
    "pair_significance_threshold", "region_qc_matching",
    "min_snps_per_block", "max_snps_per_block", "region_snp_count_matching", "required_block_snp_counts",
    "pair_snp_selection_mode", "min_pair_complete_haplotype_fraction",
    "required_complete_haplotype_count",
    "background_region_pair_p_enabled", "region_breadth_test_status",
    "region_length_tolerance_bp", "background_min_central_length_bp", "background_max_central_length_bp",
    "required_available_left_blocks", "required_available_right_blocks",
    "background_region_block_count_rejections", "region_pair_distance_matching",
    "region_block_count_tolerance", "region_block_count_allow_swap",
    "background_region_entropy_rejections",
) + REGION_FALLBACK_CONTROL_FIELDS + REGION_FALLBACK_SUMMARY_FIELDS + TEST_CONTROL_FIELDS + PAIR_CONTROL_FIELDS + ENTROPY_CONTROL_FIELDS + ENTROPY_PROFILE_FIELDS + VCF_QC_FIELDS + (
    "status",
)
BACKGROUND_REGION_FIELDS = (
    "target_region_id", "background_region_id", "max_dprime", "significant_pair_fraction",
) + tuple(field for field in STAT_FIELDS if field not in ("max_dprime", "significant_pair_fraction")) + (
    "chrom", "region_start", "region_end", "envelope_start", "envelope_end",
    "n_actual_qc_valid_pairs", "n_compared_pairs", "target_pair_mask",
    "required_block_snp_counts",
    "background_region_pair_p_enabled",
    "pair_snp_selection_mode", "min_pair_complete_haplotype_fraction",
    "required_complete_haplotype_count",
    "central_length_bp", "target_central_length_bp", "central_length_delta_bp",
    "required_available_left_blocks", "required_available_right_blocks", "background_valid_pair_labels",
    "region_block_count_tolerance", "region_block_count_allow_swap", "block_count_match_orientation",
    "matched_left_block_count_delta", "matched_right_block_count_delta",
) + REGION_FALLBACK_CONTROL_FIELDS + REGION_FALLBACK_MATCH_FIELDS + PAIR_CONTROL_FIELDS + ENTROPY_CONTROL_FIELDS + ENTROPY_MATCH_FIELDS + ENTROPY_PROFILE_FIELDS + VCF_QC_FIELDS


def existing_bool(row: Mapping[str, str], field: str, path: Path) -> bool:
    value = row.get(field, "").lower()
    if value not in ("true", "false"):
        raise ValueError(f"existing file {path} has invalid boolean {field}: {row.get(field)!r}")
    return value == "true"


def load_existing_pair_backgrounds(
    path: Path, args: argparse.Namespace,
) -> Tuple[Dict[PairBackgroundKey, List[float]], Dict[str, object]]:
    rows = read_existing_tsv(path, (
        "distance_stratum_bp", "snp_number_stratum_block1", "snp_number_stratum_block2",
        "multiallelic_dprime", "status", "left_window_start", "left_window_end",
        "right_window_start", "right_window_end", "pair_distance_min_bp",
        "pair_distance_max_bp", "pair_distance_tolerance", "pair_snp_count_matching",
        "pair_snp_selection_mode", "min_complete_haplotype_fraction",
    ))
    match_counts = getattr(args, "pair_match_snp_counts", False)
    expected_mode = "exact_ordered_selected_counts" if match_counts else "none"
    expected_selection = PAIR_SNP_SELECTION_MODE
    expected_complete_fraction = getattr(
        args, "min_pair_complete_haplotype_fraction", 0.20,
    )
    requested = args.background_pairs_per_distance
    grouped: Dict[PairBackgroundKey, List[float]] = defaultdict(list)
    for row in rows:
        if row["status"] != "ok":
            raise ValueError(f"existing pair background {path} contains non-ok row")
        if row["pair_snp_count_matching"] != expected_mode:
            raise ValueError(
                f"existing pair background {path} uses pair_snp_count_matching="
                f"{row['pair_snp_count_matching']!r}; current mode is {expected_mode!r}"
            )
        if row["pair_snp_selection_mode"] != expected_selection:
            raise ValueError(
                f"existing pair background {path} uses an incompatible SNP-selection mode"
            )
        complete_fraction = required_finite_float(
            row, "min_complete_haplotype_fraction", path,
        )
        if abs(complete_fraction - expected_complete_fraction) > STATISTIC_TIE_TOLERANCE:
            raise ValueError(
                f"existing pair background {path} uses a different complete-haplotype fraction"
            )
        tolerance = required_finite_float(row, "pair_distance_tolerance", path)
        if abs(tolerance - args.pair_distance_tolerance) > STATISTIC_TIE_TOLERANCE:
            raise ValueError(f"existing pair background {path} has a different distance tolerance")
        try:
            distance = int(row["distance_stratum_bp"])
            left_start, left_end = int(row["left_window_start"]), int(row["left_window_end"])
            right_start, right_end = int(row["right_window_start"]), int(row["right_window_end"])
            if match_counts:
                left_count = int(row["snp_number_stratum_block1"])
                right_count = int(row["snp_number_stratum_block2"])
            else:
                left_count = right_count = 0
        except ValueError as error:
            raise ValueError(f"existing pair background {path} contains invalid integer fields") from error
        if left_end - left_start != args.window_bp or right_end - right_start != args.window_bp:
            raise ValueError(f"existing pair background {path} was calculated with a different window width")
        lower, upper = pair_distance_bounds(distance, args.pair_distance_tolerance, args.window_bp)
        if int(row["pair_distance_min_bp"]) != lower or int(row["pair_distance_max_bp"]) != upper:
            raise ValueError(f"existing pair background {path} has incompatible distance bounds")
        dprime = required_finite_float(row, "multiallelic_dprime", path)
        if dprime < -STATISTIC_TIE_TOLERANCE or dprime > 1 + STATISTIC_TIE_TOLERANCE:
            raise ValueError(f"existing pair background {path} has Dprime outside [0,1]")
        grouped[(distance, left_count, right_count)].append(dprime)
    if not grouped:
        raise ValueError(f"existing pair background file is empty: {path}")
    wrong = {key: len(values) for key, values in grouped.items() if len(values) != requested}
    if wrong:
        raise ValueError(
            f"existing pair background {path} must contain exactly {requested} rows per stratum; "
            f"found {wrong}"
        )
    loaded = {key: sorted(values) for key, values in grouped.items()}
    info = {
        (f"{key[0]}:{key[1]}:{key[2]}" if match_counts else str(key[0])): {
            "distance_bp": key[0],
            "block1_selected_snp_count": key[1] if match_counts else None,
            "block2_selected_snp_count": key[2] if match_counts else None,
            "accepted": len(values), "attempts": None, "source": "reused_existing_tsv",
            "reused": True,
        }
        for key, values in loaded.items()
    }
    return loaded, info


def load_existing_region_backgrounds(
    path: Path, args: argparse.Namespace, targets: Sequence[Region], masks: Sequence[Mask],
    required_counts: Sequence[Tuple[int, int]],
    target_entropy_profiles: Sequence[Optional[RegionEntropyProfile]],
    background_pair_p_enabled: bool, prior_summary_path: Path,
) -> Dict[str, List[Dict[str, str]]]:
    summary_rows = read_existing_tsv(prior_summary_path, (
        "region_id", "chrom", "region_start", "region_end",
        "required_available_left_blocks", "required_available_right_blocks",
        "pair_significance_threshold", *ENTROPY_PROFILE_FIELDS,
    ))
    summaries_by_id = {row["region_id"]: row for row in summary_rows}
    if len(summaries_by_id) != len(summary_rows):
        raise ValueError(f"existing summary {prior_summary_path} has duplicate target IDs")
    for target, mask, counts, entropy_profile in zip(
            targets, masks, required_counts, target_entropy_profiles):
        if not mask:
            continue
        summary = summaries_by_id.get(target.region_id)
        if summary is None:
            raise ValueError(
                f"existing summary {prior_summary_path} lacks target {target.region_id!r}"
            )
        try:
            summary_start, summary_end = int(summary["region_start"]), int(summary["region_end"])
            summary_counts = (
                int(summary["required_available_left_blocks"]),
                int(summary["required_available_right_blocks"]),
            )
        except ValueError as error:
            raise ValueError(f"existing summary {prior_summary_path} has invalid target controls") from error
        if (summary["chrom"], summary_start, summary_end) != (
                target.display_chrom, target.start, target.end):
            raise ValueError(
                f"existing region background target coordinates differ for {target.region_id!r}"
            )
        if summary_counts != counts:
            raise ValueError(
                f"existing region background target block counts differ for {target.region_id!r}"
            )
        alpha = required_finite_float(summary, "pair_significance_threshold", prior_summary_path)
        if abs(alpha - args.pair_significance_threshold) > STATISTIC_TIE_TOLERANCE:
            raise ValueError(f"existing region background uses a different pair significance threshold")
        expected_entropy_fields = entropy_profile_fields(entropy_profile)
        if (entropy_profile is not None and
                any(summary[field] != expected_entropy_fields[field]
                    for field in ENTROPY_PROFILE_FIELDS)):
            raise ValueError(
                f"existing region background target entropy profile differs for {target.region_id!r}"
            )
    rows = read_existing_tsv(path, (
        "target_region_id", "background_region_id", "max_dprime", "significant_pair_fraction",
        "chrom", "region_start", "region_end", "central_length_bp", "target_central_length_bp",
        "central_length_delta_bp", "required_available_left_blocks",
        "required_available_right_blocks", "region_block_count_tolerance",
        "region_block_count_allow_swap", "region_entropy_matching_enabled",
        "region_entropy_max_wasserstein", "background_region_pair_p_enabled",
        "pair_distance_tolerance", "pair_snp_count_matching", "region_match_selection",
        "region_match_score", "region_match_strict",
    ))
    expected_ids = {target.region_id for target, mask in zip(targets, masks) if mask}
    grouped: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    seen_ids = set()
    target_by_id = {target.region_id: target for target in targets}
    counts_by_id = {target.region_id: counts for target, counts in zip(targets, required_counts)}
    expected_pair_mode = (
        "exact_ordered_selected_counts" if getattr(args, "pair_match_snp_counts", False) else "none"
    )
    entropy_enabled = not getattr(args, "skip_region_entropy_matching", False)
    for row in rows:
        target_id = row["target_region_id"]
        if target_id not in expected_ids:
            raise ValueError(f"existing region background {path} has unexpected target {target_id!r}")
        background_id = row["background_region_id"]
        if background_id in seen_ids:
            raise ValueError(f"existing region background {path} has duplicate ID {background_id!r}")
        seen_ids.add(background_id)
        target = target_by_id[target_id]
        expected_counts = counts_by_id[target_id]
        try:
            central_length = int(row["central_length_bp"])
            target_length = int(row["target_central_length_bp"])
            length_delta = int(row["central_length_delta_bp"])
            left_count = int(row["required_available_left_blocks"])
            right_count = int(row["required_available_right_blocks"])
            count_tolerance = int(row["region_block_count_tolerance"])
        except ValueError as error:
            raise ValueError(f"existing region background {path} contains invalid integer controls") from error
        if target_length != target.end - target.start or central_length - target_length != length_delta:
            raise ValueError(f"existing region background {path} has incompatible target/central lengths")
        if abs(length_delta) > getattr(args, "region_length_tolerance_bp", 10_000):
            raise ValueError(f"existing region background {path} exceeds the current length tolerance")
        if (left_count, right_count) != expected_counts:
            raise ValueError(f"existing region background {path} has incompatible target block counts")
        if count_tolerance != getattr(args, "region_block_count_tolerance", 1):
            raise ValueError(f"existing region background {path} has a different block-count tolerance")
        if existing_bool(row, "region_block_count_allow_swap", path) != getattr(
                args, "region_block_count_allow_swap", True):
            raise ValueError(f"existing region background {path} has a different swap setting")
        if existing_bool(row, "region_entropy_matching_enabled", path) != entropy_enabled:
            raise ValueError(f"existing region background {path} has a different entropy matching setting")
        entropy_threshold = required_finite_float(row, "region_entropy_max_wasserstein", path)
        if abs(entropy_threshold - getattr(args, "region_entropy_max_wasserstein", 0.05)) > STATISTIC_TIE_TOLERANCE:
            raise ValueError(f"existing region background {path} has a different entropy threshold")
        if existing_bool(row, "background_region_pair_p_enabled", path) != background_pair_p_enabled:
            raise ValueError(f"existing region background {path} has a different pair-p setting")
        if row["pair_snp_count_matching"] != expected_pair_mode:
            raise ValueError(f"existing region background {path} has a different pair SNP matching mode")
        pair_tolerance = required_finite_float(row, "pair_distance_tolerance", path)
        if abs(pair_tolerance - getattr(args, "pair_distance_tolerance", 0.10)) > STATISTIC_TIE_TOLERANCE:
            raise ValueError(f"existing region background {path} has a different pair distance tolerance")
        dprime = required_finite_float(row, "max_dprime", path)
        if dprime < -STATISTIC_TIE_TOLERANCE or dprime > 1 + STATISTIC_TIE_TOLERANCE:
            raise ValueError(f"existing region background {path} has max Dprime outside [0,1]")
        if background_pair_p_enabled:
            fraction = required_finite_float(row, "significant_pair_fraction", path)
            if not 0 <= fraction <= 1:
                raise ValueError(f"existing region background {path} has invalid significant-pair fraction")
        grouped[target_id].append(row)
    missing = expected_ids - set(grouped)
    wrong = {target_id: len(grouped.get(target_id, ())) for target_id in expected_ids
             if len(grouped.get(target_id, ())) != args.background_regions}
    if missing or wrong:
        raise ValueError(
            f"existing region background {path} must contain exactly {args.background_regions} rows "
            f"for every eligible target; counts={wrong}, missing={sorted(missing)}"
        )
    return dict(grouped)


def serialized(row: Mapping[str, object]) -> Dict[str, object]:
    return {
        key: "NA" if value is None else precise_float(value) if isinstance(value, float) else value
        for key, value in row.items()
    }


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", restval="NA")
        writer.writeheader()
        writer.writerows(serialized(row) for row in rows)


def progress(message: str) -> None:
    print(f"[{datetime.now().astimezone().isoformat(timespec='seconds')}] {message}",
          file=sys.stderr, flush=True)


class RunLogger:
    """Flushed stderr + file messages; logging never uses the sampling RNG."""

    def __init__(self, path: Path, interval_seconds: float = 30.0,
                 every_attempts: int = 1000, append: bool = False):
        self.started = time.monotonic()
        self.interval_seconds = interval_seconds
        self.every_attempts = every_attempts
        self.last_updates: Dict[str, float] = {}
        self.path = path
        self.handle = path.open("a" if append else "x", encoding="utf-8")
        self._lock = threading.RLock()

    def log(self, message: str, level: str = "INFO") -> None:
        with self._lock:
            stamp = datetime.now().astimezone().isoformat(timespec="seconds")
            line = (f"[{stamp}] [elapsed={time.monotonic() - self.started:.1f}s] "
                    f"[{level}] {message}")
            print(line, file=sys.stderr, flush=True)
            print(line, file=self.handle, flush=True)

    def periodic(self, key: str, message: str, force: bool = False) -> None:
        with self._lock:
            now = time.monotonic()
            if force or now - self.last_updates.get(key, self.started) >= self.interval_seconds:
                self.log(message)
                self.last_updates[key] = now

    def vcf_progress(self, counts: Mapping[str, int]) -> None:
        finished = bool(counts["finished"])
        skipped = bool(counts.get("maf_missing_qc_skipped", 0))
        excluded = {"finished", "maf_missing_qc_skipped"}
        if skipped:
            excluded.add("filtered_missing")
        self.periodic("vcf", "VCF " + ("DONE " if finished else "PROGRESS ") +
                      " ".join(f"{key}={value}" for key, value in counts.items() if key not in excluded) +
                      (" filtered_missing=NA maf_missing_qc=skipped" if skipped else ""),
                      force=finished)

    def close(self) -> None:
        with self._lock:
            self.handle.close()


class SamplingProgress:
    """Report completed attempts, including runs with no accepted candidates."""

    def __init__(self, logger: RunLogger, label: str, requested: int,
                 legal_starts: int, attempt_limit: int, legal_unit: str = "legal_starts"):
        self.logger = logger
        self.label = label
        self.requested = requested
        self.attempt_limit = min(legal_starts, attempt_limit)
        self.started = self.last_report = time.monotonic()
        self.last_attempts = 0
        logger.log(f"{label} START requested={requested} {legal_unit}={legal_starts} "
                   f"attempt_limit={self.attempt_limit}")

    def update(self, accepted: int, attempts: int, reasons: Mapping[str, int],
               force: bool = False) -> None:
        now = time.monotonic()
        if not (force or accepted == self.requested or
                attempts - self.last_attempts >= self.logger.every_attempts or
                now - self.last_report >= self.logger.interval_seconds):
            return
        elapsed = max(0.0, now - self.started)
        eta = f"{(self.requested - accepted) * elapsed / accepted:.1f}s" if accepted else "NA"
        status = "DONE" if accepted == self.requested else "EXHAUSTED" if force else "PROGRESS"
        count_rejected = sum(count for reason, count in reasons.items()
                             if reason == "selected_snp_count_mismatch" or "selected SNP count matching" in reason)
        block_rejected = sum(count for reason, count in reasons.items()
                             if "available block count matching" in reason)
        entropy_rejected = sum(count for reason, count in reasons.items()
                               if "block haplotype entropy distribution matching" in reason)
        rejected = sum(reasons.values())
        self.logger.log(
            f"{self.label} {status} accepted={accepted}/{self.requested} "
            f"percent={100 * accepted / self.requested:.1f}% attempts={attempts}/{self.attempt_limit} "
            f"rejected={rejected} qc_rejected={rejected - count_rejected - block_rejected - entropy_rejected} "
            f"snp_count_rejected={count_rejected} block_count_rejected={block_rejected} "
            f"entropy_rejected={entropy_rejected} "
            f"acceptance_rate={100 * accepted / attempts if attempts else 0:.2f}% "
            f"sampling_elapsed={elapsed:.1f}s eta={eta} "
            f"reasons={json.dumps(dict(reasons), sort_keys=True, ensure_ascii=False)}"
        )
        self.last_report = now
        self.last_attempts = attempts


def output_files(args: argparse.Namespace) -> Dict[str, Path]:
    paths = {
        name: Path(f"{args.output_prefix}.{suffix}") for name, suffix in (
            ("target_pairs", "target_pairs.tsv"), ("summary", "summary.tsv"),
            ("background_pairs", "background_pairs.tsv"), ("background_regions", "background_regions.tsv"),
            ("background_regions_bed", "background_regions.bed"), ("run", "run.json"),
        )
    }
    paths["log"] = getattr(args, "log_file", None) or Path(f"{args.output_prefix}.log")
    resolved = [path.resolve() for path in paths.values()]
    inputs = {path.resolve() for path in (args.regions_bed, args.background_bed, args.vcf)}
    if len(set(resolved)) != len(resolved) or any(path in inputs for path in resolved):
        raise ValueError("output/log paths must be distinct and cannot be input files")
    existing = {name for name, path in paths.items() if path.exists()}
    reusable = {"background_pairs", "background_regions"} & existing
    if existing and not reusable:
        if existing == {"log"}:
            raise ValueError("log file already exists; choose a new --log-file")
        raise ValueError(
            "output prefix already exists but contains no reusable background_pairs.tsv or "
            "background_regions.tsv; choose a new --output-prefix"
        )
    return paths


PAIR_REUSE_SIGNATURE_FIELDS = (
    "analysis_semantics_version",
    "regions_bed", "background_bed", "vcf", "window_bp", "window_step_bp",
    "block_count", "min_snps_per_block", "max_snps_per_block",
    "min_complete_haplotypes", "min_pair_complete_haplotype_fraction",
    "max_exact_entropy_combinations", "min_maf",
    "max_site_missing_rate", "max_window_missing_filtered_fraction",
    "vcf_prefiltered", "include_filtered", "allow_unphased_missing",
    "include_target_in_background", "pair_distance_tolerance",
    "pair_match_snp_counts", "background_pairs_per_distance", "skip_pair_test",
)
REGION_REUSE_SIGNATURE_FIELDS = PAIR_REUSE_SIGNATURE_FIELDS + (
    "background_regions", "region_length_tolerance_bp", "region_block_count_tolerance",
    "region_block_count_allow_swap", "region_entropy_max_wasserstein",
    "skip_region_entropy_matching", "region_match_fallback", "region_fallback_max_score",
    "skip_background_region_pair_p", "skip_region_test",
    "pair_significance_threshold",
)


def reuse_signature_from_parameters(parameters: Mapping[str, object]) -> Dict[str, object]:
    """Canonical controls that determine already-calculated background values."""
    fixed = parameters.get("snps_per_block")
    if fixed is not None:
        min_snps = max_snps = int(fixed)
    else:
        min_snps = int(parameters.get("min_snps_per_block") or 10)
        max_snps = int(parameters.get("max_snps_per_block") or 10)
    window_bp = int(parameters.get("window_bp", 10_000))

    def resolved_path(name: str) -> str:
        return str(Path(str(parameters[name])).resolve())

    return {
        "analysis_semantics_version": ANALYSIS_SEMANTICS_VERSION,
        "regions_bed": resolved_path("regions_bed"),
        "background_bed": resolved_path("background_bed"),
        "vcf": resolved_path("vcf"),
        "window_bp": window_bp,
        "window_step_bp": int(parameters.get("window_step_bp") or window_bp),
        "block_count": int(parameters.get("block_count", 10)),
        "min_snps_per_block": min_snps,
        "max_snps_per_block": max_snps,
        "min_complete_haplotypes": int(parameters.get("min_complete_haplotypes", 2)),
        "min_pair_complete_haplotype_fraction": float(
            parameters.get("min_pair_complete_haplotype_fraction", 0.20)
        ),
        "max_exact_entropy_combinations": int(parameters.get("max_exact_entropy_combinations", 100_000)),
        "min_maf": float(parameters.get("min_maf", 0.1)),
        "max_site_missing_rate": float(parameters.get("max_site_missing_rate", 0.05)),
        "max_window_missing_filtered_fraction": float(
            parameters.get("max_window_missing_filtered_fraction", 0.10)
        ),
        "vcf_prefiltered": bool(parameters.get("vcf_prefiltered", False)),
        "include_filtered": bool(parameters.get("include_filtered", False)),
        "allow_unphased_missing": bool(parameters.get("allow_unphased_missing", False)),
        "include_target_in_background": bool(parameters.get("include_target_in_background", False)),
        "pair_distance_tolerance": float(parameters.get("pair_distance_tolerance", 0.10)),
        "pair_match_snp_counts": bool(parameters.get("pair_match_snp_counts", False)),
        "background_pairs_per_distance": int(parameters.get("background_pairs_per_distance", 5000)),
        "background_regions": int(parameters.get("background_regions", 1000)),
        "region_length_tolerance_bp": int(parameters.get("region_length_tolerance_bp", 10_000)),
        "region_block_count_tolerance": int(parameters.get("region_block_count_tolerance", 1)),
        "region_block_count_allow_swap": bool(parameters.get("region_block_count_allow_swap", True)),
        "region_entropy_max_wasserstein": float(parameters.get("region_entropy_max_wasserstein", 0.05)),
        "skip_region_entropy_matching": bool(parameters.get("skip_region_entropy_matching", False)),
        "region_match_fallback": str(parameters.get("region_match_fallback", "none")),
        "region_fallback_max_score": float(parameters.get("region_fallback_max_score", 2.0)),
        "skip_background_region_pair_p": bool(parameters.get("skip_background_region_pair_p", False)),
        "skip_pair_test": bool(parameters.get("skip_pair_test", False)),
        "skip_region_test": bool(parameters.get("skip_region_test", False)),
        "pair_significance_threshold": float(parameters.get("pair_significance_threshold", 0.01)),
    }


def reuse_signature_from_args(args: argparse.Namespace) -> Dict[str, object]:
    return reuse_signature_from_parameters(vars(args))


def load_prior_run_metadata(
    path: Path, args: argparse.Namespace, fields: Sequence[str], required: bool,
) -> Optional[Dict[str, object]]:
    if not path.exists():
        if required:
            raise ValueError(
                f"cannot safely reuse background regions without completed metadata: missing {path}"
            )
        return None
    try:
        with path.open(encoding="utf-8") as handle:
            metadata = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read existing reuse metadata {path}: {error}") from error
    if required and metadata.get("status") != "complete":
        raise ValueError(f"cannot reuse background regions from an incomplete run: {path}")
    prior = metadata.get("background_reuse_signature")
    if not isinstance(prior, dict):
        parameters = metadata.get("parameters")
        if not isinstance(parameters, dict):
            raise ValueError(f"existing run metadata lacks parameters required for reuse: {path}")
        prior = reuse_signature_from_parameters(parameters)
    current = reuse_signature_from_args(args)
    mismatches = [
        f"{name}: old={prior.get(name)!r}, current={current.get(name)!r}"
        for name in fields if prior.get(name) != current.get(name)
    ]
    if mismatches:
        raise ValueError("existing background is incompatible with current controls: " + "; ".join(mismatches))
    return metadata


def read_existing_tsv(path: Path, required_fields: Sequence[str]) -> List[Dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = set(reader.fieldnames or ())
            missing = [field for field in required_fields if field not in fields]
            if missing:
                raise ValueError(f"existing file {path} is missing columns: {','.join(missing)}")
            return list(reader)
    except OSError as error:
        raise ValueError(f"cannot read existing background file {path}: {error}") from error


def required_finite_float(row: Mapping[str, str], field: str, path: Path) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"existing file {path} has invalid {field}: {row.get(field)!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"existing file {path} has non-finite {field}: {row.get(field)!r}")
    return value


def run(args: argparse.Namespace) -> Dict[str, int]:
    interval = getattr(args, "log_interval_seconds", 30.0)
    every_attempts = getattr(args, "log_every_attempts", 1000)
    if not math.isfinite(interval) or interval <= 0 or every_attempts < 1:
        raise ValueError("--log-interval-seconds must be finite and > 0; --log-every-attempts must be >= 1")
    output_paths = output_files(args)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_paths["log"].parent.mkdir(parents=True, exist_ok=True)
    reuse_enabled = {
        "background_pairs": not getattr(args, "skip_pair_test", False),
        "background_regions": not getattr(args, "skip_region_test", False),
    }
    existing_background_paths = [
        output_paths[name] for name in ("background_pairs", "background_regions")
        if output_paths[name].exists()
    ]
    reused_paths = [output_paths[name] for name in ("background_pairs", "background_regions")
                    if reuse_enabled[name] and output_paths[name].exists()]
    logger = RunLogger(
        output_paths["log"], interval, every_attempts,
        append=bool(existing_background_paths and output_paths["log"].exists()),
    )
    try:
        logger.log(f"RUN START version={VERSION} log_file={output_paths['log']} seed={args.seed}")
        if reused_paths:
            logger.log("REUSE DETECTED files=" + ",".join(str(path) for path in reused_paths))
        if (not getattr(args, "skip_pair_test", False) and
                1 / (args.background_pairs_per_distance + 1) >= args.pair_significance_threshold):
            logger.log("pair background is too small to produce p < the significance threshold", "WARNING")
        summary = _run(args, logger, output_paths)
        logger.log("RUN COMPLETE " + " ".join(f"{key}={value}" for key, value in summary.items()))
        return summary
    except Exception as error:
        logger.log(f"RUN FAILED {type(error).__name__}: {error}", "ERROR")
        raise
    except KeyboardInterrupt:
        logger.log("RUN INTERRUPTED by user", "WARNING")
        raise
    finally:
        logger.close()


def _run(args: argparse.Namespace, logger: RunLogger,
         output_paths: Mapping[str, Path]) -> Dict[str, int]:
    min_snps, max_snps = resolve_snp_limits(args)
    length_tolerance = getattr(args, "region_length_tolerance_bp", 10_000)
    block_count_tolerance = getattr(args, "region_block_count_tolerance", 1)
    block_count_allow_swap = getattr(args, "region_block_count_allow_swap", True)
    entropy_max_wasserstein = getattr(args, "region_entropy_max_wasserstein", 0.05)
    entropy_matching_enabled = not getattr(args, "skip_region_entropy_matching", False)
    pair_snp_selection_mode = PAIR_SNP_SELECTION_MODE
    pair_distance_tolerance = getattr(args, "pair_distance_tolerance", 0.10)
    pair_match_snp_counts = getattr(args, "pair_match_snp_counts", False)
    pair_batch_size = getattr(args, "pair_batch_size", 1000)
    region_batch_size = getattr(args, "region_batch_size", 1000)
    min_pair_complete_fraction = getattr(
        args, "min_pair_complete_haplotype_fraction", 0.20,
    )
    pair_snp_matching = "exact_ordered_selected_counts" if pair_match_snp_counts else "none"
    region_match_fallback = getattr(args, "region_match_fallback", "none")
    region_fallback_max_score = getattr(args, "region_fallback_max_score", 2.0)
    configured_region_attempts = getattr(args, "background_region_attempts", None)
    vcf_prefiltered = getattr(args, "vcf_prefiltered", False)
    vcf_access = getattr(args, "vcf_access", "bcftools")
    pair_test_enabled = not getattr(args, "skip_pair_test", False)
    region_test_enabled = not getattr(args, "skip_region_test", False)
    reuse_pair_backgrounds = pair_test_enabled and output_paths["background_pairs"].exists()
    reuse_region_backgrounds = region_test_enabled and output_paths["background_regions"].exists()
    prior_run_metadata: Optional[Dict[str, object]] = None
    if reuse_region_backgrounds:
        prior_run_metadata = load_prior_run_metadata(
            output_paths["run"], args, REGION_REUSE_SIGNATURE_FIELDS, required=True,
        )
        if (pair_test_enabled and not getattr(args, "skip_background_region_pair_p", False) and
                not reuse_pair_backgrounds):
            raise ValueError(
                "reusing background regions with pair-p/breadth enabled also requires the existing "
                f"pair background file {output_paths['background_pairs']}"
            )
    elif reuse_pair_backgrounds and output_paths["run"].exists():
        prior_run_metadata = load_prior_run_metadata(
            output_paths["run"], args, PAIR_REUSE_SIGNATURE_FIELDS, required=False,
        )
    if vcf_access not in ("bcftools", "preload"):
        raise ValueError("--vcf-access must be bcftools or preload")
    window_missing_qc_status = "not_assessable_prefiltered_vcf" if vcf_prefiltered else "enabled"
    if length_tolerance < 0:
        raise ValueError("--region-length-tolerance-bp must be >= 0")
    if block_count_tolerance < 0:
        raise ValueError("--region-block-count-tolerance must be >= 0")
    if not math.isfinite(entropy_max_wasserstein) or entropy_max_wasserstein < 0:
        raise ValueError("--region-entropy-max-wasserstein must be finite and >= 0")
    if not math.isfinite(pair_distance_tolerance) or not 0 <= pair_distance_tolerance < 1:
        raise ValueError("--pair-distance-tolerance must be finite and in [0, 1)")
    if pair_batch_size < 1 or region_batch_size < 1:
        raise ValueError("--pair-batch-size and --region-batch-size must be >= 1")
    if (not math.isfinite(min_pair_complete_fraction) or
            not 0 <= min_pair_complete_fraction <= 1):
        raise ValueError("--min-pair-complete-haplotype-fraction must be finite and in [0, 1]")
    if region_match_fallback not in ("none", "nearest"):
        raise ValueError("--region-match-fallback must be none or nearest")
    if not math.isfinite(region_fallback_max_score) or region_fallback_max_score < 1:
        raise ValueError("--region-fallback-max-score must be finite and >= 1")
    if (region_test_enabled and configured_region_attempts is not None and
            configured_region_attempts < args.background_regions):
        raise ValueError("--background-region-attempts must be >= --background-regions")
    logger.log(f"TEST CONTROLS pair_test_enabled={str(pair_test_enabled).lower()} "
               f"region_test_enabled={str(region_test_enabled).lower()}")
    logger.log(f"PAIR CONTROLS distance_tolerance={pair_distance_tolerance} "
               f"snp_count_matching={pair_snp_matching} distance_units=physical_bp "
               f"batch_size={pair_batch_size}")
    logger.log(f"REGION ENTROPY matching_enabled={str(region_test_enabled and entropy_matching_enabled).lower()} "
               f"max_wasserstein={entropy_max_wasserstein} units=nats weighting=equal_per_available_block")
    logger.log(
        "PAIR COMPLETE COHORT "
        f"minimum_fraction={min_pair_complete_fraction} "
        f"minimum_absolute={getattr(args, 'min_complete_haplotypes', 2)} "
        "selection_mode=pair_specific_deferred_to_dprime"
    )
    background_pair_p_enabled = (
        pair_test_enabled and region_test_enabled and
        not getattr(args, "skip_background_region_pair_p", False)
    )
    if not background_pair_p_enabled:
        logger.log("Background region pair p-values disabled; region breadth test disabled.")
    elif length_tolerance:
        logger.log("Variable central lengths with background region pair p-values may require many "
                   "additional reference-distance banks. These are sampled on demand; "
                   "use --skip-background-region-pair-p if only the region max test is needed.", "WARNING")
    logger.log("BED START reading target regions and preparing allowed background intervals")
    targets = read_regions(args.regions_bed)
    if len({target.region_id for target in targets}) != len(targets):
        raise ValueError("target BED region IDs must be unique")
    step = args.window_step_bp or args.window_bp
    span = args.window_bp + (args.block_count - 1) * step
    if any(target.start < span for target in targets):
        raise ValueError("each target must have a complete left flank; target start is below flank span")
    background = merge_background_intervals(read_regions(args.background_bed))
    envelopes = [replace(target, start=target.start - span, end=target.end + span) for target in targets]
    if not args.include_target_in_background:
        background = subtract_exclusions(background, envelopes)
    if not background:
        raise ValueError("background BED is empty after excluding target regions and their flanks")
    logger.log(f"BED DONE targets={len(targets)} background_intervals={len(background)} "
               f"background_bp={sum(item.end - item.start for item in background)} "
               f"window_bp={args.window_bp} step_bp={step} blocks_per_flank={args.block_count} "
               f"flank_span_bp={span} snp_count_range={min_snps}-{max_snps}")
    if vcf_prefiltered:
        logger.log("VCF PREFILTERED: trusting external site MAF/missing QC; skipping rate calculations "
                   "and filtering. Genotype reading and SNP/FILTER/phasing checks remain.")
        logger.log("Window missing-filtered site fraction cannot be reconstructed from a prefiltered VCF; "
                   "window missing QC disabled and missing-filtered counts/fractions reported as NA. "
                   "SNP-count QC, complete-haplotype checks and entropy matching remain enabled.", "WARNING")
    source = None
    if vcf_access == "bcftools":
        logger.log("VCF INDEX START reading header/index only; windows retrieved on demand")
        source = BCFtoolsVCFSource(args, logger)
        data = source.data
    else:
        load_intervals = background + [
            flank for target in targets for flank in (
                replace(target, start=target.start - span, end=target.start),
                replace(target, start=target.end, end=target.end + span),
            )
        ]
        logger.log("VCF START preload mode: scanning entire VCF and retaining target/background genotypes")
        data = read_vcf(
            args.vcf, load_intervals, flank_bp=0, min_maf=args.min_maf,
            min_call_rate=1.0 - args.max_site_missing_rate,
            include_filtered=args.include_filtered,
            allow_unphased_missing=args.allow_unphased_missing, include_region_interiors=True,
            progress_callback=logger.vcf_progress,
            skip_maf_missing_qc=vcf_prefiltered,
        )
    store = BlockStore(data, args, logger, source=source)
    required_complete_haplotypes = max(
        2,
        getattr(args, "min_complete_haplotypes", 2),
        math.ceil(2 * data.sample_count * min_pair_complete_fraction),
    )
    target_rows = []
    for index, target in enumerate(targets, 1):
        logger.log(f"TARGET {index}/{len(targets)} id={target.region_id} START "
                   f"windows={2 * args.block_count} cross_pairs={args.block_count ** 2}")
        rows = store.region_pairs(target)
        target_rows.append(rows)
        logger.log(f"TARGET {index}/{len(targets)} id={target.region_id} DONE "
                   f"pairs={len(rows)} valid_pairs={sum(valid_row(row) for row in rows)}")
    masks = [region_valid_pairs(rows) for rows in target_rows]
    count_requirements = [available_block_counts(rows) for rows in target_rows]
    entropy_requirements: List[Optional[RegionEntropyProfile]] = (
        [store.region_entropy_profile(target, rows) for target, rows in zip(targets, target_rows)]
        if region_test_enabled and entropy_matching_enabled else [None] * len(targets)
    )
    strata = sorted({pair_background_key(row, pair_match_snp_counts)
                     for rows in target_rows for row in rows if valid_row(row)})
    distances = sorted({key[0] for key in strata})
    if pair_test_enabled:
        logger.log(f"PAIR BACKGROUND PLAN strata={len(strata)} distances={len(distances)} "
                   f"requested_per_stratum={args.background_pairs_per_distance} "
                   f"requested_total={len(strata) * args.background_pairs_per_distance}")
    else:
        logger.log("PAIR TEST SKIPPED background_sampling=false empirical_p=false")
    if not distances:
        # Preserve QC details even when no inference is possible.
        for rows in target_rows:
            attach_pair_p(rows, {}, args.pair_significance_threshold)
            for row in rows:
                row["pair_in_region_test"] = "false"
        write_tsv(output_paths["target_pairs"], TARGET_PAIR_FIELDS, [row for rows in target_rows for row in rows])
        raise ValueError(f"no target cross-flank pair passes QC; see {output_paths['target_pairs']}")

    # Fail impossible BED geometries before spending time on Monte Carlo LD.
    pair_domains = ({distance: pair_distance_sampling_domains(
        background, args.window_bp, distance, pair_distance_tolerance)
                     for distance in distances} if pair_test_enabled else {})
    for distance, domains in pair_domains.items():
        available = sum(domain.size for domain in domains)
        if available < args.background_pairs_per_distance:
            raise ValueError(f"reference distance {distance}: only {available} legal background pair configurations; requested {args.background_pairs_per_distance}")
    region_domains: Dict[int, List[RegionLengthDomain]] = {}
    for target, mask in zip(targets, masks):
        if not region_test_enabled or not mask:
            continue
        central_bp = target.end - target.start
        if central_bp not in region_domains:
            region_domains[central_bp] = region_length_sampling_domains(background, central_bp, span, length_tolerance)
        available = sum(domain.size for domain in region_domains[central_bp])
        if available < args.background_regions:
            raise ValueError(f"target {target.region_id}: only {available} complete background region configurations; requested {args.background_regions}. Central interval AND both flanks must be inside background BED")

    if reuse_pair_backgrounds:
        pair_backgrounds, pair_sampling_info = load_existing_pair_backgrounds(
            output_paths["background_pairs"], args,
        )
        missing_initial_strata = sorted(set(strata) - set(pair_backgrounds))
        if missing_initial_strata:
            raise ValueError(
                f"existing pair background {output_paths['background_pairs']} lacks required target "
                f"strata: {missing_initial_strata}"
            )
        logger.log(
            f"PAIR BACKGROUND REUSE strata={len(pair_backgrounds)} "
            f"values={sum(len(values) for values in pair_backgrounds.values())} "
            f"file={output_paths['background_pairs']} sampling_skipped=true"
        )
    else:
        pair_backgrounds = {}
        pair_sampling_info = {}
    background_threads = getattr(args, "threads", 1)
    if background_threads < 1:
        raise ValueError("--threads must be >= 1")
    background_slots = threading.BoundedSemaphore(background_threads)

    def run_background_task(function: Callable[[T], U], item: T) -> U:
        with background_slots:
            return function(item)

    pair_domains_lock = threading.RLock()
    logger.log(
        f"BACKGROUND THREADS workers={background_threads} "
        f"pair_banks={('parallel' if pair_test_enabled else 'skipped')} "
        f"region_candidates={('parallel' if region_test_enabled else 'skipped')} "
        "output_order=deterministic"
    )
    logger.log(
        f"BACKGROUND BATCHING pair_candidates={pair_batch_size} "
        f"region_candidates={region_batch_size} "
        f"indexed_query_mode={'merged_batch' if source is not None else 'preloaded_memory'}"
    )
    if region_test_enabled:
        logger.log("REGION SCREEN two_stage=true candidate_phase=dprime_free "
                   "pair_eligibility=window_qc+snp_count+fixed_pair_complete_cohort+bilateral_polymorphism "
                   "pair_specific_entropy=false "
                   f"marginal_entropy={('enabled' if entropy_matching_enabled else 'skipped')} "
                   "snp_entropy_selection=deferred_to_dprime "
                   "selected_region_phase=full_dprime")

    def sample_pair_background(
        work: Tuple[PairBackgroundKey, str, str],
    ) -> Tuple[PairBackgroundKey, List[float], List[Dict[str, object]], str, Dict[str, object]]:
        key, label, bank_source = work
        distance, left_snp_count, right_snp_count = key
        lower, upper = pair_distance_bounds(distance, pair_distance_tolerance, args.window_bp)
        snp_label = f"{left_snp_count},{right_snp_count}" if pair_match_snp_counts else "unmatched"
        key_label = f"{distance}:{left_snp_count}:{right_snp_count}" if pair_match_snp_counts else str(distance)
        pair_id_prefix = f"pair_d{distance}" + (f"_n{left_snp_count}_{right_snp_count}"
                                               if pair_match_snp_counts else "")
        with pair_domains_lock:
            domains = pair_domains.get(distance)
        if domains is None:
            calculated = pair_distance_sampling_domains(
                background, args.window_bp, distance, pair_distance_tolerance)
            with pair_domains_lock:
                domains = pair_domains.setdefault(distance, calculated)
        available = sum(domain.size for domain in domains)
        requested = args.background_pairs_per_distance
        if available < requested:
            raise ValueError(f"reference distance {distance}: only {available} legal background pair configurations; requested {requested}")
        accepted: List[float] = []
        output_rows: List[Dict[str, object]] = []
        rejected: Counter = Counter()
        attempts = 0
        sampler = SamplingProgress(
            logger, f"PAIR stratum={label} distance_bp={distance} "
            f"distance_range_bp={lower},{upper} snps={snp_label}",
            requested, available, requested * args.max_attempt_factor,
        )
        rng = random.Random(stable_subseed(args.seed, "pair_background", key))
        candidates = iter(random_pair_windows(
            domains, requested * args.max_attempt_factor, rng,
        ))
        query_batches = 0
        while len(accepted) < requested:
            batch = list(islice(candidates, pair_batch_size))
            if not batch:
                break
            query_batches += 1
            batch_blocks = store.load_window_batch([
                key
                for domain, start, second_start in batch
                for key in ((domain.chrom, start, "left"),
                            (domain.chrom, second_start, "right"))
            ])
            logger.periodic(
                "pair_batch",
                f"PAIR BATCH stratum={label} batch={query_batches} candidates={len(batch)} "
                f"unique_windows={len(batch_blocks)}",
            )
            for domain, start, second_start in batch:
                attempts += 1
                region = Region(
                    domain.chrom, domain.display_chrom, start + args.window_bp, second_start,
                    f"{pair_id_prefix}_{len(accepted) + 1}", 0,
                )
                first = batch_blocks[(domain.chrom, start, "left")]
                second = batch_blocks[(domain.chrom, second_start, "right")]
                if not first.passes_missing_qc or not second.passes_missing_qc:
                    rejected["window_excess_missing_filtered_sites"] += 1
                    sampler.update(len(accepted), attempts, rejected)
                    continue
                if len(first.variants) < min_snps or len(second.variants) < min_snps:
                    rejected["insufficient_snps"] += 1
                    sampler.update(len(accepted), attempts, rejected)
                    continue
                if (pair_match_snp_counts and
                        (store.final_snp_count(first), store.final_snp_count(second)) !=
                        (left_snp_count, right_snp_count)):
                    rejected["selected_snp_count_mismatch"] += 1
                    sampler.update(len(accepted), attempts, rejected)
                    continue
                row = store.summarize(region, first, second)
                if not valid_row(row):
                    rejected[row["status"]] += 1
                    sampler.update(len(accepted), attempts, rejected)
                    continue
                accepted.append(float(row["multiallelic_dprime"]))
                output_rows.append(dict(
                    row, background_pair_id=region.region_id, distance_stratum_bp=distance,
                    snp_number_stratum_block1=left_snp_count if pair_match_snp_counts else "NA",
                    snp_number_stratum_block2=right_snp_count if pair_match_snp_counts else "NA",
                    pair_distance_min_bp=lower, pair_distance_max_bp=upper,
                    pair_distance_tolerance=precise_float(pair_distance_tolerance),
                    pair_snp_count_matching=pair_snp_matching,
                ))
                sampler.update(len(accepted), attempts, rejected)
                if len(accepted) == requested:
                    break
        if len(accepted) != requested:
            sampler.update(len(accepted), attempts, rejected, force=True)
            raise ValueError(f"reference distance {distance}, distance range {lower}..{upper}, SNP counts {snp_label}: accepted {len(accepted)}/{requested} background pairs after {attempts} attempts; rejections={dict(rejected)}; increase --max-attempt-factor, reduce requested count, or revise BED/QC")
        info = {
            "distance_bp": distance, "distance_min_bp": lower, "distance_max_bp": upper,
            "distance_tolerance": pair_distance_tolerance, "snp_count_matching": pair_snp_matching,
            "block1_selected_snp_count": left_snp_count if pair_match_snp_counts else None,
            "block2_selected_snp_count": right_snp_count if pair_match_snp_counts else None,
            "accepted": len(accepted), "attempts": attempts,
            "query_batches": query_batches, "candidate_batch_size": pair_batch_size,
            "source": bank_source, "legal_start_count": available,
            "legal_pair_configuration_count": available, "rejections": dict(rejected),
            "qc_rejections": {reason: count for reason, count in rejected.items() if reason != "selected_snp_count_mismatch"},
            "selected_snp_count_rejections": rejected["selected_snp_count_mismatch"],
        }
        return key, sorted(accepted), output_rows, key_label, info

    def sample_pair_backgrounds(
        keys: Iterable[PairBackgroundKey], writer: csv.DictWriter, bank_source: str,
    ) -> None:
        missing = [key for key in sorted(set(keys)) if key not in pair_backgrounds]
        base_index = len(pair_backgrounds)
        work = [(
            key,
            (f"{base_index + index}/{len(strata)}" if bank_source == "target_pairs"
             else f"{base_index + index} on_demand=true"),
            bank_source,
        ) for index, key in enumerate(missing, 1)]
        for key, accepted, output_rows, key_label, info in ordered_thread_map(
                lambda item: run_background_task(sample_pair_background, item),
                work, background_threads):
            pair_backgrounds[key] = accepted
            pair_sampling_info[key_label] = info
            writer.writerows(serialized(row) for row in output_rows)
            logger.log(f"PAIR BACKGROUND TOTAL completed_strata={len(pair_backgrounds)} "
                       f"initial_target_strata={len(strata)} "
                       f"accepted={sum(len(values) for values in pair_backgrounds.values())} "
                       f"source={bank_source}")

    if pair_test_enabled:
        with output_paths["background_pairs"].open(
                "a" if reuse_pair_backgrounds else "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=BACKGROUND_PAIR_FIELDS, delimiter="\t")
            if not reuse_pair_backgrounds:
                writer.writeheader()
            sample_pair_backgrounds(strata, writer, "target_pairs")
        logger.log("TARGET PAIR P START assigning distance-band-conditioned p-values "
                   f"snp_count_matching={pair_snp_matching}")
    elif not output_paths["background_pairs"].exists():
        write_tsv(output_paths["background_pairs"], BACKGROUND_PAIR_FIELDS, [])
    for rows, mask in zip(target_rows, masks):
        attach_pair_p(rows, pair_backgrounds if pair_test_enabled else {},
                      args.pair_significance_threshold,
                      pair_match_snp_counts, pair_distance_tolerance, args.window_bp)
        for row in rows:
            row["pair_in_region_test"] = str(
                region_test_enabled and row_key(row) in mask
            ).lower()
    write_tsv(output_paths["target_pairs"], TARGET_PAIR_FIELDS, [row for rows in target_rows for row in rows])
    logger.log(
        f"TARGET PAIR {'P DONE' if pair_test_enabled else 'DPRIME ONLY DONE'} "
        f"output={output_paths['target_pairs']}"
    )

    summaries: List[Dict[str, object]] = []
    region_sampling_info: Dict[str, object] = {}
    if region_test_enabled:
        logger.log(f"REGION BACKGROUND PLAN targets={len(targets)} eligible_targets={sum(bool(mask) for mask in masks)} "
                   f"requested_per_target={args.background_regions} "
                   f"requested_total={sum(bool(mask) for mask in masks) * args.background_regions}")
    else:
        logger.log("REGION TEST SKIPPED background_sampling=false empirical_p=false")
    if reuse_region_backgrounds:
        reused_region_rows = load_existing_region_backgrounds(
            output_paths["background_regions"], args, targets, masks, count_requirements,
            entropy_requirements, background_pair_p_enabled, output_paths["summary"],
        )
        logger.log(
            f"REGION BACKGROUND REUSE targets={len(reused_region_rows)} "
            f"regions={sum(len(rows) for rows in reused_region_rows.values())} "
            f"file={output_paths['background_regions']} sampling_skipped=true"
        )
        if not output_paths["background_regions_bed"].exists():
            with output_paths["background_regions_bed"].open("x", encoding="utf-8") as reused_bed:
                for target in targets:
                    for reused_row in reused_region_rows.get(target.region_id, ()):
                        reused_bed.write(
                            f"{reused_row['chrom']}\t{reused_row['region_start']}\t"
                            f"{reused_row['region_end']}\t{reused_row['background_region_id']}\n"
                        )
    else:
        reused_region_rows = {}
    region_handle_context = (
        nullcontext(None) if reuse_region_backgrounds or (
            not region_test_enabled and output_paths["background_regions"].exists()
        ) else
        output_paths["background_regions"].open("w", encoding="utf-8", newline="")
    )
    bed_handle_context = (
        nullcontext(None) if reuse_region_backgrounds or (
            not region_test_enabled and output_paths["background_regions_bed"].exists()
        ) else
        output_paths["background_regions_bed"].open("w", encoding="utf-8")
    )
    with region_handle_context as handle, bed_handle_context as bed_handle, \
         output_paths["background_pairs"].open("a", encoding="utf-8", newline="") as pair_handle:
        writer = (csv.DictWriter(handle, fieldnames=BACKGROUND_REGION_FIELDS, delimiter="\t")
                  if handle is not None else None)
        pair_writer = csv.DictWriter(pair_handle, fieldnames=BACKGROUND_PAIR_FIELDS, delimiter="\t")
        if writer is not None:
            writer.writeheader()
        for target_index, (target, rows, mask, required_counts, target_entropy) in enumerate(
                zip(targets, target_rows, masks, count_requirements, entropy_requirements)):
            min_length, max_length = region_length_bounds(target.end - target.start, length_tolerance)
            base: Dict[str, object] = {
                "region_id": target.region_id, "chrom": target.display_chrom,
                "region_start": target.start, "region_end": target.end,
                "central_length_bp": target.end - target.start, "flank_span_bp": span,
                "n_total_cross_pairs": len(rows), "n_target_qc_failed_pairs": len(rows) - len(mask),
                "pair_significance_threshold": args.pair_significance_threshold,
                "pair_test_enabled": str(pair_test_enabled).lower(),
                "region_test_enabled": str(region_test_enabled).lower(),
                "pair_distance_tolerance": pair_distance_tolerance,
                "pair_snp_count_matching": pair_snp_matching,
                "region_qc_matching": (
                    "not_run" if not region_test_enabled else
                    "available_block_counts_and_marginal_entropy_distributions"
                    if entropy_matching_enabled else
                    "available_block_counts_with_tolerance_and_optional_swap"
                ),
                "min_snps_per_block": min_snps, "max_snps_per_block": max_snps,
                "region_snp_count_matching": "shared_range_and_entropy_cap",
                "required_block_snp_counts": None,
                "pair_snp_selection_mode": pair_snp_selection_mode,
                "min_pair_complete_haplotype_fraction": min_pair_complete_fraction,
                "required_complete_haplotype_count": required_complete_haplotypes,
                "region_length_tolerance_bp": length_tolerance,
                "background_min_central_length_bp": min_length, "background_max_central_length_bp": max_length,
                "required_available_left_blocks": required_counts[0], "required_available_right_blocks": required_counts[1],
                "region_block_count_tolerance": block_count_tolerance,
                "region_block_count_allow_swap": str(block_count_allow_swap).lower(),
                "region_pair_distance_matching": "none",
                "region_entropy_matching_enabled": str(
                    region_test_enabled and entropy_matching_enabled
                ).lower(),
                "region_entropy_max_wasserstein": entropy_max_wasserstein,
                "region_match_fallback": region_match_fallback,
                "region_fallback_max_score": region_fallback_max_score,
                "background_region_attempt_limit": (
                    configured_region_attempts if configured_region_attempts is not None
                    else args.background_regions * args.max_attempt_factor
                ),
                "vcf_maf_missing_qc_skipped": str(vcf_prefiltered).lower(),
                "window_missing_qc_status": window_missing_qc_status,
                **entropy_profile_fields(target_entropy),
                "background_region_pair_p_enabled": str(background_pair_p_enabled).lower(),
                "region_breadth_test_status": (
                    "region_test_skipped" if not region_test_enabled else
                    "disabled_pair_test" if not pair_test_enabled else
                    "enabled" if background_pair_p_enabled else
                    "disabled_background_region_pair_p"
                ),
            }
            if not mask:
                summaries.append(dict(base, status="no_valid_cross_flank_pairs", n_background_regions=0))
                logger.log(f"REGION target={target.region_id} SKIP no valid target cross-flank pairs")
                continue
            observed = analyze_region(
                rows, mask, args.pair_significance_threshold,
                calculate_breadth=pair_test_enabled,
            )
            if not region_test_enabled:
                summaries.append(dict(
                    base, **observed,
                    region_max_empirical_p=None,
                    region_breadth_empirical_p=None,
                    n_background_regions=0,
                    background_region_attempts=0,
                    background_region_qc_rejections=0,
                    background_region_snp_count_rejections=0,
                    background_region_block_count_rejections=0,
                    background_region_entropy_rejections=0,
                    background_region_rejections=0,
                    region_strict_match_count=0,
                    region_fallback_candidate_count=0,
                    region_fallback_selected_count=0,
                    selected_region_match_score_max=None,
                    status="ok_region_test_skipped",
                ))
                region_sampling_info[target.region_id] = {
                    "accepted": 0, "attempts": 0, "skipped": True,
                    "reason": "--skip-region-test",
                }
                logger.log(f"REGION target={target.region_id} SKIP --skip-region-test")
                continue
            domains = region_domains[target.end - target.start]
            available = sum(domain.size for domain in domains)
            requested = args.background_regions
            attempt_limit = min(
                available,
                configured_region_attempts if configured_region_attempts is not None
                else requested * args.max_attempt_factor,
            )
            if available < requested:
                raise ValueError(f"target {target.region_id}: only {available} complete background region configurations; requested {requested}. Central interval AND both flanks must be inside background BED")
            if reuse_region_backgrounds:
                reused_rows = reused_region_rows[target.region_id]
                maxima = sorted(required_finite_float(
                    reused_row, "max_dprime", output_paths["background_regions"]
                ) for reused_row in reused_rows)
                fractions = sorted(required_finite_float(
                    reused_row, "significant_pair_fraction", output_paths["background_regions"]
                ) for reused_row in reused_rows) if background_pair_p_enabled else []
                selected_scores = [required_finite_float(
                    reused_row, "region_match_score", output_paths["background_regions"]
                ) for reused_row in reused_rows]
                prior_counts = (
                    prior_run_metadata.get("region_sampling_counts", {})
                    if isinstance(prior_run_metadata, dict) else {}
                )
                prior_info_value = (
                    prior_counts.get(target.region_id, {})
                    if isinstance(prior_counts, dict) else {}
                )
                prior_info = dict(prior_info_value) if isinstance(prior_info_value, dict) else {}
                attempts = prior_info.get("attempts")
                rejected = int(prior_info.get("rejections", 0))
                block_count_rejected = int(prior_info.get("block_count_rejections", 0))
                entropy_rejected = int(prior_info.get("entropy_rejections", 0))
                qc_rejected = int(prior_info.get(
                    "qc_rejections", rejected - block_count_rejected - entropy_rejected,
                ))
                strict_match_count = int(prior_info.get(
                    "strict_match_count",
                    sum(existing_bool(reused_row, "region_match_strict",
                                      output_paths["background_regions"])
                        for reused_row in reused_rows),
                ))
                fallback_selected_count = int(prior_info.get(
                    "fallback_selected_count",
                    sum(reused_row["region_match_selection"] == "nearest_fallback"
                        for reused_row in reused_rows),
                ))
                fallback_candidate_count = int(prior_info.get("fallback_candidate_count", 0))
                summaries.append(dict(
                    base, **observed,
                    region_max_empirical_p=empirical_tail_p(
                        float(observed["max_dprime"]), maxima,
                    ),
                    region_breadth_empirical_p=(
                        empirical_tail_p(float(observed["significant_pair_fraction"]), fractions)
                        if background_pair_p_enabled else None
                    ),
                    n_background_regions=len(maxima), background_region_attempts=attempts,
                    background_region_qc_rejections=qc_rejected,
                    background_region_snp_count_rejections=0,
                    background_region_block_count_rejections=block_count_rejected,
                    background_region_entropy_rejections=entropy_rejected,
                    background_region_rejections=rejected,
                    region_match_fallback=region_match_fallback,
                    region_fallback_max_score=region_fallback_max_score,
                    background_region_attempt_limit=prior_info.get("attempt_limit"),
                    region_strict_match_count=strict_match_count,
                    region_fallback_candidate_count=fallback_candidate_count,
                    region_fallback_selected_count=fallback_selected_count,
                    selected_region_match_score_max=max(selected_scores),
                    status="ok",
                ))
                prior_info.update({
                    "accepted": len(maxima), "reused": True,
                    "source": str(output_paths["background_regions"]),
                })
                region_sampling_info[target.region_id] = prior_info
                logger.log(
                    f"REGION target={target.region_id} REUSED accepted={len(maxima)}/{requested} "
                    "sampling_skipped=true"
                )
                continue
            maxima: List[float] = []
            fractions: List[float] = []
            attempts = 0
            rejected = 0
            block_count_rejected = 0
            entropy_rejected = 0
            accepted_entropy_distances: List[float] = []
            rejection_reasons: Counter = Counter()
            accepted_lengths: Counter = Counter()
            accepted_orientations: Counter = Counter()
            accepted_block_counts: Counter = Counter()
            strict_match_count = 0
            fallback_candidate_count = 0
            fallback_selected_count = 0
            selected_match_scores: List[float] = []
            fallback_heap: List[Tuple[float, int, Tuple[object, ...]]] = []
            rng = random.Random(args.seed + 1_000_003 * (target_index + 1))
            sampler = SamplingProgress(
                logger, f"REGION target={target.region_id} target_index={target_index + 1}/{len(targets)}",
                requested, available, attempt_limit, legal_unit="legal_regions",
            )
            logger.log(f"REGION CONTROLS target={target.region_id} central_length_range={min_length}-{max_length} "
                       f"target_left_blocks={required_counts[0]} target_right_blocks={required_counts[1]} "
                       f"block_count_tolerance={block_count_tolerance} "
                       f"block_count_allow_swap={str(block_count_allow_swap).lower()} "
                       "rank_pair_mask_matching=false pair_distance_matching=false "
                       f"entropy_matching_enabled={str(entropy_matching_enabled).lower()} "
                       f"entropy_max_wasserstein={entropy_max_wasserstein} "
                       f"fallback={region_match_fallback} fallback_max_score={region_fallback_max_score} "
                       f"attempt_limit={attempt_limit} batch_size={region_batch_size} "
                       f"left_right_symmetric={str(block_count_allow_swap).lower()} "
                       f"target_entropy_profile={json.dumps(entropy_profile_fields(target_entropy), sort_keys=True)}")

            def accept_region_candidate(
                domain: RegionLengthDomain, region: Region, candidate_screen: RegionCandidateScreen,
                candidate_match: RegionCandidateMatch,
                selection: str,
            ) -> None:
                nonlocal fallback_selected_count
                assert writer is not None and bed_handle is not None
                count_orientation = candidate_match.orientation
                entropy_match = candidate_match.entropy
                central_length = region.end - region.start
                candidate_rows = store.region_pairs_from_blocks(
                    region, candidate_screen.left_blocks, candidate_screen.right_blocks,
                )
                candidate_mask = region_valid_pairs(candidate_rows)
                if candidate_mask != candidate_screen.eligible_pairs:
                    raise ValueError(
                        "internal error: D'-free region screen disagrees with full pair validity"
                    )
                candidate_entropy = candidate_screen.entropy
                region = replace(region, region_id=f"{target.region_id}|bg_{len(maxima) + 1}")
                selected_match_scores.append(candidate_match.score)
                if selection == "nearest_fallback":
                    fallback_selected_count += 1
                if background_pair_p_enabled:
                    sample_pair_backgrounds(
                        {pair_background_key(row, pair_match_snp_counts)
                         for row in candidate_rows if valid_row(row)},
                        pair_writer, "background_region_pairs",
                    )
                    attach_pair_p(candidate_rows, pair_backgrounds, args.pair_significance_threshold,
                                  pair_match_snp_counts, pair_distance_tolerance, args.window_bp)
                stats = analyze_region(
                    candidate_rows, candidate_mask, args.pair_significance_threshold,
                    calculate_breadth=background_pair_p_enabled,
                )
                maxima.append(float(stats["max_dprime"]))
                if entropy_match.max_wasserstein is not None:
                    accepted_entropy_distances.append(entropy_match.max_wasserstein)
                accepted_lengths[central_length] += 1
                actual_counts = available_block_counts(candidate_rows)
                aligned_counts = actual_counts if count_orientation == "direct" else actual_counts[::-1]
                accepted_orientations[count_orientation] += 1
                accepted_block_counts[f"{actual_counts[0]},{actual_counts[1]}"] += 1
                if background_pair_p_enabled:
                    fractions.append(float(stats["significant_pair_fraction"]))
                writer.writerow(serialized(dict(stats,
                    target_region_id=target.region_id, background_region_id=region.region_id,
                    chrom=domain.display_chrom, region_start=region.start, region_end=region.end,
                    envelope_start=region.start - span, envelope_end=region.end + span,
                    n_actual_qc_valid_pairs=sum(valid_row(row) for row in candidate_rows),
                    n_compared_pairs=len(candidate_mask), target_pair_mask=None,
                    required_block_snp_counts=None,
                    central_length_bp=central_length, target_central_length_bp=target.end - target.start,
                    central_length_delta_bp=central_length - (target.end - target.start),
                    required_available_left_blocks=required_counts[0], required_available_right_blocks=required_counts[1],
                    region_block_count_tolerance=block_count_tolerance,
                    region_block_count_allow_swap=str(block_count_allow_swap).lower(),
                    block_count_match_orientation=count_orientation,
                    matched_left_block_count_delta=aligned_counts[0] - required_counts[0],
                    matched_right_block_count_delta=aligned_counts[1] - required_counts[1],
                    region_match_fallback=region_match_fallback,
                    region_fallback_max_score=region_fallback_max_score,
                    background_region_attempt_limit=attempt_limit,
                    region_match_selection=selection,
                    region_match_score=candidate_match.score,
                    region_match_strict=str(candidate_match.strict).lower(),
                    region_entropy_matching_enabled=str(entropy_matching_enabled).lower(),
                    vcf_maf_missing_qc_skipped=str(vcf_prefiltered).lower(),
                    window_missing_qc_status=window_missing_qc_status,
                    region_entropy_max_wasserstein=entropy_max_wasserstein,
                    matched_left_entropy_wasserstein=entropy_match.left_wasserstein,
                    matched_right_entropy_wasserstein=entropy_match.right_wasserstein,
                    region_entropy_wasserstein=entropy_match.max_wasserstein,
                    direct_entropy_max_wasserstein=entropy_match.direct_max_wasserstein,
                    swapped_entropy_max_wasserstein=entropy_match.swapped_max_wasserstein,
                    **entropy_profile_fields(candidate_entropy),
                    background_valid_pair_labels=",".join(f"L{left}-R{right}" for left, right in candidate_mask),
                    background_region_pair_p_enabled=str(background_pair_p_enabled).lower(),
                    pair_snp_selection_mode=pair_snp_selection_mode,
                    min_pair_complete_haplotype_fraction=min_pair_complete_fraction,
                    required_complete_haplotype_count=required_complete_haplotypes,
                    pair_distance_tolerance=pair_distance_tolerance,
                    pair_snp_count_matching=pair_snp_matching,
                )))
                bed_handle.write(f"{domain.display_chrom}\t{region.start}\t{region.end}\t{region.region_id}\n")

            def evaluate_region_candidate(
                work: Tuple[int, RegionLengthDomain, Region, int],
            ) -> Tuple[int, RegionLengthDomain, Region, Optional[RegionCandidateScreen],
                       Optional[RegionCandidateMatch], Optional[str]]:
                attempt_number, domain, region, central_length = work
                try:
                    candidate_screen = store.screen_region(region)
                    candidate_match = score_region_candidate_counts(
                        candidate_screen.available_counts, required_counts,
                        candidate_screen.entropy, target_entropy,
                        central_length - (target.end - target.start), length_tolerance,
                        block_count_tolerance, block_count_allow_swap,
                        entropy_max_wasserstein, entropy_matching_enabled,
                    )
                except ValueError as error:
                    return attempt_number, domain, region, None, None, str(error)
                return attempt_number, domain, region, candidate_screen, candidate_match, None

            sampled_candidates = iter(enumerate(
                random_region_intervals(domains, attempt_limit, rng), 1,
            ))
            region_query_batches = 0
            stop_sampling = False
            while not stop_sampling:
                raw_batch = list(islice(sampled_candidates, region_batch_size))
                if not raw_batch:
                    break
                candidate_work = [
                    (
                        attempt_number,
                        domain,
                        Region(
                            domain.chrom, domain.display_chrom, start, start + central_length,
                            f"{target.region_id}|bg_candidate_{attempt_number}", 0,
                        ),
                        central_length,
                    )
                    for attempt_number, (domain, start, central_length) in raw_batch
                ]
                region_query_batches += 1
                batch_blocks = store.load_window_batch([
                    key
                    for _attempt, _domain, candidate_region, _length in candidate_work
                    for key in store.region_window_keys(candidate_region)
                ])
                logger.log(
                    f"REGION BATCH target={target.region_id} batch={region_query_batches} "
                    f"candidates={len(candidate_work)} unique_windows={len(batch_blocks)} "
                    f"attempt_range={candidate_work[0][0]}-{candidate_work[-1][0]}"
                )
                with store._cache_lock:
                    if store._active_window_batch is not None:
                        raise RuntimeError("internal error: nested active window batch")
                    store._active_window_batch = batch_blocks
                try:
                    candidate_results = ordered_thread_map(
                        lambda item: run_background_task(evaluate_region_candidate, item),
                        candidate_work, background_threads,
                    )
                    for (attempts, domain, region, candidate_screen,
                         candidate_match, error_message) in candidate_results:
                        if error_message is not None:
                            rejected += 1
                            if "available block count matching" in error_message:
                                block_count_rejected += 1
                            elif "block haplotype entropy distribution matching" in error_message:
                                entropy_rejected += 1
                            rejection_reasons[error_message] += 1
                            sampler.update(len(maxima), attempts, rejection_reasons)
                            continue
                        assert candidate_screen is not None and candidate_match is not None
                        if not candidate_match.strict:
                            rejected += 1
                            fallback_candidate_count += 1
                            count_failed = max(
                                candidate_match.left_count_error,
                                candidate_match.right_count_error,
                            ) > block_count_tolerance
                            entropy_failed = (
                                entropy_matching_enabled and
                                candidate_match.entropy.max_wasserstein is not None and
                                candidate_match.entropy.max_wasserstein >
                                entropy_max_wasserstein + STATISTIC_TIE_TOLERANCE
                            )
                            if count_failed:
                                block_count_rejected += 1
                            elif entropy_failed:
                                entropy_rejected += 1
                            if count_failed:
                                reason = "background region fails target available block count matching"
                            elif entropy_failed:
                                reason = "background region fails block haplotype entropy distribution matching"
                            else:
                                reason = "background region fails strict covariate matching"
                            rejection_reasons[reason] += 1
                            if (region_match_fallback == "nearest" and
                                    candidate_match.score <=
                                    region_fallback_max_score + STATISTIC_TIE_TOLERANCE):
                                payload = (domain, region, candidate_screen, candidate_match)
                                entry = (-candidate_match.score, -attempts, payload)
                                if len(fallback_heap) < requested:
                                    heapq.heappush(fallback_heap, entry)
                                elif entry > fallback_heap[0]:
                                    heapq.heapreplace(fallback_heap, entry)
                            sampler.update(len(maxima), attempts, rejection_reasons)
                            continue
                        strict_match_count += 1
                        accept_region_candidate(
                            domain, region, candidate_screen, candidate_match, "strict",
                        )
                        sampler.update(len(maxima), attempts, rejection_reasons)
                        if len(maxima) == requested:
                            stop_sampling = True
                            break
                    candidate_results.close()
                finally:
                    with store._cache_lock:
                        store._active_window_batch = None
            if len(maxima) < requested and region_match_fallback == "nearest":
                ordered_fallback = sorted(fallback_heap, key=lambda item: (-item[0], -item[1]))
                needed = requested - len(maxima)
                for _, _, payload in ordered_fallback[:needed]:
                    (fallback_domain, fallback_region,
                     fallback_screen, fallback_match) = payload
                    accept_region_candidate(
                        fallback_domain, fallback_region, fallback_screen,
                        fallback_match, "nearest_fallback",
                    )
                logger.log(f"REGION FALLBACK target={target.region_id} strict_selected={len(maxima) - fallback_selected_count} "
                           f"fallback_selected={fallback_selected_count} requested={requested} "
                           f"eligible_fallback_candidates={len(fallback_heap)} "
                           f"max_selected_score={max(selected_match_scores) if selected_match_scores else 'NA'}")
            if len(maxima) != requested:
                sampler.update(len(maxima), attempts, rejection_reasons, force=True)
                fallback_note = (
                    "strict thresholds are never relaxed automatically"
                    if region_match_fallback == "none"
                    else f"nearest fallback found too few candidates with score <= {region_fallback_max_score}"
                )
                raise ValueError(f"target {target.region_id}: accepted {len(maxima)}/{requested} background regions after {attempts} attempts; block count/entropy/QC matching rejected {rejected} candidates (entropy rejected {entropy_rejected}); reasons={dict(rejection_reasons)}; {fallback_note}; increase --background-region-attempts/--max-attempt-factor, explicitly enable/expand the nearest fallback, or revise BED/QC")
            sampler.update(len(maxima), attempts, rejection_reasons, force=True)
            summaries.append(dict(base, **observed,
                region_max_empirical_p=empirical_tail_p(float(observed["max_dprime"]), sorted(maxima)),
                region_breadth_empirical_p=(
                    empirical_tail_p(float(observed["significant_pair_fraction"]), sorted(fractions))
                    if background_pair_p_enabled else None
                ),
                n_background_regions=len(maxima), background_region_attempts=attempts,
                background_region_qc_rejections=rejected - block_count_rejected - entropy_rejected,
                background_region_snp_count_rejections=0,
                background_region_block_count_rejections=block_count_rejected,
                background_region_entropy_rejections=entropy_rejected,
                background_region_rejections=rejected,
                region_match_fallback=region_match_fallback,
                region_fallback_max_score=region_fallback_max_score,
                background_region_attempt_limit=attempt_limit,
                region_strict_match_count=strict_match_count,
                region_fallback_candidate_count=fallback_candidate_count,
                region_fallback_selected_count=fallback_selected_count,
                selected_region_match_score_max=max(selected_match_scores),
                status="ok",
            ))
            region_sampling_info[target.region_id] = {"accepted": len(maxima), "attempts": attempts,
                "rejections": rejected, "qc_rejections": rejected - block_count_rejected - entropy_rejected,
                "snp_count_rejections": 0, "block_count_rejections": block_count_rejected,
                "entropy_rejections": entropy_rejected,
                "region_match_fallback": region_match_fallback,
                "region_fallback_max_score": region_fallback_max_score,
                "attempt_limit": attempt_limit,
                "query_batches": region_query_batches,
                "candidate_batch_size": region_batch_size,
                "strict_match_count": strict_match_count,
                "fallback_candidate_count": fallback_candidate_count,
                "fallback_selected_count": fallback_selected_count,
                "selected_match_score_summary": {
                    "min": min(selected_match_scores),
                    "mean": math.fsum(selected_match_scores) / len(selected_match_scores),
                    "max": max(selected_match_scores),
                },
                "region_entropy_matching_enabled": entropy_matching_enabled,
                "region_entropy_max_wasserstein": entropy_max_wasserstein,
                "target_entropy_profile": (None if target_entropy is None else {
                    "left": target_entropy.left, "right": target_entropy.right,
                    "left_complete_haplotype_counts": target_entropy.left_complete_haplotypes,
                    "right_complete_haplotype_counts": target_entropy.right_complete_haplotypes,
                    "left_ranks": target_entropy.left_ranks, "right_ranks": target_entropy.right_ranks,
                }),
                "accepted_entropy_wasserstein_summary": (
                    None if not accepted_entropy_distances else {
                        "min": min(accepted_entropy_distances),
                        "mean": math.fsum(accepted_entropy_distances) / len(accepted_entropy_distances),
                        "max": max(accepted_entropy_distances),
                    }
                ),
                "rejection_reasons": dict(rejection_reasons), "legal_region_count": available,
                "legal_start_count": sum(item.start_count for item in domains),
                "required_available_left_blocks": required_counts[0], "required_available_right_blocks": required_counts[1],
                "region_block_count_tolerance": block_count_tolerance,
                "region_block_count_allow_swap": block_count_allow_swap,
                "accepted_block_count_match_orientations": dict(accepted_orientations),
                "accepted_available_block_count_counts": dict(sorted(accepted_block_counts.items())),
                "central_length_range_bp": [min_length, max_length],
                "accepted_central_length_counts": dict(sorted(accepted_lengths.items()))}
            logger.log(f"REGION BACKGROUND TOTAL completed_targets={len(region_sampling_info)}/"
                       f"{sum(bool(item) for item in masks)} "
                       f"accepted={sum(int(info['accepted']) for info in region_sampling_info.values())}/"
                       f"{sum(bool(item) for item in masks) * args.background_regions if region_test_enabled else 0}")
    logger.log("OUTPUT START writing summary and successful-run metadata")
    if source is not None:
        logger.log("VCF QUERY DONE " + " ".join(
            f"{key}={value}" for key, value in source.statistics().items()
            if key in ("query_batches", "query_intervals", "records_returned", "retained_snp_observations")))
    write_tsv(output_paths["summary"], SUMMARY_FIELDS, summaries)
    metadata = {
        "version": VERSION, "status": "complete",
        "analysis_semantics_version": ANALYSIS_SEMANTICS_VERSION,
        "background_reuse_signature": reuse_signature_from_args(args),
        "background_reuse": {
            "pair_backgrounds_reused": reuse_pair_backgrounds,
            "region_backgrounds_reused": reuse_region_backgrounds,
            "pair_source": str(output_paths["background_pairs"]),
            "region_source": str(output_paths["background_regions"]),
        },
        "background_threads": background_threads,
        "candidate_batching": {
            "pair_batch_size": pair_batch_size,
            "region_batch_size": region_batch_size,
            "query_rule": "merge all unique block windows in a candidate batch into one indexed VCF query",
            "processing_order": "original sampled order; batch size does not change accepted results",
            "memory_rule": ("batch-local raw candidate blocks are released after filtering except selected "
                            "fallback candidates; pair-specific SNP entropy selection is deferred until "
                            "Dprime calculation"),
        },
        "pair_test_enabled": pair_test_enabled,
        "region_test_enabled": region_test_enabled,
        "parameters": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "distance_definition": "equal-width window centers, physical bp; inclusive relative distance band",
        "pair_distance_tolerance": pair_distance_tolerance,
        "pair_sampling": ("uniform unique legal (chromosome, first-window start, second-window start) configurations in each reference-distance band, same chromosome; gap need not be contained"
                          if pair_test_enabled else "skipped_by_user"),
        "region_sampling": ("uniform unique legal (chromosome, central start, central length) configurations; entire envelope contained in merged BED union"
                            if region_test_enabled else "skipped_by_user"),
        "region_candidate_screening": (
            "two-stage: candidates use no Dprime and no pair-specific entropy; "
            "marginal block entropy is also skipped when region entropy matching is disabled; "
            "pair-specific SNP entropy selection is deferred until full Dprime rows, which "
            "are calculated only for selected strict/fallback regions"
            if region_test_enabled else "skipped_by_user"
        ),
        "region_screen_pair_eligibility": ("both windows pass QC/SNP minimum; the complete cohort is "
                                           "fixed across every QC-retained candidate SNP in both windows; "
                                           "joint complete haplotypes meet the absolute and fractional "
                                           "minimum; at least two full-candidate haplotypes occur on each side"),
        "pair_snp_selection_mode": pair_snp_selection_mode,
        "min_pair_complete_haplotype_fraction": min_pair_complete_fraction,
        "required_complete_haplotype_count": required_complete_haplotypes,
        "pair_complete_haplotype_denominator": "2 * VCF sample_count (diploid haplotype chromosomes)",
        "pair_specific_selection_rule": ("fix the pair-complete haplotype cohort using every QC-retained "
                                         "candidate SNP, screen "
                                         "bilateral polymorphism on that cohort, then maximize SNP-subset "
                                         "entropy on the same cohort only when Dprime is calculated"),
        "region_qc_matching": (
            (("available_block_counts_and_marginal_entropy_distributions"
              if entropy_matching_enabled else "available_block_counts_with_tolerance_and_optional_swap")
             + "; each region uses all its own QC-valid cross-flank pairs")
            if region_test_enabled else "not_run"
        ),
        "region_entropy_matching_enabled": region_test_enabled and entropy_matching_enabled,
        "region_entropy_max_wasserstein": entropy_max_wasserstein,
        "region_match_fallback": region_match_fallback,
        "region_fallback_max_score": region_fallback_max_score,
        "background_region_attempts": configured_region_attempts,
        "region_fallback_score_rule": "min over direct/swapped of max(normalized length delta, per-side available-block count error, and when enabled per-side entropy W1); all covariates use the same orientation",
        "region_fallback_outcome_blinding": "candidate selection never uses Dprime or pair p-values",
        "region_entropy_matching_rule": "for each count-valid orientation, compute per-side W1; minimum max(left W1,right W1) <= threshold + statistic_tie_absolute_tolerance; counts and entropy use the same orientation",
        "block_entropy_definition": "-sum(p_h * ln(p_h)), nats; each available block once; selected SNPs; complete haplotypes within this block only, not partner-conditioned; empirical W1 gives every block equal weight",
        "region_entropy_orientation_rule": "when enabled minimize max W1, then count errors, then prefer direct; when disabled use count-only orientation and do not calculate region entropy/W1",
        "region_block_count_tolerance": block_count_tolerance,
        "region_block_count_allow_swap": block_count_allow_swap,
        "region_block_count_matching_rule": "each aligned side differs by <= tolerance; direct or swapped alignment allowed when enabled; both actual sides must be nonzero",
        "available_block_definition": ("participates in at least one cross-flank pair passing window/SNP "
                                       "QC, the fixed pair-complete cohort floor, and bilateral haplotype "
                                       "polymorphism"),
        "region_length_tolerance_bp": length_tolerance,
        "region_pair_distance_matching": "none",
        "region_rank_pair_mask_matching": False,
        "resolved_snp_limits": {"min_snps_per_block": min_snps, "max_snps_per_block": max_snps},
        "pair_snp_count_matching": pair_snp_matching,
        "region_snp_count_matching": "shared_range_and_entropy_cap",
        "background_region_pair_p_reference": "each pair's actual center distance defines its inclusive relative reference band; optional ordered selected SNP-count matching; missing banks sampled on demand",
        "background_region_pair_p_enabled": background_pair_p_enabled,
        "region_breadth_test_enabled": background_pair_p_enabled,
        "pair_p_interpretation": "nominal for a pair selected as the target maximum; not the haplotype-independence permutation p",
        "empirical_tail": "(1 + count(background >= observed - statistic_tie_absolute_tolerance)) / (B + 1)",
        "statistic_tie_absolute_tolerance": STATISTIC_TIE_TOLERANCE,
        "significance_rule": "pair_nominal_p < pair_significance_threshold",
        "independence_note": "pairs share blocks; sampled regions may overlap; neither is assumed independent",
        "parallelism": {
            "implementation": "bounded Python thread pools",
            "workers": background_threads,
            "parallel_units": (["pair_reference_banks"] if pair_test_enabled else []) +
                              (["background_region_candidates"] if region_test_enabled else []),
            "determinism": "task-specific stable RNG seeds and ordered result collection",
            "output_writes": "main_thread_only",
        },
        "sample_count": data.sample_count,
        "vcf_qc": {
            "maf_missing_qc_skipped": vcf_prefiltered,
            "site_maf_filter_enabled": not vcf_prefiltered,
            "site_missing_filter_enabled": not vcf_prefiltered,
            "window_missing_qc_status": window_missing_qc_status,
            "filtered_maf_count": (None if vcf_prefiltered else
                                   source.qc_counts["filtered_maf_count"] if source else data.filtered_maf_count),
            "filtered_site_missing_count": (None if vcf_prefiltered else
                                            source.qc_counts["filtered_call_rate_count"] if source else data.filtered_call_rate_count),
            "counter_scope": "queried observations, possibly repeated" if source else "preloaded BED/target scope",
            "rate_calculations": "skipped" if vcf_prefiltered else "computed_from_GT",
            "prefiltered_note": "external pre-filter SNP counts and removed sites are unavailable" if vcf_prefiltered else None,
        },
        "retained_variant_count": None if source else sum(len(items) for items in data.variants_by_chrom.values()),
        "vcf_access": dict(source.statistics() if source else dict(mode="preload"),
                           window_cache_size=store.window_cache_size,
                           cached_windows=len(store.window_cache), window_cache_hits=store.window_cache_hits),
        "allowed_background_intervals": [dict(chrom=item.display_chrom, start=item.start, end=item.end) for item in background],
        "pair_sampling_counts": pair_sampling_info, "region_sampling_counts": region_sampling_info,
        "logging": {"file": str(output_paths["log"]), "interval_seconds": logger.interval_seconds,
                    "every_attempts": logger.every_attempts,
                    "sampling_eta": "per current stratum/target, estimated from accepted samples per elapsed second"},
        "outputs": {key: str(path) for key, path in output_paths.items()},
    }
    with output_paths["run"].open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return {"region_count": len(targets),
            "ok_region_count": sum(str(row["status"]).startswith("ok") for row in summaries),
            "target_pair_count": sum(len(rows) for rows in target_rows),
            "background_pair_count": sum(len(values) for values in pair_backgrounds.values()),
            "background_region_count": sum(int(row.get("n_background_regions", 0)) for row in summaries)}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("regions_bed", type=Path, help="target BED3/BED4; IDs must be unique")
    parser.add_argument("background_bed", type=Path, help="allowed genomic background BED3/BED4")
    parser.add_argument("vcf", type=Path, help="same-cohort indexed BGZF VCF/BCF; plain VCF supported only in preload mode")
    parser.add_argument("--output-prefix", "-o", type=Path, required=True)
    parser.add_argument("--log-file", type=Path, default=None, help="stderr is also written here (default: OUTPUT_PREFIX.log; no overwrite)")
    parser.add_argument("--log-interval-seconds", type=float, default=30.0,
                        help="time interval for progress checks at processing boundaries (default: 30)")
    parser.add_argument("--log-every-attempts", type=int, default=1000,
                        help="also report sampling after this many attempts, including rejections (default: 1000)")
    parser.add_argument("--threads", "--background-threads", dest="threads", type=int, default=4,
                        help="background sampling worker threads (default: 4; 1 = sequential)")
    parser.add_argument("--vcf-access", choices=("bcftools", "preload"), default="bcftools",
                        help="default: bcftools indexed on-demand queries; preload explicitly enables the old full scan")
    parser.add_argument("--bcftools", default="bcftools", help="bcftools executable name or path")
    parser.add_argument("--vcf-window-cache-size", type=int, default=8192,
                        help="maximum cached candidate blocks; must be >= 2 * block-count (default: 8192)")
    parser.add_argument("--window-bp", type=int, default=10_000)
    parser.add_argument("--window-step-bp", type=int, default=None, help="default: window width; smaller creates sliding windows")
    parser.add_argument("--block-count", type=int, default=10, help="blocks per flank (default: 10)")
    parser.add_argument("--max-flank-bp", type=int, default=100_000)
    parser.add_argument("--snps-per-block", type=int, default=None, help="legacy fixed k; cannot combine with min/max")
    parser.add_argument("--min-snps-per-block", type=int, default=None, help="minimum eligible SNP count (default: 10)")
    parser.add_argument("--max-snps-per-block", type=int, default=None, help="upper entropy cap; keep all SNPs within range (default: 10)")
    parser.add_argument("--background-pairs-per-distance", "--background-pairs-per-stratum", dest="background_pairs_per_distance",
                        type=int, default=5000, help="accepted pairs per reference distance (or per distance/SNP-count stratum when matching is enabled)")
    parser.add_argument("--pair-batch-size", type=int, default=1000,
                        help="candidate block pairs fetched per indexed VCF batch (default: 1000)")
    parser.add_argument("--pair-distance-tolerance", type=float, default=0.10,
                        help="relative center-distance tolerance for pair backgrounds (default: 0.10 = +/-10%%; 0 = exact)")
    parser.add_argument("--pair-match-snp-counts", action="store_true",
                        help="match ordered selected SNP counts in pair backgrounds (default: no count matching; SNP-range QC still applies)")
    parser.add_argument(
        "--skip-pair-test", action="store_true",
        help="skip pair background sampling and pair empirical p-values; retain raw target pair Dprime",
    )
    parser.add_argument("--background-regions", type=int, default=1000, help="accepted backgrounds per target")
    parser.add_argument("--region-batch-size", type=int, default=1000,
                        help="candidate regions fetched per indexed VCF batch (default: 1000)")
    parser.add_argument("--background-region-attempts", type=int, default=None,
                        help="absolute legal region candidates to try per target; default: background-regions * max-attempt-factor")
    parser.add_argument("--region-match-fallback", choices=("none", "nearest"), default="nearest",
                        help="when strict matches are insufficient, fill from nearest covariate matches (default: nearest)")
    parser.add_argument("--region-fallback-max-score", type=float, default=2.0,
                        help="maximum normalized score admitted by nearest fallback (default: 2; must be >= 1)")
    parser.add_argument("--region-length-tolerance-bp", type=int, default=10_000,
                        help="maximum absolute central-length difference from target (default: 10000; 0 = exact)")
    parser.add_argument("--region-block-count-tolerance", type=int, default=1,
                        help="maximum available-block count difference on each aligned side (default: 1)")
    parser.add_argument("--no-region-block-count-swap", dest="region_block_count_allow_swap",
                        action="store_false", default=True,
                        help="disable left/right swapping for count and entropy matching (coordinates unchanged)")
    parser.add_argument("--region-entropy-max-wasserstein", type=float, default=0.05,
                        help="maximum W1 distance on each aligned side, in entropy nats (default: 0.05; 0 = exact)")
    parser.add_argument(
        "--skip-region-test", action="store_true",
        help="skip matched-region background sampling and region empirical tests",
    )
    parser.add_argument("--skip-region-entropy-matching", action="store_true",
                        help="disable region entropy filtering and entropy/W1 calculation; retain count controls")
    parser.add_argument(
        "--skip-background-region-pair-p", action="store_true",
        help=("skip pair nominal p-values in background regions and disable the region breadth test; "
              "retain target pair p-values and region max D' test"),
    )
    parser.add_argument("--pair-significance-threshold", type=float, default=0.01)
    parser.add_argument("--max-attempt-factor", type=int, default=20, help="QC rejection attempt limit / requested count")
    parser.add_argument("--min-complete-haplotypes", type=int, default=2)
    parser.add_argument(
        "--min-pair-complete-haplotype-fraction", type=float, default=0.20,
        help=("minimum fraction of the cohort's diploid haplotype chromosomes that must be "
              "complete across all candidate SNPs in a block pair (default: 0.20)"),
    )
    parser.add_argument("--max-exact-entropy-combinations", type=int, default=100_000)
    parser.add_argument("--min-maf", type=float, default=0.1)
    parser.add_argument("--max-site-missing-rate", type=float, default=0.05)
    parser.add_argument("--max-window-missing-filtered-fraction", type=float, default=0.10)
    parser.add_argument("--vcf-prefiltered", "--skip-vcf-maf-missing-qc", dest="vcf_prefiltered",
                        action="store_true", help="trust prior site MAF/missing filtering; skip recalculation/filtering "
                        "and disable unassessable window missing-filtered-fraction QC; genotype reading remains")
    parser.add_argument("--include-filtered", action="store_true")
    parser.add_argument("--allow-unphased-missing", action="store_true")
    parser.add_argument("--include-target-in-background", action="store_true", help="disable default target+flank exclusion")
    parser.add_argument("--seed", type=int, default=1)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        resolve_snp_limits(args)
    except ValueError as error:
        parser.error(str(error))
    for name in ("window_bp", "block_count", "max_flank_bp",
                 "background_pairs_per_distance", "background_regions", "max_attempt_factor",
                 "max_exact_entropy_combinations", "log_every_attempts", "threads",
                 "pair_batch_size", "region_batch_size"):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be >= 1")
    if not math.isfinite(args.log_interval_seconds) or args.log_interval_seconds <= 0:
        parser.error("--log-interval-seconds must be finite and > 0")
    if args.vcf_window_cache_size < 2 * args.block_count:
        parser.error("--vcf-window-cache-size must be >= 2 * --block-count")
    if args.region_length_tolerance_bp < 0:
        parser.error("--region-length-tolerance-bp must be >= 0")
    if args.region_block_count_tolerance < 0:
        parser.error("--region-block-count-tolerance must be >= 0")
    if (not args.skip_region_test and args.background_region_attempts is not None and
            args.background_region_attempts < args.background_regions):
        parser.error("--background-region-attempts must be >= --background-regions")
    if not math.isfinite(args.region_fallback_max_score) or args.region_fallback_max_score < 1:
        parser.error("--region-fallback-max-score must be finite and >= 1")
    if not math.isfinite(args.region_entropy_max_wasserstein) or args.region_entropy_max_wasserstein < 0:
        parser.error("--region-entropy-max-wasserstein must be finite and >= 0")
    if not math.isfinite(args.pair_distance_tolerance) or not 0 <= args.pair_distance_tolerance < 1:
        parser.error("--pair-distance-tolerance must be finite and in [0, 1)")
    if args.window_step_bp is not None and args.window_step_bp < 1:
        parser.error("--window-step-bp must be >= 1")
    span = args.window_bp + (args.block_count - 1) * (args.window_step_bp or args.window_bp)
    if span > args.max_flank_bp:
        parser.error("requested window span exceeds --max-flank-bp")
    if args.min_complete_haplotypes < 2:
        parser.error("--min-complete-haplotypes must be >= 2")
    if (not math.isfinite(args.min_pair_complete_haplotype_fraction) or
            not 0 <= args.min_pair_complete_haplotype_fraction <= 1):
        parser.error("--min-pair-complete-haplotype-fraction must be finite and in [0, 1]")
    if not 0 <= args.min_maf <= 0.5:
        parser.error("--min-maf must be in [0, 0.5]")
    for name in ("max_site_missing_rate", "max_window_missing_filtered_fraction"):
        if not 0 <= getattr(args, name) <= 1:
            parser.error(f"--{name.replace('_', '-')} must be in [0, 1]")
    if not 0 < args.pair_significance_threshold < 1:
        parser.error("--pair-significance-threshold must be in (0, 1)")
    try:
        run(args)
    except (OSError, ValueError) as error:
        progress(f"error: {error}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
