#!/usr/bin/env python3
"""Select a maximum-union-coverage compatible set in one genome BED gap.

The caller determines alignment eligibility (for example, clipped coverage
strictly greater than 10 kb) and groups candidates by population path and gap.
This module optimizes the eligible set exactly rather than greedily choosing
the longest alignment or maximizing summed, double-counted interval lengths.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from functools import cached_property
from typing import List, Optional, Sequence, Tuple

try:
    from .refine_bundle_mappings import PafBlock
except ImportError:  # Support helpers imported with scripts/ on sys.path.
    from refine_bundle_mappings import PafBlock


@dataclass(frozen=True)
class GapAlignment:
    """One raw PAF block intersected with one half-open genome BED gap.

    ``block`` retains both original coordinate axes and strand. ``start`` and
    ``end`` are exactly the block's genome interval clipped to the gap; bundle
    coordinates are deliberately not projected or changed during selection.
    """

    block: PafBlock
    gap_start: int
    gap_end: int
    start: int
    end: int

    def __post_init__(self) -> None:
        if not 0 <= self.gap_start < self.gap_end <= self.block.target_length:
            raise ValueError("gap must be a positive half-open interval within the genome path")
        expected = (
            max(self.block.target_start, self.gap_start),
            min(self.block.target_end, self.gap_end),
        )
        if (self.start, self.end) != expected or self.start >= self.end:
            raise ValueError("start/end must be the positive raw PAF genome interval clipped to its gap")

    @property
    def coverage_bp(self) -> int:
        return self.end - self.start


def _candidate_tie_key(candidate: GapAlignment) -> Tuple[object, ...]:
    block = candidate.block
    return (
        block.line_number, block.raw_query_name, block.strand,
        block.query_length, block.query_start, block.query_end,
        block.target_start, block.target_end, candidate.start, candidate.end,
    )


@dataclass
class _SelectionState:
    candidate: GapAlignment
    previous: Optional["_SelectionState"]
    covered_bp: int
    count: int
    identity_score: float
    mapq_score: int

    @cached_property
    def signature(self) -> Tuple[Tuple[object, ...], ...]:
        # Build iteratively: long compatible paths must not hit recursion limits.
        keys = []
        state: Optional[_SelectionState] = self
        while state is not None:
            keys.append(_candidate_tie_key(state.candidate))
            state = state.previous
        return tuple(reversed(keys))


def _prefer(
    first: Optional[_SelectionState],
    second: Optional[_SelectionState],
    subtract_end: bool = False,
) -> Optional[_SelectionState]:
    if first is None:
        return second
    if second is None:
        return first

    def metrics(state: _SelectionState) -> Tuple[float, ...]:
        coverage = state.covered_bp - (state.candidate.end if subtract_end else 0)
        return coverage, -state.count, state.identity_score, state.mapq_score

    first_metrics, second_metrics = metrics(first), metrics(second)
    if first_metrics != second_metrics:
        return first if first_metrics > second_metrics else second
    return first if first.signature <= second.signature else second


class _OverlapRangeMaximum:
    """Range maximum of DP coverage minus last genome endpoint."""

    def __init__(self, count: int) -> None:
        self.leaf_count = 1
        while self.leaf_count < count:
            self.leaf_count *= 2
        self.values: List[Optional[_SelectionState]] = [None] * (2 * self.leaf_count)

    def update(self, index: int, state: _SelectionState) -> None:
        index += self.leaf_count
        self.values[index] = state
        while index > 1:
            index //= 2
            self.values[index] = _prefer(
                self.values[2 * index], self.values[2 * index + 1], subtract_end=True
            )

    def query(self, left: int, right: int) -> Optional[_SelectionState]:
        best = None
        left += self.leaf_count
        right += self.leaf_count
        while left < right:
            if left % 2:
                best = _prefer(best, self.values[left], subtract_end=True)
                left += 1
            if right % 2:
                right -= 1
                best = _prefer(best, self.values[right], subtract_end=True)
            left //= 2
            right //= 2
        return best


def maximum_coverage_selection(
    candidates: Sequence[GapAlignment], max_overlap_bp: int = 1000,
) -> List[GapAlignment]:
    """Return a compatible set maximizing actual genome interval union bp.

    All candidates must refer to the same population path and the same gap.
    Eligibility/minimum covered length is the caller's responsibility. Each
    clipped interval must have length >= ``max_overlap_bp``; otherwise the
    caller must reduce the threshold. This ensures any compatible end-ordered
    chain has strictly advancing starts, so checking adjacent intervals also
    guarantees compatibility of every pair, including nonadjacent intervals.

    With a positive threshold, every pair must overlap strictly less than that
    many bp: 999 bp passes the default, 1000 bp does not. A zero threshold means
    nonoverlap, allowing zero-length intersections and touching boundaries.

    DP transition from j to i adds end_i - max(start_i, end_j), not the full
    length of i when they overlap. Nonoverlapping predecessors use prefix DP
    maxima; overlapping predecessors use a range maximum of DP[j] - end_j.
    Sorting/bisect and a segment tree give O(n log n) score optimization and
    O(n) state storage. Exact score ties may require linear path-signature
    comparisons/caching, so adversarial all-tie cases can exceed that bound.

    Equal union coverage prefers fewer alignments, then larger sums of
    clipped-bp-weighted identity and MAPQ (additive total scores, not averages),
    then smaller source-line/coordinate signatures. Results are genome-end
    ordered and independent of the input order. Raw bundle endpoints and strand
    are retained; reverse bundle direction does not affect genome selection.
    """

    if not isinstance(max_overlap_bp, int) or max_overlap_bp < 0:
        raise ValueError("max_overlap_bp must be a non-negative integer")
    if not candidates:
        return []
    group = (
        candidates[0].block.target_name, candidates[0].gap_start, candidates[0].gap_end
    )
    for candidate in candidates:
        if (candidate.block.target_name, candidate.gap_start, candidate.gap_end) != group:
            raise ValueError("maximum coverage selection requires one population path and one gap")
        if candidate.coverage_bp < max_overlap_bp:
            raise ValueError("candidate clipped length must be >= max_overlap_bp; reduce the threshold")

    ordered = sorted(candidates, key=lambda item: (item.end, item.start, _candidate_tie_key(item)))
    ends = [candidate.end for candidate in ordered]
    prefix_best: List[_SelectionState] = []
    overlap_maximum = _OverlapRangeMaximum(len(ordered))
    for index, candidate in enumerate(ordered):
        length = candidate.coverage_bp
        identity_score = length * candidate.block.identity
        mapq_score = length * candidate.block.mapq
        best = _SelectionState(candidate, None, length, 1, identity_score, mapq_score)

        nonoverlap_limit = bisect_right(ends, candidate.start, 0, index)
        predecessors = []
        if nonoverlap_limit:
            predecessors.append(prefix_best[nonoverlap_limit - 1])
        if max_overlap_bp:
            overlap_limit = bisect_left(ends, candidate.start + max_overlap_bp, 0, index)
            if overlap_limit > nonoverlap_limit:
                predecessor = overlap_maximum.query(nonoverlap_limit, overlap_limit)
                if predecessor is not None:
                    predecessors.append(predecessor)

        for predecessor in predecessors:
            added_bp = candidate.end - max(candidate.start, predecessor.candidate.end)
            extension = _SelectionState(
                candidate, predecessor, predecessor.covered_bp + added_bp,
                predecessor.count + 1, predecessor.identity_score + identity_score,
                predecessor.mapq_score + mapq_score,
            )
            best = _prefer(best, extension)
        prefix_best.append(_prefer(prefix_best[-1] if prefix_best else None, best))
        overlap_maximum.update(index, best)

    result = []
    state: Optional[_SelectionState] = prefix_best[-1]
    while state is not None:
        result.append(state.candidate)
        state = state.previous
    return list(reversed(result))
