#!/usr/bin/env python3
"""Maximum cumulative cross-flank median r² over independent physical windows.

For each (i,j), use all SNPs in the i boundary-nearest left windows and j
boundary-nearest right windows. Each side must have enough eligible SNPs.
The statistic is the maximum median r² over all admissible (i,j) pairs.
Every selected background undergoes the same search.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import math
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

try:
    from .flanking_block_ld import Region, VCFData, Variant, pairwise_r2, read_regions, read_vcf
    from .indexed_vcf import BCFtoolsVCFSource
    from .median_flanking_r2 import QuietLogger, draw_candidates, extended_background_territory, flank_bounds, flank_intervals, full_flanks, merge_territory, position_index, read_for_regions, variants_in
except ImportError:
    from flanking_block_ld import Region, VCFData, Variant, pairwise_r2, read_regions, read_vcf  # type: ignore
    from indexed_vcf import BCFtoolsVCFSource  # type: ignore
    from median_flanking_r2 import QuietLogger, draw_candidates, extended_background_territory, flank_bounds, flank_intervals, full_flanks, merge_territory, position_index, read_for_regions, variants_in  # type: ignore


@dataclass(frozen=True)
class Profile:
    region: Region
    left_counts: tuple[int, ...]
    right_counts: tuple[int, ...]

    @property
    def left_total(self) -> int:
        return sum(self.left_counts)

    @property
    def right_total(self) -> int:
        return sum(self.right_counts)

    @property
    def left_nonempty(self) -> int:
        return sum(count > 0 for count in self.left_counts)

    @property
    def right_nonempty(self) -> int:
        return sum(count > 0 for count in self.right_counts)


@dataclass(frozen=True)
class Match:
    strict: bool
    score: float
    orientation: str


@dataclass(frozen=True)
class Step:
    left_rank: int
    right_rank: int
    left_flank_bp: int
    right_flank_bp: int
    left_snps: int
    right_snps: int
    valid_pairs: int
    median_r2: Optional[float]
    eligible_for_max: bool


@dataclass(frozen=True)
class Result:
    region: Region
    profile: Profile
    steps: tuple[Step, ...]
    maximum: Optional[float]
    best_left_rank: Optional[int]
    best_right_rank: Optional[int]
    threshold_met: bool


def window_rank(variant: Variant, region: Region, width: int, side: str) -> int:
    if side == "left":
        return (region.start - 1 - variant.bed_pos) // width
    return (variant.bed_pos - region.end) // width


def profile_for(region: Region, data: VCFData, positions: dict[str, tuple[int, ...]],
                flank_bp: int, window_bp: int) -> Profile:
    count = flank_bp // window_bp
    left, right = variants_in(data, positions, region, flank_bp)
    left_counts, right_counts = [0] * count, [0] * count
    for variant in left:
        left_counts[window_rank(variant, region, window_bp, "left")] += 1
    for variant in right:
        right_counts[window_rank(variant, region, window_bp, "right")] += 1
    return Profile(region, tuple(left_counts), tuple(right_counts))


def match_profile(candidate: Profile, target: Profile, tolerance: float) -> Match:
    options = []
    for orientation, left, right in (
        ("direct", candidate.left_counts, candidate.right_counts),
        ("swapped", candidate.right_counts, candidate.left_counts),
    ):
        left_total, right_total = sum(left), sum(right)
        left_active = sum(value > 0 for value in left)
        right_active = sum(value > 0 for value in right)
        strict = (left_active == target.left_nonempty and right_active == target.right_nonempty
                  and abs(left_total - target.left_total) <= math.ceil(tolerance * target.left_total)
                  and abs(right_total - target.right_total) <= math.ceil(tolerance * target.right_total))
        count_error = (abs(left_total - target.left_total) / max(1, target.left_total)
                       + abs(right_total - target.right_total) / max(1, target.right_total))
        active_error = (abs(left_active - target.left_nonempty)
                        + abs(right_active - target.right_nonempty)) / len(left)
        shape_error = (sum(abs(a - b) for a, b in zip(left, target.left_counts)) / max(1, target.left_total)
                       + sum(abs(a - b) for a, b in zip(right, target.right_counts)) / max(1, target.right_total))
        score = count_error + active_error + 0.25 * shape_error
        options.append((not strict, score, 0 if orientation == "direct" else 1, orientation))
    not_strict, score, _, orientation = min(options)
    return Match(not not_strict, score, orientation)


def calculate(region: Region, data: VCFData, positions: dict[str, tuple[int, ...]],
              flank_bp: int, window_bp: int, min_pair_haplotypes: int,
              min_pairs: int, min_snps_per_side: int) -> Result:
    profile = profile_for(region, data, positions, flank_bp, window_bp)
    left, right = variants_in(data, positions, region, flank_bp)
    n_windows = flank_bp // window_bp
    buckets: list[list[list[float]]] = [[[] for _ in range(n_windows)] for _ in range(n_windows)]
    for first in left:
        left_rank = window_rank(first, region, window_bp, "left")
        for second in right:
            right_rank = window_rank(second, region, window_bp, "right")
            keys = first.alleles.keys() & second.alleles.keys()
            if len(keys) >= min_pair_haplotypes:
                value = pairwise_r2(first, second, tuple(keys))
                if value is not None:
                    buckets[left_rank][right_rank].append(value)
    # For each left prefix, grow the right prefix using two median heaps.
    # Each SNP pair's r² is calculated only once, then reused across prefixes.
    steps = []
    best_value: Optional[float] = None
    best_left_rank: Optional[int] = None
    best_right_rank: Optional[int] = None
    descriptive_best: Optional[tuple[float, int, int]] = None
    left_count = 0
    for left_rank in range(n_windows):
        left_count += profile.left_counts[left_rank]
        right_count = valid_pairs = 0
        lower: list[float] = []
        upper: list[float] = []
        for right_rank in range(n_windows):
            right_count += profile.right_counts[right_rank]
            for left_bucket_rank in range(left_rank + 1):
                for value in buckets[left_bucket_rank][right_rank]:
                    if not lower or value <= -lower[0]:
                        heapq.heappush(lower, -value)
                    else:
                        heapq.heappush(upper, value)
                    if len(lower) > len(upper) + 1:
                        heapq.heappush(upper, -heapq.heappop(lower))
                    elif len(upper) > len(lower):
                        heapq.heappush(lower, -heapq.heappop(upper))
                    valid_pairs += 1
            median = None
            eligible = (left_count >= min_snps_per_side and right_count >= min_snps_per_side
                        and valid_pairs >= min_pairs)
            if valid_pairs:
                median = -lower[0] if len(lower) > len(upper) else (-lower[0] + upper[0]) / 2
                choice = (median, left_rank + 1, right_rank + 1)
                if (descriptive_best is None or median > descriptive_best[0] or
                        (median == descriptive_best[0] and
                         (left_rank + right_rank, left_rank) <
                         (descriptive_best[1] + descriptive_best[2] - 2,
                          descriptive_best[1] - 1))):
                    descriptive_best = choice
            if eligible:
                # In a tie, use the smaller combined flank width, then left width.
                if (best_value is None or median > best_value or
                        (median == best_value and
                         ((left_rank + right_rank, left_rank, right_rank) <
                          (best_left_rank + best_right_rank - 2,
                           best_left_rank - 1, best_right_rank - 1)))):
                    best_value, best_left_rank, best_right_rank = median, left_rank + 1, right_rank + 1
            steps.append(Step(left_rank + 1, right_rank + 1,
                              (left_rank + 1) * window_bp, (right_rank + 1) * window_bp,
                              left_count, right_count, valid_pairs, median, eligible))
    threshold_met = best_value is not None
    if not threshold_met and descriptive_best is not None:
        best_value, best_left_rank, best_right_rank = descriptive_best
    return Result(region, profile, tuple(steps), best_value, best_left_rank,
                  best_right_rank, threshold_met)


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_bed", type=Path, help="BED3/BED4 with exactly one interval")
    parser.add_argument("background_bed", type=Path, help="allowed background territory BED")
    parser.add_argument("vcf", type=Path, help="phased VCF/VCF.gz")
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--flank-bp", type=int, default=100_000)
    parser.add_argument("--window-bp", type=int, default=5_000)
    parser.add_argument("--min-maf", type=float, default=0.05)
    parser.add_argument("--max-site-missing-rate", type=float, default=0.05)
    parser.add_argument("--min-pair-haplotypes", type=int, default=4)
    parser.add_argument("--min-pairs", type=int, default=10,
                        help="minimum valid cross-side SNP pairs for a prefix to enter the maximum search (default: 10)")
    parser.add_argument("--min-snps-per-side", type=int, default=5,
                        help="minimum eligible SNPs in EACH cumulative flank (default: 5)")
    parser.add_argument("--snp-count-tolerance", type=float, default=0.10,
                        help="strict total SNP count tolerance on each side; nonempty window counts must match exactly")
    parser.add_argument("--length-tolerance-bp", type=int, default=10_000)
    parser.add_argument("--attempts", type=int, default=20_000)
    parser.add_argument("--backgrounds", type=int, default=1_000)
    parser.add_argument("--allow-background-flank-extension", action="store_true",
                        help="allow background flanks to extend beyond background BED intervals; "
                             "do not require complete background flanks for empirical P")
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--include-filtered", action="store_true")
    parser.add_argument("--vcf-access", choices=("indexed", "preload"), default="indexed")
    parser.add_argument("--bcftools", default="bcftools")
    args = parser.parse_args(argv)
    if (args.window_bp < 1 or args.flank_bp < 1 or args.flank_bp % args.window_bp
            or args.attempts < 1 or args.backgrounds < 1 or args.batch_size < 1
            or args.min_pairs < 1 or args.min_pair_haplotypes < 2 or args.min_snps_per_side < 1
            or args.length_tolerance_bp < 0 or args.snp_count_tolerance < 0
            or not 0 <= args.min_maf < 0.5 or not 0 <= args.max_site_missing_rate < 1):
        parser.error("invalid arguments: positive windows/counts, divisible flank, and valid QC thresholds required")
    targets = read_regions(args.target_bed)
    if len(targets) != 1:
        parser.error("target BED must contain exactly one interval")
    target_region = targets[0]
    target_flanks_complete = full_flanks(target_region, args.flank_bp)
    territory = merge_territory(read_regions(args.background_bed))
    candidates = draw_candidates(target_region, territory, args.flank_bp,
                                 args.length_tolerance_bp, args.attempts, args.seed,
                                 args.allow_background_flank_extension)
    source = None
    preloaded = None
    if args.vcf_access == "indexed":
        args.allow_unphased_missing = False
        args.vcf_prefiltered = False
        source = BCFtoolsVCFSource(args, QuietLogger())
    else:
        background_intervals = (extended_background_territory(territory, args.flank_bp)
                                if args.allow_background_flank_extension else territory)
        all_intervals = [*flank_intervals(target_region, args.flank_bp), *background_intervals]
        preloaded = read_vcf(args.vcf, all_intervals, 0, args.min_maf,
                             1 - args.max_site_missing_rate, args.include_filtered,
                             False, include_region_interiors=True)
    target_data = read_for_regions(args, [target_region], source) if source else preloaded
    assert target_data is not None
    preloaded_positions = position_index(preloaded) if preloaded else None
    target = calculate(target_region, target_data,
                       preloaded_positions if preloaded_positions is not None else position_index(target_data),
                       args.flank_bp, args.window_bp, args.min_pair_haplotypes, args.min_pairs,
                       args.min_snps_per_side)
    if target.maximum is None:
        parser.error("target has no valid cross-flank SNP pair after filtering")
    screened = []
    for offset in range(0, len(candidates) if target.threshold_met else 0,
                        args.batch_size):
        batch = candidates[offset:offset + args.batch_size]
        data = read_for_regions(args, batch, source) if source else preloaded
        assert data is not None
        positions = preloaded_positions if preloaded_positions is not None else position_index(data)
        for region in batch:
            profile = profile_for(region, data, positions, args.flank_bp, args.window_bp)
            if profile.left_total and profile.right_total:
                screened.append(profile)
        print(f"screened {min(offset + len(batch), len(candidates))}/{len(candidates)} candidates",
              file=sys.stderr)
    rng = random.Random(args.seed + 1)
    rng.shuffle(screened)
    matches = {p.region.region_id: match_profile(p, target.profile, args.snp_count_tolerance)
               for p in screened}
    strict = [p for p in screened if matches[p.region.region_id].strict]
    fallback = [p for p in screened if not matches[p.region.region_id].strict]
    fallback.sort(key=lambda p: (matches[p.region.region_id].score, p.region.region_id))
    ordered = [(p, "strict") for p in strict] + [(p, "fallback") for p in fallback]
    background_rows = []
    valid_maxima = []
    for offset in range(0, len(ordered), args.batch_size):
        if len(valid_maxima) >= args.backgrounds:
            break
        batch = ordered[offset:offset + args.batch_size]
        data = read_for_regions(args, [p.region for p, _ in batch], source) if source else preloaded
        assert data is not None
        positions = preloaded_positions if preloaded_positions is not None else position_index(data)
        for profile, selection in batch:
            result = calculate(profile.region, data, positions, args.flank_bp, args.window_bp,
                               args.min_pair_haplotypes, args.min_pairs, args.min_snps_per_side)
            match = matches[profile.region.region_id]
            left_bound, right_bound = flank_bounds(profile.region, args.flank_bp)
            background_rows.append(dict(chrom=profile.region.display_chrom,
                                        start=profile.region.start, end=profile.region.end,
                                        length_bp=profile.region.end - profile.region.start,
                                        left_flank_bp=profile.region.start - left_bound,
                                        right_flank_bp=right_bound - profile.region.end,
                                        full_flanks=full_flanks(profile.region, args.flank_bp),
                                        extended_beyond_background_bed=(
                                            left_bound < profile.region.territory_start or
                                            right_bound > profile.region.territory_end),
                                        left_snps=profile.left_total, right_snps=profile.right_total,
                                        left_nonempty_windows=profile.left_nonempty,
                                        right_nonempty_windows=profile.right_nonempty,
                                        max_prefix_median_r2=result.maximum if result.maximum is not None else "NA",
                                        max_search_qualified=result.threshold_met,
                                        best_left_rank=result.best_left_rank if result.best_left_rank is not None else "NA",
                                        best_right_rank=result.best_right_rank if result.best_right_rank is not None else "NA",
                                        best_left_flank_bp=result.best_left_rank * args.window_bp if result.best_left_rank is not None else "NA",
                                        best_right_flank_bp=result.best_right_rank * args.window_bp if result.best_right_rank is not None else "NA",
                                        selection=selection, match_score=match.score,
                                        match_orientation=match.orientation))
            if result.maximum is not None and result.threshold_met:
                valid_maxima.append(result.maximum)
            if len(valid_maxima) >= args.backgrounds:
                break
        print(f"calculated {len(valid_maxima)}/{args.backgrounds} valid backgrounds; "
              f"tried {len(background_rows)}", file=sys.stderr)
    output = args.output_prefix
    output.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(Path(str(output) + ".profile.tsv"),
              ["target_id", "left_rank", "right_rank", "left_flank_bp", "right_flank_bp",
               "left_snps", "right_snps", "valid_pairs", "median_r2", "eligible_for_max", "is_maximum"],
              [dict(target_id=target.region.region_id, left_rank=step.left_rank,
                    right_rank=step.right_rank, left_flank_bp=step.left_flank_bp,
                    right_flank_bp=step.right_flank_bp,
                    left_snps=step.left_snps, right_snps=step.right_snps,
                    valid_pairs=step.valid_pairs,
                    median_r2=step.median_r2 if step.median_r2 is not None else "NA",
                    eligible_for_max=step.eligible_for_max,
                    is_maximum=(step.left_rank == target.best_left_rank and
                                step.right_rank == target.best_right_rank))
               for step in target.steps])
    write_tsv(Path(str(output) + ".backgrounds.tsv"),
              ["chrom", "start", "end", "length_bp", "left_flank_bp", "right_flank_bp",
               "full_flanks", "extended_beyond_background_bed", "left_snps", "right_snps",
               "left_nonempty_windows", "right_nonempty_windows", "max_prefix_median_r2",
               "max_search_qualified",
               "best_left_rank", "best_right_rank", "best_left_flank_bp", "best_right_flank_bp",
               "selection", "match_score", "match_orientation"],
              background_rows)
    selected_partial = sum(row["full_flanks"] is False and row["max_search_qualified"] is True
                           for row in background_rows)
    background_complete = (target_flanks_complete and target.threshold_met
                           and len(valid_maxima) == args.backgrounds
                           and (args.allow_background_flank_extension or selected_partial == 0))
    empirical_p = ((1 + sum(value >= target.maximum for value in valid_maxima)) /
                   (1 + len(valid_maxima))) if background_complete else "NA"
    write_tsv(Path(str(output) + ".summary.tsv"),
              ["target_id", "chrom", "start", "end", "flank_bp",
               "left_flank_bp", "right_flank_bp", "window_bp", "min_pairs",
               "min_snps_per_side", "max_search_qualified",
               "left_snps", "right_snps", "left_nonempty_windows", "right_nonempty_windows",
               "max_prefix_median_r2", "best_left_rank", "best_right_rank",
               "best_left_flank_bp", "best_right_flank_bp", "empirical_p_upper",
               "attempts", "unique_candidates", "screened_nonempty", "strict_available",
               "selected_strict", "selected_fallback", "selected_swapped", "selected_total",
               "selected_partial_flanks",
               "valid_backgrounds", "allow_background_flank_extension", "status"],
              [dict(target_id=target_region.region_id, chrom=target_region.display_chrom,
                    start=target_region.start, end=target_region.end,
                    flank_bp=args.flank_bp,
                    left_flank_bp=target_region.start - flank_bounds(target_region, args.flank_bp)[0],
                    right_flank_bp=flank_bounds(target_region, args.flank_bp)[1] - target_region.end,
                    window_bp=args.window_bp, min_pairs=args.min_pairs,
                    min_snps_per_side=args.min_snps_per_side,
                    max_search_qualified=target.threshold_met,
                    left_snps=target.profile.left_total, right_snps=target.profile.right_total,
                    left_nonempty_windows=target.profile.left_nonempty,
                    right_nonempty_windows=target.profile.right_nonempty,
                    max_prefix_median_r2=target.maximum,
                    best_left_rank=target.best_left_rank, best_right_rank=target.best_right_rank,
                    best_left_flank_bp=target.best_left_rank * args.window_bp,
                    best_right_flank_bp=target.best_right_rank * args.window_bp,
                    empirical_p_upper=empirical_p,
                    attempts=args.attempts, unique_candidates=len(candidates),
                    screened_nonempty=len(screened), strict_available=len(strict),
                    selected_strict=sum(row["selection"] == "strict" for row in background_rows),
                    selected_fallback=sum(row["selection"] == "fallback" for row in background_rows),
                    selected_swapped=sum(row["match_orientation"] == "swapped" for row in background_rows),
                    selected_total=len(background_rows), selected_partial_flanks=selected_partial,
                    valid_backgrounds=len(valid_maxima),
                    allow_background_flank_extension=args.allow_background_flank_extension,
                    status=("ok" if background_complete else
                            "target_short_flank" if not target_flanks_complete else
                            "target_below_search_threshold" if not target.threshold_met else
                            "no_background_candidates" if not candidates else
                            "insufficient_valid_backgrounds" if len(valid_maxima) < args.backgrounds else
                            "background_short_flanks"))])
    print(f"target max prefix median r2={target.maximum:.6g} at "
          f"left={target.best_left_rank}, right={target.best_right_rank} windows; "
          f"valid backgrounds={len(valid_maxima)}/{args.backgrounds}; empirical upper P={empirical_p}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
