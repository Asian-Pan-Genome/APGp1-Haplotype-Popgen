#!/usr/bin/env python3
"""Median cross-flank SNP r² for one BED interval and matched background gaps.

Background BED intervals define allowed central-gap sampling territory.
By default flanks stop at BED boundaries; --allow-background-flank-extension
allows them to extend outside. BED coordinates are 0-based, half-open.
Genotypes must be phased; unphased heterozygous sites are excluded.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
import sys
from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

try:
    from .flanking_block_ld import Region, VCFData, Variant, pairwise_r2, read_regions, read_vcf
    from .indexed_vcf import BCFtoolsVCFSource
except ImportError:
    from flanking_block_ld import Region, VCFData, Variant, pairwise_r2, read_regions, read_vcf  # type: ignore
    from indexed_vcf import BCFtoolsVCFSource  # type: ignore


@dataclass(frozen=True)
class Candidate:
    region: Region
    left_count: int = 0
    right_count: int = 0
    median_r2: Optional[float] = None
    valid_pairs: int = 0
    selection: str = ""
    match_score: float = math.inf
    match_orientation: str = ""
    max_r2: Optional[float] = None


@dataclass(frozen=True)
class BackgroundRegion(Region):
    territory_start: int = 0
    territory_end: int = 0
    allow_flank_extension: bool = False


def flank_bounds(region: Region, width: int) -> tuple[int, int]:
    if getattr(region, "allow_flank_extension", False):
        return max(0, region.start - width), region.end + width
    left = max(0, region.start - width, getattr(region, "territory_start", 0))
    right = min(region.end + width, getattr(region, "territory_end", region.end + width))
    return left, right


def full_flanks(region: Region, width: int) -> bool:
    left, right = flank_bounds(region, width)
    return region.start - left == width and right - region.end == width


class QuietLogger:
    def log(self, message: str, level: str = "INFO") -> None:
        if level == "WARNING":
            print(f"WARNING: {message}", file=sys.stderr)

    def periodic(self, key: str, message: str) -> None:
        pass


def flank_intervals(region: Region, width: int) -> tuple[Region, Region]:
    left_bound, right_bound = flank_bounds(region, width)
    return (
        Region(region.chrom, region.display_chrom, left_bound, region.start,
               region.region_id + ":left", 0),
        Region(region.chrom, region.display_chrom, region.end, right_bound,
               region.region_id + ":right", 0),
    )


def position_index(data: VCFData) -> dict[str, tuple[int, ...]]:
    return {chrom: tuple(variant.bed_pos for variant in variants)
            for chrom, variants in data.variants_by_chrom.items()}


def variants_in(data: VCFData, positions_by_chrom: dict[str, tuple[int, ...]],
                region: Region, width: int) -> tuple[tuple[Variant, ...], tuple[Variant, ...]]:
    variants = data.variants_by_chrom.get(region.chrom, ())
    positions = positions_by_chrom.get(region.chrom, ())
    left_bound, right_bound = flank_bounds(region, width)
    left = variants[bisect_left(positions, left_bound):bisect_left(positions, region.start)]
    right = variants[bisect_left(positions, region.end):bisect_left(positions, right_bound)]
    return left, right


def measure(region: Region, data: VCFData, positions_by_chrom: dict[str, tuple[int, ...]],
            width: int, min_pair_haplotypes: int) -> Candidate:
    left, right = variants_in(data, positions_by_chrom, region, width)
    values = []
    for first in left:
        for second in right:
            keys = first.alleles.keys() & second.alleles.keys()
            if len(keys) >= min_pair_haplotypes:
                value = pairwise_r2(first, second, tuple(keys))
                if value is not None:
                    values.append(value)
    return Candidate(region, len(left), len(right),
                     statistics.median(values) if values else None,
                     len(values), max_r2=max(values) if values else None)


def merge_territory(regions: Sequence[Region]) -> list[Region]:
    result: list[Region] = []
    for region in sorted(regions, key=lambda r: (r.chrom, r.start, r.end)):
        if result and result[-1].chrom == region.chrom and region.start <= result[-1].end:
            old = result[-1]
            result[-1] = Region(old.chrom, old.display_chrom, old.start, max(old.end, region.end),
                                old.region_id, old.line_number)
        else:
            result.append(region)
    return result


def draw_candidates(target: Region, territory: Sequence[Region], width: int,
                    tolerance: int, attempts: int, seed: int,
                    allow_background_flank_extension: bool = False) -> list[Region]:
    rng = random.Random(seed)
    target_length = target.end - target.start
    minimum_span = max(1, target_length - tolerance) + 2
    domains = [(r, r.end - r.start - minimum_span + 1) for r in territory]
    domains = [(r, n) for r, n in domains if n > 0]
    if not domains:
        return []
    cumulative = []
    total = 0
    for _, count in domains:
        total += count
        cumulative.append(total)
    seen = set()
    candidates = []
    for attempt in range(attempts):
        region, _ = domains[bisect_left(cumulative, rng.randrange(total) + 1)]
        length = rng.randint(max(1, target_length - tolerance), target_length + tolerance)
        start_min = region.start + 1
        start_max = region.end - 1 - length
        if start_min > start_max:
            continue
        start = rng.randint(start_min, start_max)
        end = start + length
        # A target and a background gap must not share the target's flank span.
        if region.chrom == target.chrom and start - width < target.end + width and end + width > target.start - width:
            continue
        key = region.chrom, start, end
        if key in seen:
            continue
        seen.add(key)
        candidates.append(BackgroundRegion(region.chrom, region.display_chrom, start, end,
                                           f"background_{attempt + 1}", 0,
                                           region.start, region.end,
                                           allow_background_flank_extension))
    return candidates


def extended_background_territory(territory: Sequence[Region], width: int) -> list[Region]:
    """Preload every position that an extended background flank could use."""
    return [Region(region.chrom, region.display_chrom, max(0, region.start - width),
                   region.end + width, region.region_id, region.line_number)
            for region in territory]


def read_for_regions(args: argparse.Namespace, regions: Sequence[Region], source: BCFtoolsVCFSource) -> VCFData:
    windows = [window for region in regions for window in flank_intervals(region, args.flank_bp)]
    return source.fetch(windows)


def match_details(candidate: Candidate, target: Candidate, count_tolerance: float,
                  length_tolerance: int) -> tuple[bool, float, str]:
    """Match left/right SNP counts in the better of direct and swapped orientations."""
    length_error = abs((candidate.region.end - candidate.region.start)
                       - (target.region.end - target.region.start)) / max(1, length_tolerance)
    options = []
    for orientation, left_count, right_count in (
        ("direct", candidate.left_count, candidate.right_count),
        ("swapped", candidate.right_count, candidate.left_count),
    ):
        left_delta = abs(left_count - target.left_count)
        right_delta = abs(right_count - target.right_count)
        strict = (left_delta <= math.ceil(count_tolerance * target.left_count)
                  and right_delta <= math.ceil(count_tolerance * target.right_count))
        score = (left_delta / max(1, target.left_count)
                 + right_delta / max(1, target.right_count) + 0.01 * length_error)
        options.append((not strict, score, 0 if orientation == "direct" else 1,
                        orientation))
    non_strict, score, _, orientation = min(options)
    return not non_strict, score, orientation


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_bed", type=Path, help="BED3/BED4 with exactly one target interval")
    parser.add_argument("background_bed", type=Path, help="BED of allowed background territory")
    parser.add_argument("vcf", type=Path, help="phased VCF/VCF.gz with GT")
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--flank-bp", type=int, required=True)
    parser.add_argument("--min-maf", type=float, default=0.05)
    parser.add_argument("--max-site-missing-rate", type=float, default=0.05)
    parser.add_argument("--min-pair-haplotypes", type=int, default=4)
    parser.add_argument("--include-filtered", action="store_true")
    parser.add_argument("--length-tolerance-bp", type=int, default=10_000)
    parser.add_argument("--snp-count-tolerance", type=float, default=0.10,
                        help="strict relative tolerance per flank, rounded up; 0 means exact (default: 0.10)")
    parser.add_argument("--attempts", type=int, default=20_000)
    parser.add_argument("--backgrounds", type=int, default=1_000)
    parser.add_argument("--allow-background-flank-extension", action="store_true",
                        help="allow background flanks to extend beyond background BED intervals; "
                             "do not require complete background flanks for empirical P")
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--vcf-access", choices=("indexed", "preload"), default="indexed",
                        help="indexed uses bcftools and a .csi/.tbi index; preload scans the VCF once")
    parser.add_argument("--bcftools", default="bcftools")
    args = parser.parse_args(argv)
    if args.flank_bp < 1 or args.length_tolerance_bp < 0 or args.attempts < 1 or args.backgrounds < 1 or args.batch_size < 1:
        parser.error("flank-bp, attempts, backgrounds and batch-size must be positive; length tolerance must be nonnegative")
    if not 0 <= args.min_maf < 0.5 or not 0 <= args.max_site_missing_rate < 1 or args.snp_count_tolerance < 0:
        parser.error("invalid MAF, missingness, or SNP-count tolerance")
    if args.min_pair_haplotypes < 2:
        parser.error("--min-pair-haplotypes must be at least 2")
    targets = read_regions(args.target_bed)
    if len(targets) != 1:
        parser.error("target BED must contain exactly one interval")
    target_region = targets[0]
    target_flanks_complete = full_flanks(target_region, args.flank_bp)
    territory = merge_territory(read_regions(args.background_bed))
    regions = draw_candidates(target_region, territory, args.flank_bp,
                              args.length_tolerance_bp, args.attempts, args.seed,
                              args.allow_background_flank_extension)
    source = None
    preloaded = None
    if args.vcf_access == "indexed":
        args.allow_unphased_missing = False
        args.vcf_prefiltered = False
        source = BCFtoolsVCFSource(args, QuietLogger())
    else:
        background_windows = (extended_background_territory(territory, args.flank_bp)
                              if args.allow_background_flank_extension else territory)
        all_windows = [*flank_intervals(target_region, args.flank_bp), *background_windows]
        preloaded = read_vcf(args.vcf, all_windows, 0, args.min_maf,
                             1 - args.max_site_missing_rate, args.include_filtered,
                             False, include_region_interiors=True)
    target_data = read_for_regions(args, [target_region], source) if source else preloaded
    assert target_data is not None
    preload_positions = position_index(preloaded) if preloaded is not None else None
    target = measure(target_region, target_data,
                     preload_positions if preload_positions is not None else position_index(target_data),
                     args.flank_bp, args.min_pair_haplotypes)
    if target.median_r2 is None:
        parser.error("target has no valid cross-flank SNP pair after filtering")
    screened: list[Candidate] = []
    for offset in range(0, len(regions), args.batch_size):
        batch = regions[offset:offset + args.batch_size]
        data = read_for_regions(args, batch, source) if source else preloaded
        assert data is not None
        positions = preload_positions if preload_positions is not None else position_index(data)
        for region in batch:
            left, right = variants_in(data, positions, region, args.flank_bp)
            if left and right:
                screened.append(Candidate(region, len(left), len(right)))
        print(f"screened {min(offset + len(batch), len(regions))}/{len(regions)} candidates", file=sys.stderr)
    rng = random.Random(args.seed + 1)
    rng.shuffle(screened)
    matches = {c.region.region_id: match_details(c, target, args.snp_count_tolerance,
                                                 args.length_tolerance_bp) for c in screened}
    strict = [c for c in screened if matches[c.region.region_id][0]]
    strict_ids = {(c.region.chrom, c.region.start, c.region.end) for c in strict}
    fallback = [c for c in screened if (c.region.chrom, c.region.start, c.region.end) not in strict_ids]
    fallback.sort(key=lambda c: (matches[c.region.region_id][1], c.region.region_id))
    ordered = [(c, "strict") for c in strict] + [(c, "fallback") for c in fallback]
    measured = []
    valid_count = 0
    for offset in range(0, len(ordered), args.batch_size):
        if valid_count >= args.backgrounds:
            break
        batch = ordered[offset:offset + args.batch_size]
        data = read_for_regions(args, [c.region for c, _ in batch], source) if source else preloaded
        assert data is not None
        positions = preload_positions if preload_positions is not None else position_index(data)
        for candidate, label in batch:
            result = measure(candidate.region, data, positions, args.flank_bp, args.min_pair_haplotypes)
            measured.append(Candidate(result.region, result.left_count, result.right_count,
                                      result.median_r2, result.valid_pairs, label,
                                      matches[result.region.region_id][1],
                                      matches[result.region.region_id][2], result.max_r2))
            if result.median_r2 is not None:
                valid_count += 1
            if valid_count >= args.backgrounds:
                break
    valid = [c for c in measured if c.median_r2 is not None]
    output = args.output_prefix
    output.parent.mkdir(parents=True, exist_ok=True)
    with Path(str(output) + ".backgrounds.tsv").open("w", newline="") as handle:
        fields = ["chrom", "start", "end", "length_bp", "left_flank_bp", "right_flank_bp",
                  "full_flanks", "extended_beyond_background_bed", "left_snps", "right_snps",
                  "valid_pairs", "median_r2", "max_r2", "selection", "match_score", "match_orientation"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for c in measured:
            left_bound, right_bound = flank_bounds(c.region, args.flank_bp)
            writer.writerow(dict(chrom=c.region.display_chrom, start=c.region.start, end=c.region.end,
                                 length_bp=c.region.end - c.region.start,
                                 left_flank_bp=c.region.start - left_bound,
                                 right_flank_bp=right_bound - c.region.end,
                                 full_flanks=full_flanks(c.region, args.flank_bp),
                                 extended_beyond_background_bed=(
                                     left_bound < c.region.territory_start or
                                     right_bound > c.region.territory_end),
                                 left_snps=c.left_count,
                                 right_snps=c.right_count, valid_pairs=c.valid_pairs,
                                 median_r2="NA" if c.median_r2 is None else c.median_r2,
                                 max_r2="NA" if c.max_r2 is None else c.max_r2,
                                 selection=c.selection, match_score=c.match_score,
                                 match_orientation=c.match_orientation))
    selected_partial = sum(not full_flanks(c.region, args.flank_bp) for c in valid)
    background_complete = (target_flanks_complete and len(valid) == args.backgrounds
                           and (args.allow_background_flank_extension or selected_partial == 0))
    empirical_p = ((1 + sum(c.median_r2 >= target.median_r2 for c in valid)) /
                   (1 + len(valid))) if background_complete else math.nan
    empirical_p_max = ((1 + sum(c.max_r2 >= target.max_r2 for c in valid)) /
                       (1 + len(valid))) if background_complete else math.nan
    with Path(str(output) + ".summary.tsv").open("w", newline="") as handle:
        fields = ["target_id", "chrom", "start", "end", "flank_bp",
                  "left_flank_bp", "right_flank_bp", "min_maf",
                  "max_site_missing_rate", "left_snps", "right_snps", "valid_pairs", "median_r2", "max_r2",
                  "background_median_r2", "background_median_max_r2", "empirical_p_upper",
                  "empirical_p_upper_max_r2", "attempts", "unique_candidates",
                  "screened_nonempty", "strict_available", "selected_strict", "selected_fallback",
                  "selected_swapped", "selected_total", "selected_partial_flanks",
                  "valid_backgrounds", "allow_background_flank_extension", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        target_left_bound, target_right_bound = flank_bounds(target.region, args.flank_bp)
        writer.writerow(dict(target_id=target.region.region_id, chrom=target.region.display_chrom,
                             start=target.region.start, end=target.region.end, flank_bp=args.flank_bp,
                             left_flank_bp=target.region.start - target_left_bound,
                             right_flank_bp=target_right_bound - target.region.end,
                             min_maf=args.min_maf, max_site_missing_rate=args.max_site_missing_rate,
                             left_snps=target.left_count, right_snps=target.right_count,
                             valid_pairs=target.valid_pairs, median_r2=target.median_r2,
                             max_r2=target.max_r2,
                             background_median_r2=statistics.median(c.median_r2 for c in valid) if valid else "NA",
                             background_median_max_r2=statistics.median(c.max_r2 for c in valid) if valid else "NA",
                             empirical_p_upper=empirical_p if background_complete else "NA",
                             empirical_p_upper_max_r2=empirical_p_max if background_complete else "NA",
                             attempts=args.attempts, unique_candidates=len(regions),
                             screened_nonempty=len(screened), strict_available=len(strict),
                             selected_strict=sum(c.selection == "strict" for c in measured),
                             selected_fallback=sum(c.selection == "fallback" for c in measured),
                             selected_swapped=sum(c.match_orientation == "swapped" for c in measured),
                             selected_total=len(measured), selected_partial_flanks=selected_partial,
                             valid_backgrounds=len(valid),
                             allow_background_flank_extension=args.allow_background_flank_extension,
                             status="ok" if background_complete else
                             ("target_short_flank" if not target_flanks_complete else
                              "no_background_candidates" if not regions else
                              "insufficient_valid_backgrounds" if len(valid) < args.backgrounds else
                              "background_short_flanks")))
    print(f"target median r2={target.median_r2:.6g}; selected={len(measured)}; "
          f"valid={len(valid)}; empirical upper P="
          f"{f'{empirical_p:.6g}' if background_complete else 'NA'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
