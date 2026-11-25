#!/usr/bin/env python3
"""Calculate multiallelic D' for ranked fixed windows flanking BED regions.

For rank r (starting at one), the upstream block is the r-th fixed-width
window moving left from the BED start, and the downstream block is the r-th
window moving right from the BED end. In the default ``matched`` mode, rows
are emitted for L1-R1, L2-R2, ..., in near-to-far order. ``all`` mode emits
every upstream/downstream window combination. ``all-blocks`` orders all left
and right blocks by genomic position and emits every unique pair, including
left-left, left-right, and right-right comparisons.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

try:
    from .flanking_block_ld import (
        Haplotype,
        HaplotypeKey,
        Region,
        Variant,
        common_haplotype_keys,
        format_float,
        haplotype_entropy,
        haplotypes_for_keys,
        multiallelic_dprime,
        permutation_statistics,
        read_regions,
        read_vcf,
        select_greedy_entropy,
    )
except ImportError:  # Direct execution from the repository root.
    from flanking_block_ld import (  # type: ignore
        Haplotype,
        HaplotypeKey,
        Region,
        Variant,
        common_haplotype_keys,
        format_float,
        haplotype_entropy,
        haplotypes_for_keys,
        multiallelic_dprime,
        permutation_statistics,
        read_regions,
        read_vcf,
        select_greedy_entropy,
    )


@dataclass(frozen=True)
class FixedBlock:
    side: str
    rank: int
    start: int
    end: int
    pre_missing_snp_count: int
    missing_filtered_snp_count: int
    missing_filtered_fraction: float
    passes_missing_qc: bool
    snp_count_before_cap: int
    selection_method: str
    variants: Tuple[Variant, ...]


def entropy_score(
    variants: Sequence[Variant],
    haplotype_keys: Optional[Sequence[HaplotypeKey]] = None,
) -> Tuple[float, int]:
    """Score a SNP subset, optionally on one pre-fixed haplotype cohort."""
    keys = (tuple(haplotype_keys) if haplotype_keys is not None
            else common_haplotype_keys(variants))
    entropy = haplotype_entropy(haplotypes_for_keys(variants, keys))
    if math.isnan(entropy):
        entropy = -math.inf
    return entropy, len(keys)


def select_greedy_entropy_for_keys(
    candidates: Sequence[Variant],
    block_snps: int,
    haplotype_keys: Sequence[HaplotypeKey],
) -> Tuple[Variant, ...]:
    """Greedy entropy selection on a fixed pair-complete haplotype set."""
    remaining = list(candidates)
    chosen: List[Variant] = []
    while remaining and len(chosen) < block_snps:
        best_index = 0
        best_score = (-math.inf, -1, 0)
        for index, candidate in enumerate(remaining):
            entropy, complete = entropy_score(chosen + [candidate], haplotype_keys)
            score = (entropy, complete, -index)
            if score > best_score:
                best_score = score
                best_index = index
        chosen.append(remaining.pop(best_index))
    return tuple(chosen)


def cap_variants_by_entropy(
    variants_nearest_first: Sequence[Variant],
    max_snps: int,
    max_exact_combinations: int,
    haplotype_keys: Optional[Sequence[HaplotypeKey]] = None,
) -> Tuple[Tuple[Variant, ...], str]:
    if len(variants_nearest_first) <= max_snps:
        return tuple(variants_nearest_first), "all_within_upper_limit"

    combination_count = math.comb(len(variants_nearest_first), max_snps)
    if combination_count <= max_exact_combinations:
        best: Optional[Tuple[Variant, ...]] = None
        best_score = (-math.inf, -1)
        for combination in itertools.combinations(variants_nearest_first, max_snps):
            score = entropy_score(combination, haplotype_keys)
            if score > best_score:
                best = combination
                best_score = score
        assert best is not None
        return best, "exact_max_entropy"

    selected = (
        select_greedy_entropy(variants_nearest_first, max_snps)
        if haplotype_keys is None else
        select_greedy_entropy_for_keys(
            variants_nearest_first, max_snps, haplotype_keys,
        )
    )
    return selected, "greedy_max_entropy"


def make_fixed_blocks(
    region: Region,
    variants: Sequence[Variant],
    pre_missing_positions: Sequence[int],
    missing_filtered_positions: Sequence[int],
    window_bp: int,
    window_step_bp: int,
    block_count: int,
    max_snps_per_block: int,
    max_window_missing_filtered_fraction: float,
    max_exact_entropy_combinations: int,
) -> Tuple[Tuple[FixedBlock, ...], Tuple[FixedBlock, ...]]:
    left: List[FixedBlock] = []
    right: List[FixedBlock] = []
    for rank in range(1, block_count + 1):
        offset = (rank - 1) * window_step_bp
        left_end = region.start - offset
        left_start = max(0, left_end - window_bp)
        right_start = region.end + offset
        right_end = right_start + window_bp
        left_candidates = tuple(
            variant
            for variant in variants
            if left_start <= variant.bed_pos < left_end
        )
        right_candidates = tuple(
            variant
            for variant in variants
            if right_start <= variant.bed_pos < right_end
        )

        left_pre_missing = sum(
            left_start <= position < left_end for position in pre_missing_positions
        )
        right_pre_missing = sum(
            right_start <= position < right_end for position in pre_missing_positions
        )
        left_missing_filtered = sum(
            left_start <= position < left_end for position in missing_filtered_positions
        )
        right_missing_filtered = sum(
            right_start <= position < right_end for position in missing_filtered_positions
        )
        left_missing_fraction = (
            left_missing_filtered / left_pre_missing if left_pre_missing else 0.0
        )
        right_missing_fraction = (
            right_missing_filtered / right_pre_missing if right_pre_missing else 0.0
        )
        left_passes_missing = (
            left_missing_fraction <= max_window_missing_filtered_fraction
        )
        right_passes_missing = (
            right_missing_fraction <= max_window_missing_filtered_fraction
        )

        if left_passes_missing:
            left_variants, left_selection = cap_variants_by_entropy(
                sorted(
                    left_candidates,
                    key=lambda variant: (-variant.pos, variant.variant_id),
                ),
                max_snps=max_snps_per_block,
                max_exact_combinations=max_exact_entropy_combinations,
            )
        else:
            left_variants = ()
            left_selection = "not_selected_failed_missing_qc"
        if right_passes_missing:
            right_variants, right_selection = cap_variants_by_entropy(
                sorted(
                    right_candidates,
                    key=lambda variant: (variant.pos, variant.variant_id),
                ),
                max_snps=max_snps_per_block,
                max_exact_combinations=max_exact_entropy_combinations,
            )
        else:
            right_variants = ()
            right_selection = "not_selected_failed_missing_qc"

        left.append(
            FixedBlock(
                "left",
                rank,
                left_start,
                left_end,
                left_pre_missing,
                left_missing_filtered,
                left_missing_fraction,
                left_passes_missing,
                len(left_candidates),
                left_selection,
                left_variants,
            )
        )
        right.append(
            FixedBlock(
                "right",
                rank,
                right_start,
                right_end,
                right_pre_missing,
                right_missing_filtered,
                right_missing_fraction,
                right_passes_missing,
                len(right_candidates),
                right_selection,
                right_variants,
            )
        )
    return tuple(left), tuple(right)


def pair_blocks(
    left: Sequence[FixedBlock],
    right: Sequence[FixedBlock],
    pairing: str,
) -> List[Tuple[FixedBlock, FixedBlock]]:
    if pairing == "matched":
        return list(zip(left, right))
    if pairing == "all":
        pairs = [(left_block, right_block) for left_block in left for right_block in right]
        return sorted(
            pairs,
            key=lambda pair: (
                max(pair[0].rank, pair[1].rank),
                pair[0].rank + pair[1].rank,
                pair[0].rank,
                pair[1].rank,
            ),
        )
    if pairing == "all-blocks":
        genomic_order = tuple(reversed(left)) + tuple(right)
        return list(itertools.combinations(genomic_order, 2))
    raise ValueError(f"unsupported pairing mode: {pairing}")


def genomic_block_index(block: FixedBlock, block_count: int) -> int:
    """Return a one-based index in genomic order: Lmax...L1,R1...Rmax."""
    if block.side == "left":
        return block_count - block.rank + 1
    return block_count + block.rank


def block_label(block: FixedBlock) -> str:
    return f"{'L' if block.side == 'left' else 'R'}{block.rank}"


def format_variants(variants: Sequence[Variant]) -> str:
    return ",".join(f"{variant.variant_id}@{variant.pos}" for variant in variants) or "."


def summarize_block_pair(
    region: Region,
    sample_count: int,
    left: FixedBlock,
    right: FixedBlock,
    min_snps_per_block: int,
    permutations: int,
    rng: random.Random,
    block_count: Optional[int] = None,
    float_formatter: Callable[[Optional[float]], str] = format_float,
    pair_specific_entropy_selection: bool = False,
    max_snps_per_block: Optional[int] = None,
    max_exact_entropy_combinations: int = 100_000,
    min_complete_haplotypes: int = 2,
    min_complete_haplotype_fraction: float = 0.0,
) -> Dict[str, str]:
    left_snp_number_passes = left.snp_count_before_cap >= min_snps_per_block
    right_snp_number_passes = right.snp_count_before_cap >= min_snps_per_block
    enough_left = left.passes_missing_qc and len(left.variants) >= min_snps_per_block
    enough_right = right.passes_missing_qc and len(right.variants) >= min_snps_per_block
    left_variants = left.variants
    right_variants = right.variants
    left_selection_method = left.selection_method
    right_selection_method = right.selection_method
    keys: Tuple[HaplotypeKey, ...] = ()
    if enough_left and enough_right:
        # In pair-specific mode this is the fixed cohort for screening,
        # entropy selection, and Dprime. It is intentionally based on every
        # QC-retained candidate SNP, including SNPs later excluded by the cap.
        keys = common_haplotype_keys(left.variants + right.variants)
    total_haplotypes = 2 * sample_count
    complete_fraction = len(keys) / total_haplotypes if total_haplotypes else 0.0
    fraction_minimum = math.ceil(min_complete_haplotype_fraction * total_haplotypes)
    required_complete = max(2, min_complete_haplotypes, fraction_minimum)
    raw_left_haplotypes = haplotypes_for_keys(left.variants, keys)
    raw_right_haplotypes = haplotypes_for_keys(right.variants, keys)
    pair_screen_passes = (
        len(keys) >= required_complete and
        len(set(raw_left_haplotypes)) >= 2 and
        len(set(raw_right_haplotypes)) >= 2
    )
    if pair_specific_entropy_selection:
        if max_snps_per_block is None:
            raise ValueError("pair-specific entropy selection requires max_snps_per_block")
        if pair_screen_passes:
            left_variants, left_selection_method = cap_variants_by_entropy(
                left.variants, max_snps_per_block, max_exact_entropy_combinations,
                haplotype_keys=keys,
            )
            right_variants, right_selection_method = cap_variants_by_entropy(
                right.variants, max_snps_per_block, max_exact_entropy_combinations,
                haplotype_keys=keys,
            )
        else:
            left_selection_method = "not_selected_pair_qc_failed"
            right_selection_method = "not_selected_pair_qc_failed"
    left_haplotypes: List[Haplotype] = haplotypes_for_keys(left_variants, keys)
    right_haplotypes: List[Haplotype] = haplotypes_for_keys(right_variants, keys)

    dprime: Optional[float] = None
    left_entropy: Optional[float] = None
    right_entropy: Optional[float] = None
    permutation_mean: Optional[float] = None
    dprime_excess: Optional[float] = None
    empirical_p: Optional[float] = None
    if len(keys) >= required_complete and (
            not pair_specific_entropy_selection or pair_screen_passes):
        dprime = multiallelic_dprime(left_haplotypes, right_haplotypes)
        left_entropy = haplotype_entropy(left_haplotypes)
        right_entropy = haplotype_entropy(right_haplotypes)
        permutation_mean, dprime_excess, empirical_p = permutation_statistics(
            left_haplotypes,
            right_haplotypes,
            dprime,
            permutations,
            rng,
        )

    statuses: List[str] = []
    if not left.passes_missing_qc:
        statuses.append("left_window_excess_missing_filtered_sites")
    elif not enough_left:
        statuses.append("insufficient_left_snps")
    if not right.passes_missing_qc:
        statuses.append("right_window_excess_missing_filtered_sites")
    elif not enough_right:
        statuses.append("insufficient_right_snps")
    if enough_left and enough_right:
        if len(keys) < 2:
            statuses.append("insufficient_complete_haplotypes")
        elif len(keys) < min_complete_haplotypes:
            statuses.append("below_min_complete_haplotypes")
        elif len(keys) < fraction_minimum:
            statuses.append("below_min_complete_haplotype_fraction")
    if (enough_left and enough_right and len(keys) >= required_complete and
            len(set(left_haplotypes)) < 2):
        statuses.append("single_left_haplotype")
    if (enough_left and enough_right and len(keys) >= required_complete and
            len(set(right_haplotypes)) < 2):
        statuses.append("single_right_haplotype")
    if not statuses:
        statuses.append("ok")

    pair_rank = (
        str(left.rank)
        if left.side != right.side and left.rank == right.rank
        else "."
    )
    if block_count is None:
        block_count = max(left.rank, right.rank)
    pair_type = f"{left.side}_{right.side}"
    center_distance = (
        (right.start + right.end) / 2.0 - (left.start + left.end) / 2.0
    )
    return {
        "chrom": region.display_chrom,
        "region_start": str(region.start),
        "region_end": str(region.end),
        "region_id": region.region_id,
        "pair_rank_near_to_far": pair_rank,
        "left_block_rank": str(left.rank),
        "right_block_rank": str(right.rank),
        "left_window_start": str(left.start),
        "left_window_end": str(left.end),
        "right_window_start": str(right.start),
        "right_window_end": str(right.end),
        "block_center_distance_bp": float_formatter(center_distance),
        "sample_count": str(sample_count),
        "left_pre_missing_snp_count": str(left.pre_missing_snp_count),
        "right_pre_missing_snp_count": str(right.pre_missing_snp_count),
        "left_missing_filtered_snp_count": str(left.missing_filtered_snp_count),
        "right_missing_filtered_snp_count": str(right.missing_filtered_snp_count),
        "left_missing_filtered_fraction": float_formatter(left.missing_filtered_fraction),
        "right_missing_filtered_fraction": float_formatter(right.missing_filtered_fraction),
        "left_snp_count_before_cap": str(left.snp_count_before_cap),
        "right_snp_count_before_cap": str(right.snp_count_before_cap),
        "left_snp_selection_method": left_selection_method,
        "right_snp_selection_method": right_selection_method,
        "left_snp_count": str(len(left_variants)),
        "right_snp_count": str(len(right_variants)),
        "left_snps": format_variants(left_variants),
        "right_snps": format_variants(right_variants),
        "complete_haplotype_count": str(len(keys)),
        "complete_haplotype_fraction": float_formatter(complete_fraction),
        "required_complete_haplotype_count": str(required_complete),
        "min_complete_haplotype_fraction": float_formatter(min_complete_haplotype_fraction),
        "left_haplotype_count": str(len(set(left_haplotypes))),
        "right_haplotype_count": str(len(set(right_haplotypes))),
        "joint_haplotype_count": str(len(set(zip(left_haplotypes, right_haplotypes)))),
        "left_haplotype_entropy": float_formatter(left_entropy),
        "right_haplotype_entropy": float_formatter(right_entropy),
        "multiallelic_dprime": float_formatter(dprime),
        "block1_snp_number_pass_threshold": str(len(left_variants)),
        "block2_snp_number_pass_threshold": str(len(right_variants)),
        "both_blocks_pass_snp_number_threshold": str(
            left_snp_number_passes and right_snp_number_passes
        ).lower(),
        "block1_missed_snp_number": str(left.missing_filtered_snp_count),
        "block2_missed_snp_number": str(right.missing_filtered_snp_count),
        "total_missed_snp_number": str(
            left.missing_filtered_snp_count + right.missing_filtered_snp_count
        ),
        "dprime_permutation_mean": float_formatter(permutation_mean),
        "dprime_excess": float_formatter(dprime_excess),
        "dprime_empirical_p": float_formatter(empirical_p),
        "pair_type": pair_type,
        "block1_side": left.side,
        "block1_rank": str(left.rank),
        "block1_label": block_label(left),
        "block1_genomic_index": str(genomic_block_index(left, block_count)),
        "block2_side": right.side,
        "block2_rank": str(right.rank),
        "block2_label": block_label(right),
        "block2_genomic_index": str(genomic_block_index(right, block_count)),
        "status": ";".join(statuses),
    }


OUTPUT_FIELDS = (
    "chrom",
    "region_start",
    "region_end",
    "region_id",
    "pair_rank_near_to_far",
    "left_block_rank",
    "right_block_rank",
    "multiallelic_dprime",
    "block1_snp_number_pass_threshold",
    "block2_snp_number_pass_threshold",
    "both_blocks_pass_snp_number_threshold",
    "block1_missed_snp_number",
    "block2_missed_snp_number",
    "total_missed_snp_number",
    "dprime_excess",
    "dprime_permutation_mean",
    "dprime_empirical_p",
    "pair_type",
    "block1_side",
    "block1_rank",
    "block1_label",
    "block1_genomic_index",
    "block2_side",
    "block2_rank",
    "block2_label",
    "block2_genomic_index",
    "left_window_start",
    "left_window_end",
    "right_window_start",
    "right_window_end",
    "block_center_distance_bp",
    "sample_count",
    "left_pre_missing_snp_count",
    "right_pre_missing_snp_count",
    "left_missing_filtered_snp_count",
    "right_missing_filtered_snp_count",
    "left_missing_filtered_fraction",
    "right_missing_filtered_fraction",
    "left_snp_count_before_cap",
    "right_snp_count_before_cap",
    "left_snp_selection_method",
    "right_snp_selection_method",
    "left_snp_count",
    "right_snp_count",
    "left_snps",
    "right_snps",
    "complete_haplotype_count",
    "complete_haplotype_fraction",
    "required_complete_haplotype_count",
    "min_complete_haplotype_fraction",
    "left_haplotype_count",
    "right_haplotype_count",
    "joint_haplotype_count",
    "left_haplotype_entropy",
    "right_haplotype_entropy",
    "status",
)


def run(args: argparse.Namespace) -> Dict[str, int]:
    regions = read_regions(args.regions_bed)
    window_step_bp = getattr(args, "window_step_bp", None)
    if window_step_bp is None:
        window_step_bp = args.window_bp
    flank_span = args.window_bp + (args.block_count - 1) * window_step_bp
    vcf_data = read_vcf(
        args.vcf,
        regions,
        flank_bp=flank_span,
        min_maf=args.min_maf,
        min_call_rate=1.0 - args.max_site_missing_rate,
        include_filtered=args.include_filtered,
        allow_unphased_missing=args.allow_unphased_missing,
    )

    rows: List[Dict[str, str]] = []
    pair_index = 0
    for region_index, region in enumerate(regions):
        chromosome_variants = vcf_data.variants_by_chrom.get(region.chrom, ())
        left, right = make_fixed_blocks(
            region,
            chromosome_variants,
            vcf_data.pre_missing_snp_positions_by_chrom.get(region.chrom, ()),
            vcf_data.missing_filtered_snp_positions_by_chrom.get(region.chrom, ()),
            window_bp=args.window_bp,
            window_step_bp=window_step_bp,
            block_count=args.block_count,
            max_snps_per_block=args.max_snps_per_block,
            max_window_missing_filtered_fraction=(
                args.max_window_missing_filtered_fraction
            ),
            max_exact_entropy_combinations=args.max_exact_entropy_combinations,
        )
        for left_block, right_block in pair_blocks(left, right, args.pairing):
            rows.append(
                summarize_block_pair(
                    region,
                    sample_count=vcf_data.sample_count,
                    left=left_block,
                    right=right_block,
                    min_snps_per_block=args.min_snps_per_block,
                    permutations=args.permutations,
                    rng=random.Random(args.seed + region_index * 1_000_003 + pair_index),
                    block_count=args.block_count,
                )
            )
            pair_index += 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    return {
        "region_count": len(regions),
        "block_pair_count": len(rows),
        "ok_block_pair_count": sum(row["status"] == "ok" for row in rows),
        "sample_count": vcf_data.sample_count,
        "candidate_variant_count": vcf_data.candidate_variant_count,
        "retained_variant_count": sum(
            len(variants) for variants in vcf_data.variants_by_chrom.values()
        ),
        "pre_missing_snp_count": sum(
            len(positions)
            for positions in vcf_data.pre_missing_snp_positions_by_chrom.values()
        ),
        "missing_filtered_snp_count": sum(
            len(positions)
            for positions in vcf_data.missing_filtered_snp_positions_by_chrom.values()
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate multiallelic D' for fixed-width SNP blocks moving "
            "outward from BED regions."
        )
    )
    parser.add_argument("regions_bed", type=Path, help="BED3/BED4 regions")
    parser.add_argument("vcf", type=Path, help="phased VCF, VCF.gz, or VCF.bgz")
    parser.add_argument("--output", "-o", type=Path, required=True, help="output TSV")
    parser.add_argument(
        "--window-bp",
        type=int,
        default=10_000,
        help="width of each fixed block (default: 10000)",
    )
    parser.add_argument(
        "--window-step-bp",
        type=int,
        default=None,
        help=(
            "distance between consecutive window starts; values smaller than "
            "--window-bp create overlapping sliding windows (default: window-bp)"
        ),
    )
    parser.add_argument(
        "--block-count",
        type=int,
        default=10,
        help="number of ranked blocks on each side (default: 10)",
    )
    parser.add_argument(
        "--max-flank-bp",
        type=int,
        default=100_000,
        help="maximum allowed span on each side (default: 100000)",
    )
    parser.add_argument(
        "--pairing",
        choices=("matched", "all", "all-blocks"),
        default="matched",
        help=(
            "matched emits L1-R1...L10-R10; all emits every L-R pair; "
            "all-blocks emits all unique pairs among both flanks"
        ),
    )
    parser.add_argument(
        "--min-snps-per-block",
        type=int,
        default=1,
        help="minimum retained SNP count required in each window (default: 1)",
    )
    parser.add_argument(
        "--max-snps-per-block",
        type=int,
        default=10,
        help=(
            "maximum SNP count per window; excess SNPs are selected by maximum "
            "haplotype entropy (default: 10)"
        ),
    )
    parser.add_argument(
        "--max-exact-entropy-combinations",
        type=int,
        default=100_000,
        help=(
            "maximum subset count for exact entropy maximization; larger searches "
            "use greedy selection (default: 100000)"
        ),
    )
    parser.add_argument("--min-maf", type=float, default=0.1)
    parser.add_argument(
        "--max-site-missing-rate",
        type=float,
        default=0.05,
        help="discard a SNP when its missing-sample fraction exceeds this value (default: 0.05)",
    )
    parser.add_argument(
        "--max-window-missing-filtered-fraction",
        type=float,
        default=0.10,
        help=(
            "discard a window when the fraction of pre-missing-QC SNPs removed "
            "for site missingness exceeds this value (default: 0.10)"
        ),
    )
    parser.add_argument("--permutations", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--include-filtered", action="store_true")
    parser.add_argument("--allow-unphased-missing", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.window_bp < 1:
        parser.error("--window-bp must be >= 1")
    if args.window_step_bp is not None and args.window_step_bp < 1:
        parser.error("--window-step-bp must be >= 1")
    if args.block_count < 1:
        parser.error("--block-count must be >= 1")
    if args.max_flank_bp < 1:
        parser.error("--max-flank-bp must be >= 1")
    window_step_bp = (
        args.window_bp if args.window_step_bp is None else args.window_step_bp
    )
    flank_span = args.window_bp + (args.block_count - 1) * window_step_bp
    if flank_span > args.max_flank_bp:
        parser.error(
            "the requested window span exceeds --max-flank-bp: "
            "--window-bp + (--block-count - 1) * --window-step-bp "
            f"({flank_span} > {args.max_flank_bp})"
        )
    if args.min_snps_per_block < 1:
        parser.error("--min-snps-per-block must be >= 1")
    if args.max_snps_per_block < args.min_snps_per_block:
        parser.error("--max-snps-per-block must be >= --min-snps-per-block")
    if args.max_exact_entropy_combinations < 1:
        parser.error("--max-exact-entropy-combinations must be >= 1")
    if not 0.0 <= args.min_maf <= 0.5:
        parser.error("--min-maf must be in [0, 0.5]")
    if not 0.0 <= args.max_site_missing_rate <= 1.0:
        parser.error("--max-site-missing-rate must be in [0, 1]")
    if not 0.0 <= args.max_window_missing_filtered_fraction <= 1.0:
        parser.error("--max-window-missing-filtered-fraction must be in [0, 1]")
    if args.permutations < 0:
        parser.error("--permutations must be >= 0")
    try:
        summary = run(args)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(
        "regions={region_count} block_pairs={block_pair_count} "
        "ok={ok_block_pair_count} samples={sample_count} "
        "candidate_variants={candidate_variant_count} "
        "pre_missing_snps={pre_missing_snp_count} "
        "missing_filtered_snps={missing_filtered_snp_count} "
        "retained_variants={retained_variant_count}".format(**summary),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
