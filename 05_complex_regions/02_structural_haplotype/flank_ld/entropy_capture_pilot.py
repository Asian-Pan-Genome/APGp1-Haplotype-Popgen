#!/usr/bin/env python3
"""Pilot fixed-window/fixed-SNP designs using random genomic windows.

The script tiles a user-supplied BED sampling universe, randomly samples a
requested number of fixed-width windows, and evaluates one or more fixed SNP
counts. For each window and SNP count k, it reports the entropy retained by a
maximum-entropy k-SNP subset relative to all eligible SNPs in that window. It
also summarizes the fraction of windows that would be removed for containing
fewer than k eligible SNPs.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
import sys
from bisect import bisect_left
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    from .fixed_window_block_ld import cap_variants_by_entropy, format_variants
    from .flanking_block_ld import (
        Region,
        Variant,
        common_haplotype_keys,
        format_float,
        haplotype_entropy,
        haplotypes_for_keys,
        read_regions,
        read_vcf,
    )
except ImportError:  # Direct execution from the repository root.
    from fixed_window_block_ld import (  # type: ignore
        cap_variants_by_entropy,
        format_variants,
    )
    from flanking_block_ld import (  # type: ignore
        Region,
        Variant,
        common_haplotype_keys,
        format_float,
        haplotype_entropy,
        haplotypes_for_keys,
        read_regions,
        read_vcf,
    )


DETAIL_FIELDS = (
    "chrom",
    "window_start",
    "window_end",
    "window_id",
    "fixed_snp_number",
    "entropy_capture",
    "entropy_loss",
    "full_entropy",
    "selected_entropy",
    "eligible_snp_number",
    "selected_snp_number",
    "window_pass_snp_number_threshold",
    "pre_missing_snp_number",
    "missed_snp_number",
    "missing_filtered_fraction",
    "window_pass_missing_qc",
    "complete_haplotype_count",
    "full_haplotype_count",
    "selected_haplotype_count",
    "selection_method",
    "selected_snps",
    "status",
)


SUMMARY_FIELDS = (
    "fixed_snp_number",
    "sampled_window_number",
    "missing_qc_fail_window_number",
    "missing_qc_fail_window_fraction",
    "sparse_window_number_all",
    "sparse_window_fraction_all",
    "missing_qc_pass_window_number",
    "sparse_window_number_among_missing_qc_pass",
    "sparse_window_fraction_among_missing_qc_pass",
    "entropy_evaluable_window_number",
    "entropy_evaluable_window_fraction",
    "entropy_capture_mean",
    "entropy_capture_median",
    "entropy_capture_q10",
    "entropy_capture_q25",
    "entropy_capture_q75",
    "entropy_capture_q90",
    "entropy_loss_mean",
    "entropy_loss_median",
)


def parse_snp_counts(value: object) -> Tuple[int, ...]:
    if isinstance(value, (tuple, list)):
        counts = [int(item) for item in value]
    else:
        tokens = str(value).replace(",", " ").split()
        try:
            counts = [int(token) for token in tokens]
        except ValueError as error:
            raise ValueError("--snp-counts must contain positive integers") from error
    if not counts or any(count < 1 for count in counts):
        raise ValueError("--snp-counts must contain positive integers")
    return tuple(sorted(set(counts)))


def tile_sampling_universe(
    regions: Sequence[Region],
    window_bp: int,
    sampling_step_bp: int,
) -> List[Region]:
    windows: List[Region] = []
    seen = set()
    for source in regions:
        if source.end - source.start < window_bp:
            continue
        window_number = 0
        for start in range(
            source.start,
            source.end - window_bp + 1,
            sampling_step_bp,
        ):
            end = start + window_bp
            key = (source.chrom, start, end)
            if key in seen:
                continue
            seen.add(key)
            window_number += 1
            windows.append(
                Region(
                    chrom=source.chrom,
                    display_chrom=source.display_chrom,
                    start=start,
                    end=end,
                    region_id=(
                        f"{source.region_id}|window_{window_number}|"
                        f"{source.display_chrom}:{start}-{end}"
                    ),
                    line_number=source.line_number,
                )
            )
    return windows


def sample_windows(
    candidates: Sequence[Region],
    requested: int,
    seed: int,
) -> List[Region]:
    if not candidates:
        raise ValueError("sampling BED contains no complete candidate windows")
    if requested < 0:
        raise ValueError("--sample-windows must be >= 0")
    if requested == 0:
        sampled = list(candidates)
    else:
        if requested > len(candidates):
            raise ValueError(
                f"requested {requested} windows but only {len(candidates)} are available"
            )
        sampled = random.Random(seed).sample(list(candidates), requested)
    return sorted(sampled, key=lambda region: (region.chrom, region.start, region.end))


def slice_variants(
    variants: Sequence[Variant],
    positions: Sequence[int],
    start: int,
    end: int,
) -> Tuple[Variant, ...]:
    left = bisect_left(positions, start)
    right = bisect_left(positions, end)
    return tuple(variants[left:right])


def count_positions(positions: Sequence[int], start: int, end: int) -> int:
    return bisect_left(positions, end) - bisect_left(positions, start)


def empirical_quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    if len(ordered) == 1:
        return ordered[0]
    index = probability * (len(ordered) - 1)
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def summarize_k(rows: Sequence[Dict[str, str]], k: int) -> Dict[str, str]:
    k_rows = [row for row in rows if int(row["fixed_snp_number"]) == k]
    total = len(k_rows)
    missing_fail = sum(row["window_pass_missing_qc"] == "false" for row in k_rows)
    missing_pass_rows = [
        row for row in k_rows if row["window_pass_missing_qc"] == "true"
    ]
    sparse_all = sum(
        int(row["eligible_snp_number"]) < k for row in k_rows
    )
    sparse_missing_pass = sum(
        int(row["eligible_snp_number"]) < k for row in missing_pass_rows
    )
    evaluable = [row for row in k_rows if row["status"] == "ok"]
    captures = [float(row["entropy_capture"]) for row in evaluable]
    losses = [float(row["entropy_loss"]) for row in evaluable]

    def fraction(numerator: int, denominator: int) -> Optional[float]:
        return numerator / denominator if denominator else None

    return {
        "fixed_snp_number": str(k),
        "sampled_window_number": str(total),
        "missing_qc_fail_window_number": str(missing_fail),
        "missing_qc_fail_window_fraction": format_float(
            fraction(missing_fail, total)
        ),
        "sparse_window_number_all": str(sparse_all),
        "sparse_window_fraction_all": format_float(fraction(sparse_all, total)),
        "missing_qc_pass_window_number": str(len(missing_pass_rows)),
        "sparse_window_number_among_missing_qc_pass": str(sparse_missing_pass),
        "sparse_window_fraction_among_missing_qc_pass": format_float(
            fraction(sparse_missing_pass, len(missing_pass_rows))
        ),
        "entropy_evaluable_window_number": str(len(evaluable)),
        "entropy_evaluable_window_fraction": format_float(
            fraction(len(evaluable), total)
        ),
        "entropy_capture_mean": format_float(
            statistics.fmean(captures) if captures else None
        ),
        "entropy_capture_median": format_float(
            statistics.median(captures) if captures else None
        ),
        "entropy_capture_q10": format_float(
            empirical_quantile(captures, 0.10) if captures else None
        ),
        "entropy_capture_q25": format_float(
            empirical_quantile(captures, 0.25) if captures else None
        ),
        "entropy_capture_q75": format_float(
            empirical_quantile(captures, 0.75) if captures else None
        ),
        "entropy_capture_q90": format_float(
            empirical_quantile(captures, 0.90) if captures else None
        ),
        "entropy_loss_mean": format_float(
            statistics.fmean(losses) if losses else None
        ),
        "entropy_loss_median": format_float(
            statistics.median(losses) if losses else None
        ),
    }


def evaluate_window(
    window: Region,
    variants: Sequence[Variant],
    pre_missing_positions: Sequence[int],
    missing_filtered_positions: Sequence[int],
    snp_counts: Sequence[int],
    max_window_missing_filtered_fraction: float,
    min_complete_haplotypes: int,
    max_exact_entropy_combinations: int,
) -> List[Dict[str, str]]:
    variant_positions = tuple(variant.bed_pos for variant in variants)
    eligible = slice_variants(
        variants,
        variant_positions,
        window.start,
        window.end,
    )
    pre_missing_count = count_positions(
        pre_missing_positions,
        window.start,
        window.end,
    )
    missed_count = count_positions(
        missing_filtered_positions,
        window.start,
        window.end,
    )
    missing_fraction = (
        missed_count / pre_missing_count if pre_missing_count else 0.0
    )
    passes_missing_qc = (
        missing_fraction <= max_window_missing_filtered_fraction
    )

    full_keys = common_haplotype_keys(eligible)
    full_haplotypes = haplotypes_for_keys(eligible, full_keys)
    full_entropy = haplotype_entropy(full_haplotypes)
    full_haplotype_count = len(set(full_haplotypes))

    rows: List[Dict[str, str]] = []
    for k in snp_counts:
        passes_snp_number = len(eligible) >= k
        selected: Tuple[Variant, ...] = ()
        selection_method = "not_selected"
        selected_entropy: Optional[float] = None
        selected_haplotype_count = 0
        capture: Optional[float] = None
        entropy_loss: Optional[float] = None

        statuses: List[str] = []
        if not passes_missing_qc:
            statuses.append("window_excess_missing_filtered_sites")
        if not passes_snp_number:
            statuses.append("insufficient_snps")
        if len(full_keys) < min_complete_haplotypes:
            statuses.append("insufficient_complete_haplotypes")
        if math.isnan(full_entropy) or full_entropy <= 0.0:
            statuses.append("zero_or_undefined_full_entropy")

        if not statuses:
            selected, selection_method = cap_variants_by_entropy(
                eligible,
                max_snps=k,
                max_exact_combinations=max_exact_entropy_combinations,
            )
            selected_haplotypes = haplotypes_for_keys(selected, full_keys)
            selected_entropy = haplotype_entropy(selected_haplotypes)
            selected_haplotype_count = len(set(selected_haplotypes))
            capture = selected_entropy / full_entropy
            capture = min(1.0, max(0.0, capture))
            entropy_loss = full_entropy - selected_entropy
            statuses.append("ok")

        rows.append(
            {
                "chrom": window.display_chrom,
                "window_start": str(window.start),
                "window_end": str(window.end),
                "window_id": window.region_id,
                "fixed_snp_number": str(k),
                "entropy_capture": format_float(capture),
                "entropy_loss": format_float(entropy_loss),
                "full_entropy": format_float(full_entropy),
                "selected_entropy": format_float(selected_entropy),
                "eligible_snp_number": str(len(eligible)),
                "selected_snp_number": str(len(selected)),
                "window_pass_snp_number_threshold": str(
                    passes_snp_number
                ).lower(),
                "pre_missing_snp_number": str(pre_missing_count),
                "missed_snp_number": str(missed_count),
                "missing_filtered_fraction": format_float(missing_fraction),
                "window_pass_missing_qc": str(passes_missing_qc).lower(),
                "complete_haplotype_count": str(len(full_keys)),
                "full_haplotype_count": str(full_haplotype_count),
                "selected_haplotype_count": str(selected_haplotype_count),
                "selection_method": selection_method,
                "selected_snps": format_variants(selected),
                "status": ";".join(statuses),
            }
        )
    return rows


def run(args: argparse.Namespace) -> Dict[str, object]:
    snp_counts = parse_snp_counts(args.snp_counts)
    sampling_step_bp = (
        args.window_bp
        if args.sampling_step_bp is None
        else args.sampling_step_bp
    )
    universe = read_regions(args.sampling_bed)
    candidate_windows = tile_sampling_universe(
        universe,
        window_bp=args.window_bp,
        sampling_step_bp=sampling_step_bp,
    )
    windows = sample_windows(candidate_windows, args.sample_windows, args.seed)

    vcf_data = read_vcf(
        args.vcf,
        windows,
        flank_bp=0,
        min_maf=args.min_maf,
        min_call_rate=1.0 - args.max_site_missing_rate,
        include_filtered=args.include_filtered,
        allow_unphased_missing=args.allow_unphased_missing,
        include_region_interiors=True,
    )

    detail_rows: List[Dict[str, str]] = []
    for window in windows:
        detail_rows.extend(
            evaluate_window(
                window,
                vcf_data.variants_by_chrom.get(window.chrom, ()),
                vcf_data.pre_missing_snp_positions_by_chrom.get(window.chrom, ()),
                vcf_data.missing_filtered_snp_positions_by_chrom.get(
                    window.chrom, ()
                ),
                snp_counts=snp_counts,
                max_window_missing_filtered_fraction=(
                    args.max_window_missing_filtered_fraction
                ),
                min_complete_haplotypes=args.min_complete_haplotypes,
                max_exact_entropy_combinations=(
                    args.max_exact_entropy_combinations
                ),
            )
        )

    summary_rows = [summarize_k(detail_rows, k) for k in snp_counts]
    detail_path = Path(f"{args.output_prefix}.windows.tsv")
    summary_path = Path(f"{args.output_prefix}.summary.tsv")
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    with detail_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DETAIL_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(detail_rows)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    return {
        "candidate_window_count": len(candidate_windows),
        "sampled_window_count": len(windows),
        "snp_counts": snp_counts,
        "detail_path": detail_path,
        "summary_path": summary_path,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Randomly sample fixed genomic windows and estimate entropy capture "
            "and SNP-sparsity filtering rates for fixed SNP counts."
        )
    )
    parser.add_argument(
        "sampling_bed",
        type=Path,
        help="BED intervals defining the genomic sampling universe",
    )
    parser.add_argument("vcf", type=Path, help="phased VCF, VCF.gz, or VCF.bgz")
    parser.add_argument(
        "--output-prefix",
        type=Path,
        required=True,
        help="prefix for .windows.tsv and .summary.tsv outputs",
    )
    parser.add_argument(
        "--sample-windows",
        type=int,
        default=1000,
        help="number of candidate windows to sample; 0 uses all (default: 1000)",
    )
    parser.add_argument(
        "--window-bp",
        type=int,
        default=10_000,
        help="physical width of each sampled window (default: 10000)",
    )
    parser.add_argument(
        "--sampling-step-bp",
        type=int,
        default=None,
        help="step used to tile the sampling BED (default: window-bp)",
    )
    parser.add_argument(
        "--snp-counts",
        default="3,5,8,10",
        help="comma- or space-separated fixed SNP counts (default: 3,5,8,10)",
    )
    parser.add_argument(
        "--max-exact-entropy-combinations",
        type=int,
        default=100_000,
        help=(
            "maximum subset count for exact entropy maximization; larger "
            "searches use the greedy approximation (default: 100000)"
        ),
    )
    parser.add_argument(
        "--min-complete-haplotypes",
        type=int,
        default=2,
        help="minimum complete haplotypes needed for entropy capture (default: 2)",
    )
    parser.add_argument("--min-maf", type=float, default=0.1)
    parser.add_argument(
        "--max-site-missing-rate",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--max-window-missing-filtered-fraction",
        type=float,
        default=0.10,
    )
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--include-filtered", action="store_true")
    parser.add_argument("--allow-unphased-missing", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.sample_windows < 0:
        parser.error("--sample-windows must be >= 0")
    if args.window_bp < 1:
        parser.error("--window-bp must be >= 1")
    if args.sampling_step_bp is not None and args.sampling_step_bp < 1:
        parser.error("--sampling-step-bp must be >= 1")
    if args.max_exact_entropy_combinations < 1:
        parser.error("--max-exact-entropy-combinations must be >= 1")
    if args.min_complete_haplotypes < 2:
        parser.error("--min-complete-haplotypes must be >= 2")
    if not 0.0 <= args.min_maf <= 0.5:
        parser.error("--min-maf must be in [0, 0.5]")
    if not 0.0 <= args.max_site_missing_rate <= 1.0:
        parser.error("--max-site-missing-rate must be in [0, 1]")
    if not 0.0 <= args.max_window_missing_filtered_fraction <= 1.0:
        parser.error("--max-window-missing-filtered-fraction must be in [0, 1]")
    try:
        parse_snp_counts(args.snp_counts)
        result = run(args)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(
        "candidate_windows={candidate_window_count} "
        "sampled_windows={sampled_window_count} snp_counts={snp_counts} "
        "details={detail_path} summary={summary_path}".format(**result),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
