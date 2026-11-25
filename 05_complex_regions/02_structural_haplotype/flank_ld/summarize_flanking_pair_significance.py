#!/usr/bin/env python3
"""Summarize significant-pair fractions and boundary-anchored prefix maxima.

Input: background_block_ld.py's complete *.target_pairs.tsv, including QC-failed
pairs. No VCF access or new LD/background calculations are performed. The maximum
fraction is descriptive: it is not a scan-adjusted empirical significance test.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple


REQUIRED_FIELDS = (
    "region_id", "block1_side", "block2_side", "block1_rank", "block2_rank",
    "status", "multiallelic_dprime",
)
METADATA_FIELDS = ("region_id", "chrom", "region_start", "region_end")
SUMMARY_FIELDS = (
    "overall_significant_pair_fraction", "max_prefix_significant_pair_fraction",
) + METADATA_FIELDS + (
    "n_eligible_pairs", "n_significant_pairs",
    "best_left_max_rank", "best_right_max_rank",
    "best_n_eligible_pairs", "best_n_significant_pairs", "n_maximizing_prefixes",
    "left_max_rank", "right_max_rank", "n_input_cross_pairs", "n_ineligible_pairs",
    "n_unreported_pairs", "n_ignored_same_side_pairs", "n_prefixes", "n_testable_prefixes",
    "expansion_mode", "p_column", "p_threshold", "min_eligible_pairs", "status",
)
EXPANSION_FIELDS = ("significant_pair_fraction",) + METADATA_FIELDS + (
    "left_max_rank", "right_max_rank", "n_eligible_pairs", "n_significant_pairs",
    "meets_min_eligible_pairs", "is_maximum", "is_selected_maximum",
    "expansion_mode", "p_column", "p_threshold", "min_eligible_pairs",
)
MISSING = {"", "NA", "N/A", ".", "NAN", "NONE", "NULL"}


@dataclass
class RegionPairs:
    metadata: Dict[str, str]
    # Key = physical near-to-far ranks, value = (eligible, significant).
    pairs: Dict[Tuple[int, int], Tuple[bool, bool]] = field(default_factory=dict)
    ignored_same_side: int = 0


def finite_number(value: str, label: str) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"{label}: expected a finite number, got {value!r}") from error
    if not math.isfinite(number):
        raise ValueError(f"{label}: expected a finite number, got {value!r}")
    return number


def positive_rank(value: str, label: str) -> int:
    try:
        rank = int(value)
    except ValueError as error:
        raise ValueError(f"{label}: rank must be a positive integer, got {value!r}") from error
    if rank < 1:
        raise ValueError(f"{label}: rank must be >= 1")
    return rank


def read_pairs(path: Path, p_column: str, threshold: float) -> Dict[str, RegionPairs]:
    regions: Dict[str, RegionPairs] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = set(REQUIRED_FIELDS + (p_column,)) - set(fields)
        if missing:
            raise ValueError(
                f"{path}: missing columns {', '.join(sorted(missing))}. "
                "Use background_block_ld.py's .target_pairs.tsv with pair tests enabled; "
                "the region summary and raw background-pair table are not suitable inputs."
            )
        if len(fields) != len(set(fields)):
            raise ValueError(f"{path}: duplicate header columns")
        for line, raw in enumerate(reader, 2):
            if None in raw or any(value is None for value in raw.values()):
                raise ValueError(f"{path}:{line}: row width does not match the header")
            row = {key: value.strip() for key, value in raw.items()}
            region_id = row["region_id"]
            label = f"{path}:{line} region={region_id}"
            if region_id.upper() in MISSING:
                raise ValueError(f"{label}: missing region_id")
            metadata = {name: row.get(name, "NA") for name in METADATA_FIELDS}
            region = regions.setdefault(region_id, RegionPairs(metadata))
            if region.metadata != metadata:
                raise ValueError(f"{label}: inconsistent coordinates for the same region_id")
            side1, side2 = row["block1_side"], row["block2_side"]
            if side1 not in ("left", "right") or side2 not in ("left", "right"):
                raise ValueError(f"{label}: block sides must be 'left' or 'right'")
            if side1 == side2:
                region.ignored_same_side += 1
                continue
            rank1 = positive_rank(row["block1_rank"], label)
            rank2 = positive_rank(row["block2_rank"], label)
            ranks = (rank1, rank2) if side1 == "left" else (rank2, rank1)
            if ranks in region.pairs:
                raise ValueError(f"{label}: duplicate cross-flank pair L{ranks[0]}-R{ranks[1]}")
            if not row["status"]:
                raise ValueError(f"{label}: missing pair QC status")
            eligible = row["status"] == "ok"
            significant = False
            if eligible:
                dprime = finite_number(row["multiallelic_dprime"], f"{label} Dprime")
                if not -1e-12 <= dprime <= 1 + 1e-12:
                    raise ValueError(f"{label}: Dprime must be in [0, 1]")
                if row[p_column].upper() in MISSING:
                    raise ValueError(
                        f"{label}: eligible pair L{ranks[0]}-R{ranks[1]} has no {p_column}. "
                        "Run background_block_ld.py without --skip-pair-test. "
                        "Missing P values cannot be counted as nonsignificant or silently dropped."
                    )
                p = finite_number(row[p_column], f"{label} {p_column}")
                if not 0 <= p <= 1:
                    raise ValueError(f"{label}: {p_column} must be in [0, 1]")
                significant = p < threshold
            region.pairs[ranks] = (eligible, significant)
    if not regions:
        raise ValueError(f"{path}: no data rows")
    return regions


def summarize_region(
    region: RegionPairs, expansion_mode: str = "symmetric", min_eligible_pairs: int = 1,
    threshold: float = 0.01, p_column: str = "pair_nominal_p",
    allow_incomplete_grid: bool = False,
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    """Count all valid pairs inside prefixes anchored at rank 1 on BOTH sides."""
    if expansion_mode not in ("symmetric", "independent"):
        raise ValueError("expansion_mode must be symmetric or independent")
    if min_eligible_pairs < 1:
        raise ValueError("min_eligible_pairs must be >= 1")
    if not region.pairs:
        raise ValueError(f"region={region.metadata['region_id']}: no cross-flank pairs")
    max_left = max(left for left, _ in region.pairs)
    max_right = max(right for _, right in region.pairs)
    unreported = max_left * max_right - len(region.pairs)
    if unreported and not allow_incomplete_grid:
        raise ValueError(
            f"region={region.metadata['region_id']}: incomplete cross-flank grid "
            f"({len(region.pairs)}/{max_left * max_right} rows). Supply the complete "
            ".target_pairs.tsv, including QC-failed rows; matched-rank pairs alone are insufficient. "
            "Use --allow-incomplete-grid only if omitted pairs are known to be ineligible."
        )

    # Two-dimensional cumulative sums; ineligible/unreported pairs add zero.
    eligible = [[0] * (max_right + 1) for _ in range(max_left + 1)]
    significant = [[0] * (max_right + 1) for _ in range(max_left + 1)]
    for left in range(1, max_left + 1):
        for right in range(1, max_right + 1):
            is_eligible, is_significant = region.pairs.get((left, right), (False, False))
            for matrix, value in ((eligible, is_eligible), (significant, is_significant)):
                matrix[left][right] = (
                    int(value) + matrix[left - 1][right] + matrix[left][right - 1]
                    - matrix[left - 1][right - 1]
                )
    if expansion_mode == "symmetric":
        bounds = [(min(k, max_left), min(k, max_right))
                  for k in range(1, max(max_left, max_right) + 1)]
    else:
        bounds = sorted(
            ((left, right) for left in range(1, max_left + 1)
             for right in range(1, max_right + 1)),
            key=lambda pair: (max(pair), sum(pair), pair),
        )
    common: Dict[str, object] = dict(
        region.metadata, expansion_mode=expansion_mode, p_column=p_column,
        p_threshold=threshold, min_eligible_pairs=min_eligible_pairs,
    )
    steps: List[Dict[str, object]] = []
    best_fraction: Optional[Fraction] = None
    best_indices: List[int] = []
    for left, right in bounds:
        n, s = eligible[left][right], significant[left][right]
        allowed = n >= min_eligible_pairs
        if allowed:
            fraction = Fraction(s, n)
            if best_fraction is None or fraction > best_fraction:
                best_fraction, best_indices = fraction, [len(steps)]
            elif fraction == best_fraction:
                best_indices.append(len(steps))
        steps.append(dict(
            common, significant_pair_fraction=s / n if n else None,
            left_max_rank=left, right_max_rank=right,
            n_eligible_pairs=n, n_significant_pairs=s,
            meets_min_eligible_pairs=allowed, is_maximum=False, is_selected_maximum=False,
        ))
    selected = None
    if best_indices:
        for index in best_indices:
            steps[index]["is_maximum"] = True
        # Prefer a better-supported maximum; break remaining ties near-to-far.
        selected = min(
            (steps[index] for index in best_indices),
            key=lambda row: (-row["n_eligible_pairs"],
                             row["left_max_rank"] + row["right_max_rank"],
                             max(row["left_max_rank"], row["right_max_rank"]),
                             row["left_max_rank"], row["right_max_rank"]),
        )
        selected["is_selected_maximum"] = True
    n, s = eligible[max_left][max_right], significant[max_left][max_right]
    status = ("no_eligible_pairs" if not n else
              "no_prefix_meets_min_eligible_pairs" if selected is None else "ok")
    summary = dict(
        common, overall_significant_pair_fraction=s / n if n else None,
        max_prefix_significant_pair_fraction=float(best_fraction) if best_fraction is not None else None,
        n_eligible_pairs=n, n_significant_pairs=s,
        best_left_max_rank=selected["left_max_rank"] if selected else None,
        best_right_max_rank=selected["right_max_rank"] if selected else None,
        best_n_eligible_pairs=selected["n_eligible_pairs"] if selected else None,
        best_n_significant_pairs=selected["n_significant_pairs"] if selected else None,
        n_maximizing_prefixes=len(best_indices), left_max_rank=max_left, right_max_rank=max_right,
        n_input_cross_pairs=len(region.pairs), n_ineligible_pairs=len(region.pairs) - n,
        n_unreported_pairs=unreported, n_ignored_same_side_pairs=region.ignored_same_side,
        n_prefixes=len(steps),
        n_testable_prefixes=sum(row["meets_min_eligible_pairs"] for row in steps), status=status,
    )
    return summary, steps


def output_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, ".17g")
    return str(value)


def write_rows(path: Path, fields: Tuple[str, ...], rows: List[Mapping[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: output_value(row[name]) for name in fields})


def run(args: argparse.Namespace) -> Dict[str, object]:
    if not math.isfinite(args.p_threshold) or not 0 < args.p_threshold <= 1:
        raise ValueError("--p-threshold must be finite and in (0, 1]")
    if args.min_eligible_pairs < 1:
        raise ValueError("--min-eligible-pairs must be >= 1")
    summary_path = Path(str(args.output_prefix) + ".summary.tsv")
    expansion_path = Path(str(args.output_prefix) + ".expansion.tsv")
    for path in (summary_path, expansion_path):
        if path.exists():
            raise ValueError(f"output already exists: {path}; choose a new --output-prefix")
    regions = read_pairs(args.pair_tsv, args.p_column, args.p_threshold)
    summaries: List[Dict[str, object]] = []
    expansions: List[Dict[str, object]] = []
    for region in regions.values():
        summary, steps = summarize_region(
            region, args.expansion_mode, args.min_eligible_pairs,
            args.p_threshold, args.p_column, args.allow_incomplete_grid,
        )
        summaries.append(summary)
        expansions.extend(steps)
    # Validate the entire input before creating either output.
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    write_rows(summary_path, SUMMARY_FIELDS, summaries)
    write_rows(expansion_path, EXPANSION_FIELDS, expansions)
    return dict(regions=len(summaries), ok_regions=sum(row["status"] == "ok" for row in summaries),
                summary=summary_path, expansion=expansion_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pair_tsv", type=Path, help="complete target_pairs.tsv (or .tsv.gz), one or more SCRs")
    parser.add_argument("--output-prefix", "-o", type=Path, required=True)
    parser.add_argument("--p-threshold", type=float, default=0.01, help="strict P < threshold (default: 0.01)")
    parser.add_argument("--p-column", default="pair_nominal_p", help="default: pair_nominal_p; never auto-use permutation P")
    parser.add_argument("--expansion-mode", choices=("symmetric", "independent"), default="symmetric",
                        help="symmetric: L1..Lk x R1..Rk; independent: L1..Li x R1..Rj (default: symmetric)")
    parser.add_argument("--min-eligible-pairs", type=int, default=1,
                        help="minimum denominator for a prefix to enter the maximum; overall fraction unchanged (default: 1)")
    parser.add_argument("--allow-incomplete-grid", action="store_true",
                        help="explicitly treat absent rank pairs as ineligible; do NOT use for matched-only/significant-only tables")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result = run(args)
    except (OSError, ValueError, csv.Error) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(" ".join(f"{key}={value}" for key, value in result.items()), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
