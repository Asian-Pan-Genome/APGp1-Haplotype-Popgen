#!/usr/bin/env python3
"""Calculate LD between phased-VCF SNP blocks flanking BED intervals.

The primary statistic is the multiallelic D' used by Bansal et al. (2007):

    D'_AB = sum_i sum_j p_i q_j |D'_ij|

where the alleles of each block are its observed multi-SNP haplotypes. The
script also reports cross-block pairwise r^2 summaries and the RV coefficient
so results can be compared with earlier pairwise-LD analyses.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import random
import statistics
import sys
from bisect import bisect_right
from collections import Counter, defaultdict
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Set, TextIO, Tuple


HaplotypeKey = Tuple[str, int]
Haplotype = Tuple[int, ...]


@dataclass(frozen=True)
class Region:
    chrom: str
    display_chrom: str
    start: int
    end: int
    region_id: str
    line_number: int


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    variant_id: str
    ref: str
    alt: str
    maf: float
    call_rate: float
    alleles: Mapping[HaplotypeKey, int]
    line_number: int

    @property
    def bed_pos(self) -> int:
        return self.pos - 1


@dataclass(frozen=True)
class SelectedRegion:
    region: Region
    left: Tuple[Variant, ...]
    right: Tuple[Variant, ...]


@dataclass(frozen=True)
class VCFData:
    sample_count: int
    variants_by_chrom: Mapping[str, Tuple[Variant, ...]]
    pre_missing_snp_positions_by_chrom: Mapping[str, Tuple[int, ...]]
    missing_filtered_snp_positions_by_chrom: Mapping[str, Tuple[int, ...]]
    candidate_variant_count: int
    filtered_non_snp_count: int
    filtered_vcf_filter_count: int
    filtered_maf_count: int
    filtered_call_rate_count: int
    filtered_unphased_count: int
    maf_missing_qc_skipped: bool = False


def normalize_chrom(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value.upper()


def open_text(path: Path) -> TextIO:
    if path.suffix in (".gz", ".bgz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_regions(path: Path) -> List[Region]:
    regions: List[Region] = []
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{line_number}: expected at least three BED columns"
                )
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: BED start/end must be integers"
                ) from error
            if start < 0 or end <= start:
                raise ValueError(
                    f"{path}:{line_number}: invalid BED interval {start}-{end}"
                )
            region_id = fields[3] if len(fields) >= 4 else f"{fields[0]}:{start}-{end}"
            regions.append(
                Region(
                    normalize_chrom(fields[0]),
                    fields[0],
                    start,
                    end,
                    region_id,
                    line_number,
                )
            )
    if not regions:
        raise ValueError(f"{path}: no BED intervals found")
    return regions


def merge_candidate_intervals(
    regions: Sequence[Region],
    flank_bp: int,
    include_region_interiors: bool = False,
) -> Dict[str, Tuple[Tuple[int, int], ...]]:
    raw: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for region in regions:
        raw[region.chrom].append((max(0, region.start - flank_bp), region.start))
        if include_region_interiors:
            raw[region.chrom].append((region.start, region.end))
        raw[region.chrom].append((region.end, region.end + flank_bp))

    merged_by_chrom: Dict[str, Tuple[Tuple[int, int], ...]] = {}
    for chrom, intervals in raw.items():
        merged: List[List[int]] = []
        for start, end in sorted(intervals):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        merged_by_chrom[chrom] = tuple((start, end) for start, end in merged)
    return merged_by_chrom


def position_in_intervals(position: int, intervals: Sequence[Tuple[int, int]]) -> bool:
    if not intervals:
        return False
    starts = [interval[0] for interval in intervals]
    index = bisect_right(starts, position) - 1
    return index >= 0 and position < intervals[index][1]


def parse_gt(gt: str) -> Tuple[Tuple[Optional[int], ...], bool]:
    if gt in ("", "."):
        return (None,), True
    if "|" in gt:
        tokens = gt.split("|")
        phased = True
    elif "/" in gt:
        tokens = gt.split("/")
        phased = False
    else:
        tokens = [gt]
        phased = True

    alleles: List[Optional[int]] = []
    for token in tokens:
        if token == ".":
            alleles.append(None)
            continue
        try:
            allele = int(token)
        except ValueError as error:
            raise ValueError(f"invalid GT allele {token!r}") from error
        if allele not in (0, 1):
            raise ValueError(f"biallelic SNP has unsupported GT allele {allele}")
        alleles.append(allele)
    return tuple(alleles), phased


def read_vcf(
    path: Path,
    regions: Sequence[Region],
    flank_bp: int,
    min_maf: float,
    min_call_rate: float,
    include_filtered: bool,
    allow_unphased_missing: bool,
    include_region_interiors: bool = False,
    progress_callback: Optional[Callable[[Mapping[str, int]], None]] = None,
    skip_maf_missing_qc: bool = False,
    vcf_lines: Optional[Iterable[str]] = None,
) -> VCFData:
    """Read genotypes; optionally trust prior site MAF/missing filtering.

    Prefiltered mode does not compute MAF/call rate (Variant values are NaN),
    nor reconstruct pre-filter site counts. SNP/FILTER/phasing checks remain.
    """
    intervals_by_chrom = merge_candidate_intervals(
        regions,
        flank_bp,
        include_region_interiors=include_region_interiors,
    )
    variants_by_chrom: Dict[str, List[Variant]] = defaultdict(list)
    pre_missing_snp_positions_by_chrom: Dict[str, List[int]] = defaultdict(list)
    missing_filtered_snp_positions_by_chrom: Dict[str, List[int]] = defaultdict(list)
    interval_starts = {
        chrom: tuple(start for start, _end in intervals)
        for chrom, intervals in intervals_by_chrom.items()
    }
    sample_names: Optional[List[str]] = None
    candidate_variant_count = 0
    filtered_non_snp_count = 0
    filtered_vcf_filter_count = 0
    filtered_maf_count = 0
    filtered_call_rate_count = 0
    filtered_unphased_count = 0
    records_scanned = retained_variant_count = 0

    def report_progress(finished: bool = False) -> None:
        if progress_callback is not None:
            counts = dict(records_scanned=records_scanned,
                          candidate_variants=candidate_variant_count,
                          retained_snps=retained_variant_count,
                          filtered_missing=filtered_call_rate_count,
                          finished=int(finished))
            if skip_maf_missing_qc:
                counts["maf_missing_qc_skipped"] = 1
            progress_callback(counts)

    with (open_text(path) if vcf_lines is None else nullcontext(vcf_lines)) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.split("\t")
                if len(header) < 10:
                    raise ValueError(f"{path}:{line_number}: VCF has no sample columns")
                sample_names = header[9:]
                if len(sample_names) != len(set(sample_names)):
                    raise ValueError(f"{path}:{line_number}: duplicate VCF sample names")
                continue
            if line.startswith("#"):
                continue
            if sample_names is None:
                raise ValueError(f"{path}:{line_number}: data encountered before #CHROM header")

            # At this point all preceding records (including out-of-BED ones)
            # have been processed. The caller controls time-based throttling.
            if records_scanned and records_scanned % 1000 == 0:
                report_progress()
            records_scanned += 1

            fields = line.split("\t")
            if len(fields) != 9 + len(sample_names):
                raise ValueError(
                    f"{path}:{line_number}: expected {9 + len(sample_names)} VCF columns, "
                    f"found {len(fields)}"
                )
            chrom = normalize_chrom(fields[0])
            if chrom not in intervals_by_chrom:
                continue
            try:
                pos = int(fields[1])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid VCF POS") from error
            if pos < 1:
                raise ValueError(f"{path}:{line_number}: VCF POS must be >= 1")
            starts = interval_starts[chrom]
            interval_index = bisect_right(starts, pos - 1) - 1
            if (
                interval_index < 0
                or pos - 1 >= intervals_by_chrom[chrom][interval_index][1]
            ):
                continue

            candidate_variant_count += 1
            ref = fields[3]
            alts = fields[4].split(",")
            if len(ref) != 1 or len(alts) != 1 or len(alts[0]) != 1 or alts[0] == ".":
                filtered_non_snp_count += 1
                continue
            if not include_filtered and fields[6] not in ("PASS", "."):
                filtered_vcf_filter_count += 1
                continue

            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise ValueError(f"{path}:{line_number}: FORMAT has no GT field")
            gt_index = format_fields.index("GT")
            allele_counts = Counter()
            allele_by_haplotype: Dict[HaplotypeKey, int] = {}
            fully_called_samples = 0
            unphased_heterozygotes = 0

            for sample_name, sample_field in zip(sample_names, fields[9:]):
                sample_values = sample_field.split(":")
                gt = sample_values[gt_index] if gt_index < len(sample_values) else "."
                try:
                    alleles, phased = parse_gt(gt)
                except ValueError as error:
                    raise ValueError(f"{path}:{line_number}: {sample_name}: {error}") from error
                called = [allele for allele in alleles if allele is not None]
                if not skip_maf_missing_qc:
                    allele_counts.update(called)
                    if called and len(called) == len(alleles):
                        fully_called_samples += 1

                is_unphased_heterozygote = (
                    not phased
                    and len(called) == len(alleles)
                    and len(set(called)) > 1
                )
                if is_unphased_heterozygote:
                    unphased_heterozygotes += 1
                    continue

                for haplotype_index, allele in enumerate(alleles):
                    if allele is not None:
                        allele_by_haplotype[(sample_name, haplotype_index)] = allele

            if skip_maf_missing_qc:
                maf = call_rate = math.nan
            else:
                pre_missing_snp_positions_by_chrom[chrom].append(pos - 1)
                total_called = sum(allele_counts.values())
                if total_called == 0:
                    filtered_call_rate_count += 1
                    missing_filtered_snp_positions_by_chrom[chrom].append(pos - 1)
                    continue
                call_rate = fully_called_samples / len(sample_names)
                if call_rate < min_call_rate:
                    filtered_call_rate_count += 1
                    missing_filtered_snp_positions_by_chrom[chrom].append(pos - 1)
                    continue
                alt_frequency = allele_counts[1] / total_called
                maf = min(alt_frequency, 1.0 - alt_frequency)
                if maf < min_maf or maf == 0.0:
                    filtered_maf_count += 1
                    continue
            if unphased_heterozygotes and not allow_unphased_missing:
                filtered_unphased_count += 1
                continue

            variant_id = fields[2]
            if variant_id == ".":
                variant_id = f"{fields[0]}:{pos}:{ref}:{alts[0]}"
            variants_by_chrom[chrom].append(
                Variant(
                    chrom=chrom,
                    pos=pos,
                    variant_id=variant_id,
                    ref=ref,
                    alt=alts[0],
                    maf=maf,
                    call_rate=call_rate,
                    alleles=allele_by_haplotype,
                    line_number=line_number,
                )
            )
            retained_variant_count += 1

    if sample_names is None:
        raise ValueError(f"{path}: no #CHROM header found")
    sorted_variants = {
        chrom: tuple(sorted(variants, key=lambda variant: (variant.pos, variant.variant_id)))
        for chrom, variants in variants_by_chrom.items()
    }
    report_progress(finished=True)
    return VCFData(
        sample_count=len(sample_names),
        variants_by_chrom=sorted_variants,
        pre_missing_snp_positions_by_chrom={
            chrom: tuple(sorted(positions))
            for chrom, positions in pre_missing_snp_positions_by_chrom.items()
        },
        missing_filtered_snp_positions_by_chrom={
            chrom: tuple(sorted(positions))
            for chrom, positions in missing_filtered_snp_positions_by_chrom.items()
        },
        candidate_variant_count=candidate_variant_count,
        filtered_non_snp_count=filtered_non_snp_count,
        filtered_vcf_filter_count=filtered_vcf_filter_count,
        filtered_maf_count=filtered_maf_count,
        filtered_call_rate_count=filtered_call_rate_count,
        filtered_unphased_count=filtered_unphased_count,
        maf_missing_qc_skipped=skip_maf_missing_qc,
    )


def common_haplotype_keys(variants: Sequence[Variant]) -> Tuple[HaplotypeKey, ...]:
    if not variants:
        return ()
    common: Set[HaplotypeKey] = set(variants[0].alleles)
    for variant in variants[1:]:
        common.intersection_update(variant.alleles)
    return tuple(sorted(common))


def haplotypes_for_keys(
    variants: Sequence[Variant], keys: Sequence[HaplotypeKey]
) -> List[Haplotype]:
    return [tuple(variant.alleles[key] for variant in variants) for key in keys]


def haplotype_entropy(haplotypes: Sequence[Haplotype]) -> float:
    if not haplotypes:
        return math.nan
    counts = Counter(haplotypes)
    total = len(haplotypes)
    return -sum((count / total) * math.log(count / total) for count in counts.values())


def select_greedy_entropy(
    candidates: Sequence[Variant], block_snps: int
) -> Tuple[Variant, ...]:
    """Greedily maximize haplotype entropy; ties retain boundary proximity."""

    remaining = list(candidates)
    chosen: List[Variant] = []
    while remaining and len(chosen) < block_snps:
        best_index = 0
        best_score = (-math.inf, -1, 0)
        for index, candidate in enumerate(remaining):
            trial = chosen + [candidate]
            keys = common_haplotype_keys(trial)
            entropy = haplotype_entropy(haplotypes_for_keys(trial, keys))
            if math.isnan(entropy):
                entropy = -math.inf
            score = (entropy, len(keys), -index)
            if score > best_score:
                best_score = score
                best_index = index
        chosen.append(remaining.pop(best_index))
    return tuple(chosen)


def choose_block(
    candidates_nearest_first: Sequence[Variant],
    block_snps: int,
    selection_method: str,
) -> Tuple[Variant, ...]:
    if selection_method == "nearest":
        return tuple(candidates_nearest_first[:block_snps])
    if selection_method == "greedy-entropy":
        return select_greedy_entropy(candidates_nearest_first, block_snps)
    raise ValueError(f"unsupported selection method: {selection_method}")


def select_flanking_variants(
    regions: Sequence[Region],
    variants_by_chrom: Mapping[str, Sequence[Variant]],
    flank_bp: int,
    block_snps: int,
    selection_method: str,
) -> List[SelectedRegion]:
    selected: List[SelectedRegion] = []
    for region in regions:
        chromosome_variants = variants_by_chrom.get(region.chrom, ())
        left_candidates = sorted(
            (
                variant
                for variant in chromosome_variants
                if region.start - flank_bp <= variant.bed_pos < region.start
            ),
            key=lambda variant: (-variant.pos, variant.variant_id),
        )
        right_candidates = sorted(
            (
                variant
                for variant in chromosome_variants
                if region.end <= variant.bed_pos < region.end + flank_bp
            ),
            key=lambda variant: (variant.pos, variant.variant_id),
        )
        selected.append(
            SelectedRegion(
                region,
                choose_block(left_candidates, block_snps, selection_method),
                choose_block(right_candidates, block_snps, selection_method),
            )
        )
    return selected


def multiallelic_dprime(
    left_haplotypes: Sequence[Haplotype],
    right_haplotypes: Sequence[Haplotype],
) -> float:
    if len(left_haplotypes) != len(right_haplotypes) or not left_haplotypes:
        return math.nan
    total = len(left_haplotypes)
    left_counts = Counter(left_haplotypes)
    right_counts = Counter(right_haplotypes)
    joint_counts = Counter(zip(left_haplotypes, right_haplotypes))
    overall = 0.0

    for left_haplotype, left_count in left_counts.items():
        p_i = left_count / total
        for right_haplotype, right_count in right_counts.items():
            q_j = right_count / total
            h_ij = joint_counts[(left_haplotype, right_haplotype)] / total
            d_ij = h_ij - p_i * q_j
            if d_ij < 0.0:
                d_max = min(p_i * q_j, (1.0 - p_i) * (1.0 - q_j))
            else:
                d_max = min(p_i * (1.0 - q_j), (1.0 - p_i) * q_j)
            if d_max > 0.0:
                overall += p_i * q_j * abs(d_ij / d_max)
    return overall


def pairwise_r2(
    first: Variant,
    second: Variant,
    haplotype_keys: Sequence[HaplotypeKey],
) -> Optional[float]:
    if not haplotype_keys:
        return None
    total = len(haplotype_keys)
    p_a = sum(first.alleles[key] for key in haplotype_keys) / total
    p_b = sum(second.alleles[key] for key in haplotype_keys) / total
    p_ab = sum(
        first.alleles[key] * second.alleles[key] for key in haplotype_keys
    ) / total
    denominator = p_a * (1.0 - p_a) * p_b * (1.0 - p_b)
    if denominator <= 0.0:
        return None
    return ((p_ab - p_a * p_b) ** 2) / denominator


def pairwise_values(
    first: Sequence[Variant],
    second: Sequence[Variant],
    haplotype_keys: Sequence[HaplotypeKey],
    same_block: bool = False,
) -> Tuple[List[float], int]:
    values: List[float] = []
    if same_block:
        expected = len(first) * (len(first) - 1) // 2
        for index, variant_a in enumerate(first):
            for variant_b in first[index + 1 :]:
                value = pairwise_r2(variant_a, variant_b, haplotype_keys)
                if value is not None:
                    values.append(value)
        return values, expected

    expected = len(first) * len(second)
    for variant_a in first:
        for variant_b in second:
            value = pairwise_r2(variant_a, variant_b, haplotype_keys)
            if value is not None:
                values.append(value)
    return values, expected


def quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = probability * (len(ordered) - 1)
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def trimmed_mean(values: Sequence[float], proportion: float = 0.1) -> float:
    ordered = sorted(values)
    trim = math.floor(len(ordered) * proportion)
    if trim and 2 * trim < len(ordered):
        ordered = ordered[trim:-trim]
    return statistics.fmean(ordered)


def permutation_statistics(
    left_haplotypes: Sequence[Haplotype],
    right_haplotypes: Sequence[Haplotype],
    observed: float,
    permutations: int,
    rng: random.Random,
) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if permutations <= 0:
        return None, None, None
    shuffled = list(right_haplotypes)
    null_values: List[float] = []
    exceedances = 0
    for _ in range(permutations):
        rng.shuffle(shuffled)
        value = multiallelic_dprime(left_haplotypes, shuffled)
        null_values.append(value)
        if value >= observed:
            exceedances += 1
    null_mean = statistics.fmean(null_values)
    empirical_p = (exceedances + 1) / (permutations + 1)
    return null_mean, observed - null_mean, empirical_p


def format_float(value: Optional[float]) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.8g}"


def format_variants(variants: Sequence[Variant]) -> str:
    return ",".join(f"{variant.variant_id}@{variant.pos}" for variant in variants) or "."


def summarize_selection(
    selection: SelectedRegion,
    sample_count: int,
    block_snps: int,
    allow_partial_blocks: bool,
    selection_method: str,
    permutations: int,
    rng: random.Random,
) -> Dict[str, str]:
    left = selection.left
    right = selection.right
    complete_blocks = len(left) == block_snps and len(right) == block_snps
    usable_blocks = complete_blocks or (
        allow_partial_blocks and len(left) > 0 and len(right) > 0
    )
    all_variants = left + right
    keys = common_haplotype_keys(all_variants) if usable_blocks else ()
    left_haplotypes = haplotypes_for_keys(left, keys)
    right_haplotypes = haplotypes_for_keys(right, keys)

    dprime: Optional[float] = None
    left_entropy: Optional[float] = None
    right_entropy: Optional[float] = None
    perm_mean: Optional[float] = None
    dprime_excess: Optional[float] = None
    empirical_p: Optional[float] = None
    if keys:
        dprime = multiallelic_dprime(left_haplotypes, right_haplotypes)
        left_entropy = haplotype_entropy(left_haplotypes)
        right_entropy = haplotype_entropy(right_haplotypes)
        perm_mean, dprime_excess, empirical_p = permutation_statistics(
            left_haplotypes,
            right_haplotypes,
            dprime,
            permutations,
            rng,
        )

    cross_values, cross_expected = pairwise_values(left, right, keys)
    left_values, left_expected = pairwise_values(left, left, keys, same_block=True)
    right_values, right_expected = pairwise_values(right, right, keys, same_block=True)

    cross_mean: Optional[float] = None
    cross_median: Optional[float] = None
    cross_trimmed_mean: Optional[float] = None
    cross_q90: Optional[float] = None
    cross_max: Optional[float] = None
    if cross_values:
        cross_mean = statistics.fmean(cross_values)
        cross_median = statistics.median(cross_values)
        cross_trimmed_mean = trimmed_mean(cross_values)
        cross_q90 = quantile(cross_values, 0.9)
        cross_max = max(cross_values)

    rv_r2: Optional[float] = None
    if (
        len(cross_values) == cross_expected
        and len(left_values) == left_expected
        and len(right_values) == right_expected
        and cross_expected > 0
    ):
        left_frobenius_sq = len(left) + 2.0 * sum(left_values)
        right_frobenius_sq = len(right) + 2.0 * sum(right_values)
        denominator = math.sqrt(left_frobenius_sq * right_frobenius_sq)
        if denominator > 0.0:
            rv_r2 = sum(cross_values) / denominator

    statuses: List[str] = []
    if len(left) < block_snps:
        statuses.append("insufficient_left_snps")
    if len(right) < block_snps:
        statuses.append("insufficient_right_snps")
    if usable_blocks and not keys:
        statuses.append("no_complete_phased_haplotypes")
    if keys and len(cross_values) < cross_expected:
        statuses.append("undefined_pairwise_r2")
    if not statuses:
        statuses.append("ok")

    region = selection.region
    return {
        "chrom": region.display_chrom,
        "start": str(region.start),
        "end": str(region.end),
        "region_id": region.region_id,
        "selection_method": selection_method,
        "sample_count": str(sample_count),
        "left_snp_count": str(len(left)),
        "right_snp_count": str(len(right)),
        "left_snps_nearest_or_selected_order": format_variants(left),
        "right_snps_nearest_or_selected_order": format_variants(right),
        "complete_haplotype_count": str(len(keys)),
        "left_haplotype_count": str(len(set(left_haplotypes))),
        "right_haplotype_count": str(len(set(right_haplotypes))),
        "joint_haplotype_count": str(len(set(zip(left_haplotypes, right_haplotypes)))),
        "left_haplotype_entropy": format_float(left_entropy),
        "right_haplotype_entropy": format_float(right_entropy),
        "multiallelic_dprime": format_float(dprime),
        "dprime_permutation_mean": format_float(perm_mean),
        "dprime_excess": format_float(dprime_excess),
        "dprime_empirical_p": format_float(empirical_p),
        "cross_r2_pairs_observed": str(len(cross_values)),
        "cross_r2_pairs_expected": str(cross_expected),
        "cross_mean_r2": format_float(cross_mean),
        "cross_median_r2": format_float(cross_median),
        "cross_trimmed_mean_r2": format_float(cross_trimmed_mean),
        "cross_q90_r2": format_float(cross_q90),
        "cross_max_r2": format_float(cross_max),
        "block_rv_r2": format_float(rv_r2),
        "status": ";".join(statuses),
    }


OUTPUT_FIELDS = (
    "chrom",
    "start",
    "end",
    "region_id",
    "selection_method",
    "sample_count",
    "left_snp_count",
    "right_snp_count",
    "left_snps_nearest_or_selected_order",
    "right_snps_nearest_or_selected_order",
    "complete_haplotype_count",
    "left_haplotype_count",
    "right_haplotype_count",
    "joint_haplotype_count",
    "left_haplotype_entropy",
    "right_haplotype_entropy",
    "multiallelic_dprime",
    "dprime_permutation_mean",
    "dprime_excess",
    "dprime_empirical_p",
    "cross_r2_pairs_observed",
    "cross_r2_pairs_expected",
    "cross_mean_r2",
    "cross_median_r2",
    "cross_trimmed_mean_r2",
    "cross_q90_r2",
    "cross_max_r2",
    "block_rv_r2",
    "status",
)


def run(args: argparse.Namespace) -> Dict[str, int]:
    regions = read_regions(args.regions_bed)
    vcf_data = read_vcf(
        args.vcf,
        regions,
        flank_bp=args.flank_bp,
        min_maf=args.min_maf,
        min_call_rate=args.min_call_rate,
        include_filtered=args.include_filtered,
        allow_unphased_missing=args.allow_unphased_missing,
    )
    selected = select_flanking_variants(
        regions,
        vcf_data.variants_by_chrom,
        flank_bp=args.flank_bp,
        block_snps=args.block_snps,
        selection_method=args.selection,
    )
    summaries = [
        summarize_selection(
            selection,
            sample_count=vcf_data.sample_count,
            block_snps=args.block_snps,
            allow_partial_blocks=args.allow_partial_blocks,
            selection_method=args.selection,
            permutations=args.permutations,
            rng=random.Random(args.seed + index),
        )
        for index, selection in enumerate(selected)
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    return {
        "region_count": len(regions),
        "ok_region_count": sum(row["status"] == "ok" for row in summaries),
        "sample_count": vcf_data.sample_count,
        "candidate_variant_count": vcf_data.candidate_variant_count,
        "retained_variant_count": sum(
            len(variants) for variants in vcf_data.variants_by_chrom.values()
        ),
        "filtered_non_snp_count": vcf_data.filtered_non_snp_count,
        "filtered_vcf_filter_count": vcf_data.filtered_vcf_filter_count,
        "filtered_maf_count": vcf_data.filtered_maf_count,
        "filtered_call_rate_count": vcf_data.filtered_call_rate_count,
        "filtered_unphased_count": vcf_data.filtered_unphased_count,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate Bansal et al. multiallelic D' between phased-VCF SNP "
            "blocks flanking BED intervals."
        )
    )
    parser.add_argument("regions_bed", type=Path, help="BED3/BED4 regions (0-based, half-open)")
    parser.add_argument("vcf", type=Path, help="phased VCF, optionally gzip-compressed")
    parser.add_argument("--output", "-o", type=Path, required=True, help="output TSV")
    parser.add_argument(
        "--block-snps",
        type=int,
        default=10,
        help="SNPs selected in each flanking block (default: 10)",
    )
    parser.add_argument(
        "--flank-bp",
        type=int,
        default=100_000,
        help="maximum flank length on each side (default: 100000)",
    )
    parser.add_argument(
        "--min-maf",
        type=float,
        default=0.1,
        help="minimum SNP minor-allele frequency (default: 0.1, as in the paper)",
    )
    parser.add_argument(
        "--min-call-rate",
        type=float,
        default=0.9,
        help="minimum fully called sample fraction per SNP (default: 0.9)",
    )
    parser.add_argument(
        "--selection",
        choices=("nearest", "greedy-entropy"),
        default="nearest",
        help=(
            "choose boundary-nearest SNPs or greedily maximize block haplotype "
            "entropy (default: nearest)"
        ),
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=0,
        help="right-block haplotype permutations for finite-sample calibration (default: 0)",
    )
    parser.add_argument("--seed", type=int, default=1, help="permutation seed (default: 1)")
    parser.add_argument(
        "--include-filtered",
        action="store_true",
        help="include variants whose VCF FILTER is neither PASS nor '.'",
    )
    parser.add_argument(
        "--allow-unphased-missing",
        action="store_true",
        help=(
            "retain SNPs containing unphased heterozygotes and exclude those "
            "sample haplotypes from complete-block calculations"
        ),
    )
    parser.add_argument(
        "--allow-partial-blocks",
        action="store_true",
        help="calculate LD when a flank contains fewer than --block-snps SNPs",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.block_snps < 1:
        parser.error("--block-snps must be >= 1")
    if args.flank_bp < 1:
        parser.error("--flank-bp must be >= 1")
    if not 0.0 <= args.min_maf <= 0.5:
        parser.error("--min-maf must be in [0, 0.5]")
    if not 0.0 <= args.min_call_rate <= 1.0:
        parser.error("--min-call-rate must be in [0, 1]")
    if args.permutations < 0:
        parser.error("--permutations must be >= 0")
    try:
        summary = run(args)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(
        "regions={region_count} ok={ok_region_count} samples={sample_count} "
        "candidate_variants={candidate_variant_count} retained_variants={retained_variant_count} "
        "filtered_non_snp={filtered_non_snp_count} filtered_vcf_filter={filtered_vcf_filter_count} "
        "filtered_maf={filtered_maf_count} filtered_call_rate={filtered_call_rate_count} "
        "filtered_unphased={filtered_unphased_count}".format(**summary),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
