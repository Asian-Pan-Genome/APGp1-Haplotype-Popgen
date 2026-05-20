#!/usr/bin/env python3
"""
Summarize Locityper leave-one-out genotyping performance at target loci.

This script implements the locus-level categorization used to compare an
expanded assembly-derived panel (APGp1+) against a baseline panel (APGp1-):

1. Filter evaluation records to haplotype-level calls.
2. Assign QV bins and no-call status.
3. Pair haplotypes between APGp1+ and APGp1- leave-one-out evaluations.
4. Compute Delta QV and fold error reduction (FER).
5. Classify loci into SR-infeasible, ceiling-effected, or panel/ancestry-dependent.
6. Export clean tables and optional publication-style summary plots.

Expected input columns in each Locityper evaluation table:
    sample, locus, avail_qv, avail_div
Optional columns used when present:
    query_type, genotype, qv, div

The optional qv/div columns are useful for samples that are not represented in a
specific leave-one-out panel and should therefore use direct genotyping accuracy
instead of haplotype availability. This behavior is controlled by
--panel-minus-use-call-prefixes.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


QV_BIN_EDGES = [-np.inf, 20, 25, 30, 35, 40, 45, np.inf]
QV_BIN_LABELS_LOW_TO_HIGH = ["<20", "20-25", "25-30", "30-35", "35-40", "40-45", ">45"]
QV_BIN_LABELS_HIGH_TO_LOW = [">45", "40-45", "35-40", "30-35", "25-30", "20-25", "<20", "No call"]

CATEGORY_ORDER = ["SR-infeasible", "Ceiling-effected", "Panel/ancestry-dependent"]
CATEGORY_COLORS = {
    "SR-infeasible": "grey",
    "Ceiling-effected": "orange",
    "Panel/ancestry-dependent": "tab:blue",
}


class InputError(ValueError):
    """Raised when an input file is missing required columns or valid records."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Clean Locityper locus-level performance analysis for comparing "
            "APGp1+ and APGp1- leave-one-out panels."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--panel-plus-eval", required=True, help="Locityper evaluation table for the APGp1+ LOO panel.")
    parser.add_argument("--panel-minus-eval", required=True, help="Locityper evaluation table for the APGp1- LOO panel.")
    parser.add_argument("--sample-pop", required=True, help="Two-column sample-to-population table: sample<TAB>population.")
    parser.add_argument("--outdir", default="locityper_performance", help="Output directory.")
    parser.add_argument("--prefix", default="locityper", help="Prefix for output files.")

    parser.add_argument("--sep", default="\t", help="Field separator for Locityper evaluation tables.")
    parser.add_argument("--sample-pop-sep", default="\t", help="Field separator for the sample population file.")
    parser.add_argument("--sample-pop-has-header", action="store_true", help="Set if --sample-pop includes a header row.")

    parser.add_argument("--plus-label", default="APGp1+ (LOO)", help="Display label for the expanded panel.")
    parser.add_argument("--minus-label", default="APGp1- (LOO)", help="Display label for the baseline panel.")
    parser.add_argument("--eas-label", default="EAS", help="Population label used to define EAS samples.")
    parser.add_argument("--non-eas-label", default="Non-EAS", help="Display label for non-EAS samples.")

    parser.add_argument(
        "--panel-minus-use-call-prefixes",
        default="",
        help=(
            "Comma-separated sample-name prefixes for which the APGp1- table should use qv/div "
            "instead of avail_qv/avail_div. Example: C,K. Leave empty to disable."
        ),
    )
    parser.add_argument(
        "--panel-plus-use-call-prefixes",
        default="",
        help="Same as --panel-minus-use-call-prefixes, but applied to the APGp1+ table.",
    )

    parser.add_argument("--qv-cap", type=float, default=45.0, help="QV cap used for Delta QV and FER calculations.")
    parser.add_argument("--no-call-threshold", type=float, default=0.05, help="No-call fraction threshold for SR-infeasible loci.")
    parser.add_argument("--delta-qv-threshold", type=float, default=1.0, help="Median EAS Delta QV threshold separating ceiling-effected and panel-dependent loci.")
    parser.add_argument("--random-seed", type=int, default=42, help="Random seed used when selecting EAS samples.")
    parser.add_argument("--eas-samples", default=None, help="Optional one-column file with preselected EAS samples. If omitted, EAS samples are sampled to match the non-EAS sample count.")
    parser.add_argument(
        "--use-all-eas-for-qv-summary",
        action="store_true",
        help="Use all EAS samples for QV-bin summaries instead of the matched EAS subset.",
    )

    parser.add_argument("--make-plots", action="store_true", help="Generate scatter and Delta QV/FER summary plots.")
    parser.add_argument("--make-qv-stack-plot", action="store_true", help="Generate a large per-locus stacked QV-bin plot.")
    parser.add_argument("--plot-format", default="pdf", choices=["pdf", "png", "svg"], help="Plot output format.")
    parser.add_argument("--stack-ncols", type=int, default=10, help="Number of columns for the stacked QV-bin locus grid.")
    parser.add_argument("--max-label-loci", type=int, default=40, help="Maximum number of loci to label in the scatter plot.")

    return parser.parse_args()


def parse_prefixes(prefix_text: str) -> tuple[str, ...]:
    return tuple(prefix.strip() for prefix in prefix_text.split(",") if prefix.strip())


def mkdir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def require_columns(df: pd.DataFrame, required: Sequence[str], source: str) -> None:
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise InputError(f"{source} is missing required column(s): {', '.join(missing)}")


def read_eval_table(path: str | Path, sep: str, panel_label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep, comment="#", header=0)
    require_columns(df, ["sample", "locus"], str(path))

    if "query_type" in df.columns:
        df = df[df["query_type"] != "gt"].copy()
    if "genotype" in df.columns:
        df = df[df["genotype"] != "*"].copy()

    if df.empty:
        raise InputError(f"{path} has no records after query_type/genotype filtering.")

    df["panel"] = panel_label
    df["hap_index"] = df.groupby(["sample", "locus"], sort=False).cumcount()
    return df


def read_sample_pop(path: str | Path, sep: str, has_header: bool) -> pd.DataFrame:
    if has_header:
        df = pd.read_csv(path, sep=sep)
        require_columns(df, ["sample", "population"], str(path))
        df = df[["sample", "population"]].copy()
    else:
        df = pd.read_csv(path, sep=sep, header=None, names=["sample", "population"], usecols=[0, 1])

    if df["sample"].duplicated().any():
        duplicated = df.loc[df["sample"].duplicated(), "sample"].head(5).tolist()
        raise InputError(f"Sample population file contains duplicated samples, for example: {duplicated}")
    return df


def starts_with_any(series: pd.Series, prefixes: Sequence[str]) -> pd.Series:
    if not prefixes:
        return pd.Series(False, index=series.index)
    return series.astype(str).apply(lambda value: value.startswith(tuple(prefixes)))


def choose_eval_metrics(
    df: pd.DataFrame,
    fallback_prefixes: Sequence[str],
    source_name: str,
) -> pd.DataFrame:
    """Create eval_qv/eval_div columns from availability or direct genotyping metrics."""
    out = df.copy()

    if "avail_qv" in out.columns:
        out["eval_qv"] = out["avail_qv"]
    elif "qv" in out.columns:
        out["eval_qv"] = out["qv"]
    else:
        raise InputError(f"{source_name} must contain avail_qv or qv.")

    if "avail_div" in out.columns:
        out["eval_div"] = out["avail_div"]
    elif "div" in out.columns:
        out["eval_div"] = out["div"]
    else:
        raise InputError(f"{source_name} must contain avail_div or div.")

    use_direct = starts_with_any(out["sample"], fallback_prefixes)
    if use_direct.any():
        if "qv" not in out.columns or "div" not in out.columns:
            raise InputError(
                f"{source_name} was asked to use qv/div for prefixes {fallback_prefixes}, "
                "but qv and/or div are not available."
            )
        out.loc[use_direct, "eval_qv"] = out.loc[use_direct, "qv"]
        out.loc[use_direct, "eval_div"] = out.loc[use_direct, "div"]

    out["eval_qv"] = pd.to_numeric(out["eval_qv"], errors="coerce")
    out["eval_div"] = pd.to_numeric(out["eval_div"], errors="coerce")
    return out


def add_qv_fields(df: pd.DataFrame, qv_cap: float) -> pd.DataFrame:
    out = df.copy()
    qv = out["eval_qv"].replace([np.inf, -np.inf], np.nan)

    out["is_no_call"] = qv.isna()
    out["qv_for_bin"] = out["eval_qv"].replace(np.inf, qv_cap + 1)
    out["qv_bin"] = pd.cut(
        out["qv_for_bin"],
        bins=QV_BIN_EDGES,
        labels=QV_BIN_LABELS_LOW_TO_HIGH,
        right=False,
    )
    out["qv_bin"] = out["qv_bin"].cat.add_categories(["No call"]).fillna("No call")

    out["qv_capped"] = out["eval_qv"].replace(np.inf, qv_cap).clip(upper=qv_cap)
    out["qv_capped_for_delta"] = out["qv_capped"].fillna(0.0)

    min_divergence = 10 ** (-qv_cap / 10.0)
    out["div_capped"] = out["eval_div"].clip(lower=min_divergence)
    return out


def prepare_panel(
    path: str | Path,
    sep: str,
    panel_label: str,
    fallback_prefixes: Sequence[str],
    qv_cap: float,
) -> pd.DataFrame:
    df = read_eval_table(path, sep=sep, panel_label=panel_label)
    df = choose_eval_metrics(df, fallback_prefixes=fallback_prefixes, source_name=str(path))
    df = add_qv_fields(df, qv_cap=qv_cap)
    return df


def add_population(df: pd.DataFrame, sample_pop: pd.DataFrame) -> pd.DataFrame:
    out = df.merge(sample_pop, on="sample", how="left")
    missing = out.loc[out["population"].isna(), "sample"].drop_duplicates().head(10).tolist()
    if missing:
        raise InputError(f"Population labels are missing for sample(s), for example: {missing}")
    return out


def select_eas_samples(
    merged: pd.DataFrame,
    eas_label: str,
    random_seed: int,
    eas_samples_file: str | None,
    out_path: Path,
) -> list[str]:
    non_eas_samples = sorted(merged.loc[merged["population"] != eas_label, "sample"].dropna().unique())
    eas_samples = sorted(merged.loc[merged["population"] == eas_label, "sample"].dropna().unique())

    if not non_eas_samples:
        raise InputError("No non-EAS samples were found after merging panel results.")
    if not eas_samples:
        raise InputError(f"No samples with population label '{eas_label}' were found after merging panel results.")

    target_n = len(non_eas_samples)

    if eas_samples_file:
        selected = pd.read_csv(eas_samples_file, header=None)[0].astype(str).tolist()
        selected = [sample for sample in selected if sample in set(eas_samples)]
        if not selected:
            raise InputError(f"No valid EAS samples from {eas_samples_file} were present in the merged table.")
    else:
        rng = np.random.default_rng(random_seed)
        replace = len(eas_samples) < target_n
        selected = rng.choice(eas_samples, size=target_n, replace=replace).tolist()

    pd.Series(selected).to_csv(out_path, index=False, header=False)
    return selected


def merge_panels(panel_plus: pd.DataFrame, panel_minus: pd.DataFrame, sample_pop: pd.DataFrame) -> pd.DataFrame:
    keep_cols = [
        "sample",
        "locus",
        "hap_index",
        "panel",
        "eval_qv",
        "eval_div",
        "qv_bin",
        "is_no_call",
        "qv_capped",
        "qv_capped_for_delta",
        "div_capped",
    ]
    plus = panel_plus[keep_cols].copy()
    minus = panel_minus[keep_cols].copy()

    merged = plus.merge(
        minus,
        on=["sample", "locus", "hap_index"],
        suffixes=("_plus", "_minus"),
        how="inner",
        validate="one_to_one",
    )
    if merged.empty:
        raise InputError("No haplotypes could be paired between APGp1+ and APGp1- tables.")

    merged = add_population(merged, sample_pop)
    merged["delta_qv"] = merged["qv_capped_for_delta_plus"] - merged["qv_capped_for_delta_minus"]
    merged["fold_error_reduction"] = merged["div_capped_minus"] / merged["div_capped_plus"]
    merged.loc[~np.isfinite(merged["fold_error_reduction"]), "fold_error_reduction"] = np.nan
    return merged


def classify_loci(
    panel_plus: pd.DataFrame,
    panel_minus: pd.DataFrame,
    merged: pd.DataFrame,
    selected_eas_samples: Sequence[str],
    plus_label: str,
    minus_label: str,
    no_call_threshold: float,
    delta_qv_threshold: float,
    eas_label: str,
) -> pd.DataFrame:
    combined = pd.concat([panel_plus, panel_minus], ignore_index=True)
    no_call = (
        combined.groupby(["locus", "panel"], observed=True)["is_no_call"]
        .mean()
        .unstack("panel")
        .rename(columns={plus_label: "no_call_frac_plus", minus_label: "no_call_frac_minus"})
    )

    for col in ["no_call_frac_plus", "no_call_frac_minus"]:
        if col not in no_call.columns:
            no_call[col] = np.nan

    no_call["sr_infeasible"] = (
        (no_call["no_call_frac_plus"] >= no_call_threshold)
        & (no_call["no_call_frac_minus"] >= no_call_threshold)
    )

    eas_subset = merged[merged["sample"].isin(set(selected_eas_samples))].copy()
    median_delta_eas = eas_subset.groupby("locus", observed=True)["delta_qv"].median().rename("median_delta_qv_eas")

    all_loci = sorted(set(panel_plus["locus"]).union(panel_minus["locus"]))
    out = pd.DataFrame(index=all_loci)
    out = out.join(no_call[["no_call_frac_plus", "no_call_frac_minus", "sr_infeasible"]], how="left")
    out = out.join(median_delta_eas, how="left")

    out["category"] = "Panel/ancestry-dependent"
    out.loc[out["sr_infeasible"].fillna(False), "category"] = "SR-infeasible"
    evaluable = ~out["sr_infeasible"].fillna(False)
    out.loc[evaluable & (out["median_delta_qv_eas"] <= delta_qv_threshold), "category"] = "Ceiling-effected"
    out.loc[evaluable & (out["median_delta_qv_eas"] > delta_qv_threshold), "category"] = "Panel/ancestry-dependent"

    out = out.reset_index(names="locus")
    out["category"] = pd.Categorical(out["category"], categories=CATEGORY_ORDER, ordered=True)
    return out


def summarize_loci(merged: pd.DataFrame, locus_categories: pd.DataFrame, eas_label: str) -> pd.DataFrame:
    def safe_median(series: pd.Series) -> float:
        return float(series.median()) if series.notna().any() else np.nan

    summaries = []
    for locus, group in merged.groupby("locus", observed=True):
        eas_group = group[group["population"] == eas_label]
        non_eas_group = group[group["population"] != eas_label]
        summaries.append(
            {
                "locus": locus,
                "n_haplotypes_paired": len(group),
                "n_samples_paired": group["sample"].nunique(),
                "median_qv_plus": safe_median(group["qv_capped_plus"]),
                "median_qv_minus": safe_median(group["qv_capped_minus"]),
                "median_delta_qv_all": safe_median(group["delta_qv"]),
                "median_delta_qv_eas_all": safe_median(eas_group["delta_qv"]),
                "median_delta_qv_non_eas": safe_median(non_eas_group["delta_qv"]),
                "median_fer_all": safe_median(group["fold_error_reduction"]),
                "median_fer_eas_all": safe_median(eas_group["fold_error_reduction"]),
                "median_fer_non_eas": safe_median(non_eas_group["fold_error_reduction"]),
            }
        )
    out = pd.DataFrame(summaries)
    category_cols = [
        "locus",
        "no_call_frac_plus",
        "no_call_frac_minus",
        "sr_infeasible",
        "median_delta_qv_eas",
        "category",
    ]
    out = out.merge(locus_categories[category_cols], on="locus", how="left")
    return out.sort_values(["category", "median_delta_qv_eas", "median_fer_all"], ascending=[True, False, False])


def make_qv_summary(
    panel_plus: pd.DataFrame,
    panel_minus: pd.DataFrame,
    sample_pop: pd.DataFrame,
    selected_eas_samples: Sequence[str],
    use_all_eas: bool,
    plus_label: str,
    minus_label: str,
    eas_label: str,
    non_eas_label: str,
) -> pd.DataFrame:
    combined = pd.concat([panel_plus, panel_minus], ignore_index=True)
    combined = add_population(combined, sample_pop)

    rows = []
    groups = [("All", combined)]
    if use_all_eas:
        groups.append((eas_label, combined[combined["population"] == eas_label]))
    else:
        groups.append((eas_label, combined[combined["sample"].isin(set(selected_eas_samples))]))
    groups.append((non_eas_label, combined[combined["population"] != eas_label]))

    for population_group, group_df in groups:
        if group_df.empty:
            continue
        for (locus, panel), sub in group_df.groupby(["locus", "panel"], observed=True):
            row = {
                "locus": locus,
                "population_group": population_group,
                "panel": panel,
                "n_haplotypes": len(sub),
                "n_samples": sub["sample"].nunique(),
            }
            for qv_bin in QV_BIN_LABELS_HIGH_TO_LOW:
                count = int((sub["qv_bin"] == qv_bin).sum())
                row[f"count_qv_{qv_bin}"] = count
                row[f"fraction_qv_{qv_bin}"] = count / len(sub) if len(sub) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def write_locus_lists(locus_categories: pd.DataFrame, outdir: Path, prefix: str) -> None:
    file_names = {
        "SR-infeasible": "sr_infeasible_loci.txt",
        "Ceiling-effected": "ceiling_effected_loci.txt",
        "Panel/ancestry-dependent": "panel_ancestry_dependent_loci.txt",
    }
    for category, file_name in file_names.items():
        loci = locus_categories.loc[locus_categories["category"].astype(str) == category, "locus"].sort_values()
        loci.to_csv(outdir / f"{prefix}.{file_name}", index=False, header=False)


def apply_plot_defaults() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["mathtext.fontset"] = "custom"
    plt.rcParams["mathtext.rm"] = "Arial"
    plt.rcParams["mathtext.it"] = "Arial:italic"
    plt.rcParams["mathtext.bf"] = "Arial:bold"
    plt.rcParams["axes.unicode_minus"] = False


def category_palette(locus_categories: pd.DataFrame) -> dict[str, str]:
    return {
        row.locus: CATEGORY_COLORS.get(str(row.category), "black")
        for row in locus_categories.itertuples(index=False)
    }


def plot_qv_comparison_scatter(summary: pd.DataFrame, outpath: Path, max_label_loci: int) -> None:
    data = summary.dropna(subset=["median_qv_plus", "median_qv_minus"]).copy()
    if data.empty:
        return

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    for category in CATEGORY_ORDER:
        sub = data[data["category"].astype(str) == category]
        if sub.empty:
            continue
        ax.scatter(
            sub["median_qv_plus"],
            sub["median_qv_minus"],
            s=22,
            alpha=0.75,
            color=CATEGORY_COLORS[category],
            edgecolors="none",
            label=category,
        )

    upper = float(np.nanmax([data["median_qv_plus"].max(), data["median_qv_minus"].max(), 45]))
    ax.plot([0, upper], [0, upper], linestyle="--", color="black", linewidth=1, alpha=0.6)

    label_data = data.sort_values(["category", "median_delta_qv_eas"], ascending=[True, False]).head(max_label_loci)
    for row in label_data.itertuples(index=False):
        if str(row.category) == "Panel/ancestry-dependent" and max_label_loci < len(data):
            continue
        ax.text(row.median_qv_plus, row.median_qv_minus, row.locus, fontsize=6, alpha=0.8)

    ax.set_xlabel("Median capped QV, APGp1+ LOO")
    ax.set_ylabel("Median capped QV, APGp1- LOO")
    ax.set_xlim(left=-1, right=upper + 1)
    ax.set_ylim(bottom=-1, top=upper + 1)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_delta_qv_fer(summary: pd.DataFrame, outpath: Path) -> None:
    data = summary.copy()
    data = data.sort_values(["category", "median_delta_qv_eas", "median_fer_all"], ascending=[True, False, False])
    if data.empty:
        return

    x = np.arange(len(data))
    colors = [CATEGORY_COLORS.get(str(cat), "black") for cat in data["category"]]

    fig, ax1 = plt.subplots(figsize=(max(10, len(data) * 0.06), 5.5))
    ax1.bar(x, data["median_fer_all"], color=colors, edgecolor="none", alpha=0.85)
    ax1.axhline(1, color="black", linestyle="--", linewidth=1)
    ax1.set_ylabel("Median fold error reduction")
    ax1.set_xlabel("Locus")
    ax1.set_xticks(x)
    ax1.set_xticklabels(data["locus"], rotation=90, fontsize=5)

    if data["median_fer_all"].dropna().max() > 50:
        ax1.set_yscale("log")

    ax2 = ax1.twinx()
    ax2.plot(x, data["median_delta_qv_eas"], marker="o", markersize=2.5, linewidth=0.8, label="EAS")
    ax2.plot(x, data["median_delta_qv_non_eas"], marker="o", markersize=2.5, linewidth=0.8, label="Non-EAS")
    ax2.axhline(1, color="black", linestyle=":", linewidth=1)
    ax2.set_ylabel("Median Delta QV")
    ax2.legend(frameon=False, fontsize=8, loc="upper right")

    legend_handles = [Patch(facecolor=CATEGORY_COLORS[c], label=c) for c in CATEGORY_ORDER]
    ax1.legend(handles=legend_handles, frameon=False, fontsize=8, loc="upper left")

    ax1.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_qv_stack_grid(
    qv_summary: pd.DataFrame,
    locus_order: Sequence[str],
    plus_label: str,
    minus_label: str,
    outpath: Path,
    ncols: int,
) -> None:
    if qv_summary.empty:
        return

    panels = [plus_label, minus_label]
    population_groups = [group for group in ["All", "EAS", "Non-EAS"] if group in set(qv_summary["population_group"])]
    x_labels = [f"{panel}\n{pop}" for pop in population_groups for panel in panels]
    x_keys = [(pop, panel) for pop in population_groups for panel in panels]

    qv_colors = plt.cm.Set1(np.linspace(0, 1, len(QV_BIN_LABELS_HIGH_TO_LOW) - 1)).tolist() + ["lightgrey"]
    n_loci = len(locus_order)
    ncols = max(1, ncols)
    nrows = math.ceil(n_loci / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(2.2 * ncols, 1.9 * nrows), squeeze=False)
    fig.subplots_adjust(left=0.04, right=0.995, top=0.98, bottom=0.04, wspace=0.25, hspace=1.15)

    qv_lookup = qv_summary.set_index(["locus", "population_group", "panel"])

    for idx, locus in enumerate(locus_order):
        ax = axes[idx // ncols, idx % ncols]
        for x_idx, (population_group, panel) in enumerate(x_keys):
            key = (locus, population_group, panel)
            if key not in qv_lookup.index:
                continue
            row = qv_lookup.loc[key]
            bottom = 0.0
            for qv_bin, color in zip(QV_BIN_LABELS_HIGH_TO_LOW, qv_colors):
                value = row.get(f"fraction_qv_{qv_bin}", 0.0)
                if pd.isna(value):
                    value = 0.0
                ax.bar(x_idx, value, bottom=bottom, color=color, edgecolor="black", linewidth=0.2, width=0.85)
                bottom += value

        ax.set_title(locus, fontsize=7, pad=2)
        ax.set_ylim(0, 1)
        ax.set_xticks(range(len(x_keys)))
        ax.set_xticklabels(x_labels, fontsize=4, rotation=90)
        ax.tick_params(axis="y", labelsize=5)
        if idx % ncols == 0:
            ax.set_ylabel("Fraction", fontsize=6)
        else:
            ax.set_yticklabels([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for idx in range(n_loci, nrows * ncols):
        axes[idx // ncols, idx % ncols].axis("off")

    legend_handles = [Patch(facecolor=color, label=qv_bin) for qv_bin, color in zip(QV_BIN_LABELS_HIGH_TO_LOW, qv_colors)]
    fig.legend(handles=legend_handles, title="QV bin", loc="upper center", ncol=8, frameon=False, fontsize=7, title_fontsize=8)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    outdir = mkdir(args.outdir)
    apply_plot_defaults()

    plus_prefixes = parse_prefixes(args.panel_plus_use_call_prefixes)
    minus_prefixes = parse_prefixes(args.panel_minus_use_call_prefixes)

    sample_pop = read_sample_pop(args.sample_pop, sep=args.sample_pop_sep, has_header=args.sample_pop_has_header)
    panel_plus = prepare_panel(args.panel_plus_eval, args.sep, args.plus_label, plus_prefixes, args.qv_cap)
    panel_minus = prepare_panel(args.panel_minus_eval, args.sep, args.minus_label, minus_prefixes, args.qv_cap)

    merged = merge_panels(panel_plus, panel_minus, sample_pop)
    selected_eas = select_eas_samples(
        merged=merged,
        eas_label=args.eas_label,
        random_seed=args.random_seed,
        eas_samples_file=args.eas_samples,
        out_path=outdir / f"{args.prefix}.selected_eas_samples.txt",
    )

    locus_categories = classify_loci(
        panel_plus=panel_plus,
        panel_minus=panel_minus,
        merged=merged,
        selected_eas_samples=selected_eas,
        plus_label=args.plus_label,
        minus_label=args.minus_label,
        no_call_threshold=args.no_call_threshold,
        delta_qv_threshold=args.delta_qv_threshold,
        eas_label=args.eas_label,
    )
    summary = summarize_loci(merged, locus_categories, eas_label=args.eas_label)
    qv_summary = make_qv_summary(
        panel_plus=panel_plus,
        panel_minus=panel_minus,
        sample_pop=sample_pop,
        selected_eas_samples=selected_eas,
        use_all_eas=args.use_all_eas_for_qv_summary,
        plus_label=args.plus_label,
        minus_label=args.minus_label,
        eas_label=args.eas_label,
        non_eas_label=args.non_eas_label,
    )

    merged.to_csv(outdir / f"{args.prefix}.paired_haplotype_metrics.tsv", sep="\t", index=False)
    locus_categories.to_csv(outdir / f"{args.prefix}.locus_categories.tsv", sep="\t", index=False)
    summary.to_csv(outdir / f"{args.prefix}.locus_performance_summary.tsv", sep="\t", index=False)
    qv_summary.to_csv(outdir / f"{args.prefix}.qv_bin_summary.tsv", sep="\t", index=False)
    write_locus_lists(locus_categories, outdir, args.prefix)

    locus_order = summary["locus"].tolist()
    pd.Series(locus_order).to_csv(outdir / f"{args.prefix}.locus_order.txt", index=False, header=False)

    if args.make_plots:
        plot_qv_comparison_scatter(
            summary,
            outdir / f"{args.prefix}.median_qv_comparison.{args.plot_format}",
            max_label_loci=args.max_label_loci,
        )
        plot_delta_qv_fer(
            summary,
            outdir / f"{args.prefix}.delta_qv_fer_summary.{args.plot_format}",
        )

    if args.make_qv_stack_plot:
        plot_qv_stack_grid(
            qv_summary=qv_summary,
            locus_order=locus_order,
            plus_label=args.plus_label,
            minus_label=args.minus_label,
            outpath=outdir / f"{args.prefix}.qv_bin_stack_by_locus.{args.plot_format}",
            ncols=args.stack_ncols,
        )

    counts = locus_categories["category"].astype(str).value_counts().reindex(CATEGORY_ORDER).fillna(0).astype(int)
    sys.stderr.write("Locus category counts:\n")
    for category, count in counts.items():
        sys.stderr.write(f"  {category}: {count}\n")
    sys.stderr.write(f"Outputs written to: {outdir}\n")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except InputError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        raise SystemExit(2)
