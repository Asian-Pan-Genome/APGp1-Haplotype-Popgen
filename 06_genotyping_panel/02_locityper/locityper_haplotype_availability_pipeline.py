#!/usr/bin/env python3
"""
Compute locus-level haplotype availability from Locityper haplotype alignments.

This pipeline keeps one full pairwise PAF per locus:

    <db-dir>/loci/<locus>/haplotypes.paf.gz

QV is then computed with the custom gt_dist.new.py logic:

    python gt_dist.new.py -G <genotypes> --loo -n 1 -i haplotypes.paf.gz
    python gt_dist.new.py -G <genotypes> --loo -n 1 -i haplotypes.paf.gz -Q samples.HPRCy1_HGSVC3.txt

The first command evaluates availability against the full panel. The second command
uses the same full PAF file but restricts available query haplotypes with -Q.
"""

from __future__ import annotations

import argparse
import math
import multiprocessing as mp
import os
from pathlib import Path
import shlex
import subprocess
import sys
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PipelineConfig:
    db_dir: Path
    genotypes: Path
    locityper: str
    gt_dist: Path
    python_exe: str
    threads_per_process: int
    align_divergence: str
    align_extra_args: tuple[str, ...]
    skip_align: bool
    skip_gt_dist: bool
    force_align: bool
    force_gt_dist: bool
    subset_name: str
    full_label: str
    subset_label: str
    apgp1_prefixes: tuple[str, ...]
    qv_bin_scheme: str
    include_no_call: bool
    input_fasta_template: str
    paf_template: str
    discarded_template: str
    subset_samples_template: str
    full_out_template: str
    subset_out_template: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run Locityper pairwise haplotype alignment and compute haplotype "
            "availability using haplotypes.paf.gz + gt_dist.new.py -Q."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    required = parser.add_argument_group("required inputs")
    required.add_argument("--loci-bed", required=True, type=Path,
                          help="BED-like file with at least four columns: chrom, start, end, locus_name.")
    required.add_argument("--samples", required=True, type=Path,
                          help="Sample list. One-column files are interpreted as sample IDs. Two-column files are interpreted as sample and source.")
    required.add_argument("--genotypes", required=True, type=Path,
                          help="Genotype/haplotype list passed to gt_dist.new.py with -G.")
    required.add_argument("--output", required=True, type=Path,
                          help="Output summary table with one row per locus, panel and QV bin.")

    paths = parser.add_argument_group("paths and commands")
    paths.add_argument("--db-dir", type=Path, default=Path("db"),
                       help="Locityper database directory containing loci/<locus>/haplotypes.fa.gz.")
    paths.add_argument("--locityper", default="locityper",
                       help="Locityper executable.")
    paths.add_argument("--gt-dist", type=Path, default=Path("gt_dist.new.py"),
                       help="Path to gt_dist.new.py. A basename is allowed if it is in the working directory.")
    paths.add_argument("--python-exe", default=sys.executable,
                       help="Python executable used to run gt_dist.new.py.")
    paths.add_argument("--raw-output", type=Path, default=None,
                       help="Optional raw per-haplotype QV table after binning.")
    paths.add_argument("--failed-output", type=Path, default=None,
                       help="Optional table listing loci that failed.")

    templates = parser.add_argument_group("per-locus path templates")
    templates.add_argument("--input-fasta-template", default="{db_dir}/loci/{locus}/haplotypes.fa.gz",
                           help="Input FASTA template for locityper align.")
    templates.add_argument("--paf-template", default="{db_dir}/loci/{locus}/haplotypes.paf.gz",
                           help="Full pairwise PAF template. This is the PAF passed to gt_dist.new.py.")
    templates.add_argument("--discarded-template", default="{db_dir}/loci/{locus}/discarded_haplotypes.txt",
                           help="Optional discarded haplotypes file template passed to gt_dist.new.py with -d when present.")
    templates.add_argument("--subset-samples-template", default="{db_dir}/loci/{locus}/samples.{subset_name}.txt",
                           help="Per-locus sample list for the -Q subset.")
    templates.add_argument("--full-out-template", default="{db_dir}/loci/{locus}/out.all.hprc_hgsvc.csv",
                           help="gt_dist.new.py output for the full panel.")
    templates.add_argument("--subset-out-template", default="{db_dir}/loci/{locus}/out.{subset_name}.hprc_hgsvc.csv",
                           help="gt_dist.new.py output for the -Q subset panel.")

    execution = parser.add_argument_group("execution")
    execution.add_argument("--num-processes", type=int, default=1,
                           help="Number of loci processed in parallel.")
    execution.add_argument("--threads-per-process", type=int, default=1,
                           help="Threads used by each locityper align process.")
    execution.add_argument("--align-divergence", default="1",
                           help="Value passed to locityper align with -D.")
    execution.add_argument("--align-extra-args", default="",
                           help="Additional arguments appended to locityper align, quoted as one string.")
    execution.add_argument("--skip-align", action="store_true",
                           help="Skip locityper align and reuse existing haplotypes.paf.gz files.")
    execution.add_argument("--skip-gt-dist", action="store_true",
                           help="Skip gt_dist.new.py and reuse existing out.*.csv files.")
    execution.add_argument("--force-align", action="store_true",
                           help="Re-run locityper align even if haplotypes.paf.gz already exists.")
    execution.add_argument("--force-gt-dist", action="store_true",
                           help="Re-run gt_dist.new.py even if output CSV files already exist.")

    sample_opts = parser.add_argument_group("sample grouping")
    sample_opts.add_argument("--subset-name", default="HPRCy1_HGSVC3",
                             help="Name used in subset output filenames and sample list filenames.")
    sample_opts.add_argument("--apgp1-prefixes", default="C,K",
                             help="Comma-separated sample prefixes treated as APGp1 when the sample file has no source column.")
    sample_opts.add_argument("--full-label", default="APGp1 + HPRCy1 + HGSVC3",
                             help="Panel label for the full-panel gt_dist.new.py run.")
    sample_opts.add_argument("--subset-label", default="HPRCy1 + HGSVC3",
                             help="Panel label for the -Q subset gt_dist.new.py run.")

    qv_opts = parser.add_argument_group("QV summarization")
    qv_opts.add_argument("--qv-bin-scheme", choices=["legacy", "locityper"], default="legacy",
                         help="legacy: <20,20-25,...,>45. locityper: <17,17-23,23-33,33-43,>=43.")
    qv_opts.add_argument("--include-no-call", action="store_true",
                         help="Include rows with missing/non-numeric QV as a No call bin in summary denominators.")

    args = parser.parse_args()
    if args.num_processes < 1:
        parser.error("--num-processes must be >= 1")
    if args.threads_per_process < 1:
        parser.error("--threads-per-process must be >= 1")
    return args


def format_path(template: str, *, cfg: PipelineConfig, locus: str) -> Path:
    return Path(template.format(
        db_dir=str(cfg.db_dir),
        locus=locus,
        subset_name=cfg.subset_name,
    ))


def run_command(cmd: list[str], *, dry_label: str) -> None:
    print(f"[CMD] {dry_label}: {' '.join(shlex.quote(x) for x in cmd)}", flush=True)
    subprocess.run(cmd, check=True)


def read_loci(loci_bed: Path) -> list[str]:
    df = pd.read_csv(loci_bed, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        raise ValueError("--loci-bed must contain at least four columns: chrom, start, end, locus_name")
    loci = df.iloc[:, 3].astype(str).tolist()
    if not loci:
        raise ValueError("No loci found in --loci-bed")
    return loci


def read_samples(samples_path: Path, apgp1_prefixes: tuple[str, ...]) -> pd.DataFrame:
    df = pd.read_csv(samples_path, sep="\t", header=None, comment="#", dtype=str)
    if df.empty:
        raise ValueError("No samples found in --samples")

    if df.shape[1] == 1:
        out = pd.DataFrame({"sample": df.iloc[:, 0].astype(str)})
        out["source"] = np.where(
            out["sample"].str.startswith(apgp1_prefixes),
            "APGp1",
            "HPRCy1 + HGSVC3",
        )
        return out

    out = pd.DataFrame({"sample": df.iloc[:, 0].astype(str), "source": df.iloc[:, 1].astype(str)})
    return out


def write_subset_samples(samples: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    subset = samples.loc[samples["source"] != "APGp1", "sample"].drop_duplicates()
    if subset.empty:
        raise ValueError("No non-APGp1 samples were found for the -Q subset")
    subset.to_csv(path, sep="\t", index=False, header=False)


def qv_bins(scheme: str) -> tuple[list[float], list[str]]:
    if scheme == "legacy":
        return [-math.inf, 20, 25, 30, 35, 40, 45, math.inf], ["<20", "20-25", "25-30", "30-35", "35-40", "40-45", ">45"]
    if scheme == "locityper":
        return [-math.inf, 17, 23, 33, 43, math.inf], ["<17", "17-23", "23-33", "33-43", ">=43"]
    raise ValueError(f"Unsupported QV bin scheme: {scheme}")


def load_and_bin_qv(path: Path, *, locus: str, panel: str, scheme: str, include_no_call: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(path, sep="\t", comment="#", header=0)
    df.insert(0, "panel", panel)
    df.insert(0, "locus", locus)

    if "qv" not in df.columns:
        raise ValueError(f"Missing qv column in {path}")

    qv = pd.to_numeric(df["qv"], errors="coerce")
    qv = qv.replace(np.inf, 100.0)
    qv = qv.replace(-np.inf, np.nan)
    df["qv_numeric"] = qv

    bins, labels = qv_bins(scheme)
    df["qv_bin"] = pd.cut(df["qv_numeric"], bins=bins, labels=labels, right=False)

    if include_no_call:
        df["qv_bin"] = df["qv_bin"].astype("object").where(df["qv_bin"].notna(), "No call")
        summary_labels = labels + ["No call"]
        denom = len(df)
    else:
        df = df.loc[df["qv_bin"].notna()].copy()
        summary_labels = labels
        denom = len(df)

    rows: list[dict[str, object]] = []
    for label in summary_labels:
        n = int((df["qv_bin"] == label).sum()) if denom else 0
        rows.append({
            "locus": locus,
            "panel": panel,
            "qv_bin": label,
            "n_haps": n,
            "total_haps": int(denom),
            "frac_haps": (n / denom) if denom else np.nan,
        })
    return pd.DataFrame(rows), df


def run_align(locus: str, cfg: PipelineConfig) -> None:
    input_fasta = format_path(cfg.input_fasta_template, cfg=cfg, locus=locus)
    paf = format_path(cfg.paf_template, cfg=cfg, locus=locus)

    if not input_fasta.exists():
        raise FileNotFoundError(f"Missing input FASTA: {input_fasta}")
    if paf.exists() and not cfg.force_align:
        print(f"[SKIP] {locus}: existing PAF {paf}", flush=True)
        return

    paf.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        cfg.locityper,
        "align",
        "-i", str(input_fasta),
        "-o", str(paf),
        "-D", cfg.align_divergence,
        f"-@{cfg.threads_per_process}",
    ]
    cmd.extend(cfg.align_extra_args)
    run_command(cmd, dry_label=f"locityper align {locus}")


def run_gt_dist(locus: str, cfg: PipelineConfig, samples: pd.DataFrame) -> None:
    paf = format_path(cfg.paf_template, cfg=cfg, locus=locus)
    discarded = format_path(cfg.discarded_template, cfg=cfg, locus=locus)
    subset_samples = format_path(cfg.subset_samples_template, cfg=cfg, locus=locus)
    full_out = format_path(cfg.full_out_template, cfg=cfg, locus=locus)
    subset_out = format_path(cfg.subset_out_template, cfg=cfg, locus=locus)

    if not paf.exists():
        raise FileNotFoundError(f"Missing PAF for gt_dist.new.py: {paf}")
    write_subset_samples(samples, subset_samples)

    base = [
        cfg.python_exe,
        str(cfg.gt_dist),
        "-G", str(cfg.genotypes),
        "--loo",
        "-n", "1",
        "-i", str(paf),
    ]
    if discarded.exists():
        base.extend(["-d", str(discarded)])

    full_out.parent.mkdir(parents=True, exist_ok=True)
    subset_out.parent.mkdir(parents=True, exist_ok=True)

    if cfg.force_gt_dist or not full_out.exists():
        run_command(base + ["-o", str(full_out)], dry_label=f"gt_dist full panel {locus}")
    else:
        print(f"[SKIP] {locus}: existing full-panel output {full_out}", flush=True)

    if cfg.force_gt_dist or not subset_out.exists():
        run_command(base + ["-Q", str(subset_samples), "-o", str(subset_out)], dry_label=f"gt_dist -Q subset {locus}")
    else:
        print(f"[SKIP] {locus}: existing subset output {subset_out}", flush=True)


def summarize_outputs(locus: str, cfg: PipelineConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    full_out = format_path(cfg.full_out_template, cfg=cfg, locus=locus)
    subset_out = format_path(cfg.subset_out_template, cfg=cfg, locus=locus)

    full_summary, full_raw = load_and_bin_qv(
        full_out,
        locus=locus,
        panel=cfg.full_label,
        scheme=cfg.qv_bin_scheme,
        include_no_call=cfg.include_no_call,
    )
    subset_summary, subset_raw = load_and_bin_qv(
        subset_out,
        locus=locus,
        panel=cfg.subset_label,
        scheme=cfg.qv_bin_scheme,
        include_no_call=cfg.include_no_call,
    )
    return pd.concat([full_summary, subset_summary], ignore_index=True), pd.concat([full_raw, subset_raw], ignore_index=True)


def process_locus(locus: str, cfg: PipelineConfig, samples: pd.DataFrame) -> dict[str, object]:
    try:
        print(f"[START] {locus}", flush=True)
        if not cfg.skip_align:
            run_align(locus, cfg)
        if not cfg.skip_gt_dist:
            run_gt_dist(locus, cfg, samples)
        summary, raw = summarize_outputs(locus, cfg)
        print(f"[DONE] {locus}", flush=True)
        return {"locus": locus, "ok": True, "summary": summary, "raw": raw, "error": ""}
    except Exception as exc:
        print(f"[ERROR] {locus}: {exc}", file=sys.stderr, flush=True)
        return {"locus": locus, "ok": False, "summary": None, "raw": None, "error": str(exc)}


def main() -> None:
    args = parse_args()
    cfg = PipelineConfig(
        db_dir=args.db_dir,
        genotypes=args.genotypes,
        locityper=args.locityper,
        gt_dist=args.gt_dist,
        python_exe=args.python_exe,
        threads_per_process=args.threads_per_process,
        align_divergence=args.align_divergence,
        align_extra_args=tuple(shlex.split(args.align_extra_args)),
        skip_align=args.skip_align,
        skip_gt_dist=args.skip_gt_dist,
        force_align=args.force_align,
        force_gt_dist=args.force_gt_dist,
        subset_name=args.subset_name,
        full_label=args.full_label,
        subset_label=args.subset_label,
        apgp1_prefixes=tuple(x for x in args.apgp1_prefixes.split(",") if x),
        qv_bin_scheme=args.qv_bin_scheme,
        include_no_call=args.include_no_call,
        input_fasta_template=args.input_fasta_template,
        paf_template=args.paf_template,
        discarded_template=args.discarded_template,
        subset_samples_template=args.subset_samples_template,
        full_out_template=args.full_out_template,
        subset_out_template=args.subset_out_template,
    )

    loci = read_loci(args.loci_bed)
    samples = read_samples(args.samples, cfg.apgp1_prefixes)

    print("=" * 80, flush=True)
    print("Locityper haplotype availability pipeline", flush=True)
    print(f"Loci:                 {len(loci)}", flush=True)
    print(f"Samples:              {samples.shape[0]}", flush=True)
    print(f"Database directory:   {cfg.db_dir}", flush=True)
    print(f"Genotypes file:       {cfg.genotypes}", flush=True)
    print(f"PAF template:         {cfg.paf_template}", flush=True)
    print(f"gt_dist script:       {cfg.gt_dist}", flush=True)
    print(f"Processes / threads:  {args.num_processes} / {cfg.threads_per_process}", flush=True)
    print(f"QV bin scheme:        {cfg.qv_bin_scheme}", flush=True)
    print("=" * 80, flush=True)

    worker_args = [(locus, cfg, samples) for locus in loci]
    if args.num_processes == 1:
        results = [process_locus(*x) for x in worker_args]
    else:
        with mp.Pool(processes=args.num_processes) as pool:
            results = pool.starmap(process_locus, worker_args)

    ok_results = [x for x in results if x["ok"]]
    failed_results = [x for x in results if not x["ok"]]

    if not ok_results:
        raise RuntimeError("No locus finished successfully. Check error messages and failed-output table.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary_df = pd.concat([x["summary"] for x in ok_results], ignore_index=True)
    summary_df.to_csv(args.output, sep="\t", index=False)
    print(f"[WRITE] summary: {args.output}", flush=True)

    if args.raw_output is not None:
        args.raw_output.parent.mkdir(parents=True, exist_ok=True)
        raw_df = pd.concat([x["raw"] for x in ok_results], ignore_index=True)
        raw_df.to_csv(args.raw_output, sep="\t", index=False)
        print(f"[WRITE] raw: {args.raw_output}", flush=True)

    if args.failed_output is not None:
        args.failed_output.parent.mkdir(parents=True, exist_ok=True)
        failed_df = pd.DataFrame(
            [{"locus": x["locus"], "error": x["error"]} for x in failed_results],
            columns=["locus", "error"],
        )
        failed_df.to_csv(args.failed_output, sep="\t", index=False)
        print(f"[WRITE] failed loci: {args.failed_output}", flush=True)

    if failed_results:
        print(f"[WARN] {len(failed_results)} loci failed; {len(ok_results)} loci succeeded.", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
