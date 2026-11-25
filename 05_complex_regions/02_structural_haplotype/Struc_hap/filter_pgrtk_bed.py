#!/usr/bin/env python3
"""Coverage- and haplotype-aware filtering of a final PGR-TK-like BED.

With only BED and FAI inputs, bundle length in bp is defined as the median
genomic span of all BED occurrences of that bundle ID. The output preserves
the original PGR annotation and any extra BED columns verbatim. An optional
TRF BED contributes immutable, interval-unioned sequence coverage.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence

try:
    from .final_bundle_filter import filter_final_bundle_classes, write_final_filter_audits
    from .merge_pgr_bundles import read_pgr_bed, read_trf_bed
except ImportError:
    from final_bundle_filter import filter_final_bundle_classes, write_final_filter_audits
    from merge_pgr_bundles import read_pgr_bed, read_trf_bed


@dataclass(frozen=True)
class Assignment:
    numeric_bundle_id: str


@dataclass(frozen=True)
class BedMapping:
    target_name: str
    target_start: int
    target_end: int
    bundle_name: str
    strand: str
    raw_query_name: str
    query_length: int
    line_number: int


def read_fai(path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected at least two FAI columns")
            sequence_id = fields[0]
            if sequence_id in lengths:
                raise ValueError(f"{path}:{line_number}: duplicate sequence ID {sequence_id}")
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: sequence length must be an integer") from exc
            if length <= 0:
                raise ValueError(f"{path}:{line_number}: sequence length must be positive")
            lengths[sequence_id] = length
    if not lengths:
        raise ValueError(f"{path}: no FAI records found")
    return lengths


def default_output_path(input_bed: Path) -> Path:
    return input_bed.with_name(f"{input_bed.stem}.filtered.bed")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Filter <=500-bp bundle classes from a final PGR-TK-like BED only when "
            "deletion merges presence/absence haplotypes, while protecting coverage, "
            "unchanged haplotype partitions, and signed structural differences."
        )
    )
    parser.add_argument("pgr_bed", type=Path)
    parser.add_argument("fai", type=Path, help="samtools .fa.fai for the complete population FASTA")
    parser.add_argument("--trf-bed", type=Path,
                        help="optional population TRF BED; passing repeat intervals contribute immutable coverage")
    parser.add_argument(
        "--trf-min-length", type=int, default=0,
        help="minimum TRF interval length (default: 0; no length filter)",
    )
    parser.add_argument("--trf-min-identity", type=float, default=80.0)
    parser.add_argument("-o", "--output-bed", type=Path, help="default: INPUT_STEM.filtered.bed")
    parser.add_argument("--max-bundle-bp", type=int, default=500,
                        help="inclusive maximum median occurrence span eligible for deletion (default: 500)")
    parser.add_argument("--min-sequence-coverage", type=float, default=0.95)
    parser.add_argument("--max-low-coverage-fraction", type=float, default=0.05)
    parser.add_argument("--force", action="store_true", help="overwrite existing output and audit files")
    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if not args.pgr_bed.is_file():
        parser.error(f"PGR BED does not exist or is not a file: {args.pgr_bed}")
    if not args.fai.is_file():
        parser.error(f"FAI does not exist or is not a file: {args.fai}")
    if args.trf_bed is not None and not args.trf_bed.is_file():
        parser.error(f"TRF BED does not exist or is not a file: {args.trf_bed}")
    if args.max_bundle_bp < 0:
        parser.error("--max-bundle-bp must be non-negative")
    if args.trf_min_length < 0:
        parser.error("--trf-min-length must be non-negative")
    if not math.isfinite(args.trf_min_identity) or not 0 <= args.trf_min_identity <= 100:
        parser.error("--trf-min-identity must be finite and in [0, 100]")
    for name in ("min_sequence_coverage", "max_low_coverage_fraction"):
        if not 0 <= getattr(args, name) <= 1:
            parser.error(f"--{name.replace('_', '-')} must be in [0, 1]")
    output_bed = args.output_bed or default_output_path(args.pgr_bed)
    protected_inputs = {args.pgr_bed.resolve(), args.fai.resolve()}
    if args.trf_bed is not None:
        protected_inputs.add(args.trf_bed.resolve())
    if output_bed.resolve() in protected_inputs:
        parser.error("output BED must differ from every input file")


def output_paths(output_bed: Path):
    prefix = output_bed.name[:-4] if output_bed.name.endswith(".bed") else output_bed.name
    audit_prefix = f"{prefix}.filter"
    audits = {
        "decisions": output_bed.parent / f"{audit_prefix}_bundle_filter.tsv",
        "merge_events": output_bed.parent / f"{audit_prefix}_bundle_haplotype_merges.tsv",
        "haplotypes_before": output_bed.parent / f"{audit_prefix}_haplotypes.before.tsv",
        "haplotypes_after": output_bed.parent / f"{audit_prefix}_haplotypes.after.tsv",
        "tandem_bundle_classes": output_bed.parent / f"{audit_prefix}_tandem_bundle_classes.tsv",
        "sequence_coverage": output_bed.parent / f"{audit_prefix}_sequence_coverage.tsv",
    }
    summary = output_bed.parent / f"{audit_prefix}_summary.json"
    return audit_prefix, audits, summary


def write_filtered_bed(input_bed: Path, output_bed: Path, retained_lines) -> None:
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_bed.name}.", suffix=".tmp", dir=output_bed.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as target, input_bed.open(
            "r", encoding="utf-8"
        ) as source:
            for line_number, line in enumerate(source, 1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#") or line_number in retained_lines:
                    target.write(line)
        os.replace(temporary_name, output_bed)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    output_bed = args.output_bed or default_output_path(args.pgr_bed)
    audit_prefix, expected_audits, summary_path = output_paths(output_bed)
    generated = [output_bed, *expected_audits.values(), summary_path]
    input_paths = {args.pgr_bed.resolve(), args.fai.resolve()}
    if args.trf_bed is not None:
        input_paths.add(args.trf_bed.resolve())
    collisions = [path for path in generated if path.resolve() in input_paths]
    if collisions:
        parser.error(f"generated output would overwrite an input file: {collisions[0]}")
    unsafe = [path for path in generated
              if path.is_symlink() or (path.exists() and not path.is_file())]
    if unsafe:
        parser.error(f"refusing unsafe output path: {unsafe[0]}")
    existing = [path for path in generated if path.exists() or path.is_symlink()]
    if existing and not args.force:
        parser.error(f"output exists; use --force to overwrite: {existing[0]}")

    sequence_lengths = read_fai(args.fai)
    trf_intervals = {}
    trf_read_summary = {
        "trf_input_row_count": 0,
        "trf_selected_row_count": 0,
        "trf_union_interval_count": 0,
    }
    if args.trf_bed is not None:
        selected_intervals, selected_records, trf_read_summary = read_trf_bed(
            args.trf_bed, args.trf_min_length, args.trf_min_identity
        )
        for path_id, records in selected_records.items():
            if path_id not in sequence_lengths:
                raise ValueError(f"TRF sequence {path_id} is missing from FAI")
            if any(not math.isfinite(record.identity) or not 0 <= record.identity <= 100
                   for record in records):
                raise ValueError(f"TRF identity must be finite and in [0, 100] for {path_id}")
        for path_id, intervals in selected_intervals.items():
            if any(interval.end > sequence_lengths[path_id] for interval in intervals):
                raise ValueError(f"TRF coordinates exceed FAI length for {path_id}")
            trf_intervals[path_id] = [(interval.start, interval.end) for interval in intervals]
    by_path = read_pgr_bed(args.pgr_bed)
    mappings = []
    spans_by_bundle = {}
    for path_id, records in by_path.items():
        if path_id not in sequence_lengths:
            raise ValueError(f"BED sequence {path_id} is missing from FAI")
        for record in records:
            if record.end > sequence_lengths[path_id]:
                raise ValueError(f"BED coordinates exceed FAI length for {path_id}")
            mappings.append(BedMapping(
                path_id, record.start, record.end, record.bundle_id, record.direction,
                record.bundle_id, record.length, record.line_number,
            ))
            spans_by_bundle.setdefault(record.bundle_id, []).append(record.length)
    class_lengths = {bundle_id: statistics.median(spans)
                     for bundle_id, spans in spans_by_bundle.items()}
    assignments = {bundle_id: Assignment(bundle_id) for bundle_id in class_lengths}
    retained, result = filter_final_bundle_classes(
        mappings, assignments, sequence_lengths,
        args.min_sequence_coverage, args.max_low_coverage_fraction,
        args.max_bundle_bp, class_lengths=class_lengths,
        fixed_coverage_intervals=trf_intervals,
    )
    retained_lines = {mapping.line_number for mapping in retained}
    write_filtered_bed(args.pgr_bed, output_bed, retained_lines)
    actual_audits = write_final_filter_audits(output_bed.parent, result, audit_prefix)
    if {Path(path).resolve() for path in actual_audits.values()} != {
        path.resolve() for path in expected_audits.values()
    }:
        raise RuntimeError("internal audit output path mismatch")
    summary = {
        "inputs": {
            "pgr_bed": str(args.pgr_bed.resolve()),
            "fai": str(args.fai.resolve()),
            "trf_bed": str(args.trf_bed.resolve()) if args.trf_bed is not None else None,
        },
        "output_bed": str(output_bed.resolve()),
        "parameters": {
            "max_bundle_bp": args.max_bundle_bp,
            "min_sequence_coverage": args.min_sequence_coverage,
            "max_low_coverage_fraction": args.max_low_coverage_fraction,
            "bundle_length_definition": "median BED occurrence span per bundle ID",
            "coverage_denominator": "full FAI sequence length for every FAI record",
            "coverage_numerator": "union of retained bundle BED and filtered TRF intervals",
            "trf_min_length": args.trf_min_length,
            "trf_min_identity": args.trf_min_identity,
            "max_bundle_bp_is_inclusive": True,
            "require_haplotype_partition_change": True,
            "copy_number_protection": (
                "only bundle classes occurring consecutively at least twice in an original "
                "signed haplotype order"
            ),
            "deleted_bundle_position_protection": False,
        },
        "counts": {**result["counts"], **trf_read_summary},
        "outputs": actual_audits,
    }
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    json.dump(summary, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
