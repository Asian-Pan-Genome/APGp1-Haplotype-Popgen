#!/usr/bin/env python3
"""Exclude every occurrence of a repeated population FASTA path ID.

The first whitespace-delimited token of each FASTA header is the path ID,
matching the readers used by the ranked-library and remapping stages.
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Callable


def fasta_id(header: str) -> str:
    fields = header[1:].split()
    if not fields:
        raise ValueError("FASTA contains an empty record ID")
    return fields[0]


def count_ids(path: Path) -> Counter:
    counts: Counter = Counter()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                counts[fasta_id(line)] += 1
    if not counts:
        raise ValueError(f"{path}: no FASTA records")
    return counts


def atomic_filter(source: Path, target: Path, keep: Callable[[str], bool]) -> tuple[int, int]:
    """Copy selected lines to a file replaced only after a complete write."""
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    kept = dropped = 0
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=target.parent,
            prefix=f".{target.name}.", delete=False,
        ) as output:
            temporary = Path(output.name)
            with source.open("r", encoding="utf-8", newline="") as input_handle:
                for line in input_handle:
                    if keep(line):
                        output.write(line)
                        kept += 1
                    else:
                        dropped += 1
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return kept, dropped


def filter_fasta(source: Path, target: Path, duplicates: set[str]) -> tuple[int, int]:
    keep_record = False
    records_kept = records_dropped = 0

    def keep(line: str) -> bool:
        nonlocal keep_record, records_kept, records_dropped
        if line.startswith(">"):
            keep_record = fasta_id(line) not in duplicates
            if keep_record:
                records_kept += 1
            else:
                records_dropped += 1
        return keep_record

    atomic_filter(source, target, keep)
    return records_kept, records_dropped


def filter_table(source: Path, target: Path, duplicates: set[str]) -> int:
    def keep(line: str) -> bool:
        stripped = line.strip()
        return not stripped or stripped.startswith("#") or stripped.split()[0] not in duplicates

    _kept, dropped = atomic_filter(source, target, keep)
    return dropped


def filter_duplicate_inputs(fasta: Path, bed: Path, trf_bed: Path | None,
                            sample_map: Path | None, output_dir: Path) -> dict:
    """Write effective inputs excluding every record of each repeated FASTA ID."""
    for path in (fasta, bed, trf_bed, sample_map):
        if path is not None and not path.is_file():
            raise ValueError(f"input file not found: {path}")
    counts = count_ids(fasta)
    duplicates = {path_id for path_id, count in counts.items() if count > 1}
    output_dir.mkdir(parents=True, exist_ok=True)
    report = output_dir / "duplicate_fasta_ids.tsv"
    atomic_filter_report = "path_id\tfasta_record_count\n" + "".join(
        f"{path_id}\t{counts[path_id]}\n" for path_id in sorted(duplicates)
    )
    # The report is small, but use the same atomic-write convention as the inputs.
    report_tmp = output_dir / ".duplicate_fasta_ids.tsv.tmp"
    report_tmp.write_text(atomic_filter_report, encoding="utf-8")
    os.replace(report_tmp, report)

    result = {
        "duplicate_id_count": len(duplicates),
        "excluded_fasta_record_count": sum(counts[path_id] for path_id in duplicates),
        "kept_fasta_record_count": sum(count for path_id, count in counts.items() if path_id not in duplicates),
        "input_fasta": str(fasta.resolve()),
        "input_bed": str(bed.resolve()),
        "report": str(report.resolve()),
        "filtered_inputs": {},
        "excluded_table_rows": {},
    }
    if duplicates:
        if result["kept_fasta_record_count"] == 0:
            raise ValueError("all FASTA IDs are duplicated; no population records remain")
        targets = {
            "fasta": (fasta, output_dir / "population.unique_ids.fa"),
            "bed": (bed, output_dir / "pgrtk.unique_ids.bed"),
            "trf_bed": (trf_bed, output_dir / "trf.unique_ids.bed"),
            "sample_map": (sample_map, output_dir / "sample_map.unique_ids.tsv"),
        }
        filter_fasta(*targets["fasta"], duplicates)
        for name in ("bed", "trf_bed", "sample_map"):
            source, target = targets[name]
            if source is not None:
                result["excluded_table_rows"][name] = filter_table(source, target, duplicates)
        if not any(line.strip() and not line.startswith("#") for line in targets["bed"][1].open()):
            raise ValueError("no BED entries remain after excluding duplicated FASTA IDs")
        result["filtered_inputs"] = {
            name: str(target.resolve()) for name, (source, target) in targets.items()
            if source is not None
        }
    summary = output_dir / "duplicate_fasta_id_filter.json"
    summary_tmp = output_dir / ".duplicate_fasta_id_filter.json.tmp"
    summary_tmp.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(summary_tmp, summary)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--trf-bed", type=Path)
    parser.add_argument("--sample-map", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = filter_duplicate_inputs(args.fasta, args.bed, args.trf_bed,
                                         args.sample_map, args.output_dir)
        print(result["duplicate_id_count"])
        return 0
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
