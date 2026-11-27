#!/usr/bin/env python3
"""Remove final bundle classes recurrently covered by TRF by carrier sequence.

For one mapped bundle occurrence, only the largest coverage supplied by one
selected merged TRF interval is considered.  Disconnected repeat intervals
are never summed to meet the coverage threshold.  A population sequence
qualifies when at least one occurrence on it reaches that threshold.  Multiple
qualifying copies on the same sequence never increase either the numerator or
denominator.  The denominator is the number of distinct population sequences
carrying the bundle class.
"""

import argparse
import bisect
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

try:
    from .merge_pgr_bundles import read_trf_bed
except ImportError:
    from merge_pgr_bundles import read_trf_bed


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path):
    records = []
    name = None
    description = None
    chunks = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, description, "".join(chunks)))
                description = line[1:]
                name = description.split(None, 1)[0] if description else ""
                if not name or name in seen:
                    raise ValueError(f"{path}:{line_number}: invalid or duplicate FASTA name {name!r}")
                seen.add(name)
                chunks = []
            elif line.strip():
                if name is None:
                    raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
                chunks.append(line.strip())
    if name is not None:
        records.append((name, description, "".join(chunks)))
    if not records:
        raise ValueError(f"{path}: no FASTA records found")
    return records


def write_fasta(path: Path, records) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for _, description, sequence in records:
            handle.write(f">{description}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")


def normalize_numeric_id(value: str) -> str:
    value = value.strip()
    return str(int(value)) if value.isdigit() else value


def read_class_map(path: Path):
    required = {"bundle_class", "numeric_bundle_id", "raw_consensus_names"}
    by_class = {}
    by_numeric = defaultdict(set)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: missing class-map columns {sorted(required)}")
        for line_number, row in enumerate(reader, 2):
            bundle_class = row["bundle_class"].strip()
            numeric_id = normalize_numeric_id(row["numeric_bundle_id"])
            if not bundle_class or not numeric_id:
                raise ValueError(f"{path}:{line_number}: empty bundle class or numeric ID")
            if bundle_class in by_class:
                raise ValueError(f"{path}:{line_number}: duplicate bundle class {bundle_class}")
            names = {value for value in row["raw_consensus_names"].split(",") if value}
            by_class[bundle_class] = (numeric_id, names)
            by_numeric[numeric_id].update(names)
    return by_class, by_numeric


def parse_footprints(row, path: Path, line_number: int):
    value = (row.get("genome_footprint_intervals") or "").strip()
    pieces = []
    if value:
        for token in value.split(","):
            token = token.strip()
            if not token:
                continue
            try:
                start_text, end_text = token.split("-", 1)
                start, end = int(start_text), int(end_text)
            except (TypeError, ValueError):
                raise ValueError(f"{path}:{line_number}: invalid footprint interval {token!r}")
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid footprint interval {token!r}")
            pieces.append((start, end))
    if not pieces:
        try:
            start, end = int(row["genome_start"]), int(row["genome_end"])
        except (KeyError, TypeError, ValueError):
            raise ValueError(f"{path}:{line_number}: invalid genome start/end")
        if start < 0 or end <= start:
            raise ValueError(f"{path}:{line_number}: invalid genome start/end {start}-{end}")
        pieces.append((start, end))
    return pieces


def merge_pieces(pieces):
    merged = []
    for start, end in sorted(pieces):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def single_merged_interval_coverage(pieces, intervals, starts):
    """Return footprint bp, union overlap, and best single-interval coverage."""
    pieces = merge_pieces(pieces)
    footprint_bp = sum(end - start for start, end in pieces)
    if not pieces or not intervals or footprint_bp <= 0:
        return footprint_bp, 0, 0, 0.0
    query_start, query_end = pieces[0][0], pieces[-1][1]
    index = max(0, bisect.bisect_left(starts, query_start) - 1)
    total_overlap_bp = 0
    max_single_overlap_bp = 0
    while index < len(intervals) and intervals[index].start < query_end:
        interval = intervals[index]
        covered = sum(
            max(0, min(end, interval.end) - max(start, interval.start))
            for start, end in pieces
        )
        total_overlap_bp += covered
        max_single_overlap_bp = max(max_single_overlap_bp, covered)
        index += 1
    return (
        footprint_bp,
        total_overlap_bp,
        max_single_overlap_bp,
        max_single_overlap_bp / footprint_bp,
    )


def filter_final_library(
    bundle_fasta: Path,
    metadata_bed: Path,
    coverage_mappings: Path,
    class_map: Path,
    trf_bed: Path,
    output_dir: Path,
    min_length: int = 0,
    min_identity: float = 80.0,
    full_coverage: float = 0.95,
    removal_frequency: float = 0.50,
):
    bundle_fasta = Path(bundle_fasta)
    metadata_bed = Path(metadata_bed)
    coverage_mappings = Path(coverage_mappings)
    class_map = Path(class_map)
    trf_bed = Path(trf_bed)
    output_dir = Path(output_dir)
    if not isinstance(min_length, int) or min_length < 0:
        raise ValueError("TRF minimum length must be a nonnegative integer")
    if not math.isfinite(min_identity) or not 0 <= min_identity <= 100:
        raise ValueError("TRF minimum identity must be finite and in [0, 100]")
    if not math.isfinite(full_coverage) or not 0 < full_coverage <= 1:
        raise ValueError("TRF full coverage must be finite and in (0, 1]")
    if not math.isfinite(removal_frequency) or not 0 < removal_frequency <= 1:
        raise ValueError("TRF removal frequency must be finite and in (0, 1]")
    for path in (bundle_fasta, metadata_bed, coverage_mappings, class_map, trf_bed):
        if not path.is_file():
            raise ValueError(f"input file does not exist: {path}")

    records = read_fasta(bundle_fasta)
    by_class, by_numeric = read_class_map(class_map)
    trf_intervals, _, trf_summary = read_trf_bed(trf_bed, min_length, min_identity)
    trf_starts = {
        path_id: [interval.start for interval in intervals]
        for path_id, intervals in trf_intervals.items()
    }

    required = {"bundle_name", "raw_consensus_name", "population_path",
                "genome_start", "genome_end"}
    carriers = defaultdict(set)
    any_overlapping = defaultdict(set)
    qualifying = defaultdict(set)
    mapping_counts = defaultdict(lambda: defaultdict(int))
    any_overlap_counts = defaultdict(lambda: defaultdict(int))
    qualifying_mapping_counts = defaultdict(lambda: defaultdict(int))
    overlap_bases = defaultdict(lambda: defaultdict(int))
    max_single_overlap_bases = defaultdict(lambda: defaultdict(int))
    max_single_coverages = defaultdict(lambda: defaultdict(float))
    raw_names = defaultdict(set)
    with coverage_mappings.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{coverage_mappings}: missing mapping columns {sorted(required)}")
        for line_number, row in enumerate(reader, 2):
            bundle_class = row["bundle_name"].strip()
            path_id = row["population_path"].strip()
            raw_name = row["raw_consensus_name"].strip()
            if not bundle_class or not path_id or not raw_name:
                raise ValueError(f"{coverage_mappings}:{line_number}: empty class, path, or consensus")
            carriers[bundle_class].add(path_id)
            raw_names[bundle_class].add(raw_name)
            mapping_counts[bundle_class][path_id] += 1
            intervals = trf_intervals.get(path_id, ())
            starts = trf_starts.get(path_id, ())
            _, total_bp, single_bp, coverage = single_merged_interval_coverage(
                parse_footprints(row, coverage_mappings, line_number), intervals, starts,
            )
            if total_bp > 0:
                any_overlapping[bundle_class].add(path_id)
                any_overlap_counts[bundle_class][path_id] += 1
                overlap_bases[bundle_class][path_id] += total_bp
            max_single_overlap_bases[bundle_class][path_id] = max(
                max_single_overlap_bases[bundle_class][path_id], single_bp,
            )
            max_single_coverages[bundle_class][path_id] = max(
                max_single_coverages[bundle_class][path_id], coverage,
            )
            if coverage >= full_coverage:
                qualifying[bundle_class].add(path_id)
                qualifying_mapping_counts[bundle_class][path_id] += 1

    decisions = []
    removed_classes = set()
    removed_numeric_ids = set()
    removed_names = set()
    all_classes = sorted(carriers, key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value))
    for bundle_class in all_classes:
        carrier_count = len(carriers[bundle_class])
        any_overlap_count = len(any_overlapping[bundle_class])
        qualifying_count = len(qualifying[bundle_class])
        fraction = qualifying_count / carrier_count if carrier_count else 0.0
        remove = fraction >= removal_frequency
        numeric_id, class_names = by_class.get(
            bundle_class, (normalize_numeric_id(bundle_class), set()),
        )
        class_names = set(class_names) | raw_names[bundle_class]
        if remove:
            removed_classes.add(bundle_class)
            removed_numeric_ids.add(numeric_id)
            removed_names.update(class_names)
            removed_names.update(by_numeric.get(numeric_id, ()))
        decisions.append((bundle_class, numeric_id, carrier_count, any_overlap_count,
                          qualifying_count,
                          f"{fraction:.9f}", int(remove),
                          "remove_all_class_consensuses" if remove else "retain"))

    retained_records = [record for record in records if record[0] not in removed_names]
    unknown_removed_names = sorted(removed_names - {record[0] for record in records})
    if removed_names and len(retained_records) == len(records):
        raise ValueError("TRF-selected classes did not match any final FASTA consensus names")

    output_dir.mkdir(parents=True, exist_ok=True)
    filtered_fasta = output_dir / "bundle.fa"
    filtered_metadata = output_dir / "bundle_metadata.bed"
    decision_path = output_dir / "bundle_summary.tsv"
    sequence_path = output_dir / "sequence_overlap.tsv"
    summary_path = output_dir / "final_trf_filter_summary.json"
    write_fasta(filtered_fasta, retained_records)
    with metadata_bed.open(encoding="utf-8", newline="") as source, \
            filtered_metadata.open("w", encoding="utf-8", newline="") as target:
        for line in source:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 4 and not line.lstrip().startswith("#"):
                numeric_id = normalize_numeric_id(fields[3].split(":", 1)[0])
                if numeric_id in removed_numeric_ids:
                    continue
            target.write(line)
    with decision_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("bundle_class", "numeric_bundle_id", "carrier_sequence_count",
                         "trf_any_overlap_sequence_count", "trf_qualifying_sequence_count",
                         "trf_qualifying_sequence_fraction",
                         "removed", "action"))
        writer.writerows(decisions)
    with sequence_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("bundle_class", "numeric_bundle_id", "population_path",
                         "mapping_count", "trf_any_overlap_mapping_count",
                         "trf_qualifying_mapping_count", "trf_total_overlap_bp",
                         "max_single_merged_trf_overlap_bp",
                         "max_single_merged_trf_coverage",
                         "trf_any_overlap_sequence", "trf_qualifying_sequence"))
        for bundle_class in all_classes:
            numeric_id = by_class.get(
                bundle_class, (normalize_numeric_id(bundle_class), set()),
            )[0]
            for path_id in sorted(carriers[bundle_class]):
                writer.writerow((bundle_class, numeric_id, path_id,
                                 mapping_counts[bundle_class][path_id],
                                 any_overlap_counts[bundle_class][path_id],
                                 qualifying_mapping_counts[bundle_class][path_id],
                                 overlap_bases[bundle_class][path_id],
                                 max_single_overlap_bases[bundle_class][path_id],
                                 f"{max_single_coverages[bundle_class][path_id]:.9f}",
                                 int(path_id in any_overlapping[bundle_class]),
                                 int(path_id in qualifying[bundle_class])))

    summary = {
        "version": 2,
        "completed": True,
        "semantics": (
            "per bundle class, each carrying population sequence contributes once; "
            "a sequence qualifies only when one selected merged TRF interval covers "
            "at least full_coverage of one mapped occurrence footprint; disconnected "
            "TRF intervals are not summed"
        ),
        "inputs": {
            "bundle_fasta": str(bundle_fasta.resolve()),
            "bundle_fasta_sha256": file_sha256(bundle_fasta),
            "bundle_metadata_bed": str(metadata_bed.resolve()),
            "coverage_mappings": str(coverage_mappings.resolve()),
            "class_map": str(class_map.resolve()),
            "trf_bed": str(trf_bed.resolve()),
            "trf_bed_sha256": file_sha256(trf_bed),
        },
        "parameters": {
            "trf_min_length": min_length,
            "trf_min_identity": min_identity,
            "trf_full_coverage": full_coverage,
            "removal_frequency": removal_frequency,
            "denominator": "distinct population sequences carrying the bundle class",
        },
        "trf": trf_summary,
        "counts": {
            "input_consensus_count": len(records),
            "retained_consensus_count": len(retained_records),
            "removed_consensus_count": len(records) - len(retained_records),
            "assessed_bundle_class_count": len(all_classes),
            "removed_bundle_class_count": len(removed_classes),
        },
        "removed_bundle_classes": sorted(removed_classes),
        "removed_numeric_bundle_ids": sorted(removed_numeric_ids, key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value)),
        "removed_consensus_names": sorted(removed_names),
        "removed_names_absent_from_fasta": unknown_removed_names,
        "outputs": {
            "bundle_fasta": str(filtered_fasta.resolve()),
            "bundle_metadata_bed": str(filtered_metadata.resolve()),
            "bundle_summary": str(decision_path.resolve()),
            "sequence_overlap": str(sequence_path.resolve()),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle_fasta", type=Path)
    parser.add_argument("bundle_metadata_bed", type=Path)
    parser.add_argument("coverage_mappings", type=Path)
    parser.add_argument("class_map", type=Path)
    parser.add_argument("trf_bed", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--trf-min-length", type=int, default=0)
    parser.add_argument("--trf-min-identity", type=float, default=80.0)
    parser.add_argument("--full-coverage", type=float, default=0.95)
    parser.add_argument("--removal-frequency", type=float, default=0.50)
    args = parser.parse_args(argv)
    try:
        summary = filter_final_library(
            args.bundle_fasta, args.bundle_metadata_bed, args.coverage_mappings,
            args.class_map, args.trf_bed, args.output_dir, args.trf_min_length,
            args.trf_min_identity, args.full_coverage, args.removal_frequency,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
