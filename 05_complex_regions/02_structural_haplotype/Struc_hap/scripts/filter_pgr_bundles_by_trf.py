#!/usr/bin/env python3
"""Filter original PGR-TK bundle IDs by population single-interval TRF coverage.

This adapter deliberately writes surviving original BED lines verbatim. It
does not alter genomic coordinates, minimizer metadata, directions or copy
flags, and does not infer repeat coverage from synthetic bundle coordinates.
"""

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Dict, Optional

try:
    from .merge_pgr_bundles import (
        filter_input_pgr_by_trf, read_pgr_bed, read_trf_bed, write_tsv,
    )
except ImportError:
    from merge_pgr_bundles import (
        filter_input_pgr_by_trf, read_pgr_bed, read_trf_bed, write_tsv,
    )


REPORT_HEADER = (
    "path_id", "bundle_start", "bundle_end", "bundle_length", "bundle_id",
    "direction", "source", "action", "repeat_bp", "repeat_coverage",
    "repeat_union_interval_count", "repeat_union_intervals",
    "internal_repeat_count", "internal_repeat_bp", "max_internal_repeat_bp",
    "internal_repeat_intervals", "future_break_candidates", "raw_trf_count",
    "raw_trf_records",
)
BUNDLE_SUMMARY_HEADER = (
    "bundle_id", "occurrence_count", "repeat_occurrence_count",
    "repeat_occurrence_fraction", "internal_repeat_occurrence_count",
    "full_repeat_occurrence_count", "dropped_occurrence_count",
    "remaining_occurrence_count", "median_repeat_bp", "max_repeat_bp",
    "median_repeat_coverage", "max_repeat_coverage", "action",
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def filter_pgr_bed_by_trf(
    bed: Path,
    trf_bed: Optional[Path],
    output_dir: Path,
    min_length: int = 0,
    min_identity: float = 80.0,
    full_coverage: float = 0.90,
    removal_frequency: float = 0.50,
) -> Dict[str, object]:
    """Remove all occurrences of IDs with recurrent single-interval coverage.

    An occurrence qualifies only when one merged TRF interval covers at least
    ``full_coverage`` of it. Disconnected intervals are not summed. Missing or
    empty annotations, no selected TRF rows, or selected TRF paths not matching
    any bundle path are safe exact-copy no-ops.
    """
    bed = Path(bed)
    trf_bed = Path(trf_bed) if trf_bed is not None else None
    output_dir = Path(output_dir)
    if not isinstance(min_length, int) or min_length < 0:
        raise ValueError("TRF minimum length must be a nonnegative integer")
    for name, value, upper in (
        ("minimum identity", min_identity, 100),
        ("full coverage", full_coverage, 1),
        ("removal frequency", removal_frequency, 1),
    ):
        if not math.isfinite(value) or not 0 <= value <= upper:
            raise ValueError("TRF %s must be finite and between 0 and %s" % (name, upper))
    output_bed = output_dir / "filtered_bundles.pgr.bed"
    report_path = output_dir / "pregroup_report.tsv"
    bundle_summary_path = output_dir / "pregroup_bundle_summary.tsv"
    summary_path = output_dir / "filter_summary.json"
    protected = {bed.resolve()}
    if trf_bed is not None:
        protected.add(trf_bed.resolve())
    for target in (output_bed, report_path, bundle_summary_path, summary_path):
        if target.resolve() in protected:
            raise ValueError("TRF filter output would overwrite an input: %s" % target)

    # Keep line-number identities from the original parser for verbatim output.
    with bed.open("r", encoding="utf-8", newline="") as handle:
        original_lines = handle.readlines()
    has_records = any(line.strip() and not line.lstrip().startswith("#")
                      for line in original_lines)
    by_path = read_pgr_bed(bed) if has_records else {}
    intervals, records = {}, {}
    read_summary = {
        "trf_input_row_count": 0,
        "trf_selected_row_count": 0,
        "trf_union_interval_count": 0,
    }
    if trf_bed is not None:
        intervals, records, read_summary = read_trf_bed(
            trf_bed, min_length, min_identity,
        )
        # The shared legacy reader converts identity with float(), which also
        # accepts NaN/infinity. Such values must never count as selected repeat
        # evidence: NaN even bypasses its identity-threshold comparison.
        for path_id, items in records.items():
            for record in items:
                if not math.isfinite(record.identity) or not 0 <= record.identity <= 100:
                    raise ValueError(
                        "%s: selected TRF identity for %s must be finite and in [0, 100]"
                        % (trf_bed, path_id)
                    )
    matched_paths = sorted(set(by_path) & set(intervals))
    enabled = bool(matched_paths)
    if trf_bed is None:
        no_op_reason = "no_trf_annotation"
    elif not read_summary["trf_input_row_count"]:
        no_op_reason = "empty_trf_annotation"
    elif not read_summary["trf_selected_row_count"]:
        no_op_reason = "no_trf_rows_pass_filters"
    elif not matched_paths:
        no_op_reason = "no_matching_bundle_paths"
    else:
        no_op_reason = None

    # A no-op must stay a no-op even when the user sets removal_frequency=0.
    filtered, report_rows, summary_rows, action_summary = filter_input_pgr_by_trf(
        by_path, intervals, records, full_coverage,
        removal_frequency if enabled else 1.0,
    )
    kept_lines = {item.line_number for items in filtered.values() for item in items}
    removed_ids = sorted(str(row[0]) for row in summary_rows
                         if row[12] == "removed_all_occurrences")
    output_dir.mkdir(parents=True, exist_ok=True)
    if not removed_ids:
        shutil.copyfile(bed, output_bed)
    else:
        with output_bed.open("w", encoding="utf-8", newline="") as handle:
            for line_number, line in enumerate(original_lines, 1):
                if not line.strip() or line.lstrip().startswith("#") or line_number in kept_lines:
                    handle.write(line)
                else:
                    # Preserve original occurrence row IDs used by ranking and
                    # sequence/TRF provenance, without retaining a BED record.
                    handle.write("# TRF_removed: " + line)
    write_tsv(report_path, REPORT_HEADER, report_rows)
    write_tsv(bundle_summary_path, BUNDLE_SUMMARY_HEADER, summary_rows)
    fingerprint_inputs = {
        "bed": str(bed.resolve()),
        "bed_sha256": file_sha256(bed),
        "trf_bed": str(trf_bed.resolve()) if trf_bed is not None else None,
        "trf_bed_sha256": file_sha256(trf_bed) if trf_bed is not None else None,
        "min_length": min_length,
        "min_identity": min_identity,
        "full_coverage": full_coverage,
        "removal_frequency": removal_frequency,
    }
    summary = {
        "version": 2,
        "enabled": enabled,
        "trf_annotation_enabled": bool(read_summary["trf_input_row_count"]),
        "no_op_reason": no_op_reason,
        "fingerprint_inputs": fingerprint_inputs,
        "fingerprint": hashlib.sha256(json.dumps(
            fingerprint_inputs, sort_keys=True, separators=(",", ":"),
        ).encode("utf-8")).hexdigest(),
        "input_bed_sha256": fingerprint_inputs["bed_sha256"],
        "trf_bed_sha256": fingerprint_inputs["trf_bed_sha256"],
        "input_path_ids": sorted(by_path),
        "retained_path_ids": sorted(path for path, items in filtered.items() if items),
        "empty_path_ids": sorted(path for path, items in filtered.items() if not items),
        "matched_trf_path_ids": matched_paths,
        "removed_bundle_ids": removed_ids,
        "input_occurrence_count": sum(map(len, by_path.values())),
        "retained_occurrence_count": len(kept_lines),
        "filtered_bed": str(output_bed.resolve()),
        "report": str(report_path.resolve()),
        "bundle_summary": str(bundle_summary_path.resolve()),
        "semantics": (
            "one merged TRF interval must meet full_coverage per occurrence; "
            "remove all occurrences when qualifying frequency >= threshold; "
            "disconnected intervals are not summed"
        ),
        **read_summary,
        **action_summary,
    }
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return {
        "bed": output_bed, "summary": summary, "summary_path": summary_path,
        "intervals": intervals, "records": records,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bed", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--trf-bed", type=Path)
    parser.add_argument(
        "--trf-min-length", type=int, default=0,
        help="minimum TRF interval length (default: 0; no length filter)",
    )
    parser.add_argument("--trf-min-identity", type=float, default=80.0)
    parser.add_argument(
        "--trf-full-coverage", type=float, default=0.90,
        help="coverage required from one merged TRF interval (default: 0.90)",
    )
    parser.add_argument(
        "--trf-bundle-removal-frequency", type=float, default=0.50,
        help="fraction of occurrences required for whole-ID removal (default: 0.50)",
    )
    args = parser.parse_args(argv)
    try:
        result = filter_pgr_bed_by_trf(
            args.bed, args.trf_bed, args.output_dir, args.trf_min_length,
            args.trf_min_identity, args.trf_full_coverage,
            args.trf_bundle_removal_frequency,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(result["bed"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
