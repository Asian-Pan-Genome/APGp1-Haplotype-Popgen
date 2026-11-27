#!/usr/bin/env python3
"""Run the iterative bundle self-comparison on a final split library."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Optional, Sequence

try:
    from .refine_bundle_mappings import normalize_bundle_name
    from .split_bundles_by_common_deletions import metadata_sizes, read_fasta
    from .split_rank1_self_alignments import split_rank1_library
except ImportError:
    from refine_bundle_mappings import normalize_bundle_name
    from split_bundles_by_common_deletions import metadata_sizes, read_fasta
    from split_rank1_self_alignments import split_rank1_library


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Iteratively self-align the final post-split bundle library, reconcile "
            "overlapping consensuses, and write a replacement FASTA plus metadata."
        )
    )
    parser.add_argument("bundle_fasta", type=Path)
    parser.add_argument("bundle_metadata_bed", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-overlap-fraction", type=float, default=0.20)
    parser.add_argument("--min-overlap-bp", type=int, default=10000)
    parser.add_argument("--min-identity", type=float, default=0.0)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--merge-identity", type=float, default=0.90)
    parser.add_argument(
        "--no-overlap-reconciliation", action="store_true",
        help="use legacy non-destructive cuts instead of shared-family reconciliation",
    )
    parser.add_argument("--keep-intermediates", action="store_true")
    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    for path in (args.bundle_fasta, args.bundle_metadata_bed):
        if not path.is_file():
            parser.error(f"input file does not exist: {path}")
    if args.threads < 1:
        parser.error("--threads must be positive")
    if not math.isfinite(args.min_overlap_fraction) or not 0 <= args.min_overlap_fraction <= 1:
        parser.error("--min-overlap-fraction must be finite and in [0, 1]")
    if args.min_overlap_bp < 0:
        parser.error("--min-overlap-bp must be non-negative")
    if not math.isfinite(args.min_identity) or not 0 <= args.min_identity <= 1:
        parser.error("--min-identity must be finite and in [0, 1]")
    if args.min_mapq < 0:
        parser.error("--min-mapq must be non-negative")
    if not math.isfinite(args.merge_identity) or not 0 <= args.merge_identity <= 1:
        parser.error("--merge-identity must be finite and in [0, 1]")
    inputs = {args.bundle_fasta.resolve(), args.bundle_metadata_bed.resolve()}
    if args.output_dir.resolve() in inputs:
        parser.error("--output-dir must be a directory distinct from both input files")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    records = read_fasta(args.bundle_fasta)
    _, sizes = metadata_sizes(args.bundle_metadata_bed)
    active_ids = {normalize_bundle_name(name, {}) for name in records}
    invalid = sorted(bundle_id for bundle_id in active_ids
                     if not bundle_id.isdigit() or str(int(bundle_id)) != bundle_id)
    if invalid:
        raise ValueError(f"final self alignment requires numeric bundle IDs: {invalid[0]}")
    missing = sorted(active_ids - sizes.keys(), key=int)
    if missing:
        raise ValueError(f"bundle metadata is missing active class {missing[0]}")
    result = split_rank1_library(
        records, sizes, args.output_dir, args.minimap2, args.threads,
        args.min_overlap_fraction, args.min_overlap_bp, 0,
        args.min_identity, args.min_mapq, stage_label="post-final-split",
        reconcile_overlaps=not args.no_overlap_reconciliation,
        merge_identity=args.merge_identity,
        keep_intermediates=args.keep_intermediates,
    )
    json.dump({
        "summary": str(Path(result["summary_path"]).resolve()),
        "counts": result["summary"]["counts"],
        "outputs": result["summary"]["outputs"],
    }, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
