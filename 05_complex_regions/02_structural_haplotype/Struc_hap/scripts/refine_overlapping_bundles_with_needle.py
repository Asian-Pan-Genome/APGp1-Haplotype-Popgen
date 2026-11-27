#!/usr/bin/env python3
"""Reconcile one coverage-filtered genomic-overlap pair with EMBOSS Needle.

The coverage-only mapping table, captured before any locus deduplication, is
used to shortlist bundle pairs that overlap on the same population path.
Needle then supplies the sequence-space overlap
coordinates used by the shared bundle reconciliation engine.  Exactly one
largest actionable pair is changed per invocation; the shell workflow remaps
the resulting library and recomputes coverage before invoking this step again.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

try:
    from .align_bundles_with_mafft import reverse_complement
    from .refine_bundle_mappings import normalize_bundle_name
    from .split_bundles_by_common_deletions import metadata_sizes, read_fasta, wrap_fasta
    from .split_rank1_self_alignments import SelfAlignment, split_rank1_library
except ImportError:
    from align_bundles_with_mafft import reverse_complement
    from refine_bundle_mappings import normalize_bundle_name
    from split_bundles_by_common_deletions import metadata_sizes, read_fasta, wrap_fasta
    from split_rank1_self_alignments import SelfAlignment, split_rank1_library


@dataclass(frozen=True)
class MappingRow:
    raw_name: str
    bundle_class: str
    path: str
    start: int
    end: int
    strand: str


@dataclass
class PairEvidence:
    first_name: str
    second_name: str
    occurrence_count: int = 0
    maximum_overlap_bp: int = 0
    maximum_first_fraction: float = 0.0
    maximum_second_fraction: float = 0.0
    forward_overlap_bp: int = 0
    reverse_overlap_bp: int = 0

    @property
    def strand(self) -> str:
        return "+" if self.forward_overlap_bp >= self.reverse_overlap_bp else "-"


def read_mapping_rows(path: Path, records: Mapping[str, str]) -> List[MappingRow]:
    required = {"bundle_name", "raw_consensus_name", "population_path",
                "genome_start", "genome_end", "strand"}
    answer = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing mapping columns: {', '.join(sorted(missing))}")
        for line_number, row in enumerate(reader, 2):
            raw_name = row["raw_consensus_name"]
            if raw_name not in records:
                raise ValueError(f"{path}:{line_number}: unknown consensus {raw_name!r}")
            try:
                start, end = int(row["genome_start"]), int(row["genome_end"])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid genome coordinates") from error
            if start < 0 or end <= start or row["strand"] not in {"+", "-"}:
                raise ValueError(f"{path}:{line_number}: invalid mapping interval or strand")
            answer.append(MappingRow(raw_name, row["bundle_name"], row["population_path"],
                                     start, end, row["strand"]))
    return answer


def genomic_overlap_candidates(rows: Sequence[MappingRow], min_fraction: float,
                               min_bp: int) -> List[PairEvidence]:
    by_path: Dict[str, List[MappingRow]] = defaultdict(list)
    for row in rows:
        by_path[row.path].append(row)
    pairs: Dict[Tuple[str, str], PairEvidence] = {}
    for path_rows in by_path.values():
        ordered = sorted(path_rows, key=lambda row: (row.start, row.end, row.raw_name))
        for index, first in enumerate(ordered):
            for second in ordered[index + 1:]:
                if second.start >= first.end:
                    break
                if first.raw_name == second.raw_name or first.bundle_class == second.bundle_class:
                    continue
                overlap = min(first.end, second.end) - max(first.start, second.start)
                if overlap <= 0:
                    continue
                first_fraction = overlap / (first.end - first.start)
                second_fraction = overlap / (second.end - second.start)
                if not (first_fraction > min_fraction or second_fraction > min_fraction
                        or overlap > min_bp):
                    continue
                left, right = sorted((first.raw_name, second.raw_name))
                evidence = pairs.setdefault((left, right), PairEvidence(left, right))
                evidence.occurrence_count += 1
                evidence.maximum_overlap_bp = max(evidence.maximum_overlap_bp, overlap)
                if first.raw_name == left:
                    evidence.maximum_first_fraction = max(evidence.maximum_first_fraction,
                                                          first_fraction)
                    evidence.maximum_second_fraction = max(evidence.maximum_second_fraction,
                                                           second_fraction)
                else:
                    evidence.maximum_first_fraction = max(evidence.maximum_first_fraction,
                                                          second_fraction)
                    evidence.maximum_second_fraction = max(evidence.maximum_second_fraction,
                                                           first_fraction)
                if first.strand == second.strand:
                    evidence.forward_overlap_bp += overlap
                else:
                    evidence.reverse_overlap_bp += overlap
    return sorted(pairs.values(), key=lambda item: (-item.maximum_overlap_bp,
                                                    item.first_name, item.second_name))


def write_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{name}\n")
        handle.writelines(wrap_fasta(sequence))


def needle_alignment_coordinates(first: str, second: str) -> Tuple[int, int, int, int, int, int]:
    if len(first) != len(second) or not first:
        raise ValueError("Needle output must contain two nonempty equal-width records")
    first_residues, second_residues = first.replace("-", ""), second.replace("-", "")
    paired = [index for index, bases in enumerate(zip(first, second)) if "-" not in bases]
    if not paired:
        raise ValueError("Needle alignment has no column containing both sequences")
    left, right = paired[0], paired[-1] + 1
    first_start = sum(base != "-" for base in first[:left])
    second_start = sum(base != "-" for base in second[:left])
    first_end = first_start + sum(base != "-" for base in first[left:right])
    second_end = second_start + sum(base != "-" for base in second[left:right])
    block = [(left_base, right_base) for left_base, right_base in zip(first[left:right], second[left:right])
             if left_base != "-" or right_base != "-"]
    matches = sum(left_base == right_base and left_base != "-" for left_base, right_base in block)
    return first_start, first_end, second_start, second_end, matches, len(block)


def run_needle_pair(executable: str, directory: Path, index: int, evidence: PairEvidence,
                    records: Mapping[str, str], gapopen: float, gapextend: float) -> Tuple[SelfAlignment, dict]:
    pair_dir = directory / "needle" / f"pair_{index:04d}"
    pair_dir.mkdir(parents=True, exist_ok=True)
    first_sequence = records[evidence.first_name]
    original_second = records[evidence.second_name]
    second_sequence = (original_second if evidence.strand == "+"
                       else reverse_complement(original_second))
    first_fasta, second_fasta = pair_dir / "first.fa", pair_dir / "second.fa"
    alignment_fasta, log = pair_dir / "alignment.fa", pair_dir / "needle.log"
    write_fasta(first_fasta, "first", first_sequence)
    write_fasta(second_fasta, "second", second_sequence)
    command = [executable, "-auto", "-asequence", str(first_fasta),
               "-bsequence", str(second_fasta), "-gapopen", str(gapopen),
               "-gapextend", str(gapextend), "-endweight", "N",
               "-aformat3", "fasta", "-outfile", str(alignment_fasta)]
    with log.open("w", encoding="utf-8") as error:
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=error, check=True)
    aligned = read_fasta(alignment_fasta)
    if set(aligned) != {"first", "second"}:
        raise ValueError(f"{alignment_fasta}: Needle did not emit first/second FASTA records")
    if aligned["first"].replace("-", "") != first_sequence.upper():
        raise ValueError(f"{alignment_fasta}: Needle changed residues in the first sequence")
    if aligned["second"].replace("-", "") != second_sequence.upper():
        raise ValueError(f"{alignment_fasta}: Needle changed residues in the second sequence")
    first_start, first_end, oriented_start, oriented_end, matches, block_length = \
        needle_alignment_coordinates(aligned["first"], aligned["second"])
    if evidence.strand == "+":
        second_start, second_end = oriented_start, oriented_end
    else:
        second_start = len(original_second) - oriented_end
        second_end = len(original_second) - oriented_start
    alignment = SelfAlignment(
        evidence.first_name, len(first_sequence), first_start, first_end,
        evidence.second_name, len(original_second), second_start, second_end,
        evidence.strand, matches, block_length, 60, index,
    )
    return alignment, {
        "pair": [evidence.first_name, evidence.second_name],
        "relative_strand": evidence.strand,
        "command": command,
        "alignment": str(alignment_fasta),
        "log": str(log),
        "identity": matches / block_length if block_length else 0.0,
        "first_interval": [first_start, first_end],
        "second_interval": [second_start, second_end],
    }


def write_candidate_table(path: Path, candidates: Sequence[PairEvidence]) -> None:
    columns = ("first_name", "second_name", "occurrence_count", "maximum_overlap_bp",
               "maximum_first_fraction", "maximum_second_fraction",
               "forward_overlap_bp", "reverse_overlap_bp", "needle_strand")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for item in candidates:
            writer.writerow((item.first_name, item.second_name, item.occurrence_count,
                             item.maximum_overlap_bp, f"{item.maximum_first_fraction:.10g}",
                             f"{item.maximum_second_fraction:.10g}", item.forward_overlap_bp,
                             item.reverse_overlap_bp, item.strand))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle_fasta", type=Path)
    parser.add_argument("bundle_metadata_bed", type=Path)
    parser.add_argument("mapping_tsv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--needle", default="needle")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-overlap-fraction", type=float, default=0.20)
    parser.add_argument("--min-overlap-bp", type=int, default=10000)
    parser.add_argument("--min-identity", type=float, default=0.0)
    parser.add_argument("--merge-identity", type=float, default=0.90)
    parser.add_argument("--gapopen", type=float, default=10.0)
    parser.add_argument("--gapextend", type=float, default=0.5)
    parser.add_argument("--keep-intermediates", action="store_true")
    return parser


def validate(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    for path in (args.bundle_fasta, args.bundle_metadata_bed, args.mapping_tsv):
        if not path.is_file():
            parser.error(f"input file does not exist: {path}")
    for name in ("min_overlap_fraction", "min_identity", "merge_identity"):
        value = getattr(args, name)
        if not math.isfinite(value) or not 0 <= value <= 1:
            parser.error(f"--{name.replace('_', '-')} must be finite and in [0, 1]")
    if args.min_overlap_bp < 0 or args.threads < 1:
        parser.error("overlap bp must be non-negative and threads must be positive")
    if not math.isfinite(args.gapopen) or args.gapopen < 0:
        parser.error("--gapopen must be finite and non-negative")
    if not math.isfinite(args.gapextend) or args.gapextend < 0:
        parser.error("--gapextend must be finite and non-negative")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate(parser, args)
    records = read_fasta(args.bundle_fasta)
    _, sizes = metadata_sizes(args.bundle_metadata_bed)
    rows = read_mapping_rows(args.mapping_tsv, records)
    candidates = genomic_overlap_candidates(rows, args.min_overlap_fraction, args.min_overlap_bp)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_candidate_table(args.output_dir / "mapping_overlap_candidates.tsv", candidates)
    needle_runs: List[dict] = []

    def provider(active_records, directory, iteration):
        if iteration != 1:
            return [], {"candidate_pair_count": 0, "needle_run_count": 0,
                        "reason": "single_operation_per_coverage_filter_round"}
        active = [item for item in candidates
                  if item.first_name in active_records and item.second_name in active_records]
        alignments = []
        for index, evidence in enumerate(active, 1):
            alignment, audit = run_needle_pair(
                args.needle, directory, index, evidence, active_records,
                args.gapopen, args.gapextend,
            )
            alignments.append(alignment)
            needle_runs.append(audit)
        return alignments, {"candidate_pair_count": len(active),
                            "needle_run_count": len(alignments)}

    result = split_rank1_library(
        records, sizes, args.output_dir, args.needle, args.threads,
        args.min_overlap_fraction, args.min_overlap_bp, 0,
        args.min_identity, 0, stage_label="post-filter-needle",
        reconcile_overlaps=True, merge_identity=args.merge_identity,
        keep_intermediates=args.keep_intermediates,
        alignment_provider=provider, alignment_backend="emboss_needle",
    )
    summary = {
        "completed": True,
        "parameters": {
            "min_genomic_overlap_fraction": args.min_overlap_fraction,
            "min_genomic_overlap_bp": args.min_overlap_bp,
            "min_sequence_overlap_fraction": args.min_overlap_fraction,
            "min_sequence_overlap_bp": args.min_overlap_bp,
            "min_identity": args.min_identity,
            "merge_identity": args.merge_identity,
            "gapopen": args.gapopen,
            "gapextend": args.gapextend,
        },
        "counts": {
            "coverage_filtered_mapping_count": len(rows),
            "mapping_overlap_candidate_pair_count": len(candidates),
            "needle_alignment_count": len(needle_runs),
            "reconciliation_operation_count": result["summary"]["counts"]["reconciliation_operation_count"],
            "split_operation_count": result["summary"]["counts"]["split_operation_count"],
        },
        "needle_alignments": needle_runs,
        "outputs": {
            "bundle_fasta": str((args.output_dir / "bundle.fa").resolve()),
            "bundle_metadata_bed": str((args.output_dir / "bundle_metadata.bed").resolve()),
            "mapping_overlap_candidates": str((args.output_dir / "mapping_overlap_candidates.tsv").resolve()),
            "reconciliation_summary": str(Path(result["summary_path"]).resolve()),
            "reconciliation_events": str((args.output_dir / "self_split_events.tsv").resolve()),
        },
    }
    summary_path = args.output_dir / "post_filter_needle_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    json.dump(summary, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
