#!/usr/bin/env python3
"""Build one majority consensus sequence from each per-bundle MSA."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class AlignmentInput:
    bundle_id: str
    occurrence_count: int
    alignment_path: Path
    status: str


@dataclass(frozen=True)
class ConsensusResult:
    bundle_id: str
    occurrence_count: int
    sequence_count: int
    alignment_length: int
    consensus_length: int
    retained_column_count: int
    dropped_low_occupancy_column_count: int
    ambiguous_base_count: int
    mean_column_occupancy: float
    mean_retained_column_occupancy: float
    median_retained_column_occupancy: float
    consensus_path: Path
    consensus_sequence: str


def bundle_sort_key(bundle_id: str) -> Tuple[int, object]:
    return (0, int(bundle_id)) if bundle_id.isdigit() else (1, bundle_id)


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = stripped[1:]
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"{path}:{line_number}: sequence before first header")
                chunks.append(stripped)
    if header is not None:
        yield header, "".join(chunks).upper()


def read_alignment(path: Path) -> List[str]:
    sequences = [sequence for _, sequence in iter_fasta(path)]
    if not sequences:
        raise ValueError(f"{path}: empty alignment")
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError(f"{path}: aligned sequences do not have one common length")
    if not sequences[0]:
        raise ValueError(f"{path}: alignment sequences are empty")
    return sequences


def resolve_alignment_path(manifest_path: Path, value: str) -> Path:
    candidate = Path(value)
    if candidate.is_file():
        return candidate
    relative = manifest_path.parent / candidate
    return relative if relative.is_file() else candidate


def read_alignment_manifest(path: Path) -> List[AlignmentInput]:
    inputs: List[AlignmentInput] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"bundle_id", "occurrence_count", "alignment_fasta", "status"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"{path}: manifest must contain {', '.join(sorted(required))}"
            )
        for row_number, row in enumerate(reader, 2):
            alignment_path = resolve_alignment_path(path, row["alignment_fasta"])
            if not alignment_path.is_file():
                raise ValueError(
                    f"{path}:{row_number}: alignment not found: {row['alignment_fasta']}"
                )
            inputs.append(
                AlignmentInput(
                    bundle_id=row["bundle_id"],
                    occurrence_count=int(row["occurrence_count"]),
                    alignment_path=alignment_path,
                    status=row["status"],
                )
            )
    if not inputs:
        raise ValueError(f"{path}: alignment manifest has no data rows")
    ids = [item.bundle_id for item in inputs]
    duplicated = sorted(
        (bundle_id for bundle_id, count in Counter(ids).items() if count > 1),
        key=bundle_sort_key,
    )
    if duplicated:
        raise ValueError(f"{path}: duplicate bundle IDs: {', '.join(duplicated)}")
    return sorted(inputs, key=lambda item: bundle_sort_key(item.bundle_id))


def build_consensus(
    sequences: Sequence[str],
    min_column_occupancy: float,
    min_base_frequency: float,
) -> Tuple[str, List[float], int, float]:
    """Return consensus, retained occupancies, dropped count, and mean occupancy."""
    sequence_count = len(sequences)
    alignment_length = len(sequences[0])
    consensus: List[str] = []
    retained_occupancies: List[float] = []
    dropped = 0
    occupancy_sum = 0.0
    for column in range(alignment_length):
        counts = Counter(sequence[column] for sequence in sequences)
        canonical = {base: counts[base] for base in "ACGT" if counts[base]}
        canonical_count = sum(canonical.values())
        occupancy = canonical_count / sequence_count
        occupancy_sum += occupancy
        if canonical_count == 0 or occupancy < min_column_occupancy:
            dropped += 1
            continue
        retained_occupancies.append(occupancy)
        largest = max(canonical.values())
        winners = sorted(base for base, count in canonical.items() if count == largest)
        if len(winners) == 1 and largest / canonical_count >= min_base_frequency:
            consensus.append(winners[0])
        else:
            consensus.append("N")
    return (
        "".join(consensus),
        retained_occupancies,
        dropped,
        occupancy_sum / alignment_length,
    )


def wrap_fasta(sequence: str, width: int) -> Iterable[str]:
    for offset in range(0, len(sequence), width):
        yield sequence[offset : offset + width] + "\n"


def process_alignment(
    item: AlignmentInput,
    consensus_dir: Path,
    min_column_occupancy: float,
    min_base_frequency: float,
    line_width: int,
) -> ConsensusResult:
    sequences = read_alignment(item.alignment_path)
    consensus, retained_occupancies, dropped, mean_column_occupancy = build_consensus(
        sequences,
        min_column_occupancy,
        min_base_frequency,
    )
    if not consensus:
        raise ValueError(
            f"bundle {item.bundle_id}: no alignment column passed the occupancy threshold"
        )
    consensus_path = consensus_dir / f"bundle_{item.bundle_id}.consensus.fa"
    with consensus_path.open("w", encoding="utf-8") as handle:
        handle.write(
            f">bundle_{item.bundle_id} occurrence_count={item.occurrence_count} "
            f"sequence_count={len(sequences)} alignment_length={len(sequences[0])} "
            f"consensus_length={len(consensus)}\n"
        )
        handle.writelines(wrap_fasta(consensus, line_width))
    return ConsensusResult(
        bundle_id=item.bundle_id,
        occurrence_count=item.occurrence_count,
        sequence_count=len(sequences),
        alignment_length=len(sequences[0]),
        consensus_length=len(consensus),
        retained_column_count=len(retained_occupancies),
        dropped_low_occupancy_column_count=dropped,
        ambiguous_base_count=sum(base not in "ACGT" for base in consensus),
        mean_column_occupancy=mean_column_occupancy,
        mean_retained_column_occupancy=float(statistics.fmean(retained_occupancies)),
        median_retained_column_occupancy=float(statistics.median(retained_occupancies)),
        consensus_path=consensus_path,
        consensus_sequence=consensus,
    )


def write_combined_fasta(path: Path, results: Sequence[ConsensusResult], width: int) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for result in results:
            handle.write(
                f">bundle_{result.bundle_id} occurrence_count={result.occurrence_count} "
                f"sequence_count={result.sequence_count} "
                f"consensus_length={result.consensus_length}\n"
            )
            handle.writelines(wrap_fasta(result.consensus_sequence, width))


def write_summary(path: Path, results: Sequence[ConsensusResult]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "bundle_id", "occurrence_count", "alignment_sequence_count",
                "alignment_length", "consensus_length", "retained_column_count",
                "dropped_low_occupancy_column_count", "ambiguous_base_count",
                "mean_column_occupancy", "mean_retained_column_occupancy",
                "median_retained_column_occupancy", "consensus_fasta",
            )
        )
        for result in results:
            writer.writerow(
                (
                    result.bundle_id,
                    result.occurrence_count,
                    result.sequence_count,
                    result.alignment_length,
                    result.consensus_length,
                    result.retained_column_count,
                    result.dropped_low_occupancy_column_count,
                    result.ambiguous_base_count,
                    f"{result.mean_column_occupancy:.6f}",
                    f"{result.mean_retained_column_occupancy:.6f}",
                    f"{result.median_retained_column_occupancy:.6f}",
                    result.consensus_path,
                )
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build one gap-filtered majority consensus per bundle MSA."
    )
    parser.add_argument(
        "alignment_manifest",
        type=Path,
        help="bundle_alignments.tsv produced by align_bundles_with_mafft.py",
    )
    parser.add_argument("-o", "--output-dir", required=True, type=Path)
    parser.add_argument("--min-column-occupancy", type=float, default=0.50)
    parser.add_argument("--min-base-frequency", type=float, default=0.50)
    parser.add_argument("--line-width", type=int, default=80)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    for name in ("min_column_occupancy", "min_base_frequency"):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            parser.error(f"--{name.replace('_', '-')} must be in [0, 1]")
    if args.line_width <= 0:
        parser.error("--line-width must be positive")
    try:
        inputs = read_alignment_manifest(args.alignment_manifest)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        consensus_dir = args.output_dir / "per_bundle"
        consensus_dir.mkdir(parents=True, exist_ok=True)
        results = [
            process_alignment(
                item,
                consensus_dir,
                args.min_column_occupancy,
                args.min_base_frequency,
                args.line_width,
            )
            for item in inputs
        ]
        write_combined_fasta(
            args.output_dir / "bundle_consensus.fa", results, args.line_width
        )
        write_summary(args.output_dir / "bundle_consensus.tsv", results)
        print(
            f"Wrote {len(results)} bundle consensus sequences to "
            f"{args.output_dir / 'bundle_consensus.fa'}"
        )
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
