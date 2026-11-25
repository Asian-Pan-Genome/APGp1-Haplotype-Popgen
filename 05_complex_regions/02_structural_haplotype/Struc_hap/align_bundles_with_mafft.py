#!/usr/bin/env python3
"""Extract every final bundle occurrence and align each bundle with MAFFT.

BED coordinates are interpreted as zero-based, half-open intervals. By default,
reverse-strand occurrences are reverse-complemented before alignment so that all
sequences for a bundle use the orientation encoded by the final PGR walk.
"""

import argparse
import csv
import shlex
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple
from urllib.parse import quote

try:
    from .build_grouped_bundle_consensus import (
        group_paths_by_order,
        select_group_bundle_ids,
        write_bundle_source_groups,
    )
except ImportError:
    from build_grouped_bundle_consensus import (
        group_paths_by_order,
        select_group_bundle_ids,
        write_bundle_source_groups,
    )

try:
    from .short_bundle_alignment import SHORT_BUNDLE_MAFFT_MIN_BP, align_short_bundle_records
except ImportError:
    from short_bundle_alignment import SHORT_BUNDLE_MAFFT_MIN_BP, align_short_bundle_records


@dataclass(frozen=True)
class BundleOccurrence:
    row_number: int
    path_id: str
    start: int
    end: int
    bundle_id: str
    direction: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class BundleInput:
    bundle_id: str
    fasta_path: Path
    occurrence_count: int
    forward_count: int
    reverse_count: int
    total_bp: int
    min_bp: int
    median_bp: float
    max_bp: int


IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
    "TGCAYRKMSWVHDBNtgcayrkmswvhdbn",
)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(IUPAC_COMPLEMENT)[::-1]


def normalize_direction(value: str) -> str:
    normalized = value.strip().lower()
    if normalized in {"0", "+", "f", "forward", "fw"}:
        return "+"
    if normalized in {"1", "-", "r", "reverse", "rev", "rc"}:
        return "-"
    raise ValueError(f"unsupported bundle direction: {value!r}")


def bundle_sort_key(bundle_id: str) -> Tuple[int, object]:
    return (0, int(bundle_id)) if bundle_id.isdigit() else (1, bundle_id)


def bundle_filename(bundle_id: str) -> str:
    return f"bundle_{quote(bundle_id, safe='._-')}.fa"


def read_pgr_occurrences(
    bed_path: Path,
    selected_bundle_ids: Optional[Set[str]] = None,
) -> Tuple[Dict[str, List[BundleOccurrence]], Dict[str, List[BundleOccurrence]]]:
    """Read native four-column PGR BED and index occurrences by path and bundle."""
    by_path: DefaultDict[str, List[BundleOccurrence]] = defaultdict(list)
    by_bundle: DefaultDict[str, List[BundleOccurrence]] = defaultdict(list)
    with bed_path.open("r", encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 4:
                fields = stripped.split()
            if len(fields) < 4:
                raise ValueError(f"{bed_path}:{row_number}: expected at least four columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{bed_path}:{row_number}: start/end must be integers"
                ) from exc
            if start < 0 or end <= start:
                raise ValueError(
                    f"{bed_path}:{row_number}: invalid half-open interval {start}-{end}"
                )
            annotation = fields[3].split(":")
            if len(annotation) < 3:
                raise ValueError(
                    f"{bed_path}:{row_number}: expected native PGR bundle annotation"
                )
            bundle_id = annotation[0]
            if selected_bundle_ids is not None and bundle_id not in selected_bundle_ids:
                continue
            occurrence = BundleOccurrence(
                row_number=row_number,
                path_id=fields[0],
                start=start,
                end=end,
                bundle_id=bundle_id,
                direction=normalize_direction(annotation[2]),
            )
            by_path[occurrence.path_id].append(occurrence)
            by_bundle[bundle_id].append(occurrence)
    if not by_bundle:
        raise ValueError(f"{bed_path}: no selected PGR bundle occurrences")
    return dict(by_path), dict(by_bundle)


def signed_path_signature(
    occurrences: Sequence[BundleOccurrence],
) -> Tuple[str, ...]:
    ordered = sorted(
        occurrences, key=lambda item: (item.start, item.end, item.row_number)
    )
    return tuple(f"{item.direction}{item.bundle_id}" for item in ordered)


def largest_path_order_group(
    occurrences_by_path: Mapping[str, Sequence[BundleOccurrence]],
) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    by_signature: DefaultDict[Tuple[str, ...], List[str]] = defaultdict(list)
    for path_id, occurrences in occurrences_by_path.items():
        by_signature[signed_path_signature(occurrences)].append(path_id)
    signature, path_ids = min(
        by_signature.items(),
        key=lambda item: (-len(item[1]), item[0]),
    )
    return signature, tuple(sorted(path_ids))


def filter_occurrences(
    occurrences_by_path: Mapping[str, Sequence[BundleOccurrence]],
    selected_path_ids: Set[str],
    selected_bundle_ids: Optional[Set[str]],
) -> Tuple[Dict[str, List[BundleOccurrence]], Dict[str, List[BundleOccurrence]]]:
    by_path: DefaultDict[str, List[BundleOccurrence]] = defaultdict(list)
    by_bundle: DefaultDict[str, List[BundleOccurrence]] = defaultdict(list)
    for path_id in sorted(selected_path_ids):
        for occurrence in occurrences_by_path[path_id]:
            if (
                selected_bundle_ids is not None
                and occurrence.bundle_id not in selected_bundle_ids
            ):
                continue
            by_path[path_id].append(occurrence)
            by_bundle[occurrence.bundle_id].append(occurrence)
    if not by_bundle:
        raise ValueError("no bundle occurrences remain after path/bundle selection")
    return dict(by_path), dict(by_bundle)


def write_selected_path_group(
    path: Path,
    signature: Sequence[str],
    path_ids: Sequence[str],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("rank", "haplotype_count", "signed_bundle_order", "path_ids")
        )
        writer.writerow((1, len(path_ids), ",".join(signature), ",".join(path_ids)))


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Stream FASTA records while holding at most one source sequence in memory."""
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                header = stripped[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                name = header.split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"{path}:{line_number}: sequence before first header")
                chunks.append(stripped)
    if name is not None:
        yield name, "".join(chunks)


def wrap_fasta(sequence: str, width: int = 80) -> Iterable[str]:
    for offset in range(0, len(sequence), width):
        yield sequence[offset : offset + width] + "\n"


def extraction_header(occurrence: BundleOccurrence) -> str:
    path_id = quote(occurrence.path_id, safe="._#-")
    return (
        f">bundle={quote(occurrence.bundle_id, safe='._-')}"
        f"|path={path_id}|start={occurrence.start}|end={occurrence.end}"
        f"|strand={occurrence.direction}|bed_row={occurrence.row_number}\n"
    )


def extract_bundle_fastas(
    fasta_path: Path,
    occurrences_by_path: Mapping[str, Sequence[BundleOccurrence]],
    occurrences_by_bundle: Mapping[str, Sequence[BundleOccurrence]],
    sequence_dir: Path,
    orient_to_bundle: bool = True,
) -> List[BundleInput]:
    """Extract BED intervals into one FASTA file per final bundle ID."""
    sequence_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {
        bundle_id: sequence_dir / bundle_filename(bundle_id)
        for bundle_id in occurrences_by_bundle
    }
    for output_path in output_paths.values():
        output_path.write_text("", encoding="utf-8")

    seen_fasta_names: Set[str] = set()
    extracted_counts: Counter[str] = Counter()
    requested_paths = set(occurrences_by_path)
    for path_id, source_sequence in iter_fasta(fasta_path):
        if path_id in seen_fasta_names:
            raise ValueError(f"{fasta_path}: duplicate FASTA record name {path_id!r}")
        seen_fasta_names.add(path_id)
        requested = occurrences_by_path.get(path_id)
        if not requested:
            continue
        output_chunks: DefaultDict[str, List[str]] = defaultdict(list)
        for occurrence in requested:
            if occurrence.end > len(source_sequence):
                raise ValueError(
                    f"{fasta_path}:{path_id}: BED interval {occurrence.start}-{occurrence.end} "
                    f"exceeds sequence length {len(source_sequence)}"
                )
            fragment = source_sequence[occurrence.start : occurrence.end]
            if orient_to_bundle and occurrence.direction == "-":
                fragment = reverse_complement(fragment)
            output_chunks[occurrence.bundle_id].append(extraction_header(occurrence))
            output_chunks[occurrence.bundle_id].extend(wrap_fasta(fragment))
            extracted_counts[occurrence.bundle_id] += 1
        for bundle_id, chunks in output_chunks.items():
            with output_paths[bundle_id].open("a", encoding="utf-8") as output_handle:
                output_handle.writelines(chunks)

    missing_paths = sorted(requested_paths - seen_fasta_names)
    if missing_paths:
        example = ", ".join(missing_paths[:5])
        raise ValueError(
            f"FASTA lacks {len(missing_paths)} BED path IDs; examples: {example}"
        )
    for bundle_id, occurrences in occurrences_by_bundle.items():
        if extracted_counts[bundle_id] != len(occurrences):
            raise ValueError(
                f"bundle {bundle_id}: extracted {extracted_counts[bundle_id]} of "
                f"{len(occurrences)} requested occurrences"
            )

    bundle_inputs: List[BundleInput] = []
    for bundle_id in sorted(occurrences_by_bundle, key=bundle_sort_key):
        occurrences = occurrences_by_bundle[bundle_id]
        lengths = [occurrence.length for occurrence in occurrences]
        directions = Counter(occurrence.direction for occurrence in occurrences)
        bundle_inputs.append(
            BundleInput(
                bundle_id=bundle_id,
                fasta_path=output_paths[bundle_id],
                occurrence_count=len(occurrences),
                forward_count=directions["+"],
                reverse_count=directions["-"],
                total_bp=sum(lengths),
                min_bp=min(lengths),
                median_bp=float(statistics.median(lengths)),
                max_bp=max(lengths),
            )
        )
    return bundle_inputs


def write_manifest(path: Path, bundle_inputs: Sequence[BundleInput]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "bundle_id", "sequence_fasta", "occurrence_count", "forward_count",
                "reverse_count", "total_bp", "min_bp", "median_bp", "max_bp",
            )
        )
        for item in bundle_inputs:
            writer.writerow(
                (
                    item.bundle_id, item.fasta_path, item.occurrence_count,
                    item.forward_count, item.reverse_count, item.total_bp,
                    item.min_bp, f"{item.median_bp:.1f}", item.max_bp,
                )
            )


def resolve_executable(value: str) -> Optional[str]:
    candidate = Path(value)
    if candidate.parent != Path(".") or candidate.is_absolute():
        return str(candidate) if candidate.is_file() else None
    return shutil.which(value)


def run_mafft_alignments(
    mafft_executable: str,
    bundle_inputs: Sequence[BundleInput],
    alignment_dir: Path,
    log_dir: Path,
    threads: int,
    extra_arguments: Sequence[str],
    resume: bool,
) -> List[Tuple[str, Path, str]]:
    alignment_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    results: List[Tuple[str, Path, str]] = []
    for index, item in enumerate(bundle_inputs, 1):
        alignment_path = alignment_dir / bundle_filename(item.bundle_id).replace(
            ".fa", ".aln.fa"
        )
        log_path = log_dir / bundle_filename(item.bundle_id).replace(".fa", ".mafft.log")
        if resume and alignment_path.is_file() and alignment_path.stat().st_size > 0:
            print(
                f"[{index}/{len(bundle_inputs)}] bundle {item.bundle_id}: reuse {alignment_path}",
                file=sys.stderr,
            )
            results.append((item.bundle_id, alignment_path, "reused"))
            continue
        if item.occurrence_count == 1:
            shutil.copyfile(item.fasta_path, alignment_path)
            log_path.write_text(
                "MAFFT skipped because the bundle has one occurrence.\n",
                encoding="utf-8",
            )
            results.append((item.bundle_id, alignment_path, "copied_single_sequence"))
            continue

        if item.max_bp < SHORT_BUNDLE_MAFFT_MIN_BP:
            records = list(iter_fasta(item.fasta_path))
            aligned = align_short_bundle_records(records)
            temporary_path = alignment_path.with_name(alignment_path.name + ".tmp")
            try:
                with temporary_path.open("w", encoding="utf-8") as handle:
                    for name, sequence in aligned.items():
                        handle.write(f">{name}\n")
                        handle.writelines(wrap_fasta(sequence))
                temporary_path.replace(alignment_path)
            except Exception:
                if temporary_path.exists():
                    temporary_path.unlink()
                raise
            log_path.write_text(
                f"MAFFT skipped: longest occurrence {item.max_bp} bp is below "
                f"{SHORT_BUNDLE_MAFFT_MIN_BP} bp; used short global alignment.\n",
                encoding="utf-8",
            )
            results.append((item.bundle_id, alignment_path, "aligned_short_global"))
            continue

        temporary_path = alignment_path.with_name(alignment_path.name + ".tmp")
        command = [
            mafft_executable,
            "--auto",
            "--thread",
            str(threads),
            "--inputorder",
            *extra_arguments,
            str(item.fasta_path),
        ]
        print(
            f"[{index}/{len(bundle_inputs)}] bundle {item.bundle_id}: "
            f"{item.occurrence_count} sequences, median {item.median_bp:.0f} bp",
            file=sys.stderr,
        )
        try:
            with temporary_path.open("w", encoding="utf-8") as output_handle, log_path.open(
                "w", encoding="utf-8"
            ) as error_handle:
                subprocess.run(
                    command,
                    stdout=output_handle,
                    stderr=error_handle,
                    check=True,
                )
            temporary_path.replace(alignment_path)
        except Exception:
            if temporary_path.exists():
                temporary_path.unlink()
            raise
        results.append((item.bundle_id, alignment_path, "aligned"))
    return results


def write_alignment_summary(
    path: Path,
    bundle_inputs: Sequence[BundleInput],
    alignment_results: Sequence[Tuple[str, Path, str]],
) -> None:
    result_by_id = {
        bundle_id: (alignment_path, status)
        for bundle_id, alignment_path, status in alignment_results
    }
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("bundle_id", "occurrence_count", "alignment_fasta", "status"))
        for item in bundle_inputs:
            alignment_path, status = result_by_id[item.bundle_id]
            writer.writerow(
                (item.bundle_id, item.occurrence_count, alignment_path, status)
            )


def parse_bundle_ids(values: Optional[Sequence[str]]) -> Optional[Set[str]]:
    if not values:
        return None
    selected: Set[str] = set()
    for value in values:
        selected.update(part.strip() for part in value.split(",") if part.strip())
    return selected


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract final PGR bundle sequences and run one MAFFT alignment per bundle."
    )
    parser.add_argument("fasta", type=Path, help="population FASTA containing the BED path IDs")
    parser.add_argument("bed", type=Path, help="final native four-column PGR bundle BED")
    parser.add_argument("-o", "--output-dir", required=True, type=Path)
    parser.add_argument(
        "--bundle-id",
        action="append",
        help="optional bundle ID or comma-separated IDs; repeat to select a subset",
    )
    parser.add_argument("--mafft", default="mafft", help="MAFFT executable name or path")
    parser.add_argument("--threads", type=int, default=4, help="threads for each MAFFT process")
    parser.add_argument(
        "--mafft-extra",
        default="",
        help="additional quoted MAFFT arguments, appended after --auto/--thread/--inputorder",
    )
    parser.add_argument(
        "--keep-genomic-orientation",
        action="store_true",
        help="do not reverse-complement occurrences whose PGR direction is reverse",
    )
    parser.add_argument(
        "--extract-only",
        action="store_true",
        help="write per-bundle FASTA files and manifest without running MAFFT",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="reuse an existing non-empty alignment instead of rerunning that bundle",
    )
    parser.add_argument(
        "--largest-path-order-group-only",
        action="store_true",
        help=(
            "before extraction, retain only paths in the most frequent exact "
            "signed-bundle-order group"
        ),
    )
    parser.add_argument(
        "--supplement-missing-bundles",
        action="store_true",
        help=(
            "extract rank-1 occurrences, then add missing bundle IDs from their "
            "highest-support lower-ranked group; already-covered IDs are not added again"
        ),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.threads == 0 or args.threads < -1:
        parser.error("--threads must be -1 or a positive integer")
    mafft_executable: Optional[str] = None
    if not args.extract_only:
        mafft_executable = resolve_executable(args.mafft)
        if mafft_executable is None:
            parser.error(
                f"MAFFT executable {args.mafft!r} was not found; install MAFFT, provide "
                "--mafft /path/to/mafft, or use --extract-only"
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sequence_dir = args.output_dir / "sequences"
    alignment_dir = args.output_dir / "alignments"
    log_dir = args.output_dir / "logs"
    try:
        selected_bundle_ids = parse_bundle_ids(args.bundle_id)
        source_groups = None
        source_bundle_ids_by_group = None
        if args.supplement_missing_bundles:
            all_occurrences_by_path, _ = read_pgr_occurrences(args.bed)
            source_groups, path_to_group = group_paths_by_order(all_occurrences_by_path)
            source_bundle_ids_by_group = select_group_bundle_ids(
                source_groups, supplement_missing_bundles=True
            )
            if selected_bundle_ids is not None:
                source_bundle_ids_by_group = {
                    group_id: bundle_ids & selected_bundle_ids
                    for group_id, bundle_ids in source_bundle_ids_by_group.items()
                    if bundle_ids & selected_bundle_ids
                }
            occurrences_by_path = {}
            occurrences_by_bundle = {}
            for path_id, occurrences in all_occurrences_by_path.items():
                group_id = path_to_group[path_id].group_id
                bundle_ids = source_bundle_ids_by_group.get(group_id, set())
                selected_occurrences = [
                    item for item in occurrences if item.bundle_id in bundle_ids
                ]
                if selected_occurrences:
                    occurrences_by_path[path_id] = selected_occurrences
                    for occurrence in selected_occurrences:
                        occurrences_by_bundle.setdefault(occurrence.bundle_id, []).append(
                            occurrence
                        )
            if not occurrences_by_bundle:
                raise ValueError("no bundle occurrences remain after rank/bundle selection")
            rank_one = source_groups[0]
            write_selected_path_group(
                args.output_dir / "selected_path_order_group.tsv",
                rank_one.signature,
                rank_one.path_ids,
            )
        elif args.largest_path_order_group_only:
            all_occurrences_by_path, _ = read_pgr_occurrences(args.bed)
            selected_signature, selected_paths = largest_path_order_group(
                all_occurrences_by_path
            )
            occurrences_by_path, occurrences_by_bundle = filter_occurrences(
                all_occurrences_by_path,
                set(selected_paths),
                selected_bundle_ids,
            )
            write_selected_path_group(
                args.output_dir / "selected_path_order_group.tsv",
                selected_signature,
                selected_paths,
            )
        else:
            occurrences_by_path, occurrences_by_bundle = read_pgr_occurrences(
                args.bed,
                selected_bundle_ids,
            )
        bundle_inputs = extract_bundle_fastas(
            args.fasta,
            occurrences_by_path,
            occurrences_by_bundle,
            sequence_dir,
            orient_to_bundle=not args.keep_genomic_orientation,
        )
        if source_groups is not None and source_bundle_ids_by_group is not None:
            write_bundle_source_groups(
                args.output_dir / "bundle_source_groups.tsv",
                source_groups,
                source_bundle_ids_by_group,
            )
        write_manifest(args.output_dir / "bundle_sequences.tsv", bundle_inputs)
        if args.extract_only:
            print(
                f"Extracted {sum(item.occurrence_count for item in bundle_inputs)} occurrences "
                f"for {len(bundle_inputs)} bundles under {sequence_dir}",
                file=sys.stderr,
            )
            return 0
        alignment_results = run_mafft_alignments(
            mafft_executable,  # type: ignore[arg-type]
            bundle_inputs,
            alignment_dir,
            log_dir,
            args.threads,
            shlex.split(args.mafft_extra),
            args.resume,
        )
        write_alignment_summary(
            args.output_dir / "bundle_alignments.tsv",
            bundle_inputs,
            alignment_results,
        )
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
