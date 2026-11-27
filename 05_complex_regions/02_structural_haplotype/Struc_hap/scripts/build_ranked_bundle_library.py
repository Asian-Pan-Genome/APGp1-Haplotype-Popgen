#!/usr/bin/env python3
"""Grow a bundle library by remapping successive signed-order haplotype groups.

New sequences come from original BED occurrences minus accepted genomic coverage;
their genomic cut points are projected through a real MSA. In quick mode,
minimap2 CIGAR maps non-backbone rows to a reference for local MAFFT additions;
the final source-to-genome cuts are still taken from the resulting real MSA.
Original numeric IDs remain reserved even after retirement.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from array import array
from bisect import bisect_left
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from itertools import zip_longest
from pathlib import Path
from statistics import median
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

try:
    from .align_bundles_with_mafft import reverse_complement
    from .build_bundle_consensus import iter_fasta, wrap_fasta, bundle_sort_key
    from .build_grouped_bundle_consensus import (
        BedOccurrence, read_bundle_bed, group_paths_by_order,
        aligned_consensus_profile, aligned_profile_identity,
    )
    from .refine_bundle_mappings import (
        merge_intervals, subtract_intervals, normalize_bundle_name,
        read_existing_bundle_metadata, read_paf, main as refine_mappings_main,
    )
    from .split_bundles_by_mapping import main as split_mapping_main
    from .split_rank1_self_alignments import split_rank1_library
    from .pipeline_intermediates import cleanup_generated_files, cleanup_rank_round
    from .filter_pgr_bundles_by_trf import filter_pgr_bed_by_trf
    from .filter_duplicate_fasta_ids import filter_duplicate_inputs
    from .refine_candidate_by_trf import REPORT_COLUMNS as TRF_REPORT_COLUMNS, refine_candidate_by_trf
    from .short_bundle_alignment import SHORT_BUNDLE_MAFFT_MIN_BP, align_short_bundle_records
except ImportError:
    from align_bundles_with_mafft import reverse_complement
    from build_bundle_consensus import iter_fasta, wrap_fasta, bundle_sort_key
    from build_grouped_bundle_consensus import (
        BedOccurrence, read_bundle_bed, group_paths_by_order,
        aligned_consensus_profile, aligned_profile_identity,
    )
    from refine_bundle_mappings import (
        merge_intervals, subtract_intervals, normalize_bundle_name,
        read_existing_bundle_metadata, read_paf, main as refine_mappings_main,
    )
    from split_bundles_by_mapping import main as split_mapping_main
    from split_rank1_self_alignments import split_rank1_library
    from pipeline_intermediates import cleanup_generated_files, cleanup_rank_round
    from filter_pgr_bundles_by_trf import filter_pgr_bed_by_trf
    from filter_duplicate_fasta_ids import filter_duplicate_inputs
    from refine_candidate_by_trf import REPORT_COLUMNS as TRF_REPORT_COLUMNS, refine_candidate_by_trf
    from short_bundle_alignment import SHORT_BUNDLE_MAFFT_MIN_BP, align_short_bundle_records


VERSION = 22
RANK_CHECKPOINT = "rank_checkpoint.json.gz"


def write_rank_checkpoint(path, value):
    """Replace the last complete state only after the new gzip is closed."""
    if path.is_symlink() or (path.exists() and not path.is_file()):
        raise ValueError(f"unsafe rank checkpoint path: {path}")
    descriptor, temporary = tempfile.mkstemp(prefix=".rank_checkpoint.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as raw:
            with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=1, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8") as handle:
                    json.dump(value, handle, separators=(",", ":"))
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def load_rank_checkpoint(path, run_fingerprint):
    if not path.is_file() or path.is_symlink():
        raise ValueError(
            f"--resume-ranks requires {path}; old runs without a checkpoint cannot be "
            "restored exactly from PAF/additions alone. Keep the old output and start "
            "a new run in a separate output directory to enable rank checkpoints."
        )
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            checkpoint = json.load(handle)
    except (OSError, EOFError, ValueError) as error:
        raise ValueError(f"invalid or incomplete rank checkpoint: {path}: {error}") from error
    if not isinstance(checkpoint, dict) or checkpoint.get("checkpoint_version") != 1:
        raise ValueError(f"unsupported rank checkpoint: {path}")
    if checkpoint.get("run_fingerprint") != run_fingerprint:
        raise ValueError("--resume-ranks inputs, parameters, or pipeline version changed; "
                         "use the original settings or a separate output directory")
    rank = checkpoint.get("completed_rank")
    required = {"groups", "entries", "sizes", "next_id", "comparison_index", "lineage",
                "iterations", "cut_templates", "retired_original_ids", "round_artifacts"}
    if not required.issubset(checkpoint) or not isinstance(rank, int) or rank < 0:
        raise ValueError(f"missing state in rank checkpoint: {path}")
    kinds = {"groups": list, "entries": list, "sizes": dict, "next_id": int,
             "comparison_index": int, "lineage": list, "iterations": list,
             "cut_templates": dict, "retired_original_ids": list, "round_artifacts": list}
    if any(not isinstance(checkpoint[key], kind) for key, kind in kinds.items()):
        raise ValueError(f"invalid state types in rank checkpoint: {path}")
    if (rank > len(checkpoint["groups"]) or len(checkpoint["iterations"]) != rank
            or len(checkpoint["round_artifacts"]) != rank
            or [row[0] for row in checkpoint["iterations"]] != list(range(1, rank + 1))):
        raise ValueError(f"inconsistent rank history in checkpoint: {path}")
    return checkpoint


@dataclass(frozen=True)
class RankedGroup:
    group_id: str
    rank: int
    signature: Tuple[str, ...]
    path_ids: Tuple[str, ...]
    sample_ids: Tuple[str, ...]

    @property
    def sample_count(self):
        return len(self.sample_ids)

    @property
    def path_count(self):
        return len(self.path_ids)


@dataclass
class LibraryEntry:
    name: str
    sequence: str
    source_bundle_id: str
    group_rank: int
    members: List[str] = field(default_factory=list)
    source_reverse: bool = False


@dataclass(frozen=True)
class RankRemapFiles:
    directory: Path
    population_fasta: Path
    before_fasta: Path
    before_metadata: Path
    before_paf: Path
    trimmed_fasta: Path
    trimmed_metadata: Path


MAPPING_OPTIONS = (
    "min_bundle_coverage", "min_identity", "min_mapq", "max_bundle_gap", "max_genome_gap",
    "max_overlap_fraction", "fixed_bundle_min_coverage", "max_fixed_overlap_fraction",
    "max_fixed_overlap_bp",
)


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_sample_map(path: Optional[Path], path_ids: Sequence[str]) -> Dict[str, str]:
    if path is None:
        return {name: name for name in path_ids}
    mapping = {}
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields or line.lstrip().startswith("#"):
                continue
            if len(fields) != 2:
                raise ValueError(f"{path}:{number}: expected path_id sample_id")
            if fields == ["path_id", "sample_id"]:
                continue
            name, sample = fields
            if name in mapping and mapping[name] != sample:
                raise ValueError(f"{path}:{number}: conflicting sample for {name}")
            mapping[name] = sample
    missing = sorted(set(path_ids) - mapping.keys())
    if missing:
        raise ValueError(f"sample-map lacks {len(missing)} BED paths: {', '.join(missing[:5])}")
    return mapping


def rank_groups(by_path, sample_map: Optional[Mapping[str, str]] = None) -> List[RankedGroup]:
    groups, _ = group_paths_by_order(by_path)
    if not sample_map:
        sample_map = {name: name for name in by_path}
    ordered = sorted(groups, key=lambda group: (
        -len({sample_map[name] for name in group.path_ids}),
        -len(group.path_ids), group.signature,
    ))
    return [RankedGroup(f"G{rank:04d}", rank, group.signature, group.path_ids,
                        tuple(sorted({sample_map[name] for name in group.path_ids})))
            for rank, group in enumerate(ordered, 1)]


def source_residuals(occurrence: BedOccurrence, covered) -> List[Tuple[int, int]]:
    """Subtract accepted genomic spans, in exact zero-based source coordinates."""
    return list(subtract_intervals([(occurrence.start, occurrence.end)], merge_intervals(covered)))


def write_tsv(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def write_json(path, value):
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_fasta(path, records):
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            handle.writelines(wrap_fasta(sequence, 80))


def write_metadata(path, sizes):
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Reserved ID/size metadata only; NOT genomic mappings\n")
        for bundle_id in sorted(sizes, key=int):
            size = sizes[bundle_id]
            handle.write(f"__metadata__\t0\t1\t{bundle_id}:{size}:0:0:{size - 1}:U\n")


def write_ranking(directory, groups, fingerprint):
    directory.mkdir(parents=True, exist_ok=True)
    write_tsv(directory / "haplotype_groups.tsv",
              ("rank", "group_id", "sample_count", "path_count", "path_ids",
               "sample_ids", "signed_bundle_order"),
              ((g.rank, g.group_id, g.sample_count, g.path_count, ",".join(g.path_ids),
                ",".join(g.sample_ids), ",".join(g.signature)) for g in groups))
    write_json(directory / "ranking_summary.json", {
        "version": VERSION, "input_fingerprint": fingerprint, "rank_count": len(groups),
        "ranking_unit": "distinct samples within each signed-order group",
        "ties": "descending path count, then signed order",
    })


def validate_ranked_library(fasta: Path, summary: Path) -> None:
    data = json.loads(summary.read_text())
    if (data.get("version") != VERSION or not data.get("completed")
            or data.get("ranks_processed") != data.get("rank_count")):
        raise ValueError("incomplete ranked library; rerun the iterative library step")
    mapping_contract = data.get("run_fingerprint", {}).get("mapping_contract", {})
    if mapping_contract != {
            "rank_minimap_gap_bp": 500,
            "self_alignment_minimap_gap_bp": 500,
            "rank1_genome_mapping": "split_only",
    }:
        raise ValueError("ranked library uses an older mapping contract; rerun the iterative library step")
    if file_digest(fasta) != data.get("fasta_sha256"):
        raise ValueError("ranked library FASTA changed since completion; rerun the iterative library step")
    records = list(iter_fasta(fasta))
    names = [header.split()[0] for header, sequence in records]
    if (not records or any(not sequence for _, sequence in records)
            or len(names) != len(set(names)) or sorted(names) != sorted(data["active_consensus_names"])):
        raise ValueError("ranked library records do not match the completed audit")
    metadata = summary.parent / "bundle_metadata.bed"
    if file_digest(metadata) != data.get("metadata_sha256"):
        raise ValueError("ranked library metadata changed since completion; rerun the iterative library step")


def write_round_cleanup(args, directory, completed=False, self_result=None):
    """Publish deletion evidence immediately, including for interrupted rounds."""
    reports = args._round_cleanup
    audit_path = directory / "intermediate_cleanup.json"
    if audit_path.is_symlink() or (audit_path.exists() and not audit_path.is_file()):
        raise ValueError(f"unsafe intermediate cleanup audit path: {audit_path}")
    write_json(audit_path, {
        "completed": completed, "mode": "keep_all" if args.keep_intermediates else "raw_paf_and_additions_only",
        "deleted": [item for report in reports for item in report["deleted"]],
        "skipped": [item for report in reports for item in report["skipped"]],
        "deleted_count": sum(report["deleted_count"] for report in reports),
        "deleted_bytes": sum(report["deleted_bytes"] for report in reports),
        "self_alignment_cleanup": "self_alignment/intermediate_cleanup.json" if self_result else None,
        "policy": ("retain all intermediates" if args.keep_intermediates else
                   "retain raw PAFs, additions.tsv, round manifest, and cleanup audit only"),
    })


def run_mafft(records, directory, stem, args):
    source = directory / f"{stem}.fa"
    alignment = directory / f"{stem}.aln.fa"
    log = directory / f"{stem}.mafft.log"
    generated_paths = [path for path in (source, alignment)
                       if not path.exists() and not path.is_symlink()]
    write_fasta(source, records)
    short_bundle = (stem.startswith("source_") or stem.startswith("compare_")) \
        and max(len(sequence) for _, sequence in records) < SHORT_BUNDLE_MAFFT_MIN_BP
    if short_bundle:
        write_fasta(alignment, align_short_bundle_records(records).items())
    elif len(records) == 1:
        write_fasta(alignment, records)
    else:
        command = [args.mafft, *shlex.split(getattr(args, "mafft_params", "--auto")),
                   "--thread", str(args.threads), "--inputorder"]
        command += shlex.split(args.mafft_extra) + [str(source)]
        with alignment.open("w") as output, log.open("w") as error:
            subprocess.run(command, stdout=output, stderr=error, check=True)
    aligned = {header.split()[0]: sequence for header, sequence in iter_fasta(alignment)}
    if (len(aligned) != len(records) or set(aligned) != {name for name, _ in records}
            or len({len(sequence) for sequence in aligned.values()}) != 1):
        raise ValueError(f"invalid MAFFT records/width: {alignment}")
    for name, sequence in records:
        if aligned[name].replace("-", "") != sequence:
            raise ValueError(f"MAFFT changed source residues: {alignment}: {name}")
    # No caller rereads these files: validated alignment profiles now live in
    # memory. Failed/invalid alignments are deliberately left for diagnosis.
    cleanup = cleanup_generated_files(
        directory, generated_paths,
        protected_paths=getattr(args, "_cleanup_protected", ()),
        keep_intermediates=getattr(args, "keep_intermediates", True),
    )
    if hasattr(args, "_round_cleanup"):
        args._round_cleanup.append(cleanup)
        write_round_cleanup(args, directory)
    return aligned


def _quick_paf_hits(path, reference_bp, queries):
    """Choose one forward, CIGAR-backed placement per source query."""
    hits = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) < 12 or fields[0] not in queries or fields[4] != "+" or fields[5] != "reference":
            continue
        try:
            qlen, qstart, qend = map(int, fields[1:4])
            tlen, tstart, tend = map(int, fields[6:9])
            matches = int(fields[9])
        except ValueError:
            continue
        if (qlen != len(queries[fields[0]]) or tlen != reference_bp
                or not (0 <= qstart < qend <= qlen and 0 <= tstart < tend <= tlen)):
            continue
        tag = next((field[5:] for field in fields[12:] if field.startswith("cg:Z:")), None)
        if not tag:
            continue
        operations = [(int(length), op) for length, op in re.findall(r"(\d+)([MID=XN])", tag)]
        if (not operations or "".join(f"{length}{op}" for length, op in operations) != tag
                or sum(length for length, op in operations if op in "MI=X") != qend - qstart
                or sum(length for length, op in operations if op in "MD=XN") != tend - tstart):
            continue
        hit = (qstart, qend, tstart, tend, operations)
        # Prefer the placement covering most of the query; otherwise a short
        # high-identity local hit could hide a usable >=95% placement.
        score = (qend - qstart, matches, tend - tstart)
        if fields[0] not in hits or score > hits[fields[0]][0]:
            hits[fields[0]] = (score, hit)
    return {name: item[1] for name, item in hits.items()}


def _quick_query_boundary(hit, ref_boundary, ref_bp, query_bp):
    """Map a reference cut to a query cut; preserve terminal unmapped tails."""
    qstart, qend, tstart, tend, operations = hit
    if ref_boundary == 0:
        return 0
    if ref_boundary == ref_bp:
        return query_bp
    if ref_boundary < tstart:
        return 0
    if ref_boundary == tstart:
        return qstart
    if ref_boundary == tend:
        return qend
    if ref_boundary > tend:
        return query_bp
    qpos, tpos = qstart, tstart
    for length, op in operations:
        if op in "M=X":
            if tpos <= ref_boundary <= tpos + length:
                return qpos + ref_boundary - tpos
            qpos += length
            tpos += length
        elif op in "DN":
            if tpos <= ref_boundary <= tpos + length:
                return qpos
            tpos += length
        else:  # I: an insertion at a cut belongs to the right-hand piece.
            qpos += length
    raise ValueError("reference boundary is outside the validated CIGAR")


def _quick_add_to_backbone(backbone, additions, directory, stem, args, option):
    """Use MAFFT --add/--addfragments without changing existing MSA columns."""
    base_path = directory / f"{stem}.backbone.fa"
    addition_path = directory / f"{stem}.additions.fa"
    result_path = directory / f"{stem}.aln.fa"
    log_path = directory / f"{stem}.mafft.log"
    write_fasta(base_path, backbone.items())
    write_fasta(addition_path, additions.items())
    command = [args.mafft, *shlex.split(getattr(args, "mafft_params", "--auto")),
               option, str(addition_path), "--thread", str(args.threads), "--inputorder"]
    command += shlex.split(args.mafft_extra) + [str(base_path)]
    with result_path.open("w") as output, log_path.open("w") as error:
        subprocess.run(command, stdout=output, stderr=error, check=True)
    aligned = {header.split()[0]: sequence for header, sequence in iter_fasta(result_path)}
    if (set(aligned) != set(backbone) | set(additions)
            or len({len(sequence) for sequence in aligned.values()}) != 1
            or any(aligned[name].replace("-", "") != sequence.replace("-", "")
                   for name, sequence in (*backbone.items(), *additions.items()))):
        raise ValueError(f"invalid MAFFT {option} result: {result_path}")
    # MAFFT preserves the reference alignment but may remove all-gap columns.
    original_columns = (column for column in zip(*backbone.values()) if any(base != "-" for base in column))
    new_columns = (column for column in zip(*(aligned[name] for name in backbone))
                   if any(base != "-" for base in column))
    if any(old != new for old, new in zip_longest(original_columns, new_columns)):
        raise ValueError(f"MAFFT changed backbone columns: {result_path}")
    return aligned


def _quick_add_fragments(backbone, fragments, directory, stem, args):
    return _quick_add_to_backbone(backbone, fragments, directory, stem, args, "--addfragments")


def _quick_add_sequences(backbone, additions, directory, stem, args):
    return _quick_add_to_backbone(backbone, additions, directory, stem, args, "--add")


def _quick_drop_short_truncated(query_bp, reference_bp):
    """Exclude whole source records <=10% of the reference AND <1000 bp."""
    return query_bp * 10 <= reference_bp and query_bp < 1000


def run_quick_source_mafft(records, backbone_names, directory, stem, args, modes=None,
                           placement_audit=None):
    """Map all source rows, align every reference window, then stitch the MSA."""
    sequences = dict(records)
    reference_name = max((name for name in backbone_names if name.startswith("row_")),
                         key=lambda name: (len(sequences[name]), name))
    reference = sequences[reference_name]
    # Apply the gate to whole occurrences/templates before either local-window
    # additions or full-backbone --add. Do not filter individual window slices.
    dropped = {name for name, sequence in records
               if _quick_drop_short_truncated(len(sequence), len(reference))}
    retained_records = [(name, sequence) for name, sequence in records if name not in dropped]
    if max(len(sequence) for _, sequence in retained_records) < SHORT_BUNDLE_MAFFT_MIN_BP:
        aligned = run_mafft(retained_records, directory, stem, args)
        if modes is not None:
            modes[stem] = "short_global_dropped_rows" if dropped else "short_global_no_mafft"
        if placement_audit is not None:
            placement_audit.extend((stem.removeprefix("source_"), name, len(sequence),
                                    "", "", "dropped_short_truncated" if name in dropped else
                                    "short_global_no_mafft") for name, sequence in records)
        return aligned
    window_bp = getattr(args, "mafft_window_bp", 10_000)
    if window_bp < 1:
        raise ValueError("MAFFT source window size must be positive")
    work = directory / "remap" / f"quick_{stem}"
    work.mkdir(parents=True, exist_ok=True)
    query_names = [name for name in sequences if name != reference_name and name not in dropped]
    hits = {reference_name: (0, len(reference), 0, len(reference), [(len(reference), "M")])}
    if query_names:
        reference_path, query_path = work / "reference.fa", work / "queries.fa"
        paf_path, log_path = work / "source_to_reference.raw.paf", work / "minimap.log"
        write_fasta(reference_path, [("reference", reference)])
        write_fasta(query_path, [(name, sequences[name]) for name in query_names])
        command = [args.minimap2, "-x", "asm20", "-g", "500", "--secondary=no", "-c",
                   "-t", str(args.threads), str(reference_path), str(query_path)]
        with paf_path.open("w") as output, log_path.open("w") as error:
            subprocess.run(command, stdout=output, stderr=error, check=True)
        hits.update(_quick_paf_hits(paf_path, len(reference),
                                    {name: sequences[name] for name in query_names}))
    windowed_records, partial_records = [], []
    for name, sequence in records:
        hit = hits.get(name)
        aligned_bp = hit[1] - hit[0] if hit else 0
        fraction = aligned_bp / len(sequence)
        if name in dropped:
            method = "dropped_short_truncated"
        elif hit is None:
            method = "dropped_no_alignment"
        elif fraction < .95:
            partial_records.append((name, sequence))
            method = "mafft_add"
        else:
            windowed_records.append((name, sequence))
            method = "windowed_backbone" if name in backbone_names else "windowed_local"
        if placement_audit is not None:
            placement_audit.append((stem.removeprefix("source_"), name, len(sequence),
                                    "" if name in dropped else aligned_bp,
                                    "" if name in dropped else f"{fraction:.10g}", method))
    cuts = {0, len(reference)}
    cuts.update(range(window_bp, len(reference), window_bp))
    cuts.update(point for name, _ in windowed_records for point in hits[name][2:4])
    cuts = sorted(cuts)
    pieces = {name: [] for name, _ in windowed_records}
    for index, (left, right) in enumerate(zip(cuts, cuts[1:])):
        fragments = {}
        for name, sequence in windowed_records:
            qleft = _quick_query_boundary(hits[name], left, len(reference), len(sequence))
            qright = _quick_query_boundary(hits[name], right, len(reference), len(sequence))
            if qright < qleft:
                raise ValueError(f"non-monotone reference-window projection: {stem}: {name}")
            if qright > qleft:
                fragments[name] = sequence[qleft:qright]
        backbone_rows = [(name, fragments[name]) for name, _ in windowed_records
                         if name in backbone_names and name in fragments]
        # The reference spans every window, so an empty backbone is a bug.
        if not backbone_rows:
            raise ValueError(f"source window lacks reference sequence: {stem}: {left}-{right}")
        aligned = run_mafft(backbone_rows, work, f"window_{index:04d}", args)
        additions = {name: fragments[name] for name, _ in windowed_records
                     if name not in backbone_names and name in fragments}
        if additions:
            aligned = _quick_add_fragments(aligned, additions, work, f"part_{index:04d}", args)
        width = len(next(iter(aligned.values())))
        for name in pieces:
            pieces[name].append(aligned.get(name, "-" * width))
    joined = {name: "".join(parts) for name, parts in pieces.items()}
    if partial_records:
        joined = _quick_add_sequences(joined, dict(partial_records), work, "global_add", args)
    kept_records = [(name, sequence) for name, sequence in records if name in joined]
    joined = {name: joined[name] for name, _ in kept_records}
    if (len({len(row) for row in joined.values()}) != 1
            or any(joined[name].replace("-", "") != sequence
                   for name, sequence in kept_records)):
        raise ValueError(f"stitched source-window MSA changed a sequence: {stem}")
    source, alignment = directory / f"{stem}.fa", directory / f"{stem}.aln.fa"
    generated = [path for path in (source, alignment) if not path.exists() and not path.is_symlink()]
    write_fasta(source, kept_records)
    write_fasta(alignment, joined.items())
    cleanup = cleanup_generated_files(
        directory, generated, protected_paths=getattr(args, "_cleanup_protected", ()),
        keep_intermediates=getattr(args, "keep_intermediates", True))
    if hasattr(args, "_round_cleanup"):
        args._round_cleanup.append(cleanup)
        write_round_cleanup(args, directory)
    if modes is not None:
        modes[stem] = ("windowed_with_mafft_add" if partial_records else
                       "windowed_dropped_rows" if len(kept_records) != len(records) else
                       "windowed_all_rows")
    return joined


def uncovered_column_ranges(occurrence, aligned, covered, columns=None):
    """Project exact BED residual boundaries through this occurrence's real MSA."""
    if columns is None:
        columns = array("I", (column for column, base in enumerate(aligned) if base != "-"))
    residuals = source_residuals(occurrence, covered)
    result = []
    for start, end in residuals:
        if occurrence.direction == "+":
            left, right = start - occurrence.start, end - occurrence.start
        else:
            left, right = occurrence.end - end, occurrence.end - start
        result.append((columns[left], columns[right - 1] + 1, start, end))
    # For whole occurrences include terminal MSA gaps to preserve common profiles.
    if residuals == [(occurrence.start, occurrence.end)]:
        result = [(0, len(aligned), occurrence.start, occurrence.end)]
    return result


def select_source_msa_rows(group, by_path, needed_ids, quick):
    """Classify every source FASTA row using 70% of its longest occurrence."""
    selected_by_id = defaultdict(set)
    audit = []
    occurrences_by_id = defaultdict(list)
    for position, token in enumerate(group.signature):
        bundle_id = token[1:]
        if bundle_id not in needed_ids:
            continue
        occurrences_by_id[bundle_id].extend(
            (position + 1, by_path[path][position]) for path in group.path_ids)
    for bundle_id, rows in occurrences_by_id.items():
        reference_bp = max(occ.end - occ.start for _, occ in rows)
        selective = quick
        lower, upper = reference_bp * .70, reference_bp * 1.30
        in_range = {occ.row_number for _, occ in rows
                    if not selective or lower <= occ.end - occ.start <= upper}
        for position, occ in rows:
            selected = occ.row_number in in_range
            if selected:
                selected_by_id[bundle_id].add(occ.row_number)
            if quick:
                audit.append((bundle_id, position, occ.path_id, occ.row_number,
                              occ.end - occ.start, reference_bp, lower if selective else "",
                              upper if selective else "", int(selected),
                              ("dropped_short_truncated"
                               if _quick_drop_short_truncated(occ.end - occ.start, reference_bp)
                               else "backbone" if selected else "local_addition")))
    return selected_by_id, audit


def position_candidates(group, by_path, genomes, coverage, directory, args, cut_templates=None,
                        trf_context=None, trf_audit=None):
    """Yield homologous residual pieces per position in signed bundle order.

    All copies of a source ID share one group-specific MSA. Partition at every
    projected coverage boundary. A record contributes only bases outside old
    coverage; records with no residual bases in a partition do not contribute.
    Reserving accepted source pieces prevents overlapping original BED rows from
    adding the same genomic region twice in this round.
    """
    alignments, forced_columns, skipped = {}, {}, []
    cut_templates = cut_templates or {}
    rows_by_id = defaultdict(list)
    for path in group.path_ids:
        for occurrence in by_path[path]:
            rows_by_id[occurrence.bundle_id].append(occurrence)
    needed_ids = {occ.bundle_id for path in group.path_ids for occ in by_path[path]
                  if source_residuals(occ, coverage.get(path, []))}
    selected_by_id, quick_audit = select_source_msa_rows(
        group, by_path, needed_ids, getattr(args, "mafft_quick", False))
    quick_modes = {}
    quick_placement_audit = []
    for bundle_id in sorted(needed_ids, key=bundle_sort_key):
        records = []
        for occurrence in rows_by_id[bundle_id]:
            sequence = genomes[occurrence.path_id][occurrence.start:occurrence.end]
            if occurrence.direction == "-":
                sequence = reverse_complement(sequence)
            records.append((f"row_{occurrence.row_number}", sequence))
        for index, template in enumerate(cut_templates.get(bundle_id, [])):
            records.append((f"cut_template_{index}", template["sequence"]))
        backbone_names = {f"row_{row}" for row in selected_by_id[bundle_id]}
        reference_bp = max(len(sequence) for name, sequence in records
                           if name.startswith("row_"))
        backbone_names.update(name for name, sequence in records
                              if name.startswith("cut_template_")
                              and len(sequence) >= reference_bp * .70)
        alignments[bundle_id] = (
            run_quick_source_mafft(records, backbone_names, directory, f"source_{bundle_id}",
                                   args, quick_modes, quick_placement_audit)
            if getattr(args, "mafft_quick", False) else
            run_mafft(records, directory, f"source_{bundle_id}", args))
        forced_columns[bundle_id] = set()
        for index, template in enumerate(cut_templates.get(bundle_id, [])):
            aligned = alignments[bundle_id].get(f"cut_template_{index}")
            if aligned is None:
                continue
            columns = [column for column, base in enumerate(aligned) if base != "-"]
            for point in template["boundaries"]:
                forced_columns[bundle_id].add(columns[point] if point < len(columns) else columns[-1] + 1)
    if getattr(args, "mafft_quick", False):
        write_tsv(directory / "mafft_quick_placement.tsv",
                  ("bundle_id", "record_name", "query_bp", "aligned_bp",
                   "query_coverage", "method"), quick_placement_audit)
        write_tsv(directory / "mafft_quick_selection.tsv",
                  ("bundle_id", "position_1based", "path_id", "bed_row", "occurrence_bp",
                   "reference_bp", "lower_bp", "upper_bp", "backbone", "reason",
                   "alignment_mode", "included_final"),
                  (row + (quick_modes[f"source_{row[0]}"],
                          int(f"row_{row[3]}" in alignments[row[0]])) for row in quick_audit))
    for position, token in enumerate(group.signature):
        bundle_id = token[1:]
        if bundle_id not in alignments:
            continue
        records = []
        boundaries = set(forced_columns[bundle_id])
        width = len(next(iter(alignments[bundle_id].values())))
        for path in group.path_ids:
            occurrence = by_path[path][position]
            if f"row_{occurrence.row_number}" not in alignments[bundle_id]:
                continue
            aligned = alignments[bundle_id][f"row_{occurrence.row_number}"]
            columns = array("I", (column for column, base in enumerate(aligned) if base != "-"))
            ranges = uncovered_column_ranges(occurrence, aligned, coverage.get(path, []), columns)
            records.append((occurrence, aligned, ranges, columns))
            for start, end, _, _ in ranges:
                boundaries.update((start, end))
        ordered = sorted(boundaries)
        for left, right in zip(ordered, ordered[1:]):
            selected, sources = [], []
            for occurrence, aligned, ranges, columns in records:
                if not any(start <= left and right <= end for start, end, _, _ in ranges):
                    continue
                segment = aligned[left:right]
                base_start = bisect_left(columns, left)
                base_end = bisect_left(columns, right)
                length = base_end - base_start
                if length < args.addition_min_bp:
                    continue
                if occurrence.direction == "+":
                    start, end = occurrence.start + base_start, occurrence.start + base_end
                else:
                    start, end = occurrence.end - base_end, occurrence.end - base_start
                selected.append(segment)
                sources.append((occurrence.path_id, start, end, occurrence.direction,
                                occurrence.row_number))
            if not selected:
                continue
            try:
                profile, consensus, _, _, _ = aligned_consensus_profile(
                    selected, args.min_column_occupancy, args.min_base_frequency)
            except ValueError as error:
                if str(error) != "no MSA column passed the candidate consensus threshold":
                    raise
                skipped.append((position + 1, bundle_id, left, right, "no_supported_canonical_column"))
                continue
            if len(consensus) < args.addition_min_bp:
                continue
            whole = all(start == by_path[path][position].start
                        and end == by_path[path][position].end
                        for path, start, end, _, _ in sources)
            candidate = {"source_bundle_id": bundle_id, "position": position + 1,
                   "msa_start": left, "msa_end": right, "msa_width": width,
                   "profile": "-" * left + profile + "-" * (width - right),
                   "aligned_sequences": ["-" * left + row + "-" * (width - right) for row in selected],
                   "sequence": consensus, "sources": sources, "whole": whole}
            refined = [candidate]
            if trf_audit is not None:
                trf_audit["before"].append(candidate)
            if trf_context and trf_context["summary"]["enabled"]:
                refined, report = refine_candidate_by_trf(candidate, trf_context["records"], args)
                if trf_audit is not None:
                    trf_audit["reports"].extend(report)
            for fragment in refined:
                if trf_audit is not None:
                    trf_audit["after"].append(fragment)
                yield fragment
                for path, start, end, _, _ in fragment["sources"]:
                    coverage[path] = merge_intervals([*coverage.get(path, []), (start, end)])
    write_tsv(directory / "skipped_candidates.tsv",
              ("position_1based", "source_bundle_id", "msa_start", "msa_end", "reason"), skipped)


def map_library(bundle_fasta, population_fasta, paf, expected_paths, genomes, args):
    """Map one specified population subset; never silently use the full FASTA.

    Always emit CIGAR for retained PAF evidence. These decisions still use only
    outer alignment coordinates, not CIGAR operations.
    """
    command = [args.minimap2, "-x", "asm20", "-p", "0.2", "-N", "10", "-g", "500", "-c",
               "-t", str(args.threads), str(bundle_fasta), str(population_fasta)]
    with paf.open("w") as output, paf.with_suffix(".minimap.log").open("w") as log:
        subprocess.run(command, stdout=output, stderr=log, check=True)
    records = {header.split()[0]: len(sequence) for header, sequence in iter_fasta(bundle_fasta)}
    expected_paths = set(expected_paths)
    for block in read_paf(paf):
        if block.target_name not in expected_paths:
            raise ValueError(f"PAF path outside requested rank subset: {block.target_name}")
        if block.target_length != len(genomes[block.target_name]):
            raise ValueError(f"PAF population length mismatch: {block.target_name}")
        if records.get(block.raw_query_name) != block.query_length:
            raise ValueError(f"PAF/active library name or length mismatch: {block.raw_query_name}")


def refine_rank_paf(paf, metadata, directory, args, log_path=None):
    """Retain filtered BED/PGR BED and tables without another mapping pass."""
    directory.mkdir(parents=True, exist_ok=True)
    options = [str(paf), "--bundle-metadata-bed", str(metadata), "--output-dir", str(directory)]
    for name in MAPPING_OPTIONS:
        options += ["--" + name.replace("_", "-"), str(getattr(args, name))]
    with (log_path or directory / "postprocess.log").open("w") as log, contextlib.redirect_stdout(log):
        refine_mappings_main(options)


def map_old_library_to_rank(entries, sizes, group, genomes, directory, args):
    """Stage 1: the library inherited from the preceding rank -> current rank."""
    remap_dir = directory / "remap"
    remap_dir.mkdir(parents=True, exist_ok=True)
    files = RankRemapFiles(
        remap_dir, directory / "rank_population.fa", directory / "library_before.fa",
        directory / "library_before.metadata.bed", remap_dir / "bundle_to_genome.raw.paf",
        remap_dir / "final_bundle.fa", remap_dir / "final_bundle_metadata.bed",
    )
    write_fasta(files.before_fasta, ((entry.name, entry.sequence) for entry in entries.values()))
    write_metadata(files.before_metadata, sizes)
    new_population = not files.population_fasta.exists() and not files.population_fasta.is_symlink()
    write_fasta(files.population_fasta, ((path, genomes[path]) for path in group.path_ids))
    if new_population and hasattr(args, "_generated_rank_fastas"):
        args._generated_rank_fastas.append(files.population_fasta)
    map_library(files.before_fasta, files.population_fasta, files.before_paf,
                group.path_ids, genomes, args)
    # This current-rank BED is an audit only. Splitting still builds its own
    # preliminary BED from current + historical evidence, as before.
    refine_rank_paf(files.before_paf, files.before_metadata, remap_dir / "before_trim", args)
    return files


def split_old_library_before_additions(entries, sizes, group, previous_groups, genomes, files, args):
    """Stage 2: split old sequences; historical mapping is evidence-only.

    Previous groups are remapped separately to the SAME pre-split library to
    keep population support and raw consensus coordinates valid. Neither their
    PAF nor their BED is used for addition discovery on the current rank.
    """
    if args.no_mapping_split:
        shutil.copyfile(files.before_fasta, files.trimmed_fasta)
        shutil.copyfile(files.before_metadata, files.trimmed_metadata)
        return dict(entries), []
    evidence_dir = files.directory / "split_evidence"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    inputs = [files.before_paf]
    prior_paths = sorted({path for previous in previous_groups for path in previous.path_ids})
    if prior_paths:
        prior_fasta, prior_paf = evidence_dir / "previous_ranks.fa", evidence_dir / "previous_ranks.raw.paf"
        new_historical = not prior_fasta.exists() and not prior_fasta.is_symlink()
        write_fasta(prior_fasta, ((path, genomes[path]) for path in prior_paths))
        if new_historical and hasattr(args, "_generated_rank_fastas"):
            args._generated_rank_fastas.append(prior_fasta)
        map_library(files.before_fasta, prior_fasta, prior_paf, prior_paths, genomes, args)
        inputs.append(prior_paf)
    evidence_paf = evidence_dir / "population_split_evidence.paf"
    with evidence_paf.open("w") as output:
        for paf in inputs:
            with paf.open() as source:
                for line in source:
                    output.write(line if line.endswith("\n") else line + "\n")
    split_dir = files.directory / "mapping_split"
    split_args = [str(files.before_fasta), str(evidence_paf), "--bundle-metadata-bed",
                  str(files.before_metadata), "--output-dir", str(split_dir)]
    renamed = {
        "min_truncated_count": "mapping_split_min_truncated_count",
        "truncated_fraction": "mapping_split_truncated_fraction",
        "max_unmapped_bp": "mapping_split_max_unmapped_bp",
        "cut_merge_bp": "mapping_split_cut_merge_bp",
        "min_cut_support": "mapping_split_min_cut_support",
        "gap_min_mapping_bp": "mapping_split_gap_min_mapping_bp",
        "gap_max_overlap_bp": "mapping_split_gap_max_overlap_bp",
        "min_terminal_fragment_bp": "mapping_split_min_terminal_bp",
        "preliminary_min_bundle_coverage": "min_bundle_coverage",
    }
    for name, attribute in renamed.items():
        split_args += ["--" + name.replace("_", "-"), str(getattr(args, attribute))]
    for name in MAPPING_OPTIONS:
        if name != "min_bundle_coverage":
            split_args += ["--" + name.replace("_", "-"), str(getattr(args, name))]
    with (files.directory / "split.log").open("w") as log, contextlib.redirect_stdout(log):
        split_mapping_main(split_args)
    shutil.copyfile(split_dir / "bundle.fa", files.trimmed_fasta)
    shutil.copyfile(split_dir / "bundle_metadata.bed", files.trimmed_metadata)
    with (split_dir / "bundle_fragments.tsv").open() as handle:
        fragments = list(csv.DictReader(handle, delimiter="\t"))
    sequences = {header.split()[0]: sequence for header, sequence in iter_fasta(files.trimmed_fasta)}
    updated, splits = dict(entries), []
    for row in fragments:
        if row["is_split"] != "1":
            continue
        parent = entries[row["parent_consensus"]]
        updated.pop(parent.name, None)
        name = row["consensus_name"]
        updated[name] = LibraryEntry(name, sequences[name], parent.source_bundle_id,
                                     parent.group_rank, [sequences[name]], parent.source_reverse)
        splits.append(row)
    _, child_metadata = read_existing_bundle_metadata(files.trimmed_metadata)
    sizes.update({key: value.bundle_size for key, value in child_metadata.items()})
    return updated, splits


def register_mapping_splits(parents, splits, entries, sizes, group_rank, lineage,
                            cut_templates, retired_original_ids, next_id):
    """Publish mapping-split lineage and inherited source-MSA boundaries."""
    next_id = max(next_id, max(map(int, sizes)) + 1)
    for row in splits:
        parent = parents[row["parent_consensus"]]
        lineage.append(("mapping_split", group_rank, row["consensus_name"],
                        row["parent_consensus"], row["parent_start"], row["parent_end"],
                        parent.source_bundle_id, "", "", "", "", "", "", "", ""))
    for parent_name in sorted({row["parent_consensus"] for row in splits}):
        parent = parents[parent_name]
        boundaries = sorted(
            {int(row["parent_start"]) for row in splits if row["parent_consensus"] == parent_name}
            | {int(row["parent_end"]) for row in splits if row["parent_consensus"] == parent_name}
        )
        template = source_cut_template(parent, boundaries)
        if template not in cut_templates[parent.source_bundle_id]:
            cut_templates[parent.source_bundle_id].append(template)
        if normalize_bundle_name(parent.name, {}) == parent.source_bundle_id:
            retired_original_ids.add(parent.source_bundle_id)
    return entries, next_id


def split_rank1_library_by_population(entries, sizes, group, genomes, directory, args):
    """Map the initialized rank-1 library once and use it only for splitting.

    There is intentionally no post-split rank-1 remap or residual-addition pass:
    rank 1 already supplied every source candidate used to initialize the library.
    """
    files = map_old_library_to_rank(entries, sizes, group, genomes, directory, args)
    updated, splits = split_old_library_before_additions(
        entries, sizes, group, (), genomes, files, args)
    write_json(files.directory / "rank_mapping_stages.json", {
        "rank": group.rank,
        "current_rank_paths": list(group.path_ids),
        "stages": ["initialize_rank1_library", "map_initialized_library_to_rank1",
                   "split_initialized_library", "stop_without_addition_discovery"],
        "mapping_split_disabled": False,
        "split_only": True,
        "before_paf": str(files.before_paf.resolve()),
        "post_trim_paf": None,
        "post_trim_library": str(files.trimmed_fasta.resolve()),
        "post_trim_names": list(updated),
        "addition_coverage_bed": None,
        "current_rank_uncovered_bed": None,
        "historical_population_mapping_is_evidence_only": True,
    })
    return updated, splits, files


def remap_trimmed_library_to_rank(entries, group, genomes, files, args):
    """Stage 3: map the post-trim library back to CURRENT rank before discovery."""
    if args.no_mapping_split:
        effective_paf = files.before_paf
    else:
        # An explicit second current-rank pass, even when no supported cuts
        # remain. Additions can never accidentally use the pre-trim PAF/BED.
        effective_paf = files.directory / "bundle_to_genome.split.raw.paf"
        map_library(files.trimmed_fasta, files.population_fasta, effective_paf,
                    group.path_ids, genomes, args)
    refine_rank_paf(effective_paf, files.trimmed_metadata, files.directory, args,
                    files.directory / "post_trim_refine.log")
    coverage = defaultdict(list)
    bed = files.directory / "refined_bundle_mappings.bed"
    with bed.open() as handle:
        for line in handle:
            fields = line.split()
            if fields and not line.startswith("#"):
                if fields[0] not in group.path_ids:
                    raise ValueError("post-trim BED contains paths outside current rank")
                coverage[fields[0]].append((int(fields[1]), int(fields[2])))
    coverage = {path: merge_intervals(intervals) for path, intervals in coverage.items()}
    uncovered = files.directory / "current_rank_uncovered.bed"
    with uncovered.open("w") as output:
        for path in group.path_ids:
            for start, end in subtract_intervals([(0, len(genomes[path]))], coverage.get(path, ())):
                output.write(f"{path}\t{start}\t{end}\n")
    write_json(files.directory / "rank_mapping_stages.json", {
        "rank": group.rank, "current_rank_paths": list(group.path_ids),
        "stages": ["map_old_library_to_current_rank", "split_old_library_before_additions",
                   "remap_trimmed_library_to_current_rank", "discover_additions_from_post_trim_bed"],
        "mapping_split_disabled": args.no_mapping_split,
        "before_paf": str(files.before_paf.resolve()), "post_trim_paf": str(effective_paf.resolve()),
        "post_trim_library": str(files.trimmed_fasta.resolve()), "post_trim_names": list(entries),
        "addition_coverage_bed": str(bed.resolve()), "current_rank_uncovered_bed": str(uncovered.resolve()),
        "historical_population_mapping_is_evidence_only": True,
    })
    return coverage


ROUND_INDEX_COLUMNS = (
    "rank", "group_id", "round_type", "split_status", "manifest", "source_pgr_bed",
    "before_trim_paf", "before_trim_pgr_bed", "after_trim_paf", "after_trim_pgr_bed",
    "split_stats", "mapping_observations", "cut_support", "bundle_fragments",
    "potential_endpoints", "library_after_fasta", "additions", "mafft_quick_selection",
    "mafft_quick_placement",
    "rank1_self_summary", "rank1_self_events", "rank1_self_first_paf", "rank1_self_final_paf",
    "trf_candidate_report", "trf_summary", "trf_candidates_before", "trf_candidates_after", "trf_fragments",
    "library_before_self_fasta", "self_alignment_summary", "self_alignment_events",
    "self_alignment_first_paf", "self_alignment_final_paf",
    "cleanup_report",
)


def retain_trf_candidate_artifacts(directory, audit, trf_context, sizes):
    """Keep source-coordinate TRF evidence, separate from PAF-driven trimming."""
    trf_dir = directory / "trf"
    trf_dir.mkdir(parents=True, exist_ok=True)
    for key, filename in (("before", "bundle_candidates.before.fa"),
                          ("after", "bundle_candidates.after.fa")):
        write_fasta(trf_dir / filename, (
            (f"candidate_{index}|source={candidate['source_bundle_id']}|position={candidate['position']}"
             f"|msa_start={candidate['msa_start']}|msa_end={candidate['msa_end']}", candidate["sequence"])
            for index, candidate in enumerate(audit[key], 1)))
    columns = list(TRF_REPORT_COLUMNS)
    columns.extend(sorted({key for row in audit["reports"] for key in row} - set(columns)))
    write_tsv(trf_dir / "candidate_report.tsv", columns,
              ([row.get(key, "") for key in columns] for row in audit["reports"]))
    occurrence_counts = defaultdict(int)
    for candidate in audit["after"]:
        for path, *_ in candidate["sources"]:
            occurrence_counts[path, candidate["source_bundle_id"]] += 1
    with (trf_dir / "candidate_fragments.pgr.bed").open("w") as output:
        output.write("# TRF-retained source fragments; source IDs, NOT final library IDs or remapping\n")
        for index, candidate in enumerate(audit["after"], 1):
            source_id = candidate["source_bundle_id"]
            for path, start, end, direction, row_number in candidate["sources"]:
                flag = "R" if occurrence_counts[path, source_id] > 1 else "U"
                size = max(1, sizes.get(source_id, 100))
                positions = (0, size - 1) if direction == "+" else (size - 1, 0)
                annotation = (f"{source_id}:{size}:{0 if direction == '+' else 1}:"
                              f"{positions[0]}:{positions[1]}:{flag}")
                output.write(f"{path}\t{start}\t{end}\t{annotation}\n")
    summary = trf_context["summary"] if trf_context else {}
    write_json(trf_dir / "trf_refinement_summary.json", {
        "completed": True, "enabled": bool(summary.get("enabled")),
        "no_op_reason": summary.get("no_op_reason", "no_trf_context"),
        "candidate_count_before": len(audit["before"]), "fragment_count_after": len(audit["after"]),
        "report_row_count": len(audit["reports"]),
        "divergence_required_for_split": True,
        "full_repeat_coverage_semantics": "maximum overlap with one merged TRF interval",
        "coordinates": "half-open real source-MSA columns and exact oriented genomic source intervals",
        "pre_group_summary": str(trf_context["summary_path"].resolve()) if trf_context else None,
    })


def inherit_trf_cut_templates(audit, cut_templates, retired_original_ids):
    """Prevent later copies from silently restoring a TRF-split parent.

    Consensus offsets come from its supported real-MSA profile, not an affine
    genomic projection. The next group's MAFFT aligns these templates alongside
    its source occurrences, exactly as for mapping/self-alignment splits.
    """
    removed_columns = defaultdict(set)
    for row in audit["reports"]:
        if row["status"] != "repeat_removed":
            continue
        key = (row["source_bundle_id"], row["position"],
               row["candidate_msa_start"], row["candidate_msa_end"])
        removed_columns[key].update((row["msa_start"], row["msa_end"]))
    for candidate in audit["before"]:
        source_id = candidate["source_bundle_id"]
        key = (source_id, candidate["position"], candidate["msa_start"], candidate["msa_end"])
        if key not in removed_columns:
            continue
        columns = [column for column, base in enumerate(candidate["profile"]) if base != "-"]
        boundaries = sorted({0, len(columns), *(bisect_left(columns, point) for point in removed_columns[key])})
        template = {"sequence": candidate["sequence"], "boundaries": boundaries}
        if template not in cut_templates[source_id]:
            cut_templates[source_id].append(template)
        retired_original_ids.add(source_id)


def source_cut_template(parent, boundaries):
    """Keep inherited cuts in the original source-MSA orientation."""
    if parent.source_reverse:
        return {"sequence": reverse_complement(parent.sequence),
                "boundaries": sorted(len(parent.sequence) - point for point in boundaries)}
    return {"sequence": parent.sequence, "boundaries": sorted(boundaries)}


def apply_library_overlap_events(entries, self_result, group_rank, cut_templates,
                                 retired_original_ids, lineage):
    """Replay simultaneous replacements, including full-length aliases/alleles.

    A discarded shared slice points to another parent's representative. Freeze
    and remove every participating parent before publishing retained outputs,
    so renaming to an existing name cannot remove its newly created record.
    """
    entries = dict(entries)
    batches = defaultdict(list)
    for event in self_result["events"]:
        batches[event["iteration"]].append(event)
    for iteration in sorted(batches, key=int):
        events = batches[iteration]
        names = {event["parent_consensus"] for event in events}
        parents = {name: entries.pop(name) for name in names}
        boundaries = defaultdict(set)
        for event in events:
            name = event["parent_consensus"]
            parent = parents[name]
            start, end = int(event["parent_start"]), int(event["parent_end"])
            boundaries[name].update((start, end))
            retained = event["retained"] in (True, 1, "1", "True", "true")
            reverse = event.get("output_reverse", False) in (True, 1, "1", "True", "true")
            sequence = parent.sequence[start:end]
            if reverse:
                sequence = reverse_complement(sequence)
            child_name = event["consensus_name"]
            if retained:
                if child_name in entries:
                    raise ValueError(f"duplicate retained overlap output: {child_name}")
                members = list(parent.members) if start == 0 and end == len(parent.sequence) else [sequence]
                if reverse and start == 0 and end == len(parent.sequence):
                    members = [reverse_complement(member) for member in members]
                entries[child_name] = LibraryEntry(child_name, sequence, parent.source_bundle_id,
                                                  parent.group_rank, members,
                                                  parent.source_reverse != reverse)
            action = ("rank1_" if group_rank == 1 else "addition_") + "overlap_" + event["action"]
            lineage.append((action, group_rank, child_name, name, start, end,
                            parent.source_bundle_id, "", "", "", "", "-" if reverse else "+", "", "", ""))
        for name, points in boundaries.items():
            parent = parents[name]
            if any(0 < point < len(parent.sequence) for point in points):
                template = source_cut_template(parent, points)
                if template not in cut_templates[parent.source_bundle_id]:
                    cut_templates[parent.source_bundle_id].append(template)
            unchanged = (name in entries and entries[name].sequence == parent.sequence
                         and entries[name].source_bundle_id == parent.source_bundle_id)
            if not unchanged and normalize_bundle_name(name, {}) == parent.source_bundle_id:
                retired_original_ids.add(parent.source_bundle_id)
    if {name: entry.sequence for name, entry in entries.items()} != self_result["records"]:
        raise ValueError("library overlap lineage does not match final active library")
    return entries


def apply_library_self_splits(entries, self_result, group_rank, cut_templates,
                              retired_original_ids, lineage):
    """Replay iterative self cuts without losing a parent's source or rank.

    Intermediate children can themselves be split before the engine returns.
    Keep each parent's entry while replaying its chronological events, and
    propagate every accepted boundary to future real source MSAs.
    """
    if any("retained" in event for event in self_result["events"]):
        return apply_library_overlap_events(entries, self_result, group_rank, cut_templates,
                                            retired_original_ids, lineage)
    entries = dict(entries)
    parent_boundaries, parents = defaultdict(set), {}
    action = "rank1_self_split" if group_rank == 1 else "addition_self_split"
    for event in self_result["events"]:
        name = event["parent_consensus"]
        parent = entries.get(name) or parents[name]
        parents[name] = parent
        entries.pop(name, None)
        child_name = event["consensus_name"]
        sequence = parent.sequence[event["parent_start"]:event["parent_end"]]
        entries[child_name] = LibraryEntry(child_name, sequence, parent.source_bundle_id,
                                          parent.group_rank, [sequence], parent.source_reverse)
        parent_boundaries[name].update((event["parent_start"], event["parent_end"]))
        lineage.append((action, group_rank, child_name, name,
                        event["parent_start"], event["parent_end"], parent.source_bundle_id,
                        "", "", "", "", "", "", "", ""))
    for name, boundaries in parent_boundaries.items():
        parent = parents[name]
        template = source_cut_template(parent, boundaries)
        if template not in cut_templates[parent.source_bundle_id]:
            cut_templates[parent.source_bundle_id].append(template)
        if normalize_bundle_name(parent.name, {}) == parent.source_bundle_id:
            retired_original_ids.add(parent.source_bundle_id)
    if {name: entry.sequence for name, entry in entries.items()} != self_result["records"]:
        raise ValueError("library self split lineage does not match final active library")
    return entries


def retain_round_artifacts(directory, group, files, args, counts, self_result=None, effective_bed=None):
    """Index retained real artifacts; skipped stages have null paths, not fake data.

    Paths are relative to the round directory so the complete output tree can
    be moved. Existing raw PAFs/diagnostics are retained in place, not duplicated.
    """
    def relative(path, required=True):
        if path is None:
            return None
        if not path.is_file():
            if required:
                raise ValueError(f"missing round artifact: {path}")
            return None
        return str(path.relative_to(directory))

    source_bed = directory / "source_bundles.pgr.bed"
    if args.keep_intermediates:
        current_paths = set(group.path_ids)
        with (effective_bed or args.bed).open() as source, source_bed.open("w") as output:
            output.write("# Current-rank source BED after pre-group TRF filtering; NOT consensus remapping\n")
            for line in source:
                fields = line.split()
                if fields and fields[0] in current_paths:
                    output.write(line if line.endswith("\n") else line + "\n")

    artifacts = {
        "source_pgr_bed": relative(source_bed, required=False),
        "before_trim_paf": None, "before_trim_pgr_bed": None,
        "after_trim_paf": None, "after_trim_pgr_bed": None,
        "historical_evidence_paf": None, "population_split_evidence_paf": None,
        "split_stats": None, "mapping_observations": None, "cut_support": None,
        "bundle_fragments": None, "potential_endpoints": None,
        "gap_alignment_candidates": None, "gap_alignment_groups": None,
        "split_preliminary_bed": None, "split_uncovered_bed": None, "split_summary": None,
        "library_before_fasta": None, "library_after_trim_fasta": None,
        "library_after_fasta": relative(directory / "library_after.fa", required=False),
        "library_after_metadata": relative(directory / "library_after.metadata.bed", required=False),
        "additions": relative(directory / "additions.tsv"),
        "mafft_quick_selection": (relative(directory / "mafft_quick_selection.tsv", required=False)
                                  if args.mafft_quick else None),
        "mafft_quick_placement": (relative(directory / "mafft_quick_placement.tsv", required=False)
                                  if args.mafft_quick else None),
        "library_initial_fasta": (relative(directory / "library_initial.fa", required=False)
                                  if group.rank == 1 else None),
        "rank1_self_summary": None, "rank1_self_events": None,
        "rank1_self_first_paf": None, "rank1_self_final_paf": None,
        "library_before_self_fasta": None, "library_before_self_metadata": None,
        "self_alignment_summary": None, "self_alignment_events": None,
        "self_alignment_first_paf": None, "self_alignment_final_paf": None,
        "trf_candidate_report": None, "trf_summary": None, "trf_candidates_before": None,
        "trf_candidates_after": None, "trf_fragments": None,
        "cleanup_report": relative(directory / "intermediate_cleanup.json"),
    }
    for key, filename in (("library_before_self_fasta", "library_before_self.fa"),
                          ("library_before_self_metadata", "library_before_self.metadata.bed")):
        path = directory / filename
        if path.is_file():
            artifacts[key] = relative(path)
    for key, filename in (("trf_candidate_report", "candidate_report.tsv"),
                          ("trf_summary", "trf_refinement_summary.json"),
                          ("trf_candidates_before", "bundle_candidates.before.fa"),
                          ("trf_candidates_after", "bundle_candidates.after.fa"),
                          ("trf_fragments", "candidate_fragments.pgr.bed")):
        path = directory / "trf" / filename
        if path.is_file():
            artifacts[key] = relative(path)
    if self_result is not None:
        self_dir = self_result["summary_path"].parent
        artifacts.update({
            "self_alignment_summary": relative(self_result["summary_path"], required=False),
            "self_alignment_events": relative(self_dir / "self_split_events.tsv", required=False),
            "self_alignment_first_paf": relative(self_dir / self_result["summary"]["iterations"][0]["paf"]),
            "self_alignment_final_paf": relative(self_dir / self_result["summary"]["iterations"][-1]["paf"]),
        })
        if group.rank == 1:
            for suffix in ("summary", "events", "first_paf", "final_paf"):
                artifacts["rank1_self_" + suffix] = artifacts["self_alignment_" + suffix]
    split_status = "not_applicable_initialization"
    if files is not None:
        split_status = "disabled" if args.no_mapping_split else "performed"
        rank1_split_only = group.rank == 1 and not args.no_mapping_split
        post_paf = (None if rank1_split_only else
                    files.before_paf if args.no_mapping_split else
                    files.directory / "bundle_to_genome.split.raw.paf")
        artifacts.update({
            "before_trim_paf": relative(files.before_paf),
            "before_trim_pgr_bed": relative(
                files.directory / "before_trim" / "refined_bundle_mappings.pgr.bed", required=False),
            "after_trim_paf": relative(post_paf),
            "after_trim_pgr_bed": (None if rank1_split_only else
                                    relative(files.directory / "refined_bundle_mappings.pgr.bed",
                                             required=False)),
            "library_before_fasta": relative(files.before_fasta, required=False),
            "library_after_trim_fasta": relative(files.trimmed_fasta, required=False),
        })
        if not args.no_mapping_split:
            evidence = files.directory / "split_evidence"
            historical = evidence / "previous_ranks.raw.paf"
            artifacts["historical_evidence_paf"] = relative(historical) if historical.is_file() else None
            artifacts["population_split_evidence_paf"] = relative(evidence / "population_split_evidence.paf")
            split_dir = files.directory / "mapping_split"
            for key, name in (
                ("split_stats", "bundle_split_stats.tsv"),
                ("mapping_observations", "mapping_observations.tsv"),
                ("cut_support", "bundle_cut_support.tsv"),
                ("bundle_fragments", "bundle_fragments.tsv"),
                ("potential_endpoints", "potential_bundle_endpoints.tsv"),
                ("gap_alignment_candidates", "gap_alignment_candidates.tsv"),
                ("gap_alignment_groups", "gap_alignment_groups.tsv"),
                ("split_preliminary_bed", "preliminary_bundles.bed"),
                ("split_uncovered_bed", "uncovered_genome_regions.bed"),
                ("split_summary", "mapping_split_summary.json"),
            ):
                artifacts[key] = relative(split_dir / name, required=False)
    manifest = directory / "round_artifacts.json"
    write_json(manifest, {
        "version": 1, "completed": True, "rank": group.rank, "group_id": group.group_id,
        "round_type": "initialization" if group.rank == 1 else "remap",
        "split_status": split_status,
        "rank1_self_status": ("not_applicable" if group.rank != 1 else
                              "performed" if self_result is not None else
                              "disabled" if group.rank == 1 and args.no_rank1_self_split else
                              "not_performed_empty_library" if group.rank == 1 else "not_applicable"),
        "self_alignment_status": ("performed" if self_result is not None else
                                  "disabled" if (args.no_rank1_self_split if group.rank == 1
                                                 else getattr(args, "no_addition_self_split", False)) else
                                  "not_performed_empty_library" if counts["library_after"] == 0 else
                                  "not_performed_no_additions"),
        "sample_count": group.sample_count, "current_rank_paths": list(group.path_ids), "counts": counts,
        "retention_mode": "keep_all" if args.keep_intermediates else "raw_paf_and_additions_only",
        "retained_raw_pafs": sorted(str(path.relative_to(directory))
                                    for path in directory.rglob("*.paf") if path.is_file()),
        "semantics": {
            "paths_relative_to": "this round directory", "paf_axes": "raw query=genome, target=bundle",
            "mapping_beds": ("rank1 before_trim is split evidence only; no after_trim/addition BED"
                             if group.rank == 1 and files is not None else
                             "filtered current-rank old-library coverage; after_trim is BEFORE additions"),
            "truncation_scope": "current and previously processed ranks, pre-split consensus coordinates",
            "source_bed": "current-rank occurrences after pre-group TRF filtering, not a library remapping",
            "trf_candidates": "source-MSA evidence before addition; TRF split requires divergence, independent of mapping split",
            "rank1_self_alignments": "bundle-to-bundle PAFs; all iterations indexed by rank1_self_summary",
            "post_addition_self_alignments": "whole-library bundle PAFs after actual additions; all iterations indexed by self_alignment_summary",
            "self_split_snapshot": "library_before_self is after TRF refinement/additions, before self splits; library_after includes self children",
            "skipped_artifacts": "null means stage was not performed; no fabricated PAF or split evidence",
        },
        "artifacts": artifacts,
    })
    # Root index uses paths relative to the library output directory.
    row = {"rank": group.rank, "group_id": group.group_id,
           "round_type": "initialization" if group.rank == 1 else "remap",
           "split_status": split_status,
           "manifest": str(manifest.relative_to(args.output_dir))}
    for key in ROUND_INDEX_COLUMNS[5:]:
        row[key] = str((directory / artifacts[key]).relative_to(args.output_dir)) if artifacts[key] else ""
    return tuple(row[key] for key in ROUND_INDEX_COLUMNS)


REFINE_OPTIONS = (
    "min_bundle_coverage", "min_identity", "min_mapq", "max_bundle_gap", "max_genome_gap",
    "max_overlap_fraction", "fixed_bundle_min_coverage", "max_fixed_overlap_fraction",
    "max_fixed_overlap_bp", "mapping_split_min_truncated_count", "mapping_split_truncated_fraction",
    "mapping_split_max_unmapped_bp", "mapping_split_cut_merge_bp", "mapping_split_min_cut_support",
    "mapping_split_gap_min_mapping_bp", "mapping_split_gap_max_overlap_bp",
    "mapping_split_min_terminal_bp",
)


def library_run_fingerprint(args, fingerprint):
    parameters = {key: value for key, value in vars(args).items()
                  if not key.startswith("_") and key not in {
                      "prepare_only", "resume", "resume_ranks", "ranking_dir", "output_dir", "bed", "fasta",
                      "sample_map", "trf_bed", "keep_intermediates"}}
    return {
        "input": fingerprint, "parameters": parameters, "version": VERSION,
        "mapping_contract": {
            "rank_minimap_gap_bp": 500,
            "self_alignment_minimap_gap_bp": 500,
            "rank1_genome_mapping": "split_only",
        },
    }


def build_library(args, by_path, genomes, groups, sizes, fingerprint, trf_context=None):
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    summary_path = output / "library_summary.json"
    run_fingerprint = library_run_fingerprint(args, fingerprint)
    resume_ranks = getattr(args, "resume_ranks", False)
    if (args.resume or resume_ranks) and summary_path.is_file():
        try:
            previous = json.loads(summary_path.read_text())
        except ValueError:
            if not resume_ranks:
                raise
            previous = {}  # A crash may have interrupted summary publication.
        if previous.get("run_fingerprint") == run_fingerprint and previous.get("completed"):
            try:
                validate_ranked_library(output / "bundle.fa", summary_path)
            except (ValueError, OSError):
                if not resume_ranks:
                    raise
            else:
                print("Reusing completed ranked library (inputs and parameters unchanged)")
                return
        if not resume_ranks:
            raise ValueError("--resume requires a completed compatible library; use --resume-ranks "
                             "for a run interrupted after a completed rank")
    entries, lineage, iterations, cut_templates = {}, [], [], defaultdict(list)
    round_artifacts = []
    retired_original_ids = set()
    next_id = max(map(int, sizes)) + 1
    comparison_index = 0
    completed_rank = 0
    group_state = json.loads(json.dumps([asdict(group) for group in groups]))
    if resume_ranks:
        checkpoint = getattr(args, "_rank_checkpoint", None)
        if checkpoint is None:
            checkpoint = load_rank_checkpoint(output / RANK_CHECKPOINT, run_fingerprint)
        if checkpoint["groups"] != group_state:
            raise ValueError("--resume-ranks haplotype ranking changed")
        completed_rank = checkpoint["completed_rank"]
        # Completion must be a contiguous prefix, paired with a saved state.
        # A later stale/partly written directory is rerun, never assumed done.
        for group in groups[:completed_rank]:
            manifest = output / "rounds" / f"rank_{group.rank:04d}" / "round_artifacts.json"
            if not manifest.is_file():
                raise ValueError(f"--resume-ranks missing completed rank manifest: {manifest}")
            data = json.loads(manifest.read_text())
            if not data.get("completed") or data.get("rank") != group.rank or data.get("group_id") != group.group_id:
                raise ValueError(f"--resume-ranks rank manifest disagrees with checkpoint: {manifest}")
        restored = [LibraryEntry(**entry) for entry in checkpoint["entries"]]
        entries = {entry.name: entry for entry in restored}
        if (len(entries) != len(restored)
                or any(not e.sequence or not 1 <= e.group_rank <= completed_rank for e in restored)):
            raise ValueError("invalid library entries in rank checkpoint")
        sizes.clear()
        sizes.update(checkpoint["sizes"])
        next_id, comparison_index = checkpoint["next_id"], checkpoint["comparison_index"]
        if next_id <= max(map(int, sizes)):
            raise ValueError("invalid reserved bundle IDs in rank checkpoint")
        lineage, iterations = checkpoint["lineage"], checkpoint["iterations"]
        cut_templates.update(checkpoint["cut_templates"])
        retired_original_ids.update(checkpoint["retired_original_ids"])
        round_artifacts = checkpoint["round_artifacts"]
        print(f"Resuming after completed rank {completed_rank}/{len(groups)}; "
              f"next rank {completed_rank + 1}" if completed_rank < len(groups) else
              f"All {len(groups)} ranks complete; publishing final library", flush=True)

    def save_checkpoint(rank):
        write_rank_checkpoint(output / RANK_CHECKPOINT, {
            "checkpoint_version": 1, "run_fingerprint": run_fingerprint,
            "completed_rank": rank, "groups": group_state,
            "entries": [asdict(entry) for entry in entries.values()], "sizes": sizes,
            "next_id": next_id, "comparison_index": comparison_index,
            "lineage": lineage, "iterations": iterations, "cut_templates": dict(cut_templates),
            "retired_original_ids": sorted(retired_original_ids), "round_artifacts": round_artifacts,
        })

    # Invalidate completion before mutating an interrupted or rerun library.
    if not resume_ranks:
        save_checkpoint(0)
    write_json(summary_path, {"version": VERSION, "completed": False, "rank_count": len(groups),
                              "ranks_processed": completed_rank, "run_fingerprint": run_fingerprint,
                              "rank_checkpoint": RANK_CHECKPOINT})
    write_ranking(output, groups, fingerprint)
    write_tsv(output / "round_artifacts.tsv", ROUND_INDEX_COLUMNS, round_artifacts)
    for group in groups[completed_rank:]:
        print(f"[library rank {group.rank}/{len(groups)}] {group.sample_count} samples, "
              f"{group.path_count} paths", flush=True)
        directory = output / "rounds" / f"rank_{group.rank:04d}"
        directory.mkdir(parents=True, exist_ok=True)
        args._round_cleanup = []
        args._generated_rank_fastas = []
        write_round_cleanup(args, directory)
        # A rerun must not leave a stale completed manifest beside partially
        # rewritten files if this round fails before publishing its new audit.
        write_json(directory / "round_artifacts.json", {
            "version": 1, "completed": False, "rank": group.rank, "group_id": group.group_id,
        })
        before = len(entries)
        coverage, splits, files, self_result = {}, [], None, None
        if entries:
            parents = dict(entries)
            print("  [1/4] Mapping the previous rank's library to current-rank samples", flush=True)
            files = map_old_library_to_rank(entries, sizes, group, genomes, directory, args)
            print("  [2/4] Splitting old bundles before additions (historical evidence pooled separately)", flush=True)
            entries, splits = split_old_library_before_additions(
                entries, sizes, group, groups[:group.rank - 1], genomes, files, args)
            print("  [3/4] Mapping the trimmed library back to the same current-rank samples", flush=True)
            coverage = remap_trimmed_library_to_rank(entries, group, genomes, files, args)
            entries, next_id = register_mapping_splits(
                parents, splits, entries, sizes, group.rank, lineage, cut_templates,
                retired_original_ids, next_id)
        write_tsv(directory / "old_library_coverage.bed", ("#path", "start", "end"),
                  ((path, start, end) for path in group.path_ids
                   for start, end in coverage.get(path, [])))
        candidates, decisions = [], []
        print("  [4/4] Checking additions against the post-trim current-rank BED" if before
              else "  Initializing the rank-1 library in source-bundle order", flush=True)
        # Copy masks: accepted old coverage remains an immutable round audit.
        masks = {path: list(items) for path, items in coverage.items()}
        trf_audit = {"before": [], "after": [], "reports": []}
        for candidate in position_candidates(group, by_path, genomes, masks, directory, args, cut_templates,
                                             trf_context, trf_audit):
            candidates.append(candidate)
        retain_trf_candidate_artifacts(directory, trf_audit, trf_context, sizes)
        inherit_trf_cut_templates(trf_audit, cut_templates, retired_original_ids)
        # Same original ID/position copies share real MSA coordinates. Complete-link
        # >90% clustering prevents chaining unrelated copies through one bridge.
        clusters = []
        for candidate in candidates:
            destination = None
            for cluster in clusters:
                if (candidate["source_bundle_id"] == cluster[0]["source_bundle_id"]
                        and candidate["whole"] == cluster[0]["whole"]
                        and all(aligned_profile_identity(candidate["profile"], member["profile"])[0]
                                > args.merge_identity for member in cluster)):
                    destination = cluster
                    break
            if destination is None:
                clusters.append([candidate])
            else:
                destination.append(candidate)
        added = reused = 0
        for cluster in clusters:
            candidate = cluster[0]
            source_id = candidate["source_bundle_id"]
            profiles = [row for member in cluster for row in member["aligned_sequences"]]
            _, sequence, _, _, _ = aligned_consensus_profile(
                profiles, args.min_column_occupancy, args.min_base_frequency)
            reused_entry = None
            # A later full copy of an already represented source class may be a
            # divergent allele. Compare by real MAFFT, not ungapped zip/length.
            for entry in entries.values():
                if (entry.source_bundle_id != source_id
                        or (candidate["whole"] and normalize_bundle_name(entry.name, {}) != source_id)):
                    continue
                compatible = True
                for member in entry.members or [entry.sequence]:
                    for new_member in cluster:
                        comparison_index += 1
                        new_sequence = (reverse_complement(new_member["sequence"])
                                        if entry.source_reverse else new_member["sequence"])
                        aligned = run_mafft([("old", member), ("new", new_sequence)], directory,
                                            f"compare_{comparison_index}", args)
                        if aligned_profile_identity(aligned["old"], aligned["new"])[0] <= args.merge_identity:
                            compatible = False
                            break
                    if not compatible:
                        break
                if compatible:
                    reused_entry = entry
                    break
            if reused_entry is not None:
                name, action = reused_entry.name, "reuse_high_identity"
                reused_entry.members.extend(reverse_complement(member["sequence"]) if reused_entry.source_reverse
                                            else member["sequence"] for member in cluster)
                reused += 1
            else:
                if candidate["whole"] and source_id not in retired_original_ids:
                    bundle_id = source_id
                    existing_names = [entry.name for entry in entries.values()
                                      if normalize_bundle_name(entry.name, {}) == bundle_id]
                    name = f"bundle_{bundle_id}"
                    if existing_names:
                        allele = 2
                        while f"bundle_{bundle_id}.allele_{allele}" in entries:
                            allele += 1
                        name = f"bundle_{bundle_id}.allele_{allele}"
                else:
                    bundle_id = str(next_id)
                    next_id += 1
                    name = f"bundle_{bundle_id}"
                    sizes[bundle_id] = max(1, sizes.get(source_id, 100))
                entries[name] = LibraryEntry(name, sequence, source_id, group.rank,
                                             [member["sequence"] for member in cluster])
                action = "initialize" if group.rank == 1 else "add_full" if candidate["whole"] else "add_residual"
                added += 1
            for member in cluster:
                for path, start, end, direction, row_number in member["sources"]:
                    lineage.append((action, group.rank, name, "", "", "", source_id,
                                    member["position"], path, start, end, direction,
                                    row_number, member["msa_start"], member["msa_end"]))
                decisions.append((member["position"], source_id, member["msa_start"], member["msa_end"],
                                  len(member["sequence"]), int(member["whole"]), action, name,
                                  len(member["sources"])))
        if group.rank == 1:
            write_fasta(directory / "library_initial.fa", ((entry.name, entry.sequence) for entry in entries.values()))
            write_metadata(directory / "library_initial.metadata.bed", sizes)
            if entries and not args.no_mapping_split:
                parents = dict(entries)
                print("  Mapping the initialized rank-1 library back to rank-1 samples for splitting only",
                      flush=True)
                entries, rank1_splits, files = split_rank1_library_by_population(
                    entries, sizes, group, genomes, directory, args)
                splits.extend(rank1_splits)
                entries, next_id = register_mapping_splits(
                    parents, rank1_splits, entries, sizes, group.rank, lineage, cut_templates,
                    retired_original_ids, next_id)
        # Compare the updated working library, not the old pre-addition FASTA.
        # A whole-library scan includes old/new and new/new pairs and revisits
        # old/old pairs whose current fragment boundaries have changed.
        write_fasta(directory / "library_before_self.fa", ((entry.name, entry.sequence) for entry in entries.values()))
        write_metadata(directory / "library_before_self.metadata.bed", sizes)
        self_enabled = not (args.no_rank1_self_split if group.rank == 1 else args.no_addition_self_split)
        if entries and self_enabled and (group.rank == 1 or added > 0):
            self_result = split_rank1_library(
                {name: entry.sequence for name, entry in entries.items()}, sizes,
                directory / "self_alignment", args.minimap2, args.threads,
                args.rank1_self_min_overlap_fraction, args.rank1_self_min_overlap_bp,
                0, args.min_identity, args.min_mapq,
                stage_label="rank1" if group.rank == 1 else f"rank{group.rank} post-addition",
                reconcile_overlaps=not args.no_self_overlap_reconciliation,
                merge_identity=args.merge_identity,
                priorities={name: entry.group_rank for name, entry in entries.items()},
                keep_intermediates=args.keep_intermediates,
            )
            sizes.update(self_result["sizes"])
            next_id = max(next_id, max(map(int, sizes)) + 1)
            entries = apply_library_self_splits(entries, self_result, group.rank,
                                               cut_templates, retired_original_ids, lineage)
        write_tsv(directory / "additions.tsv",
                  ("position_1based", "source_bundle_id", "msa_start", "msa_end", "consensus_bp",
                   "whole_source", "action", "consensus_name", "source_occurrences"), decisions)
        write_fasta(directory / "library_after.fa", ((entry.name, entry.sequence) for entry in entries.values()))
        write_metadata(directory / "library_after.metadata.bed", sizes)
        args._round_cleanup.append(cleanup_generated_files(
            directory, args._generated_rank_fastas,
            protected_paths=getattr(args, "_cleanup_protected", ()),
            keep_intermediates=args.keep_intermediates,
        ))
        args._round_cleanup.append(cleanup_rank_round(
            directory,
            protected_paths=getattr(args, "_cleanup_protected", ()),
            keep_intermediates=args.keep_intermediates,
        ))
        write_round_cleanup(args, directory, completed=True, self_result=self_result)
        round_artifacts.append(retain_round_artifacts(directory, group, files, args, {
            "library_before": before, "split_parents": len({row["parent_consensus"] for row in splits}),
            "added_consensuses": added, "reused_consensuses": reused, "library_after": len(entries),
            "rank1_self_split_parents": self_result["summary"]["counts"]["split_parent_count"] if self_result and group.rank == 1 else 0,
            "rank1_self_iterations": self_result["summary"]["counts"]["iteration_count"] if self_result and group.rank == 1 else 0,
            "self_split_parents": self_result["summary"]["counts"]["split_parent_count"] if self_result else 0,
            "self_iterations": self_result["summary"]["counts"]["iteration_count"] if self_result else 0,
        }, self_result, trf_context["bed"] if trf_context else None))
        # Publish the index after every completed round, also for interrupted runs.
        write_tsv(output / "round_artifacts.tsv", ROUND_INDEX_COLUMNS, round_artifacts)
        iterations.append((group.rank, group.group_id, group.sample_count, group.path_count,
                           before, len({row["parent_consensus"] for row in splits}), len(candidates),
                           added, reused, len(entries),
                           self_result["summary"]["counts"]["split_parent_count"] if self_result else 0,
                           self_result["summary"]["counts"]["iteration_count"] if self_result else 0))
        save_checkpoint(group.rank)
        write_json(summary_path, {"version": VERSION, "completed": False, "rank_count": len(groups),
                                  "ranks_processed": group.rank, "run_fingerprint": run_fingerprint,
                                  "rank_checkpoint": RANK_CHECKPOINT})
    if not entries:
        raise ValueError("no bundle consensus passed the library thresholds")
    write_fasta(output / "bundle.fa", ((entry.name, entry.sequence) for entry in entries.values()))
    shutil.copyfile(output / "bundle.fa", output / "bundle_consensus.fa")
    write_metadata(output / "bundle_metadata.bed", sizes)
    write_tsv(output / "library_lineage.tsv",
              ("action", "group_rank", "consensus_name", "parent_consensus", "parent_start", "parent_end",
               "source_bundle_id", "position_1based", "path_id", "genome_start", "genome_end", "strand",
               "source_bed_row", "source_msa_start", "source_msa_end"), lineage)
    write_tsv(output / "rank_iterations.tsv",
              ("rank", "group_id", "sample_count", "path_count", "library_before", "split_parents",
               "residual_candidates", "added_consensuses", "reused_consensuses", "library_after",
               "self_split_parents", "self_iterations"), iterations)
    write_tsv(output / "bundle_consensus.tsv", ("bundle_id", "output_id", "consensus_length", "group_rank"),
              ((normalize_bundle_name(entry.name, {}), entry.name, len(entry.sequence), entry.group_rank)
               for entry in entries.values()))
    # Compatibility filename, with truthful coverage semantics: all ranks examined,
    # not all original IDs retained (covered originals intentionally disappear).
    write_tsv(output / "bundle_source_groups.tsv",
              ("bundle_id", "group_rank", "source_bundle_id", "all_ranks_examined"),
              ((normalize_bundle_name(entry.name, {}), entry.group_rank, entry.source_bundle_id, 1)
               for entry in entries.values()))
    write_json(summary_path, {
        "version": VERSION, "completed": True, "rank_count": len(groups), "ranks_processed": len(groups),
        "run_fingerprint": run_fingerprint, "active_consensus_names": list(entries),
        "fasta_sha256": file_digest(output / "bundle.fa"),
        "metadata_sha256": file_digest(output / "bundle_metadata.bed"),
        "round_artifacts_index": "round_artifacts.tsv",
        "rank_checkpoint": RANK_CHECKPOINT,
        "trf_pre_group_summary": str(trf_context["summary_path"].resolve()) if trf_context else None,
        "semantics": {
            "ranking": "distinct samples per exact signed bundle order; paths if no sample-map",
            "mafft_quick": {
                "enabled": args.mafft_quick,
                "short_bundle_no_mafft_below_bp": SHORT_BUNDLE_MAFFT_MIN_BP,
                "source_reference": "longest current-rank original-BED occurrence per source ID",
                "scope": "every source-bundle FASTA; pairwise identity checks unchanged",
                "reference_window_bp": args.mafft_window_bp,
                "backbone_minimum_length_ratio": 0.7,
                "short_truncated_exclusion": "whole source occurrence or cut template length <=10% of reference AND <1000 bp; before mapping/MSA",
                "placement_minimum_query_coverage": None,
                "mafft_add_query_coverage_below": 0.95,
                "other_rows": "minimap2 -c to reference; >=95% aligned rows join windows, partially aligned rows use full-backbone MAFFT --add",
                "unmappable_rows_and_templates": "discard records without a forward CIGAR-backed placement and audit them",
                "inherited_cut_templates": "project through windowed or MAFFT --add MSA when mapped; discard when unmapped or excluded as short truncated records",
            },
            "old_library_split_population": "current rank plus separately remapped historical ranks for evidence only",
            "rank_mapping_population": "current rank only, for both before-trim and after-trim passes",
            "round_order": "rank1 initialize -> one rank1 remap for splitting only -> bundle self comparison; later ranks map old library -> split -> remap trimmed library -> additions -> bundle self comparison",
            "rank1_population_split": "one initialized-library to rank1-sequence mapping; split only, without post-split remap or addition discovery",
            "rank1_self_split": "after initialization, largest pair coverage union first; repeat minimap2 until no accepted cuts or shared-region reconciliation",
            "post_addition_self_split": "after each later rank with actual additions, scan updated whole library; repeat minimap2 until no accepted cuts or shared-region reconciliation",
            "self_overlap_reconciliation": not args.no_self_overlap_reconciliation,
            "shared_region_identity": "strictly above merge_identity: one higher-ranked representative; otherwise unique FASTA alleles sharing a numeric bundle class",
            "short_terminal_subbundles": "genome mapping only: suppress cuts creating terminal fragments < mapping_split_min_terminal_bp",
            "self_split_fragment_lengths": "no internal minimum; suppress terminal cuts only if fragment <200 bp AND <5% of current parent; genome mapping guard is independent",
            "short_internal_subbundles": "retain as independent children; no deletion",
            "new_bundle_source": "original current-rank BED minus qualified old-library genomic spans",
            "trf_group_filter": "before signed-order grouping; one merged TRF interval must reach full_coverage; remove all ID occurrences when qualifying occurrence frequency reaches removal threshold",
            "trf_candidate_refinement": "rank1 and every new candidate before addition; exact source-to-MSA projection followed by population MSA union; long union components score population-aware TR-length divergence",
            "overlap": "only uncovered source pieces are added; overlapping source rows are reserved sequentially",
            "boundary_projection": "original genome coordinates -> oriented occurrence -> real source MSA columns",
            "new_fragment_min_bp": args.addition_min_bp,
            "later_high_identity_full_copy": "retain higher-rank consensus representative; complete-link identity > threshold",
            "retired_parent_ids": "never restored; later copies inherit mapping, self-alignment and TRF cut boundaries through real source MSA templates",
            "unmapped_original_classes": "examined every rank, including already seen IDs",
            "final_all_population_refine": "performed by caller after library completion",
            "intermediate_retention": ("keep_all" if args.keep_intermediates else
                                       "per-rank raw PAFs and additions only"),
            "cleanup_audits": "round_artifacts.tsv cleanup_report; each self summary indexes its cleanup report",
        },
    })


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ranking-dir", type=Path)
    parser.add_argument("--sample-map", type=Path)
    parser.add_argument("--trf-bed", type=Path, help="optional genomic TRF annotations; empty/comment-only skips TRF processing")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--resume-ranks", action="store_true",
                        help="restore the last checkpoint with contiguous completed rank manifests; "
                             "rerun the interrupted rank and continue")
    parser.add_argument("--keep-intermediates", action="store_true",
                        help="retain disposable successful-stage FASTAs; default keeps important intermediates only")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--mafft-extra", default="")
    parser.add_argument("--mafft-params", default="--auto",
                        help="MAFFT parameters replacing the default --auto strategy")
    parser.add_argument("--mafft-quick", action="store_true",
                        help="window-align every source FASTA; sequences at least 70 percent of the "
                             "longest occurrence form the windowed backbone, while shorter "
                             "sequences are added locally; partially mapped rows use MAFFT --add; "
                             "skip records at most 10 percent of reference length and below 1000 bp")
    parser.add_argument("--mafft-window-bp", type=int, default=10000,
                        help="reference window width for --mafft-quick (default: 10000)")
    parser.add_argument("--addition-min-bp", type=int, default=1)
    parser.add_argument("--min-column-occupancy", type=float, default=.50)
    parser.add_argument("--min-base-frequency", type=float, default=.50)
    parser.add_argument("--merge-identity", type=float, default=.90)
    parser.add_argument("--rank1-self-min-overlap-fraction", type=float, default=.20)
    parser.add_argument("--rank1-self-min-overlap-bp", type=int, default=10000)
    parser.add_argument("--no-rank1-self-split", action="store_true")
    parser.add_argument("--no-addition-self-split", action="store_true",
                        help="disable later-rank self splits after additions; independent of rank1 and mapping splits")
    parser.add_argument("--no-self-overlap-reconciliation", action="store_true",
                        help="legacy self-cut mode: retain both homologous copies under separate classes")
    trf_defaults = {"trf_min_length": 0, "trf_min_identity": 80., "trf_full_coverage": .90,
                    "trf_bundle_removal_frequency": .50, "trf_split_min_bp": 1000,
                    "trf_min_sequence_coverage": .50, "trf_min_population_presence": .50,
                    "trf_min_length_divergence": .05,
                    "trf_min_nucleotide_divergence": .05, "trf_min_fragment_sequence_bp": 1000,
                    "trf_min_fragment_count": 2, "trf_min_fragment_frequency": .02}
    for name, value in trf_defaults.items():
        parser.add_argument("--" + name.replace("_", "-"), type=type(value), default=value)
    defaults = {"min_bundle_coverage": .60, "min_identity": 0., "min_mapq": 0,
                "max_bundle_gap": 1000, "max_genome_gap": 1000, "max_overlap_fraction": .10,
                "fixed_bundle_min_coverage": .80, "max_fixed_overlap_fraction": .20,
                "max_fixed_overlap_bp": 10000, "mapping_split_min_truncated_count": 5,
                "mapping_split_truncated_fraction": .60, "mapping_split_max_unmapped_bp": 10000,
                "mapping_split_cut_merge_bp": 100, "mapping_split_min_cut_support": 1,
                "mapping_split_gap_min_mapping_bp": 10000, "mapping_split_gap_max_overlap_bp": 1000,
                "mapping_split_min_terminal_bp": 1000}
    for name, value in defaults.items():
        parser.add_argument("--" + name.replace("_", "-"), type=type(value), default=value)
    parser.add_argument("--no-mapping-split", action="store_true")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.mafft_window_bp < 1:
        parser.error("--mafft-window-bp must be a positive integer")
    args._cleanup_protected = tuple(path for path in (
        args.bed, args.fasta, args.sample_map, args.trf_bed) if path is not None)
    try:
        ranking_dir = args.ranking_dir or args.output_dir
        protected_outputs = {args.output_dir / name for name in (
            "bundle.fa", "bundle_consensus.fa", "bundle_metadata.bed", "library_summary.json",
            "library_lineage.tsv", "rank_iterations.tsv", "bundle_consensus.tsv", "bundle_source_groups.tsv",
            "haplotype_groups.tsv", "ranking_summary.json", "round_artifacts.tsv", RANK_CHECKPOINT)}
        protected_outputs.update(ranking_dir / name for name in ("haplotype_groups.tsv", "ranking_summary.json"))
        protected_outputs = {path.resolve() for path in protected_outputs}
        rounds_dir = (args.output_dir / "rounds").resolve()
        trf_dirs = {(directory / "trf_filter").resolve() for directory in (ranking_dir, args.output_dir)}
        qc_dir = (args.output_dir / "00_input_qc").resolve()
        for path in (args.bed, args.fasta, args.sample_map, args.trf_bed):
            if path is not None and (path.resolve() in protected_outputs or rounds_dir in path.resolve().parents
                                     or qc_dir == path.resolve() or qc_dir in path.resolve().parents
                                     or any(directory == path.resolve() or directory in path.resolve().parents for directory in trf_dirs)):
                raise ValueError(f"output would overwrite an input; choose a separate output directory: {path}")
        if args.threads < 1 or args.addition_min_bp < 1:
            raise ValueError("threads and addition-min-bp must be positive")
        for name in ("trf_full_coverage", "trf_bundle_removal_frequency", "trf_min_sequence_coverage",
                     "trf_min_population_presence", "trf_min_length_divergence",
                     "trf_min_nucleotide_divergence", "trf_min_fragment_frequency"):
            if not 0 <= getattr(args, name) <= 1:
                raise ValueError(f"--{name.replace('_', '-')} must be in [0, 1]")
        if not 0 <= args.trf_min_identity <= 100:
            raise ValueError("--trf-min-identity must be in [0, 100] (percent)")
        if args.trf_full_coverage <= 0 or args.trf_bundle_removal_frequency <= 0:
            raise ValueError("--trf-full-coverage and --trf-bundle-removal-frequency must be positive")
        if args.trf_min_length < 0:
            raise ValueError("--trf-min-length must be non-negative")
        for name in ("trf_split_min_bp", "trf_min_fragment_sequence_bp", "trf_min_fragment_count"):
            if getattr(args, name) < 1:
                raise ValueError(f"--{name.replace('_', '-')} must be positive")
        for name in ("min_column_occupancy", "min_base_frequency", "merge_identity", "min_bundle_coverage",
                     "min_identity", "max_overlap_fraction", "fixed_bundle_min_coverage",
                     "max_fixed_overlap_fraction", "mapping_split_truncated_fraction", "rank1_self_min_overlap_fraction"):
            if not 0 <= getattr(args, name) <= 1:
                raise ValueError(f"--{name.replace('_', '-')} must be in [0, 1]")
        if args.rank1_self_min_overlap_bp < 0:
            raise ValueError("--rank1-self-min-overlap-bp must be non-negative")
        for name in REFINE_OPTIONS:
            value = getattr(args, name)
            if isinstance(value, int) and value < (1 if name in {
                    "mapping_split_min_truncated_count", "mapping_split_min_cut_support"} else 0):
                raise ValueError(f"invalid --{name.replace('_', '-')}")
        if args.mapping_split_gap_max_overlap_bp > args.mapping_split_gap_min_mapping_bp:
            raise ValueError("mapping-split-gap-max-overlap-bp must not exceed mapping-split-gap-min-mapping-bp")
        duplicate_filter = filter_duplicate_inputs(
            args.fasta, args.bed, args.trf_bed, args.sample_map, args.output_dir / "00_input_qc")
        if duplicate_filter["duplicate_id_count"]:
            effective = duplicate_filter["filtered_inputs"]
            args.fasta = Path(effective["fasta"])
            args.bed = Path(effective["bed"])
            if args.trf_bed is not None:
                args.trf_bed = Path(effective["trf_bed"])
            if args.sample_map is not None:
                args.sample_map = Path(effective["sample_map"])
            args._cleanup_protected = tuple(path for path in (
                args.bed, args.fasta, args.sample_map, args.trf_bed) if path is not None)
            print(f"Excluded {duplicate_filter['duplicate_id_count']} duplicated FASTA IDs; "
                  f"see {duplicate_filter['report']}", flush=True)
        original_by_path = read_bundle_bed(args.bed)
        source_ids = {occ.bundle_id for occurrences in original_by_path.values() for occ in occurrences}
        if any(not name.isdigit() or name != str(int(name)) for name in source_ids):
            raise ValueError("ranked library requires canonical numeric PGR bundle IDs")
        samples = read_sample_map(args.sample_map, list(original_by_path))
        fingerprint = {"bed_sha256": file_digest(args.bed), "fasta_sha256": file_digest(args.fasta),
                       "sample_map_sha256": file_digest(args.sample_map) if args.sample_map else None,
                       "trf_bed_sha256": file_digest(args.trf_bed) if args.trf_bed else None,
                       "trf_group_semantics": "single_merged_interval_coverage_v2",
                       "trf_group_parameters": {name: getattr(args, name) for name in (
                           "trf_min_length", "trf_min_identity", "trf_full_coverage", "trf_bundle_removal_frequency")}}
        # Reject incompatible restarts before replacing their retained TRF BED
        # or reports. Preparation without --resume explicitly regenerates ranks.
        ranking_summary = ranking_dir / "ranking_summary.json"
        if not args.prepare_only and not args.resume_ranks and ranking_summary.is_file():
            if json.loads(ranking_summary.read_text()).get("input_fingerprint") != fingerprint:
                raise ValueError("ranking inputs changed; rerun ranking preparation")
        library_summary = args.output_dir / "library_summary.json"
        if not args.prepare_only and args.resume_ranks:
            args._rank_checkpoint = load_rank_checkpoint(
                args.output_dir / RANK_CHECKPOINT, library_run_fingerprint(args, fingerprint))
        if not args.prepare_only and args.resume and not args.resume_ranks and library_summary.is_file():
            if json.loads(library_summary.read_text()).get("run_fingerprint") != library_run_fingerprint(args, fingerprint):
                raise ValueError("--resume inputs/parameters changed; rerun without --resume in a separate output directory")
            validate_ranked_library(args.output_dir / "bundle.fa", library_summary)
        trf_context = filter_pgr_bed_by_trf(
            args.bed, args.trf_bed, ranking_dir / "trf_filter", args.trf_min_length, args.trf_min_identity,
            args.trf_full_coverage, args.trf_bundle_removal_frequency)
        try:
            by_path = read_bundle_bed(trf_context["bed"])
        except ValueError as error:
            if "no PGR bundle occurrences" not in str(error):
                raise
            raise ValueError("no bundles remain after pre-group TRF filtering; see trf_filter/filter_summary.json") from error
        groups = rank_groups(by_path, samples)
        genomes = {}
        seen = set()
        for header, sequence in iter_fasta(args.fasta):
            name = header.split()[0]
            if name in seen:
                raise ValueError(f"duplicate FASTA path: {name}")
            seen.add(name)
            if name in original_by_path:
                genomes[name] = sequence
        if set(original_by_path) - genomes.keys():
            raise ValueError("population FASTA lacks BED paths")
        for path, occurrences in original_by_path.items():
            if any(occ.end > len(genomes[path]) for occ in occurrences):
                raise ValueError(f"BED occurrence exceeds FASTA length: {path}")
        if args.prepare_only:
            write_ranking(ranking_dir, groups, fingerprint)
            print(f"Ranked {len(groups)} haplotype groups by sample count")
            return 0
        if not ranking_summary.is_file():
            write_ranking(ranking_dir, groups, fingerprint)
        for name in ("mafft", "minimap2"):
            if shutil.which(getattr(args, name)) is None:
                raise ValueError(f"{name} executable not found: {getattr(args, name)}")
        reserved_ids, metadata = read_existing_bundle_metadata(args.bed)
        sizes = {key: metadata[key].bundle_size if key in metadata else 100 for key in reserved_ids}
        build_library(args, by_path, genomes, groups, sizes, fingerprint, trf_context)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
