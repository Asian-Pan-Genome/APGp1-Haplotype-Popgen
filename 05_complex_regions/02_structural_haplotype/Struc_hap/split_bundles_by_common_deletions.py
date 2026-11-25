#!/usr/bin/env python3
"""Split final bundle consensuses at recurrent internal MSA indels.

The input mapping table is the preliminary, locus-deduplicated output from
``refine_bundle_mappings.py``.  Each retained genomic occurrence is extracted,
oriented to its raw bundle consensus, and aligned together with that consensus.
Only internal gap runs are deletion evidence: terminal gaps from incomplete
mapping cannot create cuts.  A qualifying deletion is retained in the bundle
library as its own child segment; the operation cuts at both breakpoints rather
than discarding the sequence between them.  A population-supported insertion
relative to every consensus of a bundle class creates one cut at its consensus
anchor; the inserted sequence is not invented as a child consensus.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shlex
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple
from urllib.parse import quote, unquote

try:
    from .align_bundles_with_mafft import reverse_complement
    from .refine_bundle_mappings import (
        natural_name_key,
        normalize_bundle_name,
        read_existing_bundle_metadata,
    )
except ImportError:
    from align_bundles_with_mafft import reverse_complement
    from refine_bundle_mappings import (
        natural_name_key,
        normalize_bundle_name,
        read_existing_bundle_metadata,
    )


@dataclass(frozen=True)
class MappingOccurrence:
    mapping_id: str
    raw_consensus: str
    path_id: str
    start: int
    end: int
    strand: str
    bundle_coverage: float


@dataclass(frozen=True)
class DeletionEvent:
    mapping_id: str
    support_unit: str
    start: int
    end: int
    alignment_start: int
    alignment_end: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class CommonDeletion:
    start: int
    end: int
    support_count: int
    population_unit_count: int
    event_count: int
    alignment_start: int
    alignment_end: int

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def support_frequency(self) -> float:
        return self.support_count / self.population_unit_count


@dataclass(frozen=True)
class InsertionEvent:
    mapping_id: str
    support_unit: str
    anchor: int
    inserted_bp: int
    alignment_start: int
    alignment_end: int


@dataclass(frozen=True)
class CommonInsertion:
    anchor: int
    inserted_bp: int
    support_count: int
    population_unit_count: int
    event_count: int
    alignment_start: int
    alignment_end: int

    @property
    def support_frequency(self) -> float:
        return self.support_count / self.population_unit_count


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                header = stripped[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                name, chunks = header.split()[0], []
            else:
                if name is None:
                    raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
                chunks.append(stripped)
    if name is not None:
        yield name, "".join(chunks).upper()


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
    for name, sequence in iter_fasta(path):
        if name in records:
            raise ValueError(f"{path}: duplicate FASTA record {name!r}")
        if not sequence:
            raise ValueError(f"{path}: empty FASTA record {name!r}")
        records[name] = sequence
    if not records:
        raise ValueError(f"{path}: no FASTA records")
    return records


def wrap_fasta(sequence: str, width: int = 80) -> Iterable[str]:
    for offset in range(0, len(sequence), width):
        yield sequence[offset : offset + width] + "\n"


def read_mapping_table(path: Path, consensus_names: Sequence[str]) -> Dict[str, List[MappingOccurrence]]:
    required = {
        "mapping_id", "raw_consensus_name", "population_path", "genome_start",
        "genome_end", "strand", "bundle_coverage",
    }
    selected = set(consensus_names)
    answer: Dict[str, List[MappingOccurrence]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing mapping columns: {', '.join(sorted(missing))}")
        seen_ids = set()
        for line_number, row in enumerate(reader, 2):
            name = row["raw_consensus_name"]
            if name not in selected:
                raise ValueError(
                    f"{path}:{line_number}: mapped consensus {name!r} is absent from bundle FASTA"
                )
            mapping_id = row["mapping_id"]
            if not mapping_id or mapping_id in seen_ids:
                raise ValueError(f"{path}:{line_number}: duplicate/empty mapping_id {mapping_id!r}")
            seen_ids.add(mapping_id)
            try:
                start, end = int(row["genome_start"]), int(row["genome_end"])
                coverage = float(row["bundle_coverage"])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinate or coverage") from error
            if start < 0 or end <= start or not 0 <= coverage <= 1:
                raise ValueError(f"{path}:{line_number}: invalid mapping interval/coverage")
            if row["strand"] not in {"+", "-"}:
                raise ValueError(f"{path}:{line_number}: invalid strand {row['strand']!r}")
            answer[name].append(
                MappingOccurrence(
                    mapping_id, name, row["population_path"], start, end,
                    row["strand"], coverage,
                )
            )
    return dict(answer)


def read_sample_map(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    answer: Dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_number}: expected path_id and sample_id")
            path_id, sample_id = fields
            previous = answer.get(path_id)
            if previous is not None and previous != sample_id:
                raise ValueError(f"{path}:{line_number}: conflicting sample for {path_id!r}")
            answer[path_id] = sample_id
    if not answer:
        raise ValueError(f"{path}: no path/sample records")
    return answer


def extract_occurrences(
    genome_fasta: Path,
    occurrences_by_consensus: Mapping[str, Sequence[MappingOccurrence]],
) -> Dict[str, List[Tuple[MappingOccurrence, str]]]:
    by_path: Dict[str, List[MappingOccurrence]] = defaultdict(list)
    for occurrences in occurrences_by_consensus.values():
        for occurrence in occurrences:
            by_path[occurrence.path_id].append(occurrence)
    requested_paths = set(by_path)
    seen_paths = set()
    extracted: Dict[str, List[Tuple[MappingOccurrence, str]]] = defaultdict(list)
    for path_id, sequence in iter_fasta(genome_fasta):
        if path_id not in requested_paths:
            continue
        if path_id in seen_paths:
            raise ValueError(f"{genome_fasta}: duplicate requested record {path_id!r}")
        seen_paths.add(path_id)
        for occurrence in by_path[path_id]:
            if occurrence.end > len(sequence):
                raise ValueError(
                    f"{genome_fasta}:{path_id}: mapping interval "
                    f"{occurrence.start}-{occurrence.end} exceeds length {len(sequence)}"
                )
            fragment = sequence[occurrence.start : occurrence.end]
            if occurrence.strand == "-":
                fragment = reverse_complement(fragment)
            extracted[occurrence.raw_consensus].append((occurrence, fragment))
    missing = sorted(requested_paths - seen_paths, key=natural_name_key)
    if missing:
        raise ValueError(
            f"{genome_fasta}: missing {len(missing)} mapped path(s), including {missing[0]!r}"
        )
    return dict(extracted)


def safe_stem(name: str) -> str:
    return quote(name, safe="._-")


def write_alignment_input(
    path: Path,
    consensuses: Sequence[Tuple[str, str]],
    occurrences: Sequence[Tuple[MappingOccurrence, str]],
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for consensus_name, consensus in consensuses:
            handle.write(f">reference|consensus={quote(consensus_name, safe='._-')}\n")
            handle.writelines(wrap_fasta(consensus))
        for occurrence, sequence in occurrences:
            handle.write(
                f">occurrence|mapping={quote(occurrence.mapping_id, safe='._-')}"
                f"|path={quote(occurrence.path_id, safe='._-')}\n"
            )
            handle.writelines(wrap_fasta(sequence))


def run_mafft(
    executable: str,
    source: Path,
    alignment: Path,
    log: Path,
    threads: int,
    extra: Sequence[str],
    params: Sequence[str] = ("--auto",),
) -> None:
    temporary = alignment.with_suffix(alignment.suffix + ".tmp")
    command = [executable, *params, "--thread", str(threads), "--inputorder", *extra, str(source)]
    try:
        with temporary.open("w", encoding="utf-8") as output, log.open("w", encoding="utf-8") as error:
            subprocess.run(command, stdout=output, stderr=error, check=True)
        temporary.replace(alignment)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise


def read_alignment(path: Path) -> List[Tuple[str, str]]:
    records = list(iter_fasta(path))
    if not records:
        raise ValueError(f"{path}: empty alignment")
    widths = {len(sequence) for _, sequence in records}
    if len(widths) != 1:
        raise ValueError(f"{path}: MAFFT records do not share one alignment width")
    return records


def parse_mapping_id(header: str) -> str:
    for field in header.split("|"):
        if field.startswith("mapping="):
            return field.split("=", 1)[1]
    raise ValueError(f"alignment occurrence lacks mapping ID: {header}")


def parse_consensus_name(header: str) -> str:
    for field in header.split("|"):
        if field.startswith("consensus="):
            return unquote(field.split("=", 1)[1])
    raise ValueError(f"alignment reference lacks consensus name: {header}")


def deletion_events(
    reference: str,
    occurrence: str,
    mapping_id: str,
    support_unit: str,
    min_bp: int,
) -> List[DeletionEvent]:
    if len(reference) != len(occurrence):
        raise ValueError("reference/occurrence alignment widths differ")
    non_gap_columns = [index for index, base in enumerate(occurrence) if base != "-"]
    if not non_gap_columns:
        return []
    first, last = non_gap_columns[0], non_gap_columns[-1]
    reference_prefix = [0]
    for base in reference:
        reference_prefix.append(reference_prefix[-1] + (base != "-"))
    answer: List[DeletionEvent] = []
    column = first + 1
    while column < last:
        if occurrence[column] != "-":
            column += 1
            continue
        start_column = column
        while column <= last and occurrence[column] == "-":
            column += 1
        end_column = column
        start = reference_prefix[start_column]
        end = reference_prefix[end_column]
        # Strictly greater matches the requested "> 1000 bp" semantics.
        if end - start > min_bp and 0 < start < end < reference_prefix[-1]:
            answer.append(
                DeletionEvent(mapping_id, support_unit, start, end, start_column, end_column)
            )
    return answer


def insertion_events(
    reference: str,
    occurrence: str,
    mapping_id: str,
    support_unit: str,
    min_bp: int,
) -> List[InsertionEvent]:
    """Return internal occurrence insertions relative to the consensus backbone.

    A reference-gap run is an insertion locus only when the occurrence has
    reference-aligned sequence on both sides.  This is the insertion analogue
    of ignoring terminal deletion gaps from incomplete mappings.
    """

    if len(reference) != len(occurrence):
        raise ValueError("reference/occurrence alignment widths differ")
    reference_prefix = [0]
    reference_supported_columns = []
    for column, base in enumerate(reference):
        reference_prefix.append(reference_prefix[-1] + (base != "-"))
        if base != "-" and occurrence[column] != "-":
            reference_supported_columns.append(column)
    if not reference_supported_columns:
        return []
    first, last = reference_supported_columns[0], reference_supported_columns[-1]
    reference_length = reference_prefix[-1]
    answer: List[InsertionEvent] = []
    column = 0
    while column < len(reference):
        if reference[column] != "-":
            column += 1
            continue
        start_column = column
        while column < len(reference) and reference[column] == "-":
            column += 1
        end_column = column
        anchor = reference_prefix[start_column]
        inserted_bp = sum(base != "-" for base in occurrence[start_column:end_column])
        if (
            inserted_bp > min_bp
            and 0 < anchor < reference_length
            and first < start_column
            and end_column <= last
        ):
            answer.append(
                InsertionEvent(
                    mapping_id, support_unit, anchor, inserted_bp,
                    start_column, end_column,
                )
            )
    return answer


def reciprocal_overlap(left: DeletionEvent, right: DeletionEvent) -> float:
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start))
    if overlap == 0:
        return 0.0
    return min(overlap / left.length, overlap / right.length)


def collapse_events_by_support_unit(
    events: Sequence[DeletionEvent],
) -> List[DeletionEvent]:
    """Give every population unit one breakpoint vote within an event cluster."""

    by_unit: Dict[str, List[DeletionEvent]] = defaultdict(list)
    for event in events:
        by_unit[event.support_unit].append(event)
    collapsed = []
    for support_unit in sorted(by_unit, key=natural_name_key):
        unit_events = by_unit[support_unit]
        collapsed.append(
            DeletionEvent(
                min(item.mapping_id for item in unit_events),
                support_unit,
                round(statistics.median(item.start for item in unit_events)),
                round(statistics.median(item.end for item in unit_events)),
                round(statistics.median(item.alignment_start for item in unit_events)),
                round(statistics.median(item.alignment_end for item in unit_events)),
            )
        )
    return collapsed


def cluster_common_deletions(
    events: Sequence[DeletionEvent],
    population_unit_count: int,
    min_bp: int,
    min_count: int,
    min_frequency: float,
    event_overlap: float,
) -> List[CommonDeletion]:
    required = max(min_count, math.ceil(min_frequency * population_unit_count))
    clusters: List[List[DeletionEvent]] = []
    for event in sorted(events, key=lambda item: (item.start, item.end, item.mapping_id)):
        best_index: Optional[int] = None
        best_overlap = -1.0
        for index, cluster in enumerate(clusters):
            unit_events = collapse_events_by_support_unit(cluster)
            starts = [item.start for item in unit_events]
            ends = [item.end for item in unit_events]
            representative = DeletionEvent(
                "cluster", "cluster", round(statistics.median(starts)),
                round(statistics.median(ends)), 0, 0
            )
            overlap = reciprocal_overlap(event, representative)
            if overlap >= event_overlap and overlap > best_overlap:
                best_index, best_overlap = index, overlap
        if best_index is None:
            clusters.append([event])
        else:
            clusters[best_index].append(event)

    candidates: List[CommonDeletion] = []
    for cluster in clusters:
        unit_events = collapse_events_by_support_unit(cluster)
        if len(unit_events) < required:
            continue
        start = round(statistics.median(item.start for item in unit_events))
        end = round(statistics.median(item.end for item in unit_events))
        if end - start <= min_bp:
            continue
        candidates.append(
            CommonDeletion(
                start, end, len(unit_events), population_unit_count, len(cluster),
                round(statistics.median(item.alignment_start for item in unit_events)),
                round(statistics.median(item.alignment_end for item in unit_events)),
            )
        )

    # Ambiguous overlapping clusters cannot contribute two contradictory pairs
    # of cuts. Prefer the most-supported, most-frequent, then longest event.
    selected: List[CommonDeletion] = []
    for candidate in sorted(
        candidates,
        key=lambda item: (-item.support_count, -item.length, item.start, item.end),
    ):
        if any(min(candidate.end, kept.end) > max(candidate.start, kept.start) for kept in selected):
            continue
        selected.append(candidate)
    return sorted(selected, key=lambda item: (item.start, item.end))


def cluster_common_insertions(
    events: Sequence[InsertionEvent],
    population_unit_count: int,
    min_count: int,
    min_frequency: float,
) -> List[CommonInsertion]:
    """Select recurrent insertions, giving each population unit one vote/locus."""

    required = max(min_count, math.ceil(min_frequency * population_unit_count))
    by_locus: Dict[Tuple[int, int, int], List[InsertionEvent]] = defaultdict(list)
    for event in events:
        by_locus[(event.anchor, event.alignment_start, event.alignment_end)].append(event)

    selected = []
    for (anchor, alignment_start, alignment_end), locus_events in sorted(by_locus.items()):
        by_unit: Dict[str, List[InsertionEvent]] = defaultdict(list)
        for event in locus_events:
            by_unit[event.support_unit].append(event)
        if len(by_unit) < required:
            continue
        unit_lengths = [
            round(statistics.median(event.inserted_bp for event in unit_events))
            for unit_events in by_unit.values()
        ]
        selected.append(
            CommonInsertion(
                anchor=anchor,
                inserted_bp=round(statistics.median(unit_lengths)),
                support_count=len(by_unit),
                population_unit_count=population_unit_count,
                event_count=len(locus_events),
                alignment_start=alignment_start,
                alignment_end=alignment_end,
            )
        )
    return selected


def validate_alignment(
    alignment: Sequence[Tuple[str, str]],
    consensuses: Mapping[str, str],
    occurrences: Mapping[str, str],
) -> Tuple[Dict[str, str], Dict[str, str]]:
    aligned_references: Dict[str, str] = {}
    aligned_occurrences: Dict[str, str] = {}
    seen_occurrence = False
    for record_header, sequence in alignment:
        if record_header.startswith("reference|"):
            if seen_occurrence:
                raise ValueError("MAFFT changed input order: reference follows an occurrence")
            name = parse_consensus_name(record_header)
            source = consensuses.get(name)
            if source is None or name in aligned_references:
                raise ValueError(f"unexpected/duplicate aligned reference {name!r}")
            if sequence.replace("-", "").upper() != source.upper():
                raise ValueError(f"MAFFT changed consensus residues for {name!r}")
            aligned_references[name] = sequence
            continue
        seen_occurrence = True
        occurrence_header = record_header
        mapping_id = parse_mapping_id(occurrence_header)
        if mapping_id in aligned_occurrences:
            raise ValueError(f"duplicate aligned occurrence {mapping_id!r}")
        source = occurrences.get(mapping_id)
        if source is None:
            raise ValueError(f"unexpected aligned occurrence {mapping_id!r}")
        if sequence.replace("-", "").upper() != source.upper():
            raise ValueError(f"MAFFT changed occurrence residues for {mapping_id!r}")
        aligned_occurrences[mapping_id] = sequence
    if set(aligned_references) != set(consensuses):
        raise ValueError("MAFFT omitted one or more consensus references")
    if set(aligned_occurrences) != set(occurrences):
        raise ValueError("MAFFT omitted one or more occurrence records")
    return aligned_references, aligned_occurrences


def metadata_sizes(path: Path) -> Tuple[set[str], Dict[str, int]]:
    ids, metadata = read_existing_bundle_metadata(path)
    return set(ids), {name: value.bundle_size for name, value in metadata.items()}


def aligned_identity(left: str, right: str) -> float:
    compared = 0
    matched = 0
    for left_base, right_base in zip(left.upper(), right.upper()):
        if left_base == "-" and right_base == "-":
            continue
        compared += 1
        if left_base == right_base and left_base != "-":
            matched += 1
    return matched / compared if compared else 1.0


def prefix_non_gaps(sequence: str) -> List[int]:
    prefix = [0]
    for base in sequence:
        prefix.append(prefix[-1] + (base != "-"))
    return prefix


def write_outputs(
    records: Mapping[str, str],
    names_by_class: Mapping[str, Sequence[str]],
    common_by_class: Mapping[str, Sequence[CommonDeletion]],
    aligned_references_by_class: Mapping[str, Mapping[str, str]],
    occurrence_counts: Mapping[str, int],
    population_unit_counts: Mapping[str, int],
    all_events: Mapping[str, Sequence[DeletionEvent]],
    statuses: Mapping[str, str],
    metadata_bed: Path,
    output_dir: Path,
    args: argparse.Namespace,
    insertion_events_by_class: Optional[Mapping[str, Sequence[InsertionEvent]]] = None,
    common_insertions_by_class: Optional[Mapping[str, Sequence[CommonInsertion]]] = None,
) -> Dict[str, object]:
    insertion_events_by_class = insertion_events_by_class or {}
    common_insertions_by_class = common_insertions_by_class or {}
    existing_ids, sizes = metadata_sizes(metadata_bed)
    used_ids = set(existing_ids)
    used_ids.update(bundle_class for bundle_class in names_by_class if bundle_class.isdigit())
    next_id = max((int(value) for value in used_ids), default=0) + 1
    fragment_rows = []
    output_records: List[Tuple[str, str, str]] = []
    output_fasta = output_dir / "bundle.fa"
    for parent_class in sorted(names_by_class, key=natural_name_key):
        parent_names = sorted(names_by_class[parent_class], key=natural_name_key)
        deletions = common_by_class.get(parent_class, ())
        insertions = common_insertions_by_class.get(parent_class, ())
        parent_size = sizes.get(parent_class)
        if not deletions and not insertions:
            for name in parent_names:
                sequence = records[name]
                output_records.append((name, sequence, f"parent_consensus={name} common_indel_split=0"))
                fragment_rows.append(
                    (name, parent_class, len(sequence), 1, 0, len(sequence), 0, len(sequence),
                     name, parent_class, parent_size or "", 0)
                )
            continue

        aligned_references = aligned_references_by_class[parent_class]
        width = len(next(iter(aligned_references.values())))
        alignment_cuts = sorted(
            {point for deletion in deletions for point in
             (deletion.alignment_start, deletion.alignment_end)}
            | {insertion.alignment_start for insertion in insertions}
        )
        boundaries = (0, *alignment_cuts, width)
        backbone = [
            any(sequence[column] != "-" for sequence in aligned_references.values())
            for column in range(width)
        ]
        backbone_total = sum(backbone)
        parent_prefixes = {
            name: prefix_non_gaps(aligned_references[name]) for name in parent_names
        }
        for part, (alignment_start, alignment_end) in enumerate(
            zip(boundaries, boundaries[1:]), 1
        ):
            candidates = []
            for name in parent_names:
                aligned_segment = aligned_references[name][alignment_start:alignment_end]
                sequence = aligned_segment.replace("-", "")
                if not sequence:
                    continue
                prefix = parent_prefixes[name]
                candidates.append(
                    (name, sequence, aligned_segment, prefix[alignment_start], prefix[alignment_end])
                )
            if not candidates:
                continue
            child_id = str(next_id)
            next_id += 1
            backbone_bp = sum(backbone[alignment_start:alignment_end])
            child_size = (
                max(1, math.ceil(parent_size * backbone_bp / backbone_total))
                if parent_size and backbone_total else 100
            )
            sizes[child_id] = child_size
            representatives: List[Tuple[str, str, str, int, int]] = []
            child_by_parent: Dict[str, str] = {}
            for candidate in candidates:
                matching_index = next(
                    (
                        index for index, representative in enumerate(representatives)
                        if aligned_identity(candidate[2], representative[2]) > args.merge_identity
                    ),
                    None,
                )
                if matching_index is None:
                    representatives.append(candidate)
                    matching_index = len(representatives) - 1
                child_name = (
                    f"bundle_{child_id}" if matching_index == 0
                    else f"bundle_{child_id}.allele_{matching_index + 1}"
                )
                child_by_parent[candidate[0]] = child_name
            for index, representative in enumerate(representatives):
                child_name = (
                    f"bundle_{child_id}" if index == 0 else f"bundle_{child_id}.allele_{index + 1}"
                )
                output_records.append(
                    (child_name, representative[1],
                     f"parent_class={parent_class} alignment_interval={alignment_start}-{alignment_end} "
                     "common_indel_split=1")
                )
            for name, _sequence, _aligned, start, end in candidates:
                fragment_rows.append(
                    (name, parent_class, len(records[name]), part, start, end,
                     alignment_start, alignment_end, child_by_parent[name], child_id,
                     child_size, 1)
                )

    if any(common_by_class.values()) or any(common_insertions_by_class.values()):
        with output_fasta.open("w", encoding="utf-8") as handle:
            for name, sequence, description in output_records:
                handle.write(f">{name} {description}\n")
                handle.writelines(wrap_fasta(sequence))
    else:
        # Byte-identical pass-through lets final cleanup safely remove an older
        # redundant mapping-split copy when this stage found no common indel.
        shutil.copyfile(args.bundle_fasta, output_fasta)

    metadata_output = output_dir / "bundle_metadata.bed"
    with metadata_output.open("w", encoding="utf-8") as handle:
        handle.write("# ID/size metadata only; synthetic coordinates, NOT genomic mappings\n")
        for bundle_id in sorted(sizes, key=int):
            size = sizes[bundle_id]
            handle.write(f"__metadata__\t0\t1\t{bundle_id}:{size}:0:0:{size - 1}:U\n")

    with (output_dir / "bundle_fragments.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("parent_consensus", "parent_class", "parent_length", "part", "parent_start",
             "parent_end", "alignment_start_column", "alignment_end_column", "consensus_name",
             "numeric_bundle_id", "bundle_size", "is_split")
        )
        writer.writerows(fragment_rows)

    with (output_dir / "deletion_events.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("bundle_class", "mapping_id", "support_unit", "deletion_start", "deletion_end", "deletion_bp",
             "alignment_start_column", "alignment_end_column")
        )
        for bundle_class in sorted(all_events, key=natural_name_key):
            for event in all_events[bundle_class]:
                writer.writerow(
                    (bundle_class, event.mapping_id, event.support_unit, event.start, event.end, event.length,
                     event.alignment_start, event.alignment_end)
                )

    with (output_dir / "common_deletions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("bundle_class", "deletion_start", "deletion_end", "deletion_bp",
             "alignment_start_column", "alignment_end_column", "supporting_unit_count",
             "population_unit_count", "support_frequency", "event_count", "selected")
        )
        for bundle_class in sorted(names_by_class, key=natural_name_key):
            for deletion in common_by_class.get(bundle_class, ()):
                writer.writerow(
                    (bundle_class, deletion.start, deletion.end, deletion.length,
                     deletion.alignment_start, deletion.alignment_end,
                     deletion.support_count, deletion.population_unit_count,
                     f"{deletion.support_frequency:.10g}", deletion.event_count, 1)
                )

    with (output_dir / "insertion_events.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("bundle_class", "mapping_id", "support_unit", "insertion_anchor",
             "inserted_bp", "alignment_start_column", "alignment_end_column")
        )
        for bundle_class in sorted(insertion_events_by_class, key=natural_name_key):
            for event in insertion_events_by_class[bundle_class]:
                writer.writerow(
                    (bundle_class, event.mapping_id, event.support_unit, event.anchor,
                     event.inserted_bp, event.alignment_start, event.alignment_end)
                )

    with (output_dir / "common_insertions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("bundle_class", "insertion_anchor", "median_inserted_bp",
             "alignment_start_column", "alignment_end_column", "supporting_unit_count",
             "population_unit_count", "support_frequency", "event_count", "selected")
        )
        for bundle_class in sorted(names_by_class, key=natural_name_key):
            for insertion in common_insertions_by_class.get(bundle_class, ()):
                writer.writerow(
                    (bundle_class, insertion.anchor, insertion.inserted_bp,
                     insertion.alignment_start, insertion.alignment_end,
                     insertion.support_count, insertion.population_unit_count,
                     f"{insertion.support_frequency:.10g}", insertion.event_count, 1)
                )

    with (output_dir / "bundle_split_stats.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("bundle_class", "raw_consensus_names", "max_bundle_length", "occurrence_count",
             "population_unit_count", "internal_deletion_event_count", "common_deletion_count",
             "insertion_event_count", "common_insertion_count",
             "cut_positions", "alignment_cut_columns", "status")
        )
        for bundle_class in sorted(names_by_class, key=natural_name_key):
            parent_names = names_by_class[bundle_class]
            deletions = common_by_class.get(bundle_class, ())
            insertions = common_insertions_by_class.get(bundle_class, ())
            cuts = sorted(
                {point for deletion in deletions for point in (deletion.start, deletion.end)}
                | {insertion.anchor for insertion in insertions}
            )
            alignment_cuts = sorted({point for deletion in deletions for point in
                                     (deletion.alignment_start, deletion.alignment_end)}
                                    | {insertion.alignment_start for insertion in insertions})
            writer.writerow(
                (bundle_class, ",".join(parent_names), max(len(records[name]) for name in parent_names),
                 occurrence_counts.get(bundle_class, 0), population_unit_counts.get(bundle_class, 0),
                 len(all_events.get(bundle_class, ())), len(deletions),
                 len(insertion_events_by_class.get(bundle_class, ())), len(insertions),
                 ",".join(map(str, cuts)),
                 ",".join(map(str, alignment_cuts)),
                 "split" if cuts else statuses.get(bundle_class, "retained"))
            )

    split_consensus_count = sum(
        bool(common_by_class.get(name) or common_insertions_by_class.get(name))
        for name in names_by_class
    )
    summary = {
        "completed": True,
        "inputs": {
            "bundle_fasta": str(args.bundle_fasta.resolve()),
            "genome_fasta": str(args.genome_fasta.resolve()),
            "mapping_tsv": str(args.mapping_tsv.resolve()),
            "bundle_metadata_bed": str(metadata_bed.resolve()),
            "sample_map": str(args.sample_map.resolve()) if args.sample_map is not None else None,
        },
        "parameters": {
            "min_deletion_bp": args.min_deletion_bp,
            "min_insertion_bp": args.min_deletion_bp,
            "length_rule": "strictly_greater_than",
            "min_deletion_count": args.min_deletion_count,
            "min_insertion_count": args.min_deletion_count,
            "min_deletion_frequency": args.min_deletion_frequency,
            "min_insertion_frequency": args.min_deletion_frequency,
            "support_unit": "sample" if args.sample_map is not None else "population_path",
            "event_overlap_fraction": args.event_overlap_fraction,
            "same_child_allele_merge_identity": args.merge_identity,
            "terminal_gap_policy": "ignored",
            "deleted_interval_policy": "retained_as_independent_child_bundle",
            "inserted_interval_policy": "not_present_in_consensus; single_cut_at_anchor",
        },
        "counts": {
            "input_consensus_count": len(records),
            "input_bundle_class_count": len(names_by_class),
            "aligned_bundle_class_count": sum(status == "aligned" for status in statuses.values()),
            "occurrence_count": sum(occurrence_counts.values()),
            "internal_deletion_event_count": sum(len(value) for value in all_events.values()),
            "population_unit_count_sum_across_classes": sum(population_unit_counts.values()),
            "common_deletion_count": sum(len(value) for value in common_by_class.values()),
            "internal_insertion_event_count": sum(
                len(value) for value in insertion_events_by_class.values()
            ),
            "common_insertion_count": sum(
                len(value) for value in common_insertions_by_class.values()
            ),
            "split_consensus_count": split_consensus_count,
            "output_consensus_count": len(output_records),
        },
        "outputs": {
            "bundle_fasta": str(output_fasta.resolve()),
            "bundle_metadata_bed": str(metadata_output.resolve()),
            "bundle_fragments": str((output_dir / "bundle_fragments.tsv").resolve()),
            "deletion_events": str((output_dir / "deletion_events.tsv").resolve()),
            "common_deletions": str((output_dir / "common_deletions.tsv").resolve()),
            "insertion_events": str((output_dir / "insertion_events.tsv").resolve()),
            "common_insertions": str((output_dir / "common_insertions.tsv").resolve()),
            "bundle_split_stats": str((output_dir / "bundle_split_stats.tsv").resolve()),
        },
    }
    (output_dir / "common_deletion_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Split final bundle consensuses at recurrent internal indels in a population MSA."
    )
    parser.add_argument("bundle_fasta", type=Path)
    parser.add_argument("genome_fasta", type=Path)
    parser.add_argument("mapping_tsv", type=Path)
    parser.add_argument("--bundle-metadata-bed", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--mafft-extra", default="")
    parser.add_argument("--mafft-params", default="--auto",
                        help="MAFFT parameters replacing the default --auto strategy")
    parser.add_argument("--sample-map", type=Path)
    parser.add_argument(
        "--min-deletion-bp", type=int, default=1000,
        help=("internal deletion or insertion length must be strictly greater than "
              "this value (default: 1000)"),
    )
    parser.add_argument("--min-deletion-count", type=int, default=5)
    parser.add_argument("--min-deletion-frequency", type=float, default=0.01)
    parser.add_argument("--event-overlap-fraction", type=float, default=0.80)
    parser.add_argument(
        "--merge-identity", type=float, default=0.90,
        help="aligned child alleles with identity strictly above this value share one sequence",
    )
    parser.add_argument("--keep-intermediates", action="store_true")
    return parser


def resolve_executable(value: str) -> Optional[str]:
    candidate = Path(value)
    if candidate.parent != Path(".") or candidate.is_absolute():
        return str(candidate) if candidate.is_file() and candidate.stat().st_mode & 0o111 else None
    return shutil.which(value)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.threads == 0 or args.threads < -1:
        parser.error("--threads must be -1 or a positive integer")
    if args.min_deletion_bp < 0:
        parser.error("--min-deletion-bp must be non-negative")
    if args.min_deletion_count < 1:
        parser.error("--min-deletion-count must be positive")
    for name in ("min_deletion_frequency", "event_overlap_fraction", "merge_identity"):
        value = getattr(args, name)
        if not 0 <= value <= 1 or (name != "merge_identity" and value == 0):
            interval = "[0, 1]" if name == "merge_identity" else "(0, 1]"
            parser.error(f"--{name.replace('_', '-')} must be in {interval}")
    for path in (args.bundle_fasta, args.genome_fasta, args.mapping_tsv, args.bundle_metadata_bed):
        if not path.is_file():
            parser.error(f"input file not found: {path}")
    if args.sample_map is not None and not args.sample_map.is_file():
        parser.error(f"input file not found: {args.sample_map}")

    try:
        records = read_fasta(args.bundle_fasta)
        occurrences_by_name = read_mapping_table(args.mapping_tsv, tuple(records))
        extracted = extract_occurrences(args.genome_fasta, occurrences_by_name)
        sample_by_path = read_sample_map(args.sample_map)
        names_by_class: Dict[str, List[str]] = defaultdict(list)
        for name in records:
            names_by_class[normalize_bundle_name(name, {})].append(name)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        source_dir = args.output_dir / "sequences"
        alignment_dir = args.output_dir / "alignments"
        log_dir = args.output_dir / "logs"
        for directory in (source_dir, alignment_dir, log_dir):
            directory.mkdir(parents=True, exist_ok=True)

        mafft: Optional[str] = None
        statuses: Dict[str, str] = {}
        all_events: Dict[str, List[DeletionEvent]] = {}
        common_by_class: Dict[str, List[CommonDeletion]] = {}
        all_insertions: Dict[str, List[InsertionEvent]] = {}
        common_insertions_by_class: Dict[str, List[CommonInsertion]] = {}
        aligned_references_by_class: Dict[str, Dict[str, str]] = {}
        occurrence_counts: Dict[str, int] = {}
        population_unit_counts: Dict[str, int] = {}
        for bundle_class in sorted(names_by_class, key=natural_name_key):
            parent_names = sorted(names_by_class[bundle_class], key=natural_name_key)
            occurrences = [
                item for name in parent_names for item in extracted.get(name, ())
            ]
            support_units = set()
            for occurrence, _sequence in occurrences:
                if args.sample_map is not None and occurrence.path_id not in sample_by_path:
                    raise ValueError(
                        f"{args.sample_map}: no sample ID for mapped path {occurrence.path_id!r}"
                    )
                support_units.add(sample_by_path.get(occurrence.path_id, occurrence.path_id))
            population_unit_count = len(support_units)
            occurrence_counts[bundle_class] = len(occurrences)
            population_unit_counts[bundle_class] = population_unit_count
            required = max(
                args.min_deletion_count,
                math.ceil(args.min_deletion_frequency * population_unit_count),
            )
            if max(len(records[name]) for name in parent_names) < 100:
                statuses[bundle_class] = "bundle_below_100bp_no_mafft"
                continue
            if max(len(records[name]) for name in parent_names) <= args.min_deletion_bp:
                statuses[bundle_class] = "bundle_not_long_enough"
                continue
            if population_unit_count < required:
                statuses[bundle_class] = "insufficient_population_support"
                continue
            if mafft is None:
                mafft = resolve_executable(args.mafft)
                if mafft is None:
                    raise ValueError(
                        f"MAFFT executable {args.mafft!r} was not found for eligible final bundles"
                    )
            stem = safe_stem(bundle_class)
            source = source_dir / f"{stem}.fa"
            alignment = alignment_dir / f"{stem}.aln.fa"
            log = log_dir / f"{stem}.mafft.log"
            parent_records = [(name, records[name]) for name in parent_names]
            write_alignment_input(source, parent_records, occurrences)
            run_mafft(mafft, source, alignment, log, args.threads,
                      shlex.split(args.mafft_extra), shlex.split(args.mafft_params))
            occurrence_sources = {item.mapping_id: fragment for item, fragment in occurrences}
            occurrence_units = {
                item.mapping_id: sample_by_path.get(item.path_id, item.path_id)
                for item, _ in occurrences
            }
            aligned_references, aligned_occurrences = validate_alignment(
                read_alignment(alignment), dict(parent_records), occurrence_sources
            )
            width = len(next(iter(aligned_references.values())))
            reference_backbone = "".join(
                "N" if any(sequence[column] != "-" for sequence in aligned_references.values())
                else "-"
                for column in range(width)
            )
            events = [
                event
                for mapping_id, aligned in aligned_occurrences.items()
                for event in deletion_events(
                    reference_backbone, aligned, mapping_id, occurrence_units[mapping_id],
                    args.min_deletion_bp
                )
            ]
            common = cluster_common_deletions(
                events, population_unit_count, args.min_deletion_bp,
                args.min_deletion_count, args.min_deletion_frequency,
                args.event_overlap_fraction,
            )
            insertions = [
                event
                for mapping_id, aligned in aligned_occurrences.items()
                for event in insertion_events(
                    reference_backbone, aligned, mapping_id, occurrence_units[mapping_id],
                    args.min_deletion_bp,
                )
            ]
            common_insertions = cluster_common_insertions(
                insertions, population_unit_count, args.min_deletion_count,
                args.min_deletion_frequency,
            )
            all_events[bundle_class] = events
            common_by_class[bundle_class] = common
            all_insertions[bundle_class] = insertions
            common_insertions_by_class[bundle_class] = common_insertions
            aligned_references_by_class[bundle_class] = aligned_references
            statuses[bundle_class] = "aligned"

        summary = write_outputs(
            records, names_by_class, common_by_class, aligned_references_by_class,
            occurrence_counts, population_unit_counts, all_events, statuses,
            args.bundle_metadata_bed, args.output_dir, args,
            insertion_events_by_class=all_insertions,
            common_insertions_by_class=common_insertions_by_class,
        )
        if not args.keep_intermediates:
            for directory in (source_dir, alignment_dir, log_dir):
                shutil.rmtree(directory)
            summary["intermediate_retention"] = {
                "kept": False,
                "removed": ["sequences", "alignments", "logs"],
                "failure_policy": "inputs, partial alignments and MAFFT logs remain when the stage fails",
            }
            (args.output_dir / "common_deletion_summary.json").write_text(
                json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
