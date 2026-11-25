#!/usr/bin/env python3
"""Build bundle consensuses after grouping haplotypes by signed bundle order."""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple
from urllib.parse import unquote

try:
    from .build_bundle_consensus import (
        bundle_sort_key,
        iter_fasta,
        read_alignment_manifest,
        wrap_fasta,
    )
except ImportError:
    from build_bundle_consensus import (
        bundle_sort_key,
        iter_fasta,
        read_alignment_manifest,
        wrap_fasta,
    )


@dataclass(frozen=True)
class BedOccurrence:
    row_number: int
    path_id: str
    start: int
    end: int
    bundle_id: str
    direction: str

    @property
    def key(self) -> Tuple[str, int, int, str, str]:
        return self.path_id, self.start, self.end, self.bundle_id, self.direction


@dataclass(frozen=True)
class PathOrderGroup:
    group_id: str
    rank: int
    signature: Tuple[str, ...]
    path_ids: Tuple[str, ...]

    @property
    def haplotype_count(self) -> int:
        return len(self.path_ids)


@dataclass
class CandidateConsensus:
    candidate_id: str
    group_id: str
    group_rank: int
    haplotype_count: int
    position_1based: int
    bundle_id: str
    direction: str
    record_headers: Tuple[str, ...]
    sequences: Tuple[str, ...]
    aligned_profile: str
    consensus_sequence: str
    retained_column_count: int
    dropped_column_count: int
    ambiguous_base_count: int
    mean_retained_occupancy: float


@dataclass(frozen=True)
class FinalConsensus:
    output_id: str
    bundle_id: str
    allele_index: int
    candidate_ids: Tuple[str, ...]
    group_ids: Tuple[str, ...]
    unique_path_count: int
    source_occurrence_count: int
    consensus_sequence: str
    ambiguous_base_count: int
    output_path: Path


def normalize_direction(value: str) -> str:
    if value in {"0", "+"}:
        return "+"
    if value in {"1", "-"}:
        return "-"
    raise ValueError(f"unsupported PGR direction: {value!r}")


def read_bundle_bed(path: Path) -> Dict[str, List[BedOccurrence]]:
    by_path: DefaultDict[str, List[BedOccurrence]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 4:
                fields = stripped.split()
            if len(fields) < 4:
                raise ValueError(f"{path}:{row_number}: expected at least four columns")
            annotation = fields[3].split(":")
            if len(annotation) < 3:
                raise ValueError(f"{path}:{row_number}: invalid PGR annotation")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{row_number}: invalid interval {start}-{end}")
            by_path[fields[0]].append(
                BedOccurrence(
                    row_number=row_number,
                    path_id=fields[0],
                    start=start,
                    end=end,
                    bundle_id=annotation[0],
                    direction=normalize_direction(annotation[2]),
                )
            )
    if not by_path:
        raise ValueError(f"{path}: no PGR bundle occurrences")
    for occurrences in by_path.values():
        occurrences.sort(key=lambda item: (item.start, item.end, item.row_number))
    return dict(by_path)


def group_paths_by_order(
    occurrences_by_path: Mapping[str, Sequence[BedOccurrence]],
) -> Tuple[List[PathOrderGroup], Dict[str, PathOrderGroup]]:
    by_signature: DefaultDict[Tuple[str, ...], List[str]] = defaultdict(list)
    for path_id, occurrences in occurrences_by_path.items():
        ordered_occurrences = sorted(
            occurrences, key=lambda item: (item.start, item.end, item.row_number)
        )
        signature = tuple(
            f"{occurrence.direction}{occurrence.bundle_id}"
            for occurrence in ordered_occurrences
        )
        by_signature[signature].append(path_id)
    ordered = sorted(
        by_signature.items(),
        key=lambda item: (-len(item[1]), item[0]),
    )
    groups: List[PathOrderGroup] = []
    path_to_group: Dict[str, PathOrderGroup] = {}
    for rank, (signature, paths) in enumerate(ordered, 1):
        group = PathOrderGroup(
            group_id=f"G{rank:04d}",
            rank=rank,
            signature=signature,
            path_ids=tuple(sorted(paths)),
        )
        groups.append(group)
        for path_id in group.path_ids:
            path_to_group[path_id] = group
    return groups, path_to_group


def select_group_bundle_ids(
    groups: Sequence[PathOrderGroup],
    all_haplotype_groups: bool = False,
    supplement_missing_bundles: bool = False,
) -> Dict[str, Set[str]]:
    """Select IDs by support rank, never importing lower-rank copies of a seen ID.

    Each ID contributed by a group retains every position in that group, so
    duplicated positions can still produce separate consensus alleles.
    """
    if all_haplotype_groups and supplement_missing_bundles:
        raise ValueError(
            "--all-haplotype-groups and --supplement-missing-bundles are mutually exclusive"
        )
    selected: Dict[str, Set[str]] = {}
    covered: Set[str] = set()
    for group in sorted(groups, key=lambda item: item.rank):
        bundle_ids = {token[1:] for token in group.signature}
        if all_haplotype_groups:
            contributed = bundle_ids
        elif supplement_missing_bundles:
            contributed = bundle_ids - covered
        else:
            contributed = bundle_ids if not selected else set()
        if contributed:
            selected[group.group_id] = contributed
            covered.update(contributed)
        if not all_haplotype_groups and not supplement_missing_bundles:
            break
    return selected


def write_bundle_source_groups(
    path: Path,
    groups: Sequence[PathOrderGroup],
    selected_bundle_ids_by_group: Mapping[str, Set[str]],
) -> None:
    """Report exactly the selected IDs/positions, not all IDs in selected groups."""
    input_bundle_ids = {token[1:] for group in groups for token in group.signature}
    selected_bundle_ids = {
        bundle_id for bundle_ids in selected_bundle_ids_by_group.values()
        for bundle_id in bundle_ids
    }
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "bundle_id", "group_id", "group_rank", "haplotype_count",
                "positions_1based", "signed_bundle_order", "path_ids",
                "selected_occurrence_count", "input_bundle_count",
                "selected_bundle_count", "covers_all_input_bundles",
            )
        )
        for group in sorted(groups, key=lambda item: item.rank):
            for bundle_id in sorted(
                selected_bundle_ids_by_group.get(group.group_id, set()),
                key=bundle_sort_key,
            ):
                positions = [
                    position for position, token in enumerate(group.signature, 1)
                    if token[1:] == bundle_id
                ]
                writer.writerow(
                    (
                        bundle_id, group.group_id, group.rank, group.haplotype_count,
                        ",".join(map(str, positions)), ",".join(group.signature),
                        ",".join(group.path_ids),
                        len(positions) * group.haplotype_count,
                        len(input_bundle_ids), len(selected_bundle_ids),
                        int(input_bundle_ids <= selected_bundle_ids),
                    )
                )


def validate_bundle_consensus_coverage(
    consensus_fasta: Path,
    bundle_bed: Optional[Path] = None,
    source_groups_tsv: Optional[Path] = None,
) -> None:
    """Reject stale rank-1-only FASTA before a resumed refinement drops missing IDs.

    Use the final BED when available, otherwise the source audit emitted by a
    successful supplemented consensus run. Exactly one source is required.
    """
    if (bundle_bed is None) == (source_groups_tsv is None):
        raise ValueError("provide exactly one of bundle_bed or source_groups_tsv")
    if bundle_bed is not None:
        expected = {
            occurrence.bundle_id
            for occurrences in read_bundle_bed(bundle_bed).values()
            for occurrence in occurrences
        }
    else:
        with source_groups_tsv.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "bundle_id" not in reader.fieldnames:
                raise ValueError(f"{source_groups_tsv}: expected bundle_id column")
            rows = list(reader)
            if (
                "covers_all_input_bundles" not in reader.fieldnames
                or any(row["covers_all_input_bundles"] != "1" for row in rows)
            ):
                raise ValueError(
                    f"{source_groups_tsv}: source audit does not confirm coverage of all "
                    "input bundle IDs; rerun extraction/alignment and consensus with "
                    "--supplement-missing-bundles before refining"
                )
            expected = {row["bundle_id"] for row in rows if row["bundle_id"]}
        if not expected:
            raise ValueError(f"{source_groups_tsv}: no bundle source rows")

    def canonical_id(bundle_id: str) -> str:
        return str(int(bundle_id)) if bundle_id.isdigit() else bundle_id

    expected = {canonical_id(bundle_id) for bundle_id in expected}
    observed: Set[str] = set()
    for header, sequence in iter_fasta(consensus_fasta):
        if not sequence:
            raise ValueError(f"{consensus_fasta}: empty consensus sequence: {header}")
        name = header.split()[0]
        if name.startswith("bundle_"):
            name = name[len("bundle_"):]
        name = re.sub(r"\.allele_\d+$", "", name)
        observed.add(canonical_id(name))
    missing = sorted(expected - observed, key=bundle_sort_key)
    if missing:
        raise ValueError(
            f"{consensus_fasta}: lacks consensus for {len(missing)} bundle ID(s): "
            f"{', '.join(missing[:10])}; rerun extraction/alignment and consensus "
            "with --supplement-missing-bundles before refining"
        )


def parse_alignment_header(header: str) -> Tuple[str, int, int, str, str]:
    values: Dict[str, str] = {}
    for field in header.split("|"):
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = unquote(value)
    required = {"bundle", "path", "start", "end", "strand"}
    missing = sorted(required - values.keys())
    if missing:
        raise ValueError(f"alignment header lacks {', '.join(missing)}: {header}")
    return (
        values["path"],
        int(values["start"]),
        int(values["end"]),
        values["bundle"],
        values["strand"],
    )


def aligned_consensus_profile(
    sequences: Sequence[str],
    min_column_occupancy: float,
    min_base_frequency: float,
) -> Tuple[str, str, int, int, float]:
    """Return fixed-width profile plus its gap-removed consensus sequence."""
    if not sequences:
        raise ValueError("cannot build a consensus from zero sequences")
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("candidate sequences do not have a common MSA length")
    profile: List[str] = []
    sequence: List[str] = []
    retained_occupancies: List[float] = []
    dropped = 0
    for column in range(len(sequences[0])):
        counts = Counter(item[column].upper() for item in sequences)
        canonical = {base: counts[base] for base in "ACGT" if counts[base]}
        canonical_count = sum(canonical.values())
        occupancy = canonical_count / len(sequences)
        if canonical_count == 0 or occupancy < min_column_occupancy:
            profile.append("-")
            dropped += 1
            continue
        retained_occupancies.append(occupancy)
        largest = max(canonical.values())
        winners = sorted(base for base, count in canonical.items() if count == largest)
        base = (
            winners[0]
            if len(winners) == 1 and largest / canonical_count >= min_base_frequency
            else "N"
        )
        profile.append(base)
        sequence.append(base)
    if not sequence:
        raise ValueError("no MSA column passed the candidate consensus threshold")
    return (
        "".join(profile),
        "".join(sequence),
        len(sequence),
        dropped,
        float(statistics.fmean(retained_occupancies)),
    )


def build_candidate_definitions(
    groups: Sequence[PathOrderGroup],
    occurrences_by_path: Mapping[str, Sequence[BedOccurrence]],
    selected_bundle_ids_by_group: Optional[Mapping[str, Set[str]]] = None,
) -> Tuple[
    Dict[Tuple[str, int], Tuple[str, PathOrderGroup, int, str, str]],
    Dict[Tuple[str, int, int, str, str], Tuple[str, int]],
]:
    definitions: Dict[
        Tuple[str, int], Tuple[str, PathOrderGroup, int, str, str]
    ] = {}
    occurrence_to_candidate: Dict[
        Tuple[str, int, int, str, str], Tuple[str, int]
    ] = {}
    serial = 0
    for group in groups:
        representative = occurrences_by_path[group.path_ids[0]]
        for position, occurrence in enumerate(representative, 1):
            if (
                selected_bundle_ids_by_group is not None
                and occurrence.bundle_id
                not in selected_bundle_ids_by_group.get(group.group_id, set())
            ):
                continue
            serial += 1
            key = group.group_id, position
            candidate_id = f"C{serial:06d}"
            definitions[key] = (
                candidate_id,
                group,
                position,
                occurrence.bundle_id,
                occurrence.direction,
            )
            for path_id in group.path_ids:
                path_occurrence = occurrences_by_path[path_id][position - 1]
                if (
                    path_occurrence.bundle_id != occurrence.bundle_id
                    or path_occurrence.direction != occurrence.direction
                ):
                    raise ValueError(
                        f"path-order group {group.group_id} is internally inconsistent"
                    )
                occurrence_to_candidate[path_occurrence.key] = key
    return definitions, occurrence_to_candidate


def build_candidates(
    bundle_bed: Path,
    alignment_manifest: Path,
    min_column_occupancy: float,
    min_base_frequency: float,
    all_haplotype_groups: bool = False,
    supplement_missing_bundles: bool = False,
) -> Tuple[List[PathOrderGroup], List[CandidateConsensus]]:
    occurrences_by_path = read_bundle_bed(bundle_bed)
    groups, _ = group_paths_by_order(occurrences_by_path)
    selected_bundle_ids_by_group = select_group_bundle_ids(
        groups, all_haplotype_groups, supplement_missing_bundles
    )
    selected_groups = [
        group for group in groups if group.group_id in selected_bundle_ids_by_group
    ]
    selected_path_ids = {
        path_id for group in selected_groups for path_id in group.path_ids
    }
    definitions, occurrence_to_candidate = build_candidate_definitions(
        selected_groups, occurrences_by_path, selected_bundle_ids_by_group
    )
    known_occurrences = {
        occurrence.key
        for occurrences in occurrences_by_path.values() for occurrence in occurrences
    }
    sequences_by_candidate: DefaultDict[Tuple[str, int], List[Tuple[str, str]]] = (
        defaultdict(list)
    )
    seen_occurrences = set()
    for alignment_input in read_alignment_manifest(alignment_manifest):
        records = list(iter_fasta(alignment_input.alignment_path))
        if not records:
            raise ValueError(f"{alignment_input.alignment_path}: empty alignment")
        lengths = {len(sequence) for _, sequence in records}
        if len(lengths) != 1:
            raise ValueError(
                f"{alignment_input.alignment_path}: sequences have unequal alignment lengths"
            )
        for header, sequence in records:
            occurrence_key = parse_alignment_header(header)
            if occurrence_key[3] != alignment_input.bundle_id:
                raise ValueError(
                    f"alignment manifest/header bundle mismatch: "
                    f"{alignment_input.bundle_id} versus {occurrence_key[3]}"
                )
            candidate_key = occurrence_to_candidate.get(occurrence_key)
            if candidate_key is None:
                # A selected lower-rank group can also carry already-covered IDs;
                # those copies must not enter the higher-rank consensus.
                if occurrence_key in known_occurrences:
                    continue
                if occurrence_key[0] not in selected_path_ids:
                    continue
                raise ValueError(
                    "an occurrence from a selected haplotype group is absent from "
                    f"the final BED candidate map: {occurrence_key}"
                )
            if occurrence_key in seen_occurrences:
                raise ValueError(f"duplicate aligned occurrence: {occurrence_key}")
            seen_occurrences.add(occurrence_key)
            sequences_by_candidate[candidate_key].append((header, sequence.upper()))
    missing = sorted(set(occurrence_to_candidate) - seen_occurrences)
    if missing:
        raise ValueError(
            f"alignment manifest lacks {len(missing)} final BED occurrences; "
            f"example: {missing[0]}; rerun extraction/alignment with "
            "--supplement-missing-bundles, without --resume for stale alignments"
        )

    candidates: List[CandidateConsensus] = []
    for key, definition in definitions.items():
        candidate_id, group, position, bundle_id, direction = definition
        records = sequences_by_candidate[key]
        if len(records) != group.haplotype_count:
            raise ValueError(
                f"{candidate_id}: expected {group.haplotype_count} sequences, "
                f"observed {len(records)}"
            )
        headers = tuple(header for header, _ in records)
        sequences = tuple(sequence for _, sequence in records)
        profile, consensus, retained, dropped, mean_occupancy = (
            aligned_consensus_profile(
                sequences,
                min_column_occupancy,
                min_base_frequency,
            )
        )
        candidates.append(
            CandidateConsensus(
                candidate_id=candidate_id,
                group_id=group.group_id,
                group_rank=group.rank,
                haplotype_count=group.haplotype_count,
                position_1based=position,
                bundle_id=bundle_id,
                direction=direction,
                record_headers=headers,
                sequences=sequences,
                aligned_profile=profile,
                consensus_sequence=consensus,
                retained_column_count=retained,
                dropped_column_count=dropped,
                ambiguous_base_count=sum(base == "N" for base in consensus),
                mean_retained_occupancy=mean_occupancy,
            )
        )
    candidates.sort(
        key=lambda item: (
            item.group_rank,
            item.position_1based,
            bundle_sort_key(item.bundle_id),
        )
    )
    return groups, candidates


def aligned_profile_identity(first: str, second: str) -> Tuple[float, int, int]:
    if len(first) != len(second):
        raise ValueError("candidate profiles from one bundle have unequal MSA lengths")
    matches = 0
    compared = 0
    for left, right in zip(first, second):
        if left == "-" and right == "-":
            continue
        compared += 1
        if left == right and left in "ACGT":
            matches += 1
    return (matches / compared if compared else 0.0), matches, compared


def cluster_candidates(
    candidates: Sequence[CandidateConsensus],
    merge_identity: float,
) -> Tuple[
    Dict[str, List[List[CandidateConsensus]]],
    Dict[Tuple[str, str], Tuple[float, int, int]],
]:
    by_bundle: DefaultDict[str, List[CandidateConsensus]] = defaultdict(list)
    for candidate in candidates:
        by_bundle[candidate.bundle_id].append(candidate)
    identities: Dict[Tuple[str, str], Tuple[float, int, int]] = {}
    clusters_by_bundle: Dict[str, List[List[CandidateConsensus]]] = {}
    for bundle_id in sorted(by_bundle, key=bundle_sort_key):
        ordered = sorted(
            by_bundle[bundle_id],
            key=lambda item: (-item.haplotype_count, item.group_rank, item.position_1based),
        )
        for index, first in enumerate(ordered):
            for second in ordered[index + 1 :]:
                key = tuple(sorted((first.candidate_id, second.candidate_id)))
                identities[key] = aligned_profile_identity(
                    first.aligned_profile, second.aligned_profile
                )
        clusters: List[List[CandidateConsensus]] = []
        for candidate in ordered:
            destination: Optional[List[CandidateConsensus]] = None
            for cluster in clusters:
                if all(
                    identities[
                        tuple(sorted((candidate.candidate_id, member.candidate_id)))
                    ][0]
                    > merge_identity
                    for member in cluster
                ):
                    destination = cluster
                    break
            if destination is None:
                clusters.append([candidate])
            else:
                destination.append(candidate)
        clusters_by_bundle[bundle_id] = clusters
    return clusters_by_bundle, identities


def final_consensuses(
    clusters_by_bundle: Mapping[str, Sequence[Sequence[CandidateConsensus]]],
    output_dir: Path,
    min_column_occupancy: float,
    min_base_frequency: float,
    line_width: int,
) -> List[FinalConsensus]:
    per_bundle_dir = output_dir / "per_bundle"
    per_bundle_dir.mkdir(parents=True, exist_ok=True)
    results: List[FinalConsensus] = []
    for bundle_id in sorted(clusters_by_bundle, key=bundle_sort_key):
        clusters = clusters_by_bundle[bundle_id]
        multiple = len(clusters) > 1
        for allele_index, cluster in enumerate(clusters, 1):
            records: Dict[str, str] = {}
            for candidate in cluster:
                for header, sequence in zip(candidate.record_headers, candidate.sequences):
                    previous = records.setdefault(header, sequence)
                    if previous != sequence:
                        raise ValueError(f"alignment header maps to two sequences: {header}")
            _, consensus, _, _, _ = aligned_consensus_profile(
                tuple(records.values()),
                min_column_occupancy,
                min_base_frequency,
            )
            output_id = (
                f"bundle_{bundle_id}.allele_{allele_index}"
                if multiple
                else f"bundle_{bundle_id}"
            )
            output_path = per_bundle_dir / f"{output_id}.consensus.fa"
            unique_paths = {
                parse_alignment_header(header)[0] for header in records
            }
            result = FinalConsensus(
                output_id=output_id,
                bundle_id=bundle_id,
                allele_index=allele_index,
                candidate_ids=tuple(candidate.candidate_id for candidate in cluster),
                group_ids=tuple(sorted({candidate.group_id for candidate in cluster})),
                unique_path_count=len(unique_paths),
                source_occurrence_count=len(records),
                consensus_sequence=consensus,
                ambiguous_base_count=sum(base == "N" for base in consensus),
                output_path=output_path,
            )
            with output_path.open("w", encoding="utf-8") as handle:
                handle.write(
                    f">{output_id} bundle_id={bundle_id} allele={allele_index} "
                    f"candidate_count={len(cluster)} unique_path_count={len(unique_paths)} "
                    f"source_occurrence_count={len(records)} "
                    f"consensus_length={len(consensus)}\n"
                )
                handle.writelines(wrap_fasta(consensus, line_width))
            results.append(result)
    return results


def write_path_groups(
    path: Path,
    groups: Sequence[PathOrderGroup],
    selected_group_ids: Sequence[str],
) -> None:
    selected = set(selected_group_ids)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "group_id", "rank", "haplotype_count", "selected_for_consensus",
                "signed_bundle_order", "path_ids",
            )
        )
        for group in groups:
            writer.writerow(
                (
                    group.group_id,
                    group.rank,
                    group.haplotype_count,
                    int(group.group_id in selected),
                    ",".join(group.signature),
                    ",".join(group.path_ids),
                )
            )


def write_candidate_outputs(
    output_dir: Path,
    candidates: Sequence[CandidateConsensus],
    line_width: int,
) -> None:
    with (output_dir / "candidate_consensus.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "candidate_id", "group_id", "group_rank", "haplotype_count",
                "position_1based", "bundle_id", "direction", "sequence_count",
                "consensus_length", "retained_column_count", "dropped_column_count",
                "ambiguous_base_count", "mean_retained_column_occupancy",
            )
        )
        for candidate in candidates:
            writer.writerow(
                (
                    candidate.candidate_id,
                    candidate.group_id,
                    candidate.group_rank,
                    candidate.haplotype_count,
                    candidate.position_1based,
                    candidate.bundle_id,
                    candidate.direction,
                    len(candidate.sequences),
                    len(candidate.consensus_sequence),
                    candidate.retained_column_count,
                    candidate.dropped_column_count,
                    candidate.ambiguous_base_count,
                    f"{candidate.mean_retained_occupancy:.6f}",
                )
            )
    with (output_dir / "candidate_consensus.fa").open("w", encoding="utf-8") as handle:
        for candidate in candidates:
            handle.write(
                f">{candidate.candidate_id} bundle_id={candidate.bundle_id} "
                f"group_id={candidate.group_id} group_rank={candidate.group_rank} "
                f"haplotype_count={candidate.haplotype_count} "
                f"position={candidate.position_1based}\n"
            )
            handle.writelines(wrap_fasta(candidate.consensus_sequence, line_width))


def write_identity_report(
    path: Path,
    candidates: Sequence[CandidateConsensus],
    identities: Mapping[Tuple[str, str], Tuple[float, int, int]],
    cluster_by_candidate: Mapping[str, str],
    merge_identity: float,
) -> None:
    by_id = {candidate.candidate_id: candidate for candidate in candidates}
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "bundle_id", "candidate_a", "candidate_b", "identity", "matches",
                "compared_alignment_columns", "identity_above_threshold",
                "same_final_cluster",
            )
        )
        for (first_id, second_id), (identity, matches, compared) in sorted(
            identities.items()
        ):
            writer.writerow(
                (
                    by_id[first_id].bundle_id,
                    first_id,
                    second_id,
                    f"{identity:.6f}",
                    matches,
                    compared,
                    int(identity > merge_identity),
                    int(cluster_by_candidate[first_id] == cluster_by_candidate[second_id]),
                )
            )


def write_final_outputs(
    output_dir: Path,
    results: Sequence[FinalConsensus],
    line_width: int,
) -> None:
    with (output_dir / "bundle_consensus.fa").open(
        "w", encoding="utf-8"
    ) as handle, (output_dir / "bundle.fa").open("w", encoding="utf-8") as bundle_handle:
        for result in results:
            header = (
                f">{result.output_id} bundle_id={result.bundle_id} "
                f"allele={result.allele_index} candidate_count={len(result.candidate_ids)} "
                f"unique_path_count={result.unique_path_count} "
                f"source_occurrence_count={result.source_occurrence_count} "
                f"consensus_length={len(result.consensus_sequence)}\n"
            )
            handle.write(header)
            bundle_handle.write(header)
            handle.writelines(wrap_fasta(result.consensus_sequence, line_width))
            bundle_handle.writelines(wrap_fasta(result.consensus_sequence, line_width))
    with (output_dir / "bundle_consensus.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "output_id", "bundle_id", "allele_index", "candidate_count",
                "candidate_ids", "group_ids", "unique_path_count",
                "source_occurrence_count", "consensus_length", "ambiguous_base_count",
                "consensus_fasta",
            )
        )
        for result in results:
            writer.writerow(
                (
                    result.output_id,
                    result.bundle_id,
                    result.allele_index,
                    len(result.candidate_ids),
                    ",".join(result.candidate_ids),
                    ",".join(result.group_ids),
                    result.unique_path_count,
                    result.source_occurrence_count,
                    len(result.consensus_sequence),
                    result.ambiguous_base_count,
                    result.output_path,
                )
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Group haplotypes by signed bundle order, build group-position candidate "
            "consensuses, and merge same-bundle candidates above an identity threshold."
        )
    )
    parser.add_argument("bundle_bed", type=Path)
    parser.add_argument("alignment_manifest", type=Path)
    parser.add_argument("-o", "--output-dir", required=True, type=Path)
    parser.add_argument("--min-column-occupancy", type=float, default=0.50)
    parser.add_argument("--min-base-frequency", type=float, default=0.50)
    parser.add_argument("--merge-identity", type=float, default=0.90)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--all-haplotype-groups",
        action="store_true",
        help=(
            "build candidates from every bundle-order class; by default only the "
            "largest (rank-1) class contributes consensus sequences"
        ),
    )
    selection.add_argument(
        "--supplement-missing-bundles",
        action="store_true",
        help=(
            "start with rank-1 bundle candidates, then visit lower-support groups "
            "in rank order and add only IDs absent from all earlier groups"
        ),
    )
    parser.add_argument("--line-width", type=int, default=80)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    for name in (
        "min_column_occupancy", "min_base_frequency", "merge_identity"
    ):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            parser.error(f"--{name.replace('_', '-')} must be in [0, 1]")
    if args.line_width <= 0:
        parser.error("--line-width must be positive")
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        groups, candidates = build_candidates(
            args.bundle_bed,
            args.alignment_manifest,
            args.min_column_occupancy,
            args.min_base_frequency,
            args.all_haplotype_groups,
            args.supplement_missing_bundles,
        )
        clusters_by_bundle, identities = cluster_candidates(
            candidates, args.merge_identity
        )
        results = final_consensuses(
            clusters_by_bundle,
            args.output_dir,
            args.min_column_occupancy,
            args.min_base_frequency,
            args.line_width,
        )
        cluster_by_candidate: Dict[str, str] = {}
        for result in results:
            for candidate_id in result.candidate_ids:
                cluster_by_candidate[candidate_id] = result.output_id
        selected_group_ids = sorted({candidate.group_id for candidate in candidates})
        write_path_groups(
            args.output_dir / "path_order_groups.tsv",
            groups,
            selected_group_ids,
        )
        selected_bundle_ids_by_group: DefaultDict[str, Set[str]] = defaultdict(set)
        for candidate in candidates:
            selected_bundle_ids_by_group[candidate.group_id].add(candidate.bundle_id)
        write_bundle_source_groups(
            args.output_dir / "bundle_source_groups.tsv",
            groups,
            selected_bundle_ids_by_group,
        )
        write_candidate_outputs(args.output_dir, candidates, args.line_width)
        write_identity_report(
            args.output_dir / "candidate_identity.tsv",
            candidates,
            identities,
            cluster_by_candidate,
            args.merge_identity,
        )
        write_final_outputs(args.output_dir, results, args.line_width)
        print(
            f"Grouped {sum(group.haplotype_count for group in groups)} haplotypes into "
            f"{len(groups)} bundle-order classes; selected {len(selected_group_ids)} class(es) "
            f"and wrote {len(results)} final consensus sequences for "
            f"{len(clusters_by_bundle)} bundle IDs"
        )
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
