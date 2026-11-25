#!/usr/bin/env python3
"""Non-destructively partition library consensuses at inter-bundle PAF boundaries.

PAF is native here: query and target are both current library records. CIGAR is
always emitted, but no CIGAR projection or consensus merging is performed.
The legacy mode only partitions records. Optional overlap reconciliation also
unifies shared classes and removes high-identity shared sequence duplicates.
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

try:
    from .build_bundle_consensus import wrap_fasta
    from .refine_bundle_mappings import read_paf, normalize_bundle_name, merge_intervals, interval_total
    from .pipeline_intermediates import cleanup_generated_files
except ImportError:
    from build_bundle_consensus import wrap_fasta
    from refine_bundle_mappings import read_paf, normalize_bundle_name, merge_intervals, interval_total
    from pipeline_intermediates import cleanup_generated_files


SELF_TERMINAL_MIN_BP = 200
SELF_TERMINAL_FRACTION_DENOMINATOR = 20
SELF_TERMINAL_MIN_FRACTION = 1 / SELF_TERMINAL_FRACTION_DENOMINATOR


@dataclass(frozen=True)
class SelfAlignment:
    first_name: str
    first_length: int
    first_start: int
    first_end: int
    second_name: str
    second_length: int
    second_start: int
    second_end: int
    strand: str
    matching_bases: int
    alignment_block_length: int
    mapq: int
    line_number: int


def read_self_paf(path, records, min_identity=0., min_mapq=0):
    """Validate native self PAF; ignore same-name rows and low-quality hits."""
    alignments = []
    # The shared parser validates PAF then swaps its semantic axes. Swap back
    # explicitly because self alignment has no population/genome axis.
    for block in read_paf(Path(path)):
        for name, length in ((block.target_name, block.target_length),
                             (block.raw_query_name, block.query_length)):
            if name not in records:
                raise ValueError(f"self PAF contains a name outside the active library: {name}")
            if len(records[name]) != length:
                raise ValueError(f"self PAF/library length mismatch: {name}")
        if block.target_name == block.raw_query_name:
            continue
        if block.identity < min_identity or block.mapq < min_mapq:
            continue
        alignments.append(SelfAlignment(
            block.target_name, block.target_length, block.target_start, block.target_end,
            block.raw_query_name, block.query_length, block.query_start, block.query_end,
            block.strand, block.matching_bases, block.alignment_block_length,
            block.mapq, block.line_number,
        ))
    return alignments


def protect_self_terminal_fragments(cuts, length):
    """Suppress terminal cuts only when their fragment is <200 bp AND <5%.

    Each cut is compared to both current parent sequence ends. Checking every
    cut also prevents suppression of an outer cut from exposing another tiny
    terminal fragment. Internal fragment lengths do not impose a minimum.
    The integer fraction test keeps the exact 5% boundary inclusive.
    """
    if length <= 0:
        raise ValueError("self-alignment parent length must be positive")
    points = sorted(set(cuts))
    if any(point < 0 or point > length for point in points):
        raise ValueError("self-alignment cut is outside the parent sequence")
    accepted, suppressed = [], []
    for point in points:
        if not 0 < point < length:
            continue
        left_bp, right_bp = point, length - point
        left_tiny = left_bp < SELF_TERMINAL_MIN_BP and left_bp * SELF_TERMINAL_FRACTION_DENOMINATOR < length
        right_tiny = right_bp < SELF_TERMINAL_MIN_BP and right_bp * SELF_TERMINAL_FRACTION_DENOMINATOR < length
        (suppressed if left_tiny or right_tiny else accepted).append(point)
    return tuple(accepted), tuple(suppressed)


def alignment_decision(alignment, min_overlap_fraction=.20, min_overlap_bp=10000,
                       min_terminal_bp=0):
    """Test overlap and apply the fixed dual-threshold terminal protection.

    ``min_terminal_bp`` is retained for API compatibility but is ignored.
    Self comparison protects only terminal pieces that are both <200 bp and
    <5% of their current parent; the genome-remapping guard is independent.
    """
    if not 0 <= min_overlap_fraction <= 1 or min_overlap_bp < 0 or min_terminal_bp < 0:
        raise ValueError("invalid self-alignment overlap/terminal thresholds")
    if alignment.alignment_block_length <= 0:
        raise ValueError("self-alignment block length must be positive")
    for length, start, end in ((alignment.first_length, alignment.first_start, alignment.first_end),
                               (alignment.second_length, alignment.second_start, alignment.second_end)):
        if not 0 <= start < end <= length:
            raise ValueError("invalid self-alignment interval or length")
    first_span = alignment.first_end - alignment.first_start
    second_span = alignment.second_end - alignment.second_start
    first_fraction = first_span / alignment.first_length
    second_fraction = second_span / alignment.second_length
    overlap_bp = min(first_span, second_span)
    threshold_pass = (first_fraction > min_overlap_fraction or second_fraction > min_overlap_fraction
                      or max(first_span, second_span) > min_overlap_bp)
    first_raw = tuple(point for point in (alignment.first_start, alignment.first_end)
                      if 0 < point < alignment.first_length)
    second_raw = tuple(point for point in (alignment.second_start, alignment.second_end)
                       if 0 < point < alignment.second_length)
    first_cuts, first_suppressed = protect_self_terminal_fragments(first_raw, alignment.first_length)
    second_cuts, second_suppressed = protect_self_terminal_fragments(second_raw, alignment.second_length)
    splittable = threshold_pass and bool(first_cuts or second_cuts)
    status = ("below_threshold" if not threshold_pass else "splittable" if splittable
              else "terminal_guard_blocked" if first_raw or second_raw
              else "full_length_no_internal_endpoint")
    return {
        "first_name": alignment.first_name, "first_length": alignment.first_length,
        "first_start": alignment.first_start, "first_end": alignment.first_end,
        "second_name": alignment.second_name, "second_length": alignment.second_length,
        "second_start": alignment.second_start, "second_end": alignment.second_end,
        "strand": alignment.strand, "first_span_bp": first_span, "second_span_bp": second_span,
        "overlap_bp": overlap_bp, "first_fraction": first_fraction, "second_fraction": second_fraction,
        "identity": alignment.matching_bases / alignment.alignment_block_length,
        "mapq": alignment.mapq, "paf_line": alignment.line_number,
        "first_cuts": first_cuts, "second_cuts": second_cuts,
        "first_suppressed_cuts": first_suppressed, "second_suppressed_cuts": second_suppressed,
        "threshold_pass": threshold_pass, "splittable": splittable, "status": status,
    }


def canonical_alignment(alignment):
    """Represent reciprocal rows with the same axes (strand remains relative)."""
    if alignment.first_name <= alignment.second_name:
        return alignment
    return SelfAlignment(
        alignment.second_name, alignment.second_length, alignment.second_start, alignment.second_end,
        alignment.first_name, alignment.first_length, alignment.first_start, alignment.first_end,
        alignment.strand, alignment.matching_bases, alignment.alignment_block_length,
        alignment.mapq, alignment.line_number,
    )


def decision_sort_key(item):
    return (-item.get("pair_overlap_bp", item["overlap_bp"]), -item["overlap_bp"],
            item["first_name"], item["second_name"],
            item["first_start"], item["first_end"], item["second_start"], item["second_end"],
            item["strand"], -item["identity"], -item["mapq"], item["paf_line"])


CANDIDATE_COLUMNS = (
    "first_name", "first_length", "first_start", "first_end", "second_name", "second_length",
    "second_start", "second_end", "strand", "first_span_bp", "second_span_bp", "overlap_bp",
    "first_fraction", "second_fraction", "identity", "mapq", "paf_line",
    "pair_first_union_bp", "pair_second_union_bp", "pair_overlap_bp",
    "pair_first_fraction", "pair_second_fraction", "alignment_threshold_pass",
    "first_cuts", "second_cuts", "first_suppressed_cuts", "second_suppressed_cuts",
    "threshold_pass", "splittable", "status", "selected",
)
EVENT_COLUMNS = (
    "iteration", "parent_consensus", "parent_bundle_id", "parent_length", "part",
    "parent_start", "parent_end", "consensus_name", "bundle_id", "bundle_size", "partner_consensus",
    "action", "role", "retained", "output_reverse", "identity",
)


def write_table(path, columns, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            values = []
            for key in columns:
                value = row.get(key, "")
                if isinstance(value, tuple):
                    value = ",".join(map(str, value))
                elif isinstance(value, bool):
                    value = int(value)
                values.append(value)
            writer.writerow(values)


def write_library(directory, records, sizes, fasta_name="bundle.fa", metadata_name="bundle_metadata.bed"):
    with (directory / fasta_name).open("w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            handle.writelines(wrap_fasta(sequence, 80))
    with (directory / metadata_name).open("w") as handle:
        handle.write("# Reserved numeric ID/size metadata; NOT genomic mappings\n")
        for bundle_id in sorted(sizes, key=int):
            size = sizes[bundle_id]
            handle.write(f"__metadata__\t0\t1\t{bundle_id}:{size}:0:0:{size - 1}:U\n")


class _SelfSnapshotRetention:
    """Retire owned snapshots only after their replacement input is written."""

    def __init__(self, output_dir, keep_intermediates):
        self.output_dir = output_dir
        self.keep_intermediates = keep_intermediates
        self.generated_paths = set()
        self.removed_intermediates = []

    def write(self, directory, records, sizes, fasta_name="bundle.fa",
              metadata_name="bundle_metadata.bed"):
        path = directory / fasta_name
        newly_created = not path.exists() and not path.is_symlink()
        write_library(directory, records, sizes, fasta_name, metadata_name)
        if newly_created:
            self.generated_paths.add(path)

    def retire_previous(self, iterations, current_input):
        if not iterations or self.keep_intermediates:
            return
        prior = iterations[-1]
        fields = ["output_fasta"]
        if prior["iteration"] != 1:
            fields.append("library_fasta")
        candidates = {
            self.output_dir / prior[field]: field
            for field in fields if prior.get(field) is not None
        }
        owned = set(candidates) & self.generated_paths
        audit = cleanup_generated_files(
            self.output_dir, sorted(owned, key=str),
            protected_paths=(self.output_dir / "iterations/iteration_0001/bundle.fa", current_input),
        )
        for removed in audit["deleted"]:
            path = self.output_dir / removed["path"]
            field = candidates[path]
            entry = {**removed, "iteration": prior["iteration"], "field": field,
                     "reason": "superseded_by_successfully_written_next_iteration_input"}
            prior[field] = None
            prior.setdefault("removed_intermediates", []).append(entry)
            self.removed_intermediates.append(entry)
            self.generated_paths.discard(path)

    def summary(self):
        audit_path = self.output_dir / "intermediate_cleanup.json"
        if audit_path.is_symlink() or (audit_path.exists() and not audit_path.is_file()):
            raise ValueError(f"unsafe self-alignment cleanup audit target: {audit_path}")
        audit = {
            "stage": "self_alignment", "keep_intermediates": self.keep_intermediates,
            "deleted": self.removed_intermediates,
            "deleted_count": len(self.removed_intermediates),
            "deleted_bytes": sum(item["bytes"] for item in self.removed_intermediates),
        }
        audit_path.write_text(json.dumps(audit, indent=2) + "\n")
        return {
            "retention": {
                "keep_intermediates": self.keep_intermediates,
                "policy": "retain initial input and final iteration input/output FASTAs; "
                          "retain every PAF, report, log and metadata file; retire only "
                          "newly generated superseded FASTAs after replacement input succeeds",
            },
            "removed_intermediates": self.removed_intermediates,
        }


def _reverse_complement(sequence):
    return sequence.translate(str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn"))[::-1]


def _priority_rank(priorities, name):
    value = priorities.get(name, priorities.get(normalize_bundle_name(name, {}), float("inf")))
    if isinstance(value, dict):
        value = value.get("group_rank", float("inf"))
    return float("inf") if value is None else float(value)


def _shared_partition(decision, side):
    """Exact outer PAF slices; suppressed short ends stay on the shared slice."""
    length = decision[f"{side}_length"]
    boundaries = (0, *decision[f"{side}_cuts"], length)
    aligned_start, aligned_end = decision[f"{side}_start"], decision[f"{side}_end"]
    return [(part, start, end, "shared" if start <= aligned_start and aligned_end <= end else "flank")
            for part, (start, end) in enumerate(zip(boundaries, boundaries[1:]), 1)]


def _self_decisions_from_alignments(alignments, records, min_overlap_fraction, min_overlap_bp):
    """Deduplicate reciprocal alignments and apply the pair-wide interval-union gate."""
    unique = {}
    for alignment in alignments:
        alignment = canonical_alignment(alignment)
        key = (alignment.first_name, alignment.second_name, alignment.strand,
               alignment.first_start, alignment.first_end, alignment.second_start, alignment.second_end)
        previous = unique.get(key)
        quality = lambda item: (item.matching_bases / item.alignment_block_length, item.mapq, -item.line_number)
        if previous is None or quality(alignment) > quality(previous):
            unique[key] = alignment
    decisions = [alignment_decision(item, min_overlap_fraction, min_overlap_bp) for item in unique.values()]
    pairs = {}
    for decision in decisions:
        pairs.setdefault((decision["first_name"], decision["second_name"]), []).append(decision)
    for pair, members in pairs.items():
        first_bp = interval_total(merge_intervals((item["first_start"], item["first_end"]) for item in members))
        second_bp = interval_total(merge_intervals((item["second_start"], item["second_end"]) for item in members))
        first_fraction, second_fraction = first_bp / len(records[pair[0]]), second_bp / len(records[pair[1]])
        passed = (first_fraction > min_overlap_fraction or second_fraction > min_overlap_fraction
                  or max(first_bp, second_bp) > min_overlap_bp)
        for decision in members:
            decision.update({"alignment_threshold_pass": decision["threshold_pass"],
                             "pair_first_union_bp": first_bp, "pair_second_union_bp": second_bp,
                             "pair_overlap_bp": min(first_bp, second_bp),
                             "pair_first_fraction": first_fraction, "pair_second_fraction": second_fraction,
                             "threshold_pass": passed})
            decision["splittable"] = passed and bool(decision["first_cuts"] or decision["second_cuts"])
    return sorted(decisions, key=decision_sort_key)


def _self_decisions(paf, records, min_identity, min_mapq, min_overlap_fraction, min_overlap_bp):
    """Read native PAF and build pair-wide self-alignment decisions."""
    return _self_decisions_from_alignments(
        read_self_paf(paf, records, min_identity, min_mapq), records,
        min_overlap_fraction, min_overlap_bp,
    )


def _write_synthetic_self_paf(path, alignments):
    """Persist non-PAF aligner results in native bundle-to-bundle PAF form."""
    with Path(path).open("w") as handle:
        for alignment in alignments:
            fields = (
                alignment.first_name, alignment.first_length,
                alignment.first_start, alignment.first_end, alignment.strand,
                alignment.second_name, alignment.second_length,
                alignment.second_start, alignment.second_end,
                alignment.matching_bases, alignment.alignment_block_length,
                alignment.mapq,
            )
            handle.write("\t".join(map(str, fields)) + "\n")


def _reconcile_rank1_library(records, sizes, output_dir, minimap2, threads,
                             min_overlap_fraction, min_overlap_bp, min_terminal_bp,
                             min_identity, min_mapq, stage_label, merge_identity, priorities,
                             keep_intermediates, alignment_provider=None,
                             alignment_backend="minimap2"):
    """One pair per iteration: partition, assign a shared family, then deduplicate.

    All coordinates refer to frozen input records. Every participating parent
    receives events, including unsplit retained aliases and suppressed shared
    copies. Replay therefore removes all iteration parents before adding only
    retained outputs; a suppressed output name identifies its representative.
    """
    if not records or any(not sequence for sequence in records.values()):
        raise ValueError("rank1 self alignment requires a nonempty library")
    if threads < 1 or not 0 <= min_identity <= 1 or min_mapq < 0 or not 0 <= merge_identity <= 1:
        raise ValueError("invalid self-alignment threads/quality thresholds")
    if not 0 <= min_overlap_fraction <= 1 or min_overlap_bp < 0 or min_terminal_bp < 0:
        raise ValueError("invalid self-alignment overlap/terminal thresholds")
    records, sizes, priorities = dict(records), dict(sizes), dict(priorities or {})
    for name in records:
        bundle_id = normalize_bundle_name(name, {})
        if not bundle_id.isdigit() or str(int(bundle_id)) != bundle_id:
            raise ValueError(f"rank1 self alignment requires numeric bundle classes: {name}")
        sizes.setdefault(bundle_id, 100)
    if any(not key.isdigit() or str(int(key)) != key or value < 1 for key, value in sizes.items()):
        raise ValueError("invalid reserved numeric bundle metadata")
    next_id = max(map(int, sizes)) + 1
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    retention = _SelfSnapshotRetention(output_dir, keep_intermediates)
    summary_path = output_dir / "self_alignment_summary.json"
    summary_path.write_text(json.dumps({"version": 4, "completed": False, **retention.summary()}) + "\n")
    initial_count, initial_bp = len(records), sum(map(len, records.values()))
    events, iterations, alias_events, aliases, states = [], [], [], {}, set()
    split_operations = merge_operations = allele_operations = 0
    iteration = 0
    write_table(output_dir / "self_split_events.tsv", EVENT_COLUMNS, [])

    def fresh_class():
        nonlocal next_id
        bundle_id = str(next_id)
        next_id += 1
        return bundle_id

    while True:
        state = tuple(sorted(records.items()))
        if state in states:
            raise ValueError("overlap reconciliation revisited an active library state")
        states.add(state)
        iteration += 1
        directory = output_dir / "iterations" / f"iteration_{iteration:04d}"
        directory.mkdir(parents=True, exist_ok=True)
        retention.write(directory, records, sizes)
        fasta, paf = directory / "bundle.fa", directory / "bundle_self.raw.paf"
        retention.retire_previous(iterations, fasta)
        summary_path.write_text(json.dumps({"version": 4, "completed": False,
                                           "iterations": iterations, **retention.summary()}, indent=2) + "\n")
        print(f"  [{stage_label} self {iteration}] Reconciling {len(records)} bundles", flush=True)
        provider_info = {}
        if alignment_provider is None:
            command = [minimap2, "-x", "asm20", "-X", "-g", "500", "-c", "-t", str(threads), str(fasta), str(fasta)]
            with paf.open("w") as output, (directory / "minimap.log").open("w") as log:
                subprocess.run(command, stdout=output, stderr=log, check=True)
            decisions = _self_decisions(
                paf, records, min_identity, min_mapq,
                min_overlap_fraction, min_overlap_bp,
            )
        else:
            alignments, provider_info = alignment_provider(records, directory, iteration)
            alignments = [item for item in alignments
                          if item.matching_bases / item.alignment_block_length >= min_identity
                          and item.mapq >= min_mapq]
            _write_synthetic_self_paf(paf, alignments)
            decisions = _self_decisions_from_alignments(
                alignments, records, min_overlap_fraction, min_overlap_bp,
            )
        selected = None
        for decision in decisions:
            same_class = (normalize_bundle_name(decision["first_name"], {})
                          == normalize_bundle_name(decision["second_name"], {}))
            high = decision["identity"] > merge_identity
            actionable = decision["threshold_pass"] and (high or not same_class or decision["splittable"])
            decision["status"] = ("below_threshold" if not decision["threshold_pass"] else
                                  "merge_overlap" if high else "unify_overlap" if actionable else "same_class_reconciled")
            if selected is None and actionable:
                selected = decision
        for decision in decisions:
            decision["selected"] = decision is selected
        write_table(directory / "alignment_candidates.tsv", CANDIDATE_COLUMNS, decisions)
        write_table(directory / "selected_split.tsv", CANDIDATE_COLUMNS, [selected] if selected else [])
        info = {"iteration": iteration, "input_bundle_count": len(records),
                "alignment_backend": alignment_backend,
                "paf": str(paf.relative_to(output_dir)), "library_fasta": str(fasta.relative_to(output_dir)),
                "metadata": str((directory / "bundle_metadata.bed").relative_to(output_dir)),
                "candidates": str((directory / "alignment_candidates.tsv").relative_to(output_dir)),
                "selected_split": str((directory / "selected_split.tsv").relative_to(output_dir)),
                "qualifying_alignment_count": sum(item["threshold_pass"] for item in decisions),
                "splittable_alignment_count": sum(item["splittable"] for item in decisions),
                "qualifying_pair_count": len({(item["first_name"], item["second_name"])
                                             for item in decisions if item["threshold_pass"]}),
                "selected_pair": [selected["first_name"], selected["second_name"]] if selected else None,
                "alignment_provider": provider_info}
        iterations.append(info)
        batch = []
        if selected is not None:
            frozen = dict(records)
            frozen_ranks = {name: _priority_rank(priorities, name) for name in frozen}
            partitions = {side: _shared_partition(selected, side) for side in ("first", "second")}
            shared = {side: next(piece for piece in pieces if piece[3] == "shared")
                      for side, pieces in partitions.items()}
            high = selected["identity"] > merge_identity

            def representative_key(side):
                name = selected[f"{side}_name"]
                piece = shared[side]
                return (_priority_rank(priorities, name), -(piece[2] - piece[1]),
                        int(normalize_bundle_name(name, {})), name)

            representative = min(("first", "second"), key=representative_key)
            whole_sides = [side for side in ("first", "second") if not selected[f"{side}_cuts"]]
            class_side = min(whole_sides, key=representative_key) if whole_sides else None
            shared_class = normalize_bundle_name(selected[f"{class_side}_name"], {}) if class_side else fresh_class()
            participating = {selected["first_name"], selected["second_name"]}
            # A whole-member class union must carry every existing allele along.
            migrated = {}
            for side in whole_sides:
                old_class = normalize_bundle_name(selected[f"{side}_name"], {})
                reverse_family = side != representative and selected["strand"] == "-"
                if old_class == shared_class and not reverse_family:
                    continue
                if old_class != shared_class:
                    aliases[old_class] = shared_class
                    alias_events.append({"iteration": iteration, "source_bundle_id": old_class,
                                         "bundle_id": shared_class, "reason": "whole_family_union"})
                for name in frozen:
                    if normalize_bundle_name(name, {}) == old_class and name not in participating:
                        migrated[name] = side
                participating.update(migrated)
            for name in participating:
                records.pop(name)
            shared_size = max(max(1, math.ceil(sizes[normalize_bundle_name(selected[f"{side}_name"], {})]
                                                   * (shared[side][2] - shared[side][1]) / len(frozen[selected[f"{side}_name"]])))
                              for side in ("first", "second"))
            sizes[shared_class] = max(sizes.get(shared_class, 1), shared_size)

            def output_name(bundle_id, preferred=None):
                if preferred and normalize_bundle_name(preferred, {}) == bundle_id and preferred not in records:
                    return preferred
                base = f"bundle_{bundle_id}"
                if base not in records:
                    return base
                allele = 2
                while f"{base}.allele_{allele}" in records:
                    allele += 1
                return f"{base}.allele_{allele}"

            def add_event(name, partner, part, start, end, child_name, child_id, size,
                          role, action, retained=True, reverse=False):
                row = {"iteration": iteration, "parent_consensus": name,
                       "parent_bundle_id": normalize_bundle_name(name, {}), "parent_length": len(frozen[name]),
                       "part": part, "parent_start": start, "parent_end": end,
                       "consensus_name": child_name, "bundle_id": child_id, "bundle_size": size,
                       "partner_consensus": partner, "action": action, "role": role,
                       "retained": retained, "output_reverse": reverse, "identity": selected["identity"]}
                batch.append(row)
                if retained:
                    sequence = frozen[name][start:end]
                    records[child_name] = _reverse_complement(sequence) if reverse else sequence
                    priorities[child_name] = frozen_ranks[name]

            # Allocate the representative first, so suppressed aliases always
            # reference the exact actual FASTA record, not just a numeric class.
            rep_name = selected[f"{representative}_name"]
            rep_piece = shared[representative]
            rep_output = output_name(shared_class, rep_name if not selected[f"{representative}_cuts"] else None)
            other_side = "second" if representative == "first" else "first"
            add_event(rep_name, selected[f"{other_side}_name"], *rep_piece[:3], rep_output,
                      shared_class, sizes[shared_class], "shared", "merge" if high else "allele")
            other_name, other_piece = selected[f"{other_side}_name"], shared[other_side]
            other_output = rep_output if high else output_name(
                shared_class, other_name if not selected[f"{other_side}_cuts"] else None)
            add_event(other_name, rep_name, *other_piece[:3], other_output, shared_class,
                      sizes[shared_class], "shared", "suppress" if high else "allele",
                      retained=not high, reverse=selected["strand"] == "-")
            for side in ("first", "second"):
                name, partner = selected[f"{side}_name"], selected[f"{'second' if side == 'first' else 'first'}_name"]
                for part, start, end, role in partitions[side]:
                    if role == "shared":
                        continue
                    child_id = fresh_class()
                    size = max(1, math.ceil(sizes[normalize_bundle_name(name, {})] * (end - start) / len(frozen[name])))
                    sizes[child_id] = size
                    add_event(name, partner, part, start, end, output_name(child_id), child_id,
                              size, "flank", "split")
            for name, side in sorted(migrated.items()):
                child_name = output_name(shared_class, name)
                add_event(name, rep_name, 1, 0, len(frozen[name]), child_name, shared_class,
                          sizes[shared_class], "shared", "family_alias",
                          reverse=side != representative and selected["strand"] == "-")
            # Keep event order chronological and each parent's pieces in
            # forward-coordinate order; output creation above is already done.
            batch.sort(key=lambda row: (row["parent_consensus"], row["part"]))
            events.extend(batch)
            split_operations += bool(selected["first_cuts"] or selected["second_cuts"])
            merge_operations += high
            allele_operations += not high
            info.update({"action": "merge" if high else "unify_alleles", "shared_bundle_id": shared_class,
                         "representative_consensus": rep_output, "identity": selected["identity"],
                         "suppressed_record_count": sum(not row["retained"] for row in batch)})
            if tuple(sorted(records.items())) == state:
                raise ValueError("selected overlap reconciliation made no library change")
            write_table(output_dir / "self_split_events.tsv", EVENT_COLUMNS, events)
        retention.write(directory, records, sizes, "bundle_after.fa", "bundle_after.metadata.bed")
        write_table(directory / "bundle_fragments.tsv", EVENT_COLUMNS, batch)
        info.update({"output_fasta": str((directory / "bundle_after.fa").relative_to(output_dir)),
                     "output_metadata": str((directory / "bundle_after.metadata.bed").relative_to(output_dir)),
                     "fragments": str((directory / "bundle_fragments.tsv").relative_to(output_dir)),
                     "output_bundle_count": len(records)})
        if selected is None:
            remaining = [item for item in decisions if item["threshold_pass"]]
            break

    for old_class in aliases:
        target, seen = aliases[old_class], {old_class}
        while target in aliases:
            if target in seen:
                raise ValueError("cycle in reconciled bundle-family aliases")
            seen.add(target)
            target = aliases[target]
        aliases[old_class] = target
    families = [{"consensus_name": name, "bundle_id": normalize_bundle_name(name, {}),
                 "sequence_length": len(sequence), "group_rank": _priority_rank(priorities, name)}
                for name, sequence in records.items()]
    suppressed = [row for row in events if not row["retained"]]
    write_library(output_dir, records, sizes)
    write_table(output_dir / "self_split_events.tsv", EVENT_COLUMNS, events)
    write_table(output_dir / "self_overlap_aliases.tsv", ("iteration", "source_bundle_id", "bundle_id", "reason"), alias_events)
    write_table(output_dir / "self_suppressed_sequences.tsv", EVENT_COLUMNS, suppressed)
    write_table(output_dir / "self_bundle_families.tsv", ("consensus_name", "bundle_id", "sequence_length", "group_rank"), families)
    final_bp = sum(map(len, records.values()))
    summary = {
        "version": 4, "completed": True, "stop_reason": "no_actionable_overlaps", "stage_label": stage_label,
        "parameters": {"reconcile_overlaps": True, "merge_identity": merge_identity,
                       "alignment_backend": alignment_backend,
                       "min_overlap_fraction": min_overlap_fraction, "min_overlap_bp": min_overlap_bp,
                       "min_terminal_bp": SELF_TERMINAL_MIN_BP, "min_terminal_fraction": SELF_TERMINAL_MIN_FRACTION,
                       "legacy_min_terminal_bp": min_terminal_bp, "min_identity": min_identity, "min_mapq": min_mapq},
        "counts": {"initial_bundle_count": initial_count, "final_bundle_count": len(records),
                   "initial_sequence_bp": initial_bp, "total_sequence_bp": final_bp,
                   "removed_sequence_bp": initial_bp - final_bp, "iteration_count": iteration,
                   "reconciliation_operation_count": iteration - 1, "split_operation_count": split_operations,
                   "merge_operation_count": merge_operations, "allele_unification_operation_count": allele_operations,
                   "split_parent_count": len({row["parent_consensus"] for row in events if row["action"] == "split"}),
                   "participating_parent_count": len({row["parent_consensus"] for row in events}),
                   "suppressed_sequence_count": len(suppressed), "final_family_count": len({row["bundle_id"] for row in families}),
                   "remaining_qualifying_alignments": len(remaining)},
        "semantics": {"paf_axes": "native query=bundle, target=bundle; forward sequence coordinates",
                      "fraction_test": "either pair-wide interval union / parent length strictly > threshold",
                      "bp_test": "max(query union bp, target union bp) strictly > threshold",
                      "coverage_basis": "union of forward outer PAF intervals; no CIGAR or linear projection",
                      "identity_test": "PAF matching_bases / alignment_block_length strictly > merge_identity",
                      "priority": "largest pair union overlap, largest block; representative: earliest group rank, longest shared-containing slice, numeric class/name",
                      "one_pair_per_iteration": True, "same_name_hits_ignored": True, "reciprocal_rows_deduplicated": True,
                      "cigar_used": False, "all_sequences_retained": False, "terminal_guard_applied": True,
                      "terminal_guard_rule": "suppress only terminal fragments strictly <200 bp AND strictly <5% of current parent",
                      "minimum_fragment_length_applied": False, "internal_fragment_minimum_applied": False,
                      "legacy_min_terminal_bp_ignored": True,
                      "low_identity_overlaps": "unique FASTA alleles share one numeric bundle class",
                      "reverse_shared_alleles": "oriented to representative; flanks retain parent orientation",
                      "event_replay": "freeze and remove all iteration parents, then create retained outputs from parent slices; output_reverse applies reverse complement",
                      "paths_relative_to": "self-alignment output directory"},
        "aliases": aliases, "iterations": iterations, "remaining_qualifying_alignments": remaining,
        **retention.summary(),
        "outputs": {"bundle_fasta": "bundle.fa", "metadata": "bundle_metadata.bed", "events": "self_split_events.tsv",
                    "aliases": "self_overlap_aliases.tsv", "suppressed_sequences": "self_suppressed_sequences.tsv",
                    "families": "self_bundle_families.tsv", "cleanup_report": "intermediate_cleanup.json"},
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    return {"records": records, "sizes": sizes, "events": events, "summary": summary,
            "summary_path": summary_path, "aliases": aliases, "suppressed": suppressed, "families": families}


def split_rank1_library(records, sizes, output_dir, minimap2, threads=8,
                        min_overlap_fraction=.20, min_overlap_bp=10000,
                        min_terminal_bp=0, min_identity=0., min_mapq=0, stage_label="rank1", *,
                        reconcile_overlaps=False, merge_identity=.90, priorities=None,
                        keep_intermediates=True, alignment_provider=None,
                        alignment_backend="minimap2"):
    """Re-align after each largest valid pair split until a fixed point.

    Both parents are partitioned when they have accepted internal boundaries;
    a full-length member of a containment pair is retained unchanged. Each
    child's sequence is an exact parent slice. Fully homologous children remain
    distinct: this operation does not grant authority to merge/drop records.
    Internal pieces have no minimum length. A terminal cut is suppressed only
    when its fragment is both <200 bp and <5% of its current parent. The legacy
    ``min_terminal_bp`` argument is ignored, independently of genome mapping.
    By default the standalone Python API retains every intermediate snapshot.
    Pipelines may disable this to retire superseded, newly generated FASTAs.
    """
    if reconcile_overlaps:
        return _reconcile_rank1_library(
            records, sizes, output_dir, minimap2, threads, min_overlap_fraction,
            min_overlap_bp, min_terminal_bp, min_identity, min_mapq, stage_label,
            merge_identity, priorities, keep_intermediates,
            alignment_provider, alignment_backend,
        )
    if not records or any(not sequence for sequence in records.values()):
        raise ValueError("rank1 self alignment requires a nonempty library")
    if threads < 1 or not 0 <= min_identity <= 1 or min_mapq < 0:
        raise ValueError("invalid self-alignment threads/quality thresholds")
    if not 0 <= min_overlap_fraction <= 1 or min_overlap_bp < 0 or min_terminal_bp < 0:
        raise ValueError("invalid self-alignment overlap/terminal thresholds")
    records, sizes = dict(records), dict(sizes)
    for name in records:
        bundle_id = normalize_bundle_name(name, {})
        if not bundle_id.isdigit() or str(int(bundle_id)) != bundle_id:
            raise ValueError(f"rank1 self alignment requires numeric bundle classes: {name}")
        sizes.setdefault(bundle_id, 100)
    if any(not key.isdigit() or str(int(key)) != key or value < 1 for key, value in sizes.items()):
        raise ValueError("invalid reserved numeric bundle metadata")
    next_id = max(map(int, sizes)) + 1
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    retention = _SelfSnapshotRetention(output_dir, keep_intermediates)
    summary_path = output_dir / "self_alignment_summary.json"
    summary_path.write_text(json.dumps({"version": 3, "completed": False, **retention.summary()}) + "\n")
    initial_count, total_bp = len(records), sum(map(len, records.values()))
    events, iterations = [], []
    write_table(output_dir / "self_split_events.tsv", EVENT_COLUMNS, events)
    iteration = 0
    while True:
        iteration += 1
        directory = output_dir / "iterations" / f"iteration_{iteration:04d}"
        directory.mkdir(parents=True, exist_ok=True)
        retention.write(directory, records, sizes)
        fasta, paf = directory / "bundle.fa", directory / "bundle_self.raw.paf"
        retention.retire_previous(iterations, fasta)
        summary_path.write_text(json.dumps({"version": 3, "completed": False,
                                           "iterations": iterations, **retention.summary()}, indent=2) + "\n")
        # All chains are requested so diagonal/primary hits do not hide other
        # bundle matches. Same-name hits are also excluded again by the parser.
        command = [minimap2, "-x", "asm20", "-X", "-g", "500", "-c", "-t", str(threads), str(fasta), str(fasta)]
        print(f"  [{stage_label} self {iteration}] Aligning {len(records)} bundles", flush=True)
        with paf.open("w") as output, (directory / "minimap.log").open("w") as log:
            subprocess.run(command, stdout=output, stderr=log, check=True)
        unique = {}
        for alignment in read_self_paf(paf, records, min_identity, min_mapq):
            alignment = canonical_alignment(alignment)
            key = (alignment.first_name, alignment.second_name, alignment.strand,
                   alignment.first_start, alignment.first_end, alignment.second_start, alignment.second_end)
            previous = unique.get(key)
            quality = lambda item: (item.matching_bases / item.alignment_block_length, item.mapq, -item.line_number)
            if previous is None or quality(alignment) > quality(previous):
                unique[key] = alignment
        decisions = [alignment_decision(item, min_overlap_fraction, min_overlap_bp)
                     for item in unique.values()]
        pairs = {}
        for decision in decisions:
            pair = (decision["first_name"], decision["second_name"])
            pairs.setdefault(pair, []).append(decision)
        for pair, members in pairs.items():
            first_bp = interval_total(merge_intervals((item["first_start"], item["first_end"]) for item in members))
            second_bp = interval_total(merge_intervals((item["second_start"], item["second_end"]) for item in members))
            first_fraction, second_fraction = first_bp / len(records[pair[0]]), second_bp / len(records[pair[1]])
            threshold_pass = (first_fraction > min_overlap_fraction or second_fraction > min_overlap_fraction
                              or max(first_bp, second_bp) > min_overlap_bp)
            for decision in members:
                decision.update({"alignment_threshold_pass": decision["threshold_pass"],
                                 "pair_first_union_bp": first_bp, "pair_second_union_bp": second_bp,
                                 "pair_overlap_bp": min(first_bp, second_bp),
                                 "pair_first_fraction": first_fraction, "pair_second_fraction": second_fraction,
                                 "threshold_pass": threshold_pass})
                decision["splittable"] = threshold_pass and bool(decision["first_cuts"] or decision["second_cuts"])
                decision["status"] = ("below_threshold" if not threshold_pass else "splittable" if decision["splittable"]
                                      else "terminal_guard_blocked" if decision["first_suppressed_cuts"] or decision["second_suppressed_cuts"]
                                      else "full_length_no_internal_endpoint")
        decisions.sort(key=decision_sort_key)
        selected = next((item for item in decisions if item["splittable"]), None)
        for decision in decisions:
            decision["selected"] = decision is selected
        write_table(directory / "alignment_candidates.tsv", CANDIDATE_COLUMNS, decisions)
        write_table(directory / "selected_split.tsv", CANDIDATE_COLUMNS, [selected] if selected else [])
        iteration_info = {
            "iteration": iteration, "input_bundle_count": len(records),
            "paf": str(paf.relative_to(output_dir)),
            "library_fasta": str(fasta.relative_to(output_dir)),
            "metadata": str((directory / "bundle_metadata.bed").relative_to(output_dir)),
            "candidates": str((directory / "alignment_candidates.tsv").relative_to(output_dir)),
            "selected_split": str((directory / "selected_split.tsv").relative_to(output_dir)),
            "qualifying_alignment_count": sum(item["threshold_pass"] for item in decisions),
            "splittable_alignment_count": sum(item["splittable"] for item in decisions),
            "qualifying_pair_count": len({(item["first_name"], item["second_name"])
                                          for item in decisions if item["threshold_pass"]}),
            "selected_pair": [selected["first_name"], selected["second_name"]] if selected else None,
        }
        iterations.append(iteration_info)
        if selected is None:
            retention.write(directory, records, sizes, "bundle_after.fa", "bundle_after.metadata.bed")
            write_table(directory / "bundle_fragments.tsv", EVENT_COLUMNS, [])
            iteration_info.update({"output_fasta": str((directory / "bundle_after.fa").relative_to(output_dir)),
                                   "output_metadata": str((directory / "bundle_after.metadata.bed").relative_to(output_dir)),
                                   "fragments": str((directory / "bundle_fragments.tsv").relative_to(output_dir))})
            remaining = [item for item in decisions if item["threshold_pass"]]
            break
        prior_count = len(records)
        first_event = len(events)
        for side, other in (("first", "second"), ("second", "first")):
            cuts = selected[f"{side}_cuts"]
            if not cuts:
                continue
            name = selected[f"{side}_name"]
            partner = selected[f"{other}_name"]
            sequence = records.pop(name)
            bundle_id = normalize_bundle_name(name, {})
            boundaries = (0, *cuts, len(sequence))
            for part, (start, end) in enumerate(zip(boundaries, boundaries[1:]), 1):
                child_id, child_name = str(next_id), f"bundle_{next_id}"
                next_id += 1
                size = max(1, math.ceil(sizes[bundle_id] * (end - start) / len(sequence)))
                sizes[child_id] = size
                records[child_name] = sequence[start:end]
                events.append({"iteration": iteration, "parent_consensus": name, "parent_bundle_id": bundle_id,
                               "parent_length": len(sequence), "part": part, "parent_start": start,
                               "parent_end": end, "consensus_name": child_name, "bundle_id": child_id,
                               "bundle_size": size, "partner_consensus": partner})
        if len(records) <= prior_count or sum(map(len, records.values())) != total_bp:
            raise ValueError("self split did not make monotonic non-destructive progress")
        retention.write(directory, records, sizes, "bundle_after.fa", "bundle_after.metadata.bed")
        write_table(directory / "bundle_fragments.tsv", EVENT_COLUMNS, events[first_event:])
        iteration_info.update({"output_fasta": str((directory / "bundle_after.fa").relative_to(output_dir)),
                               "output_metadata": str((directory / "bundle_after.metadata.bed").relative_to(output_dir)),
                               "fragments": str((directory / "bundle_fragments.tsv").relative_to(output_dir))})
        write_table(output_dir / "self_split_events.tsv", EVENT_COLUMNS, events)
    write_library(output_dir, records, sizes)
    write_table(output_dir / "self_split_events.tsv", EVENT_COLUMNS, events)
    summary = {
        "version": 3, "completed": True, "stop_reason": "no_splittable_alignments",
        "stage_label": stage_label,
        "parameters": {"min_overlap_fraction": min_overlap_fraction, "min_overlap_bp": min_overlap_bp,
                       "min_terminal_bp": SELF_TERMINAL_MIN_BP,
                       "min_terminal_fraction": SELF_TERMINAL_MIN_FRACTION,
                       "legacy_min_terminal_bp": min_terminal_bp,
                       "min_identity": min_identity, "min_mapq": min_mapq},
        "counts": {"initial_bundle_count": initial_count, "final_bundle_count": len(records),
                   "iteration_count": iteration, "split_operation_count": iteration - 1,
                   "split_parent_count": len({row["parent_consensus"] for row in events}),
                   "remaining_qualifying_alignments": len(remaining), "total_sequence_bp": total_bp},
        "semantics": {"paf_axes": "native query=bundle, target=bundle; forward sequence coordinates",
                      "fraction_test": "either pair-wide PAF interval union / its bundle length strictly > threshold",
                      "bp_test": "max(query union bp, target union bp) strictly > threshold",
                      "coverage_basis": "union of forward outer PAF intervals; no CIGAR or linear projection",
                      "priority": "largest min(pair union bp), then largest min(block spans), deterministic names/coordinates",
                      "one_pair_per_iteration": True, "same_name_hits_ignored": True,
                      "reciprocal_rows_deduplicated": True, "cigar_used": False,
                      "all_sequences_retained": True, "terminal_guard_applied": True,
                      "terminal_guard_rule": "suppress only terminal fragments strictly <200 bp AND strictly <5% of the current parent length",
                      "minimum_fragment_length_applied": False, "legacy_min_terminal_bp_ignored": True,
                      "internal_fragment_minimum_applied": False,
                      "remaining_homologs": "reported, not merged/deleted; no artificial cuts for full-length matches",
                      "paths_relative_to": "self-alignment output directory"},
        "iterations": iterations, "remaining_qualifying_alignments": remaining,
        **retention.summary(),
        "outputs": {"bundle_fasta": "bundle.fa", "metadata": "bundle_metadata.bed",
                    "events": "self_split_events.tsv", "cleanup_report": "intermediate_cleanup.json"},
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    return {"records": records, "sizes": sizes, "events": events, "summary": summary,
            "summary_path": summary_path}
