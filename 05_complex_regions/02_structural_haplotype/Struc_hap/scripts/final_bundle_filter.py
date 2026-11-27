"""Filter final mapping classes without erasing inversion/order distinctions.

The library FASTA is not modified. Coverage is evaluated on every population
FASTA sequence, including sequences without surviving mappings. Signed bundle
orders are compared before and after each atomic whole-class deletion. A class
is retained when deleting it would leave the population haplotype partition
unchanged; eligible deletions may merge haplotypes through presence/absence,
non-tandem copy-number differences, or placement of a deleted class in
different repeated units. Copy-number differences remain protected when that
class occurs in a consecutive tandem run in any original haplotype. Opposite
orientation patterns in the same surviving-bundle context are protected when
they cannot be explained by copy indels, whereas differences confined to
separate contexts may represent multiple indels.
"""

from __future__ import annotations

import statistics
from collections import Counter, defaultdict
from typing import Tuple

try:
    from .merge_pgr_bundles import read_sequence_lengths, union_intervals, write_tsv
except ImportError:
    from merge_pgr_bundles import read_sequence_lengths, union_intervals, write_tsv


Signature = Tuple[Tuple[str, str], ...]


def signed_orders(mappings, assignments, sequence_lengths):
    by_path = {path_id: [] for path_id in sequence_lengths}
    for index, mapping in enumerate(mappings):
        by_path[mapping.target_name].append((mapping.target_start, mapping.target_end, index, mapping))
    return {
        path_id: tuple((assignments[item.bundle_name].numeric_bundle_id, item.strand)
                       for _, _, _, item in sorted(records))
        for path_id, records in by_path.items()
    }


def order_text(signature: Signature) -> str:
    return ",".join(direction + bundle_id for bundle_id, direction in signature) or "."


def haplotype_rows(signatures):
    by_signature = defaultdict(list)
    for path_id, signature in signatures.items():
        by_signature[signature].append(path_id)
    return [
        (f"H{rank:04d}", len(paths), order_text(signature), ",".join(sorted(paths)))
        for rank, (signature, paths) in enumerate(
            sorted(by_signature.items(), key=lambda item: (-len(item[1]), item[0])), 1)
    ]


def tandem_repeat_rows(signatures):
    """Report classes occurring at least twice consecutively on any path."""
    observations = defaultdict(lambda: {"paths": set(), "max_run": 0})
    for path_id, signature in signatures.items():
        run_id = None
        run_length = 0
        for bundle_id, _ in signature:
            if bundle_id == run_id:
                run_length += 1
            else:
                run_id = bundle_id
                run_length = 1
            if run_length >= 2:
                observations[bundle_id]["paths"].add(path_id)
                observations[bundle_id]["max_run"] = max(
                    observations[bundle_id]["max_run"], run_length
                )
    return [
        (bundle_id, len(values["paths"]), values["max_run"], ",".join(sorted(values["paths"])))
        for bundle_id, values in sorted(
            observations.items(),
            key=lambda item: (
                not item[0].isdigit(),
                int(item[0]) if item[0].isdigit() else item[0],
                item[0],
            ),
        )
    ]


def is_subsequence(shorter, longer):
    """Whether one signed-copy order can arise by deleting copies from another."""
    remaining = iter(longer)
    return all(any(token == observed for observed in remaining) for token in shorter)


def deleted_class_gaps(signature, bundle_id, removed_bundle_ids):
    """Locate a deleted class between consecutive surviving signed tokens.

    The number of surviving tokens already traversed identifies a gap in the
    full bundle order. For two signatures that collapse after deletion, the
    surviving token orders agree, so equal gap numbers describe the same
    structural context even when deleted copies occur at several loci.
    """
    gap = 0
    by_gap = defaultdict(list)
    for current_id, direction in signature:
        if current_id not in removed_bundle_ids:
            gap += 1
        elif current_id == bundle_id:
            by_gap[gap].append(direction)
    return by_gap


def classify_collapsing_pair(left, right, tandem_bundle_ids,
                             ignored_position_bundle_ids=frozenset()):
    """Classify the structural distinction erased by a proposed deletion."""
    shared = {bundle_id for bundle_id, _ in left} & {bundle_id for bundle_id, _ in right}
    left_shared = tuple(token for token in left if token[0] in shared)
    right_shared = tuple(token for token in right if token[0] in shared)
    if left_shared == right_shared:
        return "presence_absence"

    left_counts = Counter(bundle_id for bundle_id, _ in left_shared)
    right_counts = Counter(bundle_id for bundle_id, _ in right_shared)
    count_differences = {
        bundle_id for bundle_id in shared
        if left_counts[bundle_id] != right_counts[bundle_id]
    }

    # Non-deleted classes retain their original orientation protection.
    if any(
        {direction for current_id, direction in left_shared if current_id == bundle_id}
        != {direction for current_id, direction in right_shared if current_id == bundle_id}
        for bundle_id in count_differences - set(ignored_position_bundle_ids)
    ):
        return "inversion_or_orientation"

    # Compare deleted copies at their positions in the *whole* bundle order.
    # A -B copy lost in one gap and a +B copy gained in another are two indels,
    # not an inversion of the same occurrence. Conflicting orientations in a
    # shared gap remain protected unless one signed run is a subsequence of the
    # other, in which case extra copies can be explained as indels.
    ignored_shared = shared & set(ignored_position_bundle_ids)
    for bundle_id in ignored_shared:
        left_gaps = deleted_class_gaps(left, bundle_id, ignored_position_bundle_ids)
        right_gaps = deleted_class_gaps(right, bundle_id, ignored_position_bundle_ids)
        for gap in left_gaps.keys() | right_gaps.keys():
            left_directions = tuple(left_gaps.get(gap, ()))
            right_directions = tuple(right_gaps.get(gap, ()))
            if left_directions and right_directions and not (
                is_subsequence(left_directions, right_directions)
                or is_subsequence(right_directions, left_directions)
            ):
                return "inversion_or_orientation"

    # Copy-number protection is deliberately restricted to classes observed in
    # a consecutive tandem run somewhere in the original population.
    if count_differences & tandem_bundle_ids:
        return "copy_number"

    # A deleted class may occupy different repeated units without making that
    # placement a protected structural distinction. Contextual orientation
    # differences explainable by indels are also allowed.

    # Ignore non-tandem multiplicity and the placement of classes being
    # deleted while testing the order/orientation of all remaining classes.
    stable_shared = shared - count_differences - ignored_shared
    left_comparable = tuple(token for token in left_shared if token[0] in stable_shared)
    right_comparable = tuple(token for token in right_shared if token[0] in stable_shared)
    if left_comparable != right_comparable:
        if any(
            tuple(direction for current_id, direction in left_comparable if current_id == bundle_id)
            != tuple(direction for current_id, direction in right_comparable if current_id == bundle_id)
            for bundle_id in stable_shared
        ):
            return "inversion_or_orientation"
        return "position_or_order"
    if count_differences:
        return "non_tandem_copy_number"
    if ignored_shared:
        if any(
            tuple(direction for current_id, direction in left_shared if current_id == bundle_id)
            != tuple(direction for current_id, direction in right_shared if current_id == bundle_id)
            for bundle_id in ignored_shared
        ):
            return "contextual_indels"
        return "deleted_bundle_position"
    if any(
        tuple(direction for current_id, direction in left_shared if current_id == bundle_id)
        != tuple(direction for current_id, direction in right_shared if current_id == bundle_id)
        for bundle_id in shared
    ):
        return "inversion_or_orientation"
    return "position_or_order"


def check_cumulative_collapses(original_orders, previous_removed, candidate, cache,
                               tandem_bundle_ids):
    """Do not let later deletions hide same-locus orientation distinctions.

    For any newly collapsed original pair, surviving-class order and
    orientation remain protected. Deleted copies are compared within gaps of
    that surviving order, so changes at different loci can be indels. Copy
    number remains protected for a class with an original adjacent tandem run.
    """
    removed = previous_removed | {candidate}
    ignored_positions = frozenset(removed)
    groups = defaultdict(set)
    counts = Counter(original_orders.values())
    for original in counts:
        after = tuple(token for token in original if token[0] not in removed)
        groups[after].add(original)
    events = []
    protected = False
    for after, originals in sorted(groups.items()):
        before_current = {tuple(token for token in original if token[0] not in previous_removed)
                          for original in originals}
        if len(before_current) < 2:
            continue
        ordered = sorted(originals)
        reason = "presence_absence"
        for index, left in enumerate(ordered):
            for right in ordered[index + 1:]:
                key = ignored_positions, left, right
                if key not in cache:
                    cache[key] = classify_collapsing_pair(
                        left, right, tandem_bundle_ids, ignored_positions
                    )
                if cache[key] not in {
                    "presence_absence", "non_tandem_copy_number", "deleted_bundle_position",
                    "contextual_indels",
                }:
                    reason = cache[key]
                    break
                if cache[key] in {
                    "non_tandem_copy_number", "deleted_bundle_position", "contextual_indels"
                }:
                    reason = cache[key]
            if reason not in {
                "presence_absence", "non_tandem_copy_number", "deleted_bundle_position",
                "contextual_indels",
            }:
                break
        allowed = reason in {
            "presence_absence", "non_tandem_copy_number", "deleted_bundle_position",
            "contextual_indels",
        }
        protected = protected or not allowed
        events.append((order_text(after), len(originals), sum(counts[original] for original in originals),
                       int(allowed), reason, ";".join(order_text(original) for original in ordered)))
    return not protected, events


def mapping_coverage_bp(mappings, fixed_intervals=()):
    intervals = [
        (mapping.target_start, mapping.target_end) for mapping in mappings
    ]
    intervals.extend(fixed_intervals)
    return sum(end - start for start, end in union_intervals(intervals))


def filter_final_bundle_classes(mappings, assignments, sequence_lengths,
                                min_sequence_coverage=0.95, max_low_sequence_fraction=0.05,
                                max_bundle_bp=500, class_lengths=None,
                                fixed_coverage_intervals=None):
    """Try mapped classes no longer than ``max_bundle_bp``, shortest first."""
    if not 0 <= min_sequence_coverage <= 1 or not 0 <= max_low_sequence_fraction <= 1:
        raise ValueError("final sequence coverage thresholds must be in [0, 1]")
    if max_bundle_bp < 0:
        raise ValueError("final maximum bundle length must be non-negative")
    if not sequence_lengths or any(length <= 0 for length in sequence_lengths.values()):
        raise ValueError("final bundle filtering requires positive population FASTA sequence lengths")
    fixed_coverage_intervals = fixed_coverage_intervals or {}
    unknown_fixed_paths = set(fixed_coverage_intervals) - set(sequence_lengths)
    if unknown_fixed_paths:
        raise ValueError(
            "fixed coverage sequence is missing from population lengths: "
            f"{sorted(unknown_fixed_paths)[0]}"
        )
    fixed_by_path = {path_id: [] for path_id in sequence_lengths}
    for path_id, intervals in fixed_coverage_intervals.items():
        for start, end in intervals:
            if start < 0 or end <= start or end > sequence_lengths[path_id]:
                raise ValueError(f"fixed coverage coordinates exceed sequence length for {path_id}")
        fixed_by_path[path_id] = list(union_intervals(intervals))
    by_path = {path_id: [] for path_id in sequence_lengths}
    lengths_by_class = defaultdict(dict)
    for mapping in mappings:
        if mapping.target_name not in sequence_lengths:
            raise ValueError(f"mapped sequence {mapping.target_name} is missing from population FASTA")
        if mapping.target_start < 0 or mapping.target_end > sequence_lengths[mapping.target_name]:
            raise ValueError(f"mapping coordinates exceed FASTA length for {mapping.target_name}")
        by_path[mapping.target_name].append(mapping)
        class_id = assignments[mapping.bundle_name].numeric_bundle_id
        lengths_by_class[class_id][mapping.raw_query_name] = mapping.query_length
    inferred_lengths = {class_id: statistics.median(alleles.values())
                        for class_id, alleles in lengths_by_class.items()}
    if class_lengths is None:
        lengths = inferred_lengths
    else:
        missing = set(inferred_lengths) - set(class_lengths)
        if missing:
            raise ValueError(f"bundle length is missing for classes: {', '.join(sorted(missing)[:10])}")
        lengths = {class_id: class_lengths[class_id] for class_id in inferred_lengths}
        if any(length < 0 for length in lengths.values()):
            raise ValueError("bundle class lengths must be non-negative")
    ordered_classes = sorted(lengths, key=lambda class_id: (
        lengths[class_id], not class_id.isdigit(), int(class_id) if class_id.isdigit() else class_id, class_id))
    candidates = {class_id for class_id in ordered_classes if lengths[class_id] <= max_bundle_bp}
    before_signatures = signed_orders(mappings, assignments, sequence_lengths)
    tandem_rows = tandem_repeat_rows(before_signatures)
    tandem_bundle_ids = {row[0] for row in tandem_rows}
    current_signatures = before_signatures.copy()
    coverage = {
        path_id: mapping_coverage_bp(records, fixed_by_path[path_id])
        for path_id, records in by_path.items()
    }
    baseline = coverage.copy()
    def low(path_id, bp):
        return bp / sequence_lengths[path_id] < min_sequence_coverage
    total_paths = len(sequence_lengths)
    low_count = sum(low(path_id, bp) for path_id, bp in coverage.items())
    baseline_low_count = low_count
    decisions = []
    merge_events = []
    deleted_classes = set()
    collapse_cache = {}
    coverage_protected = haplotype_protected = haplotype_unchanged = 0
    for attempt, class_id in enumerate(ordered_classes, 1):
        if class_id not in candidates:
            occurrence_count = sum(
                assignments[record.bundle_name].numeric_bundle_id == class_id
                for records in by_path.values() for record in records
            )
            affected_count = sum(
                any(assignments[record.bundle_name].numeric_bundle_id == class_id for record in records)
                for records in by_path.values()
            )
            haplotype_count = len(set(current_signatures.values()))
            decisions.append((
                attempt, class_id, lengths[class_id], occurrence_count, affected_count, 0,
                low_count, f"{low_count / total_paths:.9f}", 1,
                haplotype_count, haplotype_count, 0, 1, "retain_above_max_bundle_bp",
            ))
            continue
        trial = {}
        trial_bp = {}
        trial_low_count = low_count
        removed_count = 0
        for path_id, records in by_path.items():
            remaining = [record for record in records
                         if assignments[record.bundle_name].numeric_bundle_id != class_id]
            if len(remaining) == len(records):
                continue
            removed_count += len(records) - len(remaining)
            trial[path_id] = remaining
            trial_bp[path_id] = mapping_coverage_bp(remaining, fixed_by_path[path_id])
            trial_low_count += int(low(path_id, trial_bp[path_id])) - int(low(path_id, coverage[path_id]))
        trial_low_fraction = trial_low_count / total_paths
        haplotype_allowed, events = check_cumulative_collapses(
            before_signatures, deleted_classes, class_id, collapse_cache, tandem_bundle_ids)
        before_count = len(set(current_signatures.values()))
        trial_signatures = {
            path_id: tuple(token for token in signature if token[0] != class_id)
            for path_id, signature in current_signatures.items()
        }
        after_count = len(set(trial_signatures.values()))
        haplotype_changed = after_count < before_count
        coverage_allowed = trial_low_fraction <= max_low_sequence_fraction
        lost_bp = sum(coverage[path_id] - bp for path_id, bp in trial_bp.items())
        if not haplotype_changed:
            action = "retain_haplotype_unchanged"
            haplotype_unchanged += 1
        elif not haplotype_allowed:
            action = "retain_haplotype_protected"
            haplotype_protected += 1
        elif not coverage_allowed:
            action = "retain_coverage_protected"
            coverage_protected += 1
        else:
            action = "delete"
            by_path.update(trial)
            coverage.update(trial_bp)
            low_count = trial_low_count
            current_signatures = trial_signatures
            deleted_classes.add(class_id)
        decisions.append((attempt, class_id, lengths[class_id], removed_count, len(trial), lost_bp,
                          trial_low_count, f"{trial_low_fraction:.9f}", int(coverage_allowed),
                          before_count, after_count, int(haplotype_changed),
                          int(haplotype_allowed), action))
        merge_events.extend((attempt, class_id, action, *event) for event in events)
    retained = [mapping for mapping in mappings
                if assignments[mapping.bundle_name].numeric_bundle_id not in deleted_classes]
    coverage_rows = [
        (path_id, sequence_lengths[path_id], baseline[path_id],
         f"{baseline[path_id] / sequence_lengths[path_id]:.9f}", coverage[path_id],
         f"{coverage[path_id] / sequence_lengths[path_id]:.9f}", int(low(path_id, coverage[path_id])))
        for path_id in sorted(sequence_lengths)
    ]
    return retained, {
        "decisions": decisions,
        "merge_events": merge_events,
        "haplotypes_before": haplotype_rows(before_signatures),
        "haplotypes_after": haplotype_rows(current_signatures),
        "tandem_bundle_classes": tandem_rows,
        "sequence_coverage": coverage_rows,
        "counts": {
            "mapped_bundle_classes": len(ordered_classes),
            "candidate_bundle_classes": len(candidates),
            "above_max_bundle_bp_classes": len(ordered_classes) - len(candidates),
            "deleted_bundle_classes": len(deleted_classes),
            "coverage_protected_bundle_classes": coverage_protected,
            "haplotype_protected_bundle_classes": haplotype_protected,
            "haplotype_unchanged_bundle_classes": haplotype_unchanged,
            "tandem_repeat_bundle_classes": len(tandem_bundle_ids),
            "sequence_count": total_paths,
            "baseline_low_coverage_sequences": baseline_low_count,
            "final_low_coverage_sequences": low_count,
            "haplotypes_before": len(set(before_signatures.values())),
            "haplotypes_after": len(set(current_signatures.values())),
            "mappings_before": len(mappings),
            "mappings_after": len(retained),
            "fixed_coverage_sequence_count": sum(bool(items) for items in fixed_by_path.values()),
            "fixed_coverage_bp": sum(
                end - start for items in fixed_by_path.values() for start, end in items
            ),
        },
    }


def write_final_filter_audits(directory, result, prefix="final"):
    specs = {
        "decisions": (f"{prefix}_bundle_filter.tsv", (
            "attempt_index", "bundle_id", "consensus_length_bp", "candidate_occurrence_count",
            "affected_sequence_count", "trial_lost_coverage_bp", "trial_low_sequence_count",
            "trial_low_sequence_fraction", "coverage_allowed", "haplotypes_before",
            "trial_haplotypes_after", "haplotype_changed", "haplotype_allowed", "action")),
        "merge_events": (f"{prefix}_bundle_haplotype_merges.tsv", (
            "attempt_index", "bundle_id", "action", "trial_signed_order", "before_haplotype_count",
            "sequence_count", "haplotype_merge_allowed", "reason", "original_signed_orders")),
        "haplotypes_before": (f"{prefix}_haplotypes.before.tsv", (
            "haplotype_id", "sequence_count", "signed_bundle_order", "path_ids")),
        "haplotypes_after": (f"{prefix}_haplotypes.after.tsv", (
            "haplotype_id", "sequence_count", "signed_bundle_order", "path_ids")),
        "tandem_bundle_classes": (f"{prefix}_tandem_bundle_classes.tsv", (
            "bundle_id", "sequence_count", "maximum_consecutive_copies", "path_ids")),
        "sequence_coverage": (f"{prefix}_sequence_coverage.tsv", (
            "path_id", "sequence_length_bp", "baseline_covered_bp", "baseline_coverage",
            "final_covered_bp", "final_coverage", "below_min_coverage")),
    }
    outputs = {}
    for key, (filename, header) in specs.items():
        path = directory / filename
        write_tsv(path, header, result[key])
        outputs[key] = str(path.resolve())
    return outputs
