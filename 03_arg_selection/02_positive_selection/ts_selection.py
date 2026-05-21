import argparse
import json
import math
import sys
import tskit
from scipy.special import gammaln
from typing import Callable, List, Dict, Tuple, Any, Optional
from collections import defaultdict
from tqdm import tqdm

# Absurdly large number.
MAX_GEN_TIME = 2**63


# Compute the number of lineages at every time t, such that t is the time of some node
# in the tree. Computation starts at the given node and traverses downwards.
#
# Definition: The number of lineages is equal to the number of active edges. An edge is active
# at time t if the source (parent) node of the edge has a time t(p) >= t and the destination (child)
# node of the edge has a time t(c) < t.
def lineages_by_time(nodes_time, tree):
    lineage_decrease_times = []
    samples_below = {}

    # We could get lineages faster, but we do it recursively here to get the number of samples
    # below each mutation as well. TODO: there are probably APIs for getting mutation frequency?
    def traverse(node):
        nonlocal lineage_decrease_times
        nonlocal samples_below
        if tree.is_sample(node):
            return 1
        t = nodes_time[node]
        children = tree.children(node)
        lineage_decrease_times.append((t, len(children) - 1))
        sb = 0
        for child in children:
            sb += traverse(child)
        samples_below[node] = sb
        return sb

    for r in tree.roots:
        traverse(r)
    lineage_decrease_times.sort(reverse=True)

    # Pairs (time, num_lineages) in descending time order, ascending lineages order. The time is the time when the number of
    # lineages increases to num_lineages.
    lineages = 1
    time2lins = {MAX_GEN_TIME: lineages}
    for t, c in lineage_decrease_times:
        lineages += c
        # Collapses entries with the same time.
        time2lins[t] = lineages

    # We have samples_below each node, now we need counts for mutations, which differs only when there are nested mutations
    # of the same site (i.e., recurrent).
    mut_samples = defaultdict(int)
    for mut in tree.mutations():
        mut_sample_ct = samples_below.get(mut.node, 1)
        mut_samples[mut.id] += mut_sample_ct
        if mut.parent != -1:
            mut_samples[mut.parent] -= mut_sample_ct

    return time2lins, mut_samples


def logfact(n):
    """
    The natural log of the factorial, which is implemented using the log-gamma function from scipy, since
    Gamma(n+1) == n!
    """
    return gammaln(n + 1)


def logchoose(n, k):
    """
    The natural log of (n choose k).
    """
    return logfact(n) - (logfact(k) + logfact(n - k))


def logprob(N, fn, k, x):
    """
    Natural log of the probability of observed the N, fn, k triple:
    * N is the number of samples (haplotypes) in the dataset.
    * fn is the number of samples carrying the mutation of interest in the dataset.
    * k is the number of lineages at the time the mutation increased from 1 to 2
      lineages carrying it.

    Note: the result is undefined (NaN) when k < 3, which implies that the mutation is
    on all currently-active lineages.
    """
    return logchoose(fn - 1, x - 1) + (
        logchoose(N - fn - 1, k - x - 1) - logchoose(N - 1, k - 1)
    )


def pvalue(N, fn, k, x):
    """
    The pvalue for the given parameters.

    :param N: Number of total samples.
    :param fn: Number of samples containing the mutation of interest.
    :param k: Number of lineages active at current time.
    :param x: Number of lineages beneath the mutation at current time.
    """
    # The number of samples beneath the mutation must be less than the number of lineages
    # remaining after removing the others ("N - (k - 1 - x)")
    #
    if fn >= (N - (k - 1 - x)):
        print(f"OUT OF RANGE: {fn} >= {N} - ({k} - 1 - {x})")
        return None
    pv = 0
    for f in range(fn, N - (k - 1 - x)):
        pv += math.exp(logprob(N, f, k, x))
    return pv


def selection_tree_stat(
    ts: tskit.TreeSequence,
    pv_threshold: float,
    mapped: Callable[[int], int] = lambda l: l,
    verbose: bool = True,
    compute_aggressive: bool = False,
    region: Optional[Tuple[int, int]] = None,
    max_mut_map_ratio: float = 2.5,
) -> Tuple[List[Dict[str, Any]], List[float], List[float]]:
    """
    Selection statistic from:
    L. Speidel, M. Forest, S. Shi, and S. R. Myers. A method for genome-wide genealogy estimation
    for thousands of samples. Nature genetics, 51(9):1321–1329, 2019.

    Extended slightly to handle more than 2 lineages beneath the mutation, because the general tskit
    TreeSequence may contain non-binary trees.

    :param ts: The tree sequence.
    :type ts: tskit.TreeSequence
    :param pv_threshold: Just used to reduce the output size. All pvalues (for all mutations) are returned,
        but only mutations that exceed this threshold will be returned with full details (statistics about
        the containing tree, N/Fn/k values, position, etc.)
    :type pv_threshold: float
    :param mapped: Function that maps positions to positions. For example, if you built a tree-sequence against
        GRCh38 coordinates you can pass in a function that maps to CHM13 coordinates, if you want to compare
        between GRCh38 and CHM13 datasets.
    :type mapped: Callable[[int], int]
    :param verbose: Whether to print information to stderr.
    :type verbose: bool
    :param compute_aggressive: When True, compute both the conservative (always computed) and aggressive forms
        of the p-value. The conservative p-value only looks at the oldest mutation in each set of equivalent
        mutations (same position and allele). The aggressive p-value looks at the oldest mutation, but considers
        the times and frequency of all equivalent mutations together, which can produce a smaller p-value (under-
        estimate). Default: False.
    :type compute_aggressive: bool
    :param max_mut_map_ratio: Default: 2.5. Do not compute p-values for any tree that has a mutation map ratio
        higher than this. Given a variant (position, allele) it can be mapped to multiple tskit Mutation objects
        (each having it's own node mapping). If, on average, each "variant" in the tree is mapped to more Mutations
        than this threshold, skip the tree.
    :type max_mut_map_ratio: float
    :return: (list_of_results, list_of_pvalues). The list_of_results contains a dictionary for each mutation
        that exceeded the pv_threshold. The list_of_pvalues contains a single floating point value for every
        mutation in the tree.
    :rtype: Tuple[List[Dict[str, Any]], List[float], List[float]]
    """
    if region is None:
        region = (0, 2**64)
    mapped_elsewhere = 0
    below_thresh_list = []
    conservative_pvalues = []
    aggressive_pvalues = []
    N = ts.num_samples
    nodes_time = ts.nodes_time

    def inspan(span: Tuple[int, int], value: int) -> bool:
        return value >= span[0] and value <= span[1]

    def overlaps(span1, span2):
        return (
            inspan(span1, span2[0])
            or inspan(span1, span2[1])
            or inspan(span2, span1[0])
            or inspan(span2, span1[1])
        )

    for i, tree in tqdm(enumerate(ts.trees()), total=ts.num_trees):
        if not overlaps(tree.interval, region):
            continue

        # Collect all mutations by (site, allele) combination and then flag the oldest in each group
        collated_muts = defaultdict(list)
        for mut in tree.mutations():
            collated_muts[(mut.site, mut.derived_state)].append(mut)

        muts_per_variant = []
        for _, mut_list in collated_muts.items():
            mut_list.sort(key=lambda mut: mut.time, reverse=True)
            muts_per_variant.append(len(mut_list))

        # Skip poorly mapped mutations
        if muts_per_variant:
            mut_map_ratio = sum(muts_per_variant) / len(muts_per_variant)
            if mut_map_ratio > max_mut_map_ratio:
                continue
        else:
            mut_map_ratio = 0

        time2lins, mut_samples = lineages_by_time(nodes_time, tree)
        for mut_list in collated_muts.values():
            # Representative mutation: the oldest one in each set.
            rep_mut = mut_list[0]
            # The time that the tree below the mutation splits into more than one lineage is the time of the node below it
            num_children = len(tree.children(rep_mut.node))
            if num_children >= 2:
                oldest_mut_split_time = nodes_time[rep_mut.node]
                k = time2lins[oldest_mut_split_time]
                if k < N:
                    fn_under = mut_samples[rep_mut.id]
                    fn_over = sum(map(lambda m: mut_samples[m.id], mut_list))
                    assert fn_over <= N

                    # TODO: Improve the handling of recurrent mutations: we should have a "blockers" list,
                    # and ignore children that are blocked by another mutation at the site. This should be an
                    # extremely rare case, because we are only looking at the nodes immediately below a mutation,
                    # so it would have to have another mutation on the same site as a direct child. For now,
                    # manual inspection of the results is sufficient.

                    site = ts.site(rep_mut.site)
                    mut_pos = mapped(site.position)
                    if mut_pos != site.position:
                        mapped_elsewhere += 1

                    num_children_aggressive = num_children
                    k_aggressive = k
                    if compute_aggressive and len(mut_list) > 1:
                        next_time = mut_list[1].time
                        assert next_time <= rep_mut.time
                        # The other equivalent mutation is older than the time of our split, so it goes to
                        # two lineages at _that_ time instead.
                        if next_time > oldest_mut_split_time:
                            num_children_aggressive = 2
                            lineages = 1
                            inc_times = list(sorted(time2lins.keys()))
                            for i in range(len(inc_times)):
                                if inc_times[i] < next_time:
                                    break
                                lineages = time2lins[inc_times[i]]
                            k_aggressive = lineages

                    entry = {
                        "tree": i,
                        "tree_muts": tree.num_mutations,
                        "tree_mut_map_ratio": mut_map_ratio,
                        "N": N,
                        "k": k,
                        "k_aggressive": k_aggressive,
                        "x": num_children,
                        "x_aggressive": num_children_aggressive,
                        "tree_start_unmapped": tree.interval[0],
                        "tree_start": mapped(tree.interval[0]),
                        "mut_pos": mut_pos,
                        "mut_allele": rep_mut.derived_state,
                    }

                    # Since this statistic has to do with lineages with/without the
                    # mutation, it doesn't tell us anything if the mutation has fixed.
                    below_thresh = False
                    if fn_under < N and fn_under > 0:
                        pv = pvalue(N, fn_under, k, num_children)
                        if pv is not None:
                            if pv <= pv_threshold:
                                below_thresh = True
                                entry["pvalue_conservative"] = pv
                            conservative_pvalues.append(pv)
                    if (
                        compute_aggressive
                        and fn_over < N
                        and fn_over > 0
                        and k_aggressive < N
                    ):
                        pv = pvalue(N, fn_over, k_aggressive, num_children_aggressive)
                        if pv is not None:
                            if pv <= pv_threshold:
                                below_thresh = True
                                entry["pvalue_aggressive"] = pv
                            aggressive_pvalues.append(pv)
                    if below_thresh:
                        below_thresh_list.append(entry)

    if verbose:
        print(
            f"Mapped {mapped_elsewhere} mutations to new positions based on input map",
            file=sys.stderr,
        )
    return below_thresh_list, conservative_pvalues, aggressive_pvalues


def region_string(input: str) -> Tuple[int, int]:
    try:
        lower, upper = input.split("-")
        lower, upper = int(lower), int(upper)
    except ValueError:
        raise RuntimeError(f"Invalid region string: '{input}', must be 'lower-upper'")
    return lower, upper


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree_sequence", help="tskit tree sequence file")
    parser.add_argument(
        "pv_threshold",
        type=float,
        help="Only report _detailed_ results for mutations below this p-value",
    )
    parser.add_argument(
        "--aggressive",
        action="store_true",
        help="Compute aggressive p-values (in addition to conservative)",
    )
    parser.add_argument(
        "--region",
        type=region_string,
        help="Only compute for trees overlapping with the given region (lower, upper)",
    )
    args = parser.parse_args()

    ts = tskit.load(args.tree_sequence)

    below_thresh, conservative_pvalues, aggressive_pvalues = selection_tree_stat(
        ts,
        args.pv_threshold,
        lambda l: l,
        compute_aggressive=args.aggressive,
        region=args.region,
    )

    print(
        json.dumps(
            {
                "conservative_pvalues": conservative_pvalues,
                "aggressive_pvalues": aggressive_pvalues,
                "below_threshold": below_thresh,
            },
            indent=2,
        )
    )
