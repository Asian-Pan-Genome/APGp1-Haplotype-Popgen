import tskit
import sys
import json
import numpy as np
from tqdm import tqdm
from collections import defaultdict


def coalescences(nodes_time, tree: tskit.Tree, root: int):
    min_tmrca = 2**64
    max_tmrca = 0
    avg_tmrca = 0
    pairs = 0

    def sb(tree, node):
        nonlocal min_tmrca
        nonlocal max_tmrca
        nonlocal avg_tmrca
        nonlocal pairs
        if tree.is_sample(node):
            return 1
        below = 0
        prod = 1
        for c in tree.children(node):
            b = sb(tree, c)
            below += b
            prod *= b
        t = nodes_time[node]
        min_tmrca = min(min_tmrca, t)
        max_tmrca = max(max_tmrca, t)
        avg_tmrca += t * prod
        pairs += prod
        return below

    sb(tree, root)
    if pairs != 0:
        return pairs, min_tmrca, max_tmrca, avg_tmrca / pairs
    else:
        return pairs, None, None, None


if __name__ == "__main__":
    with open("WINDOW_SIZE") as f:
        BP_WIN_SIZE = int(f.read())
    TIME_WINDOWS = np.array([0, float("inf")])

    coal_rates = []
    coal_counts = []
    avg_pw_times = []
    window_max_times = []
    window_min_times = []
    window_avg_times = []
    window_avg_pw_times = []
    window_max_ttotal = []
    window_min_ttotal = []
    window_avg_ttotal = []
    print(
        "Filename,NumTrees,BranchDiversity,CoalRate,Ne,Muts,UniqueMuts,NumSites,SeqLen,RecurrentMuts"
    )
    prefix = sys.argv[1]
    print(f"Using prefix {prefix}", file=sys.stderr)
    for filename in sys.argv[2:]:
        parts = filename.split(":")

        # You can optionally pass in a coordinate map that maps the IGD file's positions
        # to another coordinate system.
        mapped = None
        if len(parts) > 1:
            assert len(parts) == 2, f"Unexpected ':' in filename: {filename}"
            with open(parts[1]) as f:
                mapped = {int(k): int(v) for k, v in json.load(f).items()}
            filename = parts[0]
            print(f"Using coordinate map {parts[1]} for {filename}", file=sys.stderr)

        ts = tskit.load(filename)

        ##############################################################################
        # Simple statistics

        seen_pos = set()
        unique_muts = set()
        site2alleles = defaultdict(set)
        dups = 0
        recurrent = 0
        for site in ts.sites():
            # There should not be sites with different ancestral alleles (not possible)
            assert site.position not in seen_pos
            seen_pos.add(site.position)

            for mutation in site.mutations:
                if mutation.derived_state == site.ancestral_state:
                    recurrent += 1
                else:
                    site2alleles[site.position].add(mutation.derived_state)
                    unique_muts.add((site.position, mutation.derived_state))

        branch_diversity = ts.diversity(mode="branch")
        coal_rate = ts.pair_coalescence_rates(TIME_WINDOWS)[0]
        ne = 1 / (2 * coal_rate)

        print(
            f"{filename},{ts.num_trees},{branch_diversity},{coal_rate},{ne},{ts.num_mutations},{len(unique_muts)},{ts.num_sites},{ts.sequence_length},{recurrent}"
        )

        ##############################################################################
        # Per-tree/window statistics

        # Construct added/deleted edge dictionaries for computing T_total later.
        added = defaultdict(list)
        deleted = defaultdict(list)
        edges = list(ts.tables.edges)
        for edge in edges:
            added[edge.left].append(edge)
            deleted[edge.right].append(edge)

        # Window-based TMRCA statistics.
        win_tmrca_max = defaultdict(float)
        win_tmrca_min = defaultdict(lambda: float(2**64))
        win_tmrca_avg = defaultdict(float)
        win_tmrca_pw_avg = defaultdict(float)
        # Window-based T_total statistics.
        win_ttotal_max = defaultdict(float)
        win_ttotal_min = defaultdict(lambda: float(2**64))
        win_ttotal_avg = defaultdict(float)

        total_span = defaultdict(float)
        skipped_trees = 0

        prev_t_total = 0.0

        nodes_time = ts.nodes_time
        for tree in tqdm(ts.trees(), total=ts.num_trees):
            # Info about tree location
            tree_start = tree.interval[0]
            tree_end = tree.interval[1]
            if mapped is not None:
                mapped_tree_start = mapped.get(tree_start)
            else:
                mapped_tree_start = tree_start

            # TMRCA of the root(s)
            min_tmrca = 2**64
            max_tmrca = 0
            avg_pw_tmrca = 0
            for r in tree.roots:
                r_num_pairs, r_min_tmrca, r_max_tmrca, r_avg_tmrca = coalescences(
                    nodes_time, tree, r
                )
                if r_num_pairs == 0:
                    continue

                if r_max_tmrca > max_tmrca:
                    max_tmrca = r_max_tmrca
                if r_min_tmrca < min_tmrca:
                    min_tmrca = r_min_tmrca
                avg_pw_tmrca = r_avg_tmrca
            if tree.roots:
                avg_pw_tmrca /= len(tree.roots)
            if mapped_tree_start is None:
                skipped_trees += 1
            else:
                avg_pw_times.append(
                    {
                        "filename": filename,
                        "tree_start": mapped_tree_start,
                        "avg_pw_tmrca": avg_pw_tmrca,
                    }
                )

            # T_total (compute based on changed edges only)
            current_t_total = prev_t_total
            for edge in added[tree_start]:
                edge_length = ts.nodes_time[edge.parent] - ts.nodes_time[edge.child]
                current_t_total += edge_length
            for edge in deleted[tree_start]:
                edge_length = ts.nodes_time[edge.parent] - ts.nodes_time[edge.child]
                current_t_total -= edge_length
            prev_t_total = current_t_total

            # Record window-based statistics.
            prev_position = None
            for position in range(int(tree_start), int(tree_end)):
                if mapped is not None:
                    mapped_pos = mapped.get(position)
                else:
                    mapped_pos = position
                if mapped_pos is not None:
                    window = mapped_pos // BP_WIN_SIZE
                    win_tmrca_max[window] = max(win_tmrca_max[window], max_tmrca)
                    win_tmrca_min[window] = min(win_tmrca_min[window], max_tmrca)
                    win_ttotal_max[window] = max(
                        win_ttotal_max[window], current_t_total
                    )
                    win_ttotal_min[window] = min(
                        win_ttotal_min[window], current_t_total
                    )
                    if prev_position is not None:
                        span = position - prev_position
                        assert span > 0
                        total_span[window] += span
                        win_tmrca_avg[window] += span * max_tmrca
                        win_tmrca_pw_avg[window] += span * avg_pw_tmrca
                        win_ttotal_avg[window] += span * current_t_total
                    prev_position = position

        for window in win_tmrca_max.keys():
            window_max_times.append(
                {
                    "filename": filename,
                    "window_start": window * BP_WIN_SIZE,
                    "window_size": BP_WIN_SIZE,
                    "max_tmrca": win_tmrca_max[window],
                }
            )
            window_min_times.append(
                {
                    "filename": filename,
                    "window_start": window * BP_WIN_SIZE,
                    "window_size": BP_WIN_SIZE,
                    "min_tmrca": win_tmrca_min[window],
                }
            )
            # This is the average of the root TMRCA
            window_avg_times.append(
                {
                    "filename": filename,
                    "window_start": window * BP_WIN_SIZE,
                    "window_size": BP_WIN_SIZE,
                    "avg_tmrca": win_tmrca_avg[window] / total_span[window],
                }
            )
            # This is the average of the pairwise TMRCAs
            window_avg_pw_times.append(
                {
                    "filename": filename,
                    "window_start": window * BP_WIN_SIZE,
                    "window_size": BP_WIN_SIZE,
                    "avg_tmrca": win_tmrca_pw_avg[window] / total_span[window],
                }
            )
            window_max_ttotal.append(
                {
                    "filename": filename,
                    "window_start": window * BP_WIN_SIZE,
                    "window_size": BP_WIN_SIZE,
                    "max_ttotal": win_ttotal_max[window],
                }
            )
            window_min_ttotal.append(
                {
                    "filename": filename,
                    "window_start": window * BP_WIN_SIZE,
                    "window_size": BP_WIN_SIZE,
                    "min_ttotal": win_ttotal_min[window],
                }
            )
            window_avg_ttotal.append(
                {
                    "filename": filename,
                    "window_start": window * BP_WIN_SIZE,
                    "window_size": BP_WIN_SIZE,
                    "avg_ttotal": win_ttotal_avg[window] / total_span[window],
                }
            )

        # We no longer calculation coalescence rate (using ts.pair_coalescence_rates) because
        # we need to do the calculation of windows across different coordinate systems.

        # Coalescence counts by time, across all samples, and across the entire genome.
        generations = np.array(np.arange(0, 500_000, 1_000))
        time_coal_counts = ts.pair_coalescence_counts(
            time_windows=generations, pair_normalise=True
        )
        cumu_coal_counts = np.cumsum(time_coal_counts)
        coal_counts.extend(
            [
                {
                    "filename": filename,
                    "coals": float(c),
                    "cumu_coals": float(cs),
                    "generation": int(gen),
                }
                for gen, c, cs in zip(generations, time_coal_counts, cumu_coal_counts)
            ]
        )

    with open(f"{prefix}.ts.pw_mrca.json", "w") as fout:
        fout.write(json.dumps(avg_pw_times, indent=2))
    with open(f"{prefix}.ts.winmrca.max.json", "w") as fout:
        fout.write(json.dumps(window_max_times, indent=2))
    with open(f"{prefix}.ts.winmrca.min.json", "w") as fout:
        fout.write(json.dumps(window_min_times, indent=2))
    with open(f"{prefix}.ts.winmrca.avg.json", "w") as fout:
        fout.write(json.dumps(window_avg_times, indent=2))
    with open(f"{prefix}.ts.winmrca.pwavg.json", "w") as fout:
        fout.write(json.dumps(window_avg_pw_times, indent=2))

    with open(f"{prefix}.ts.winttotal.max.json", "w") as fout:
        fout.write(json.dumps(window_max_ttotal, indent=2))
    with open(f"{prefix}.ts.winttotal.min.json", "w") as fout:
        fout.write(json.dumps(window_min_ttotal, indent=2))
    with open(f"{prefix}.ts.winttotal.avg.json", "w") as fout:
        fout.write(json.dumps(window_avg_ttotal, indent=2))

    with open(f"{prefix}.ts.coal_counts.json", "w") as fout:
        fout.write(json.dumps(coal_counts, indent=2))
