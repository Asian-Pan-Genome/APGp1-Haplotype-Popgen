# Generate figures for the tree at a position.
import itertools
import tskit
import sys
import os

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: ts_export_newick.py <tree seq> <position>", file=sys.stderr)
        exit(1)
    TS_FILE = sys.argv[1]
    POS = int(sys.argv[2])
    scale = sys.argv[3] if len(sys.argv) > 3 else "time"

    ts = tskit.load(TS_FILE)

    def get_tree_for_position(ts: tskit.TreeSequence, position: int):
        for tree in ts.trees():
            if tree.interval[0] <= position and tree.interval[1] > position:
                return tree

    def get_first_pos_for_node(ts: tskit.TreeSequence, node: int):
        position = None
        for e in ts.tables.edges:
            if e.parent == node or e.child == node:
                position = e.left
                break
        assert position is not None
        return position

    def get_muts_for_position(ts: tskit.TreeSequence, tree: tskit.Tree, position: int):
        result = []
        for mut in tree.mutations():
            site = ts.site(mut.site)
            if site.position == position:
                result.append(mut)
        return result

    if "USE_NODE_ID" in os.environ:
        position = get_first_pos_for_node(ts, POS)
    else:
        position = POS
    tree = get_tree_for_position(ts, position)

    # Lets see how long-lived the node below the mutations are.
    el = ts.edges_left
    er = ts.edges_right
    ep = ts.edges_parent
    ec = ts.edges_child

    def get_node_spans(node):
        edges = set()
        assert len(ep) == len(ec)
        for i in range(len(ep)):
            if ep[i] == node:
                edges.add(i)
            if ec[i] == node:
                edges.add(i)

        spans = []
        for e in edges:
            spans.append([el[e], er[e]])
        spans.sort()
        new_spans = []
        for s in spans:
            if new_spans and new_spans[-1][1] == s[0]:
                new_spans[-1][1] = s[1]
            else:
                new_spans.append(s)
        return list(map(tuple, new_spans))

    muts = get_muts_for_position(ts, tree, position)
    print("Mutations:")
    for m in muts:
        print(m)
        print(f"  SPANS: {get_node_spans(m.node)}")

    # Now get the smaller tree-sequence for only the tree interval.
    ts2 = ts.keep_intervals([[tree.interval[0], tree.interval[1]]])
    tree = get_tree_for_position(ts2, position)
    muts = get_muts_for_position(ts2, tree, position)
    print(f"Showing {ts2.num_trees} trees")

    roots = tree.roots
    print(f"Roots: {roots}")
    root_times = [ts2.nodes_time[r] for r in roots]
    print(f"Root times: {root_times}")
    # Mutation node times (below the mut)
    node_times = [ts2.nodes_time[m.node] for m in muts]
    print(f"Node times: {node_times}")

    # Get the node _above_ the mutation
    parent_nodes = [tree.parent(m.node) for m in muts]
    print(f"Parents: {parent_nodes}")
    parent_times = [ts2.nodes_time[p] for p in parent_nodes]
    print(f"Parent times: {parent_times}")

    i = 0
    label = "Selected"
    node_labels = {m.node: label for m in muts}
    for n in range(ts2.num_nodes):
        if tree.is_sample(n):
            node_labels[n] = f"s{i}"
            i += 1

    newick = tree.as_newick(node_labels=node_labels)
    with open(f"{os.path.basename(TS_FILE)}.{POS}.txt", "w") as fout:
        fout.write(newick)
