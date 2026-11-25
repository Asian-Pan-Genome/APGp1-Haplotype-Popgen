#!/usr/bin/env python3
"""Compare signed PGR-TK bundle haplotypes and draw a Sankey diagram.

Haplotypes are exact ordered tuples of (bundle ID, direction) on each BED path.
Coordinates and PGR-TK R/U flags do not define separate haplotypes.  Each BED
gets its own class IDs, calculated only from path IDs shared by both BEDs.
The diagram may aggregate rare classes into Other, while TSV outputs always
retain every class and every shared-path transition.
"""

from __future__ import annotations

import argparse
import colorsys
import hashlib
import html
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

try:
    from .merge_pgr_bundles import read_pgr_bed, write_tsv
except ImportError:
    from merge_pgr_bundles import read_pgr_bed, write_tsv


Signature = Tuple[Tuple[str, str], ...]


@dataclass(frozen=True)
class Haplotype:
    hap_id: str
    signature: Signature
    paths: Tuple[str, ...]


@dataclass(frozen=True)
class PlotNode:
    node_id: str
    count: int
    member_ids: Tuple[str, ...]


def signature_text(signature: Signature) -> str:
    return ",".join(f"{direction}{bundle_id}" for bundle_id, direction in signature)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_haplotypes(by_path, prefix: str, shared: set):
    groups = defaultdict(list)
    for path in sorted(shared):
        records = by_path[path]
        signature = tuple(record.token for record in records)
        groups[signature].append(path)
    ordered = sorted(
        groups,
        key=lambda signature: (-len(groups[signature]), signature_text(signature)),
    )
    haplotypes = []
    path_to_haplotype = {}
    for rank, signature in enumerate(ordered, 1):
        paths = tuple(sorted(groups[signature]))
        hap_id = f"{prefix}{rank:04d}"
        haplotypes.append(Haplotype(hap_id, signature, paths))
        path_to_haplotype.update((path, hap_id) for path in paths)
    return haplotypes, path_to_haplotype


def plot_nodes(haplotypes: Sequence[Haplotype], prefix: str, limit: int):
    used = sorted(
        haplotypes,
        key=lambda hap: (-len(hap.paths), hap.hap_id),
    )
    kept = used[:limit] if limit else used
    omitted = used[len(kept):]
    nodes = [PlotNode(hap.hap_id, len(hap.paths), (hap.hap_id,)) for hap in kept]
    mapping = {hap.hap_id: hap.hap_id for hap in kept}
    if omitted:
        other_id = f"{prefix}_OTHER"
        nodes.append(PlotNode(
            other_id,
            sum(len(hap.paths) for hap in omitted),
            tuple(hap.hap_id for hap in omitted),
        ))
        mapping.update((hap.hap_id, other_id) for hap in omitted)
    return nodes, mapping


def node_color(index: int, other: bool) -> str:
    if other:
        return "#8290a3"
    hue = (index * 0.61803398875 + 0.06) % 1.0
    rgb = colorsys.hsv_to_rgb(hue, 0.57, 0.82)
    return "#" + "".join(f"{round(channel * 255):02x}" for channel in rgb)


def sankey_svg(
    left_nodes: Sequence[PlotNode],
    right_nodes: Sequence[PlotNode],
    edges: Mapping[Tuple[str, str], int],
    total: int,
    left_label: str,
    right_label: str,
) -> str:
    width = 1180
    left_x, right_x, node_width = 270, 890, 22
    top, gap = 76, 12
    largest_side = max(len(left_nodes), len(right_nodes))
    chart_height = max(620, largest_side * 27, min(1800, total * 2))
    height = top + chart_height + 55
    scale = min(
        (chart_height - gap * (len(nodes) - 1)) / total
        for nodes in (left_nodes, right_nodes)
    )
    if scale <= 0:
        raise ValueError("too many plotted classes for the available chart height")

    def positions(nodes):
        used_height = total * scale + gap * (len(nodes) - 1)
        cursor = top + (chart_height - used_height) / 2
        result = {}
        for node in nodes:
            result[node.node_id] = (cursor, node.count * scale)
            cursor += node.count * scale + gap
        return result

    left_pos, right_pos = positions(left_nodes), positions(right_nodes)
    left_index = {node.node_id: index for index, node in enumerate(left_nodes)}
    right_index = {node.node_id: index for index, node in enumerate(right_nodes)}
    left_cursor = {node.node_id: left_pos[node.node_id][0] for node in left_nodes}
    right_cursor = {node.node_id: right_pos[node.node_id][0] for node in right_nodes}
    left_colors = {
        node.node_id: node_color(index, node.node_id.endswith("_OTHER"))
        for index, node in enumerate(left_nodes)
    }
    right_colors = {
        node.node_id: node_color(index + 3, node.node_id.endswith("_OTHER"))
        for index, node in enumerate(right_nodes)
    }
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="{width}" height="{height}" role="img" '
        f'aria-label="Haplotype correspondence Sankey diagram">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{left_x}" y="38" text-anchor="middle" class="heading">'
        f'{html.escape(left_label)}</text>',
        f'<text x="{right_x + node_width}" y="38" text-anchor="middle" class="heading">'
        f'{html.escape(right_label)}</text>',
        f'<text x="{width / 2}" y="38" text-anchor="middle" class="count">'
        f'{total} shared paths</text>',
    ]
    middle_x = (left_x + node_width + right_x) / 2
    for (left_id, right_id), count in sorted(
        edges.items(), key=lambda item: (
            left_index[item[0][0]], right_index[item[0][1]],
        ),
    ):
        thickness = count * scale
        y_left, y_right = left_cursor[left_id], right_cursor[right_id]
        left_cursor[left_id] += thickness
        right_cursor[right_id] += thickness
        path = (
            f"M {left_x + node_width:.2f} {y_left:.2f} "
            f"C {middle_x:.2f} {y_left:.2f}, {middle_x:.2f} {y_right:.2f}, "
            f"{right_x:.2f} {y_right:.2f} "
            f"L {right_x:.2f} {y_right + thickness:.2f} "
            f"C {middle_x:.2f} {y_right + thickness:.2f}, "
            f"{middle_x:.2f} {y_left + thickness:.2f}, "
            f"{left_x + node_width:.2f} {y_left + thickness:.2f} Z"
        )
        parts.append(
            f'<path d="{path}" fill="{left_colors[left_id]}" fill-opacity="0.42">'
            f'<title>{html.escape(left_id)} → {html.escape(right_id)}: {count} paths</title>'
            '</path>'
        )
    for nodes, locations, x, colors, label_x, anchor in (
        (left_nodes, left_pos, left_x, left_colors, left_x - 9, "end"),
        (right_nodes, right_pos, right_x, right_colors, right_x + node_width + 9, "start"),
    ):
        for node in nodes:
            y, node_height = locations[node.node_id]
            label = f"{node.node_id} ({node.count})"
            members = ", ".join(node.member_ids)
            parts.append(
                f'<rect x="{x}" y="{y:.2f}" width="{node_width}" '
                f'height="{node_height:.2f}" rx="2" fill="{colors[node.node_id]}">'
                f'<title>{html.escape(label)}; classes: {html.escape(members)}</title>'
                '</rect>'
            )
            if node_height >= 12:
                parts.append(
                    f'<text x="{label_x}" y="{y + node_height / 2 + 4:.2f}" '
                    f'text-anchor="{anchor}" class="node-label">'
                    f'{html.escape(label)}</text>'
                )
    parts.append(
        '<style>.heading{font:600 18px sans-serif;fill:#172238}'
        '.count{font:13px sans-serif;fill:#526176}'
        '.node-label{font:12px sans-serif;fill:#172238}</style></svg>'
    )
    return "\n".join(parts) + "\n"


def sankey_html(
    svg: str,
    left_nodes: Sequence[PlotNode],
    right_nodes: Sequence[PlotNode],
    left_label: str,
    right_label: str,
    total: int,
) -> str:
    rows = []
    for side, nodes in ((left_label, left_nodes), (right_label, right_nodes)):
        for node in nodes:
            rows.append(
                "<tr><td>" + html.escape(side) + "</td><td>" +
                html.escape(node.node_id) + "</td><td>" + str(node.count) +
                "</td><td>" + html.escape(", ".join(node.member_ids)) + "</td></tr>"
            )
    return (
        "<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
        "<title>PGR-TK haplotype correspondence</title>"
        "<style>body{font:14px system-ui,sans-serif;color:#172238;max-width:1240px;"
        "margin:24px auto;padding:0 16px}svg{max-width:100%;height:auto;border:1px solid #dde3ea}"
        "table{border-collapse:collapse;width:100%}th,td{border-bottom:1px solid #dde3ea;"
        "padding:6px 9px;text-align:left}td:last-child{overflow-wrap:anywhere}</style>"
        "</head><body><h1>PGR-TK haplotype correspondence</h1><p>"
        + str(total) + " shared paths. Hover over a ribbon or node for counts. "
        "The TSV files contain all unaggregated haplotypes and transitions.</p>"
        + svg + "<h2>Plotted classes</h2><table><thead><tr><th>BED</th><th>Node</th>"
        "<th>Shared paths</th><th>Original class IDs</th></tr></thead><tbody>"
        + "".join(rows) + "</tbody></table></body></html>\n"
    )


def compare_pgrtk_haplotypes(
    bed_a: Path,
    bed_b: Path,
    output_prefix: Path,
    label_a: Optional[str] = None,
    label_b: Optional[str] = None,
    max_plotted_haplotypes: int = 20,
):
    bed_a, bed_b, output_prefix = Path(bed_a), Path(bed_b), Path(output_prefix)
    if not bed_a.is_file() or not bed_b.is_file():
        raise ValueError("both PGR-TK BED inputs must exist")
    if max_plotted_haplotypes < 0:
        raise ValueError("--max-plotted-haplotypes must be nonnegative")
    label_a = label_a or bed_a.stem
    label_b = label_b or bed_b.stem
    output_names = (
        "haplotypes_a.tsv", "haplotypes_b.tsv", "path_assignments.tsv",
        "transitions.tsv", "plot_transitions.tsv", "unmatched_paths.tsv",
        "representative_a.tsv", "representative_b.tsv",
        "sankey.svg", "sankey.html", "summary.json",
    )
    outputs = {name: Path(f"{output_prefix}.{name}") for name in output_names}
    protected = {bed_a.resolve(), bed_b.resolve()}
    if any(path.resolve() in protected for path in outputs.values()):
        raise ValueError("an output path may not overwrite an input BED")

    by_a, by_b = read_pgr_bed(bed_a), read_pgr_bed(bed_b)
    paths_a, paths_b = set(by_a), set(by_b)
    shared = paths_a & paths_b
    if not shared:
        raise ValueError("the BED files have no shared path IDs; no correspondence can be drawn")
    hap_a, path_to_a = build_haplotypes(by_a, "A", shared)
    hap_b, path_to_b = build_haplotypes(by_b, "B", shared)
    transition_paths = defaultdict(list)
    for path in sorted(shared):
        transition_paths[(path_to_a[path], path_to_b[path])].append(path)

    left_nodes, left_map = plot_nodes(hap_a, "A", max_plotted_haplotypes)
    right_nodes, right_map = plot_nodes(hap_b, "B", max_plotted_haplotypes)
    plot_edges = Counter()
    for (left_id, right_id), paths in transition_paths.items():
        plot_edges[(left_map[left_id], right_map[right_id])] += len(paths)

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for name, haplotypes in (
        ("haplotypes_a.tsv", hap_a),
        ("haplotypes_b.tsv", hap_b),
    ):
        write_tsv(
            outputs[name],
            ("haplotype_id", "path_count", "path_frequency", "representative_path_id",
             "signed_bundle_order", "path_ids"),
            (
                (hap.hap_id, len(hap.paths), f"{len(hap.paths) / len(shared):.9f}",
                 hap.paths[0], signature_text(hap.signature), ",".join(hap.paths))
                for hap in haplotypes
            ),
        )
    for name, haplotypes, by_path in (
        ("representative_a.tsv", hap_a, by_a),
        ("representative_b.tsv", hap_b, by_b),
    ):
        write_tsv(
            outputs[name],
            ("molecule", "start", "end", "gene", "orientation"),
            (
                (f"{hap.hap_id}({len(hap.paths)})", record.start, record.end,
                 record.bundle_id, 0 if record.direction == "+" else 1)
                for hap in haplotypes
                for record in by_path[hap.paths[0]]
            ),
        )
    write_tsv(
        outputs["path_assignments.tsv"],
        ("path_id", "status", "haplotype_a", "haplotype_b"),
        (
            (path, "shared", path_to_a[path], path_to_b[path])
            for path in sorted(shared)
        ),
    )
    write_tsv(
        outputs["transitions.tsv"],
        ("haplotype_a", "haplotype_b", "path_count", "shared_path_frequency", "path_ids"),
        (
            (left_id, right_id, len(paths), f"{len(paths) / len(shared):.9f}",
             ",".join(paths))
            for (left_id, right_id), paths in sorted(transition_paths.items())
        ),
    )
    write_tsv(
        outputs["plot_transitions.tsv"],
        ("plot_haplotype_a", "plot_haplotype_b", "path_count", "shared_path_frequency"),
        (
            (left_id, right_id, count, f"{count / len(shared):.9f}")
            for (left_id, right_id), count in sorted(plot_edges.items())
        ),
    )
    write_tsv(
        outputs["unmatched_paths.tsv"],
        ("bed", "path_id"),
        [("a", path) for path in sorted(paths_a - shared)]
        + [("b", path) for path in sorted(paths_b - shared)],
    )
    svg = sankey_svg(left_nodes, right_nodes, plot_edges, len(shared), label_a, label_b)
    outputs["sankey.svg"].write_text(svg, encoding="utf-8")
    outputs["sankey.html"].write_text(
        sankey_html(svg, left_nodes, right_nodes, label_a, label_b, len(shared)),
        encoding="utf-8",
    )
    summary = {
        "inputs": {
            "bed_a": str(bed_a.resolve()), "bed_a_sha256": file_sha256(bed_a),
            "bed_b": str(bed_b.resolve()), "bed_b_sha256": file_sha256(bed_b),
        },
        "parameters": {
            "label_a": label_a, "label_b": label_b,
            "max_plotted_haplotypes_per_side": max_plotted_haplotypes,
            "haplotype_rule": "exact signed bundle order within each BED",
            "path_selection": "intersection of BED path IDs",
        },
        "counts": {
            "bed_a_path_count": len(paths_a), "bed_b_path_count": len(paths_b),
            "shared_path_count": len(shared),
            "bed_a_only_path_count": len(paths_a - shared),
            "bed_b_only_path_count": len(paths_b - shared),
            "bed_a_haplotype_count": len(hap_a), "bed_b_haplotype_count": len(hap_b),
            "transition_count": len(transition_paths),
            "plotted_left_node_count": len(left_nodes),
            "plotted_right_node_count": len(right_nodes),
            "plotted_transition_count": len(plot_edges),
        },
        "outputs": {name: str(path.resolve()) for name, path in outputs.items()},
    }
    outputs["summary.json"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bed_a", type=Path, help="first four-column PGR-TK BED")
    parser.add_argument("bed_b", type=Path, help="second four-column PGR-TK BED")
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--label-a", help="display label for the first BED")
    parser.add_argument("--label-b", help="display label for the second BED")
    parser.add_argument(
        "--max-plotted-haplotypes", type=int, default=20,
        help="top classes per BED shown individually; remaining classes become Other (default: 20; 0: all)",
    )
    args = parser.parse_args(argv)
    try:
        summary = compare_pgrtk_haplotypes(
            args.bed_a, args.bed_b, args.output_prefix,
            args.label_a, args.label_b, args.max_plotted_haplotypes,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
