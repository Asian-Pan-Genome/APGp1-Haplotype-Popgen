import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tskit
import os
import seaborn as sns

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def read_results_json(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)
    return pd.DataFrame(data["results"])

import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42

def simple_plot_genes_df(plot_df, ax,color='blue', fontsize=10):
    for x in plot_df['x']:
        ax.axvline(x=x, color='gray', linestyle='--', alpha=0.3,lw=0.5)
    for _, row in plot_df.iterrows():
        # ax.text(row['x']+1.5e5, row['y'], row['gene_id'], 
        #         fontsize=fontsize, ha='center', va='bottom')
        ax.annotate(
            row['gene_id'],
            (row['x'], row['y']),
            xytext=(3, 0),                 
            textcoords='offset points',
            fontsize=fontsize,
            ha='left',
            va='bottom'
        )
import sys
import gzip
CHROM = 6
chromosome_map = [
    "INVALID",
    "NC_060925.1", # chr1
    "NC_060926.1",
    "NC_060927.1",
    "NC_060928.1",
    "NC_060929.1",
    "NC_060930.1",
    "NC_060931.1",
    "NC_060932.1",
    "NC_060933.1",
    "NC_060934.1",
    "NC_060935.1",
    "NC_060936.1",
    "NC_060937.1",
    "NC_060938.1",
    "NC_060939.1",
    "NC_060940.1",
    "NC_060941.1",
    "NC_060942.1",
    "NC_060943.1",
    "NC_060944.1",
    "NC_060945.1",
    "NC_060946.1", # chr22
    "NC_060947.1",
    "NC_060948.1"
]
chrom_name_to_id = {name: idx for idx, name in enumerate(chromosome_map) if idx != 0}
biotype_colors = {
    "protein_coding": "green",
    "lncRNA": "#fc8803",
    "pseudogene": "blue",
    "transcribed_pseudogene": "blue",
    "miRNA": "purple",
    "other": "purple",
    "ncRNA": "purple",
    "snRNA": "purple",
    "snoRNA": "purple",
    "C_region": "purple",
    "misc_RNA": "purple",
    "C_region_pseudogene": "purple",
    "J_segment": "purple",
    "V_segment": "purple",
    "tRNA": "purple",
    "V_segment_pseudogene": "blue",
    "NA":"red"
}
# Ignored:
#NC_060947.1 chrX
#NC_060948.1 chrY

def gene_map_for_chr(map_file, chrom_num):
    map_result = []
    only_seq = chromosome_map[chrom_num]
    with gzip.open(map_file, "rt") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if line.startswith("#"):
                continue
            line = line.split("\t")
            assert len(line) == 9
            seqname, source, feature, start, end, score, strand, frame, attribute = line
            if feature == "gene" and seqname == only_seq:
                map_result.append([int(start), int(end), attribute])
    return map_result

gene_map_6 = sorted(gene_map_for_chr(os.path.join(SCRIPT_DIR, "GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"), CHROM))

base_path = f"{SCRIPT_DIR}/plot_data/"

def load_bundle(prefix):
    df = pd.read_csv(base_path + f"{prefix}_mask_manhattan_df.csv")
    gene_df = pd.read_csv(base_path + f"{prefix}_mask_gene_df.csv")
    
    with open(base_path + f"{prefix}_mask_plot_config.json") as f:
        config = json.load(f)
        
    return df, gene_df, config


def compute_midpoints(df):
    chromosomes = sorted(df["chr"].unique())
    midpoints = []
    for chr_id in chromosomes:
        df_chr = df[df["chr"] == chr_id]
        mid = df_chr["x"].min() + (df_chr["x"].max() - df_chr["x"].min()) / 2
        midpoints.append(mid)
    return chromosomes, midpoints


def plot_panel(ax, df, gene_df, config, panel_label, panel_name):

    colors = {
        0: str(config["colors"]["0"]),
        1: str(config["colors"]["1"])
    }

    ax.scatter(
        df["x"],
        df["value"],
        c=df["chr_parity"].map(colors),
        s=config["marker_size"],
        alpha=config["alpha"]
    )

    ax.axhline(
        y=config["percentile_9999"],
        color="red",
        linestyle="--",
        lw=0.5,
        label="0.01% empirical threshold (excluding HLA region)"
    )

    simple_plot_genes_df(gene_df, ax)

    yticks = list(np.arange(0, config["ymax"]+5, config["ytick_step"]))
    yticks.append(config["percentile_9999"])
    ax.set_yticks(sorted(yticks))

    ax.set_ylabel(f"Pairwise TMRCA (Mya)\n{panel_name}")

    ax.text(
        -0.04, 1.1,
        panel_label,
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="bottom"
    )



def main():
    # Load both datasets
    lr_df, lr_gene_df, lr_config = load_bundle("lr")
    sr_df, sr_gene_df, sr_config = load_bundle("sr")

    ys=10

    # 11 groups, total 55
    partitions = [8,8,8]+[8 for i in range(5)]+[9]
    start_ys = [10 for i in range(10)]+[10]
    step = 2
    ys = []

    for n, y0 in zip(partitions, start_ys):
        ys.extend([y0 + i * step for i in range(n)])

    # assert len(ys) == len(gene_df)
    lr_gene_df["y"] = ys[:len(lr_gene_df)]

    # lr_gene_df.loc[lr_gene_df.index[-1], "y"] = 30

    datasets = [
        (lr_df, lr_gene_df, lr_config),
        (sr_df, sr_gene_df, sr_config)
    ]

    panel_labels = ["a", "b"]
    panel_names = ["Assembly", "NGS"]

    # ---- Create stacked subplots ----
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    for i, ax in enumerate(axes):
        df, gene_df, config = datasets[i]
        plot_panel(ax, df, gene_df, config, panel_labels[i], panel_names[i])

    # ---- Shared X limits ----
    xmin_all = min(lr_df["x"].min(), sr_df["x"].min())
    xmax_all = max(lr_df["x"].max(), sr_df["x"].max())
    x_range = xmax_all - xmin_all
    xmax_extended = xmax_all + 0.11 * x_range

    for ax in axes:
        ax.set_xlim(-0.03 * x_range, xmax_extended + 10)

    # ---- Chromosome ticks (bottom panel only) ----
    chromosomes, midpoints = compute_midpoints(lr_df)
    axes[-1].set_xticks(midpoints, chromosomes)
    axes[-1].set_xlabel("Chromosome")

    # ---- Legend (top panel only) ----
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(
        handles,
        labels,
        loc="upper right",
        frameon=True,
        fontsize=11
    )

    plt.tight_layout(h_pad=3.0)
    # mpl.rcParams["ps.fonttype"] = 42
    plt.savefig(
        os.path.join(SCRIPT_DIR, "mask_tmrca_all_integrate.pdf"),
        bbox_inches="tight"
    )

    plt.show()

if __name__ == "__main__":
    main()