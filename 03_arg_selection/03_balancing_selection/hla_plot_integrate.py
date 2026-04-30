import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tskit
import json
import os
import seaborn as sns
import pickle


import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def simple_plot_genes_df(plot_df, ax,color='blue', fontsize=10,biotype_color=True,x_shift=0,line_color="grey",line_lw=0.5):
    for x in plot_df['x']:
        ax.axvline(x=x, color=line_color, linestyle='--', alpha=0.5,lw=line_lw)
    for _, row in plot_df.iterrows():
        # ax.text(row['x']+1.5e5, row['y'], row['gene_id'], 
        #         fontsize=fontsize, ha='center', va='bottom')
        color=row["color"]
        ax.annotate(
            row['gene_id'],
            (row['x'], row['y']),
            xytext=(x_shift, 4),                 
            textcoords='offset points',
            fontsize=fontsize,
            ha='left',
            va='bottom',
            color=color
        )
def line_plot_genes_df(plot_df, ax,color='black', fontsize=10,width=1,plot_lines=True,biotype_color=True,x_shift=0,line_color="grey"):
    if plot_lines:
        for x in plot_df['x']:
            ax.axvline(x=x, color=line_color, linestyle='--', alpha=0.5,lw=0.5)
    half_w = width / 2

    for _, row in plot_df.iterrows():
        ax.fill_between(
            [row['start'], row['end']],
            row['y'] - half_w,
            row['y'] + half_w,
            color=color,
            linewidth=0
        )
    for _, row in plot_df.iterrows():
        # ax.text(row['x']+1.5e5, row['y'], row['gene_id'], 
        #         fontsize=fontsize, ha='center', va='bottom')
        color=row["color"]
        ax.annotate(
            row['gene_id'],
            (row['x'], row['y']),
            xytext=(x_shift, 4),                 
            textcoords='offset points',
            fontsize=fontsize,
            ha='left',
            va='bottom',
            color=color
        )



all_names=["singer","tsinfer","relate"]
protein_and_selected_genes_list=[]
tmrca_df_list=[]
for i in range(3):
    name=all_names[i]
    protein_and_selected_genes_list.append(pd.read_csv(f"{SCRIPT_DIR}/plot_data/{name}_protein_and_selected_genes.csv"))
    tmrca_df_list.append(pd.read_csv(f"{SCRIPT_DIR}/plot_data/{name}_tmrca.csv"))

import matplotlib.colors as mcolors
def darken(color, factor=0.7):
    r, g, b = mcolors.to_rgb(color)
    return (r*factor, g*factor, b*factor)

palette = {
    "LR (Singer)": "skyblue", 
    "SR (Singer)":  "darkslategrey",
    "LR (tsinfer)": "skyblue", 
    "SR (tsinfer)":  "darkslategrey",
    "LR (relate)": "skyblue", 
    "SR (relate)":  "darkslategrey"
}

line_color="olive"


plot_names=["Singer","tsinfer","relate"]
plot_y_label=["SINGER","tsinfer+tsdate","Relate"]

fig, axes = plt.subplots(
    nrows=3, 
    ncols=1, 
    figsize=(12, 10), 
    sharex=True
)

for i, ax in enumerate(axes):

    # --- gene track ---
    simple_plot_genes_df(
        protein_and_selected_genes_list[i],  # assume list of 3 dfs
        ax,
        x_shift=1,
        line_color=line_color,
        line_lw=0.6
    )

    # --- TMRCA lines ---
    sns.lineplot(
        data=tmrca_df_list[i],   # assume list of 3 dfs
        x="Position (BP)",
        y="Pairwise TMRCA (Mya)",
        hue="Method",
        hue_order=[f"LR ({plot_names[i]})", f"SR ({plot_names[i]})"],
        palette=palette,
        alpha=0.9,
        lw=0.4,
        ax=ax
    )

    # threshold line
    ax.axhline(
        y=12,
        color=line_color,
        linestyle='--',
        lw=0.6,
        alpha=0.5
    )

    ax.set_ylim([0, 50])

    yticks = list(np.arange(0, 55, 5))
    if 12 not in yticks:
        yticks.append(12)
        yticks.remove(10)
    ax.set_yticks(sorted(yticks))

    ax.set_ylabel(f"Pairwise TMRCA (Mya)\n{plot_y_label[i]}")

# ----- Legend (only once, from first axis) -----
handles, labels = axes[0].get_legend_handles_labels()
label_map={}
for name in plot_names:
    label_map[f"LR ({name})"] ="Assembly"
    label_map[f"SR ({name})"] ="NGS" 
new_labels = [label_map.get(l, l) for l in labels]

leg = axes[0].legend(
    handles=handles[:2],
    labels=new_labels[:2],
    title=None
)

for line in leg.get_lines():
    line.set_linewidth(2.0)
    line.set_alpha(1.0)

# Remove legends from other panels
for ax in axes[1:]:
    ax.get_legend().remove()

plt.tight_layout(h_pad=4.0)

panel_labels = ["a", "b", "c"]

fig.canvas.draw()  # ensures positions are updated

for i, ax in enumerate(axes):
    pos = ax.get_position()
    
    fig.text(
        pos.x0 - 0.04,      # shift slightly left of axis
        pos.y1 + 0.02,      # slightly above axis
        panel_labels[i],
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="bottom"
    )

plt.savefig(
    f"{SCRIPT_DIR}/figures/hla_anno_integrate.pdf",
    bbox_inches="tight"
)

plt.show()