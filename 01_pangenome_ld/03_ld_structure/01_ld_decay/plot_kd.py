#!/usr/bin/env python3

###############################################################################
# KDE contour plot for SV–SNP LD relationships
#
# Usage:
#   python plot_kd.py \
#       --input sv.100k.ld.dis.maxr2.tsv \
#       --output sv.100k.ld.dis.maxr2.pdf
#
# Description:
#   This script generates kernel density contour plots and scatter plots
#   for SV–SNP linkage disequilibrium (LD) relationships.
#
# Input format (tab-delimited):
#   column1 : x-axis values (e.g. r²)
#   column2 : y-axis values (e.g. SV-SNP distance in kb)
#   column3 : group/class label
#
# Example input:
#   0.85    12.4    Assembly
#   0.72    35.1    NGS
#
# Output:
#   PDF figure containing:
#       - KDE contour density plots
#       - Scatter plot overlay
#
# Requirements:
#   - pandas
#   - numpy
#   - matplotlib
#   - seaborn
###############################################################################

import argparse
import pandas as pd
import seaborn as sns
import matplotlib

# Non-interactive backend for cluster/HPC environments
matplotlib.use("Agg")

import matplotlib.pyplot as plt

# ------------------------------ Arguments ---------------------------------- #

parser = argparse.ArgumentParser(
    description="Plot KDE contours for SV–SNP LD relationships"
)

parser.add_argument(
    "--input",
    required=True,
    help="Input TSV file"
)

parser.add_argument(
    "--output",
    required=True,
    help="Output PDF file"
)

args = parser.parse_args()

# ------------------------------- Load data --------------------------------- #

data = pd.read_csv(
    args.input,
    sep="\t",
    header=None,
    names=["r2", "distance_kb", "group"]
)

# ------------------------------ Plot settings ------------------------------ #

sns.set_style("whitegrid")

fig, ax = plt.subplots(figsize=(10, 8))

# ------------------------------- KDE plots --------------------------------- #

for group in data["group"].unique():

    group_data = data[data["group"] == group]

    sns.kdeplot(
        data=group_data,
        x="r2",
        y="distance_kb",
        levels=5,
        alpha=0.7,
        linewidths=1.5,
        ax=ax,
        label=group
    )

# ----------------------------- Scatter plots ------------------------------- #

sns.scatterplot(
    data=data,
    x="r2",
    y="distance_kb",
    hue="group",
    alpha=0.3,
    s=20,
    edgecolor=None,
    ax=ax
)

# ----------------------------- Figure labels ------------------------------- #

ax.set_title(
    "SV–SNP LD Density Distribution",
    fontsize=16
)

ax.set_xlabel(
    r"Maximum LD ($r^2$, MAF ≥ 0.01)",
    fontsize=12
)

ax.set_ylabel(
    "SV–SNP Distance (kb)",
    fontsize=12
)

ax.legend(
    title="Group",
    frameon=True
)

# -------------------------------- Save plot -------------------------------- #

plt.tight_layout()

plt.savefig(
    args.output,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("============================================================")
print("Plot generated successfully")
print(f"Input : {args.input}")
print(f"Output: {args.output}")
print("============================================================")
