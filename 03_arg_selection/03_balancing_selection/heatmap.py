import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tskit
import os
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext
from scipy.stats import spearmanr
from sklearn.metrics import mean_squared_error
from matplotlib.colors import LogNorm,LinearSegmentedColormap

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

mpl.rcParams["pdf.fonttype"] = 42


cmap2 = LinearSegmentedColormap.from_list(
    "white_to_skyblue",
    ["white", "#87CEFA", "darkslategrey"]  
)


step=1000
cuts = np.array(list(range(0, int(5e6) + 1, step)))

singer_lr_div=new_singer_lr.diversity(windows=cuts,mode="branch")
singer_sr_div=singer_sr.diversity(windows=cuts,mode="branch")

x = singer_sr_div
y = singer_lr_div


#HLA region, determined by the annotation
regions = [
    [29563926, 29909798],
    [31034915, 31313072],
    [32192811, 32589842],
    [32785920,33011047]
]

midpoints = 0.5 * (cuts[:-1] + cuts[1:])

mask = np.zeros_like(x, dtype=bool)

for start, end in regions:
    mask |= (midpoints+2.9e7 >= start) & (midpoints+2.9e7 <= end)

#mask
x_genes = x[mask]
y_genes = y[mask]
mask_neutral = ~mask
x_neutral = x[mask_neutral]
y_neutral = y[mask_neutral]





mpl.rcParams.update({
    "font.size": 14,        
    "axes.titlesize": 14,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.titlesize": 14
})


bins = 100
xbins = np.logspace(np.log10(1e4), np.log10(3e6), bins)
ybins = np.logspace(np.log10(1e4), np.log10(3e6), bins)


def plot_panel(ax, x, y, title=None):
    H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins])

    pcm = ax.pcolormesh(
        xedges,
        yedges,
        H.T + 1,
        cmap=cmap2,
        shading="auto",
        norm=LogNorm(vmin=1),
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlim(1e4, 3e6)
    ax.set_ylim(1e4, 3e6)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    line_min = max(xmin, ymin)
    line_max = min(xmax, ymax)

    ax.plot(
        [line_min, line_max],
        [line_min, line_max],
        linestyle="--",
        linewidth=1,
        color="black",
        alpha=0.7,
        zorder=3
    )

    ticks = [1e4, 1e5, 1e6]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.xaxis.set_major_formatter(LogFormatterMathtext())
    ax.yaxis.set_major_formatter(LogFormatterMathtext())

    rho, _ = spearmanr(x, y)
    mse_val = mean_squared_error(x, y)
    log10_bias = np.mean(np.log10(y) - np.log10(x))
    mean_ratio = np.mean(y) / np.mean(x)

    ax.text(
        0.05, 0.95,
                # f"Spearman ρ = {rho:.3f}\nMean log10 bias = {log10_bias:.3f}\nMean ratio = {mean_ratio:.3f}",
        f"Spearman ρ = {rho:.3f}\nRatio of means = {mean_ratio:.3f}",
        transform=ax.transAxes,
        verticalalignment='top',
        horizontalalignment='left',
        fontsize=12,
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
    )

    if title:
        ax.set_title(title)

    return pcm


def add_panel_label(fig, ax, label, dx=-0.02, dy=0.01, fontsize=14):
    bbox = ax.get_position()

    fig.text(
        bbox.x0 + dx,
        bbox.y1 + dy,
        label,
        fontsize=fontsize,
        fontweight="bold",
        ha="right",
        va="bottom"
    )



fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True,gridspec_kw={"wspace": 0.15})
axes[1].tick_params(axis='y', labelleft=True)

pcm1 = plot_panel(axes[0], x_genes, y_genes)
pcm2 = plot_panel(axes[1], x_neutral, y_neutral)


axes[0].set_ylabel("Assembly",labelpad=15)  # only left panel
axes[1].set_ylabel("")          

axes[0].set_xlabel("")         
fig.supxlabel("NGS", y=-0.02, fontsize=14)    


cbar_ax = fig.add_axes([0.95, 0.1, 0.02, 0.8])

cbar = fig.colorbar(pcm1, cax=cbar_ax)
cbar.set_label("Count")


add_panel_label(fig, axes[0], "a", dx=-0.02, dy=0.03)
add_panel_label(fig, axes[1], "b", dx=-0.02, dy=0.03)


plt.tight_layout(rect=[0, 0, 0.9, 1]) 
plt.subplots_adjust(wspace=0.12)

plt.savefig(
    f"{SCRIPT_DIR}/new_heat_map.pdf",
    bbox_inches="tight"
)

plt.show()