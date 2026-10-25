#!/usr/bin/env python3
"""Regenerate Figure 4C (UMAP overlays) from the Run B 80,000-cell subset.

Top-left  : all WBP5-high cells   (n = 2,345, double-high included)
Top-right : all EGFR-high cells   (n = 878,  double-high included)
Bottom-left : double-high cells   (n = 84)
Bottom-right: combined overlay, mutually exclusive colouring so that the four
              categories partition the 80,000 cells.
"""
import numpy as np, pandas as pd, matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "axes.linewidth": 0.8, "savefig.dpi": 600, "savefig.bbox": "tight",
})

GREY, TEAL, RED, NAVY = "#D9D9D9", "#2AA198", "#E8503A", "#1F3A93"

obs = pd.read_csv("obs.csv")
um = pd.read_csv("umap_coords.csv")
assert (um["barcode"].values == obs["barcode"].values).all()
obs["UMAP1"], obs["UMAP2"] = um["UMAP1"].values, um["UMAP2"].values

a, b = obs["expr_WBP5"], obs["expr_EGFR"]
tw, te = a[a > 0].quantile(.90), b[b > 0].quantile(.90)
obs["W"] = a >= tw
obs["E"] = b >= te
obs["D"] = obs["W"] & obs["E"]
obs["Wo"] = obs["W"] & ~obs["D"]
obs["Eo"] = obs["E"] & ~obs["D"]
obs["Other"] = ~obs["W"] & ~obs["E"]

nW, nE, nD, nWo, nEo, nO = (int(obs[c].sum()) for c in
                            ["W", "E", "D", "Wo", "Eo", "Other"])


def base(ax):
    ax.scatter(obs["UMAP1"], obs["UMAP2"], s=1.2, c=GREY, linewidths=0,
               rasterized=True)


def style(ax, title, sub):
    ax.set_title(title, fontsize=13, fontweight="bold", pad=16)
    ax.text(0.5, 1.015, sub, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=8, color="#666666")
    ax.set_xlabel("UMAP1", fontsize=10)
    ax.set_ylabel("UMAP2", fontsize=10)
    ax.tick_params(labelsize=9)
    ax.spines[["top", "right"]].set_visible(False)


fig, axes = plt.subplots(2, 2, figsize=(11.4, 9.4))
fig.subplots_adjust(wspace=0.24, hspace=0.42, right=0.80)

# --- top-left : all WBP5-high -------------------------------------------
ax = axes[0, 0]; base(ax)
s = obs[obs["W"]]
ax.scatter(s["UMAP1"], s["UMAP2"], s=2.4, c=TEAL, linewidths=0, rasterized=True)
style(ax, "WBP5-high overlay",
      f"Threshold \u2265 {tw:.2f} (top 10% of WBP5-positive cells)  n = {nW:,}")

# --- top-right : all EGFR-high ------------------------------------------
ax = axes[0, 1]; base(ax)
s = obs[obs["E"]]
ax.scatter(s["UMAP1"], s["UMAP2"], s=2.4, c=RED, linewidths=0, rasterized=True)
style(ax, "EGFR-high overlay",
      f"Threshold \u2265 {te:.2f} (top 10% of EGFR-positive cells)  n = {nE:,}")

# --- bottom-left : double-high ------------------------------------------
ax = axes[1, 0]; base(ax)
s = obs[obs["D"]]
ax.scatter(s["UMAP1"], s["UMAP2"], s=6.0, c=NAVY, linewidths=0)
style(ax, "WBP5-high / EGFR-high overlay", f"Double-high  n = {nD:,}")

# --- bottom-right : combined (mutually exclusive) ------------------------
ax = axes[1, 1]; base(ax)
for flag, col, sz in [("Eo", RED, 2.4), ("Wo", TEAL, 2.4), ("D", NAVY, 6.0)]:
    s = obs[obs[flag]]
    ax.scatter(s["UMAP1"], s["UMAP2"], s=sz, c=col, linewidths=0,
               rasterized=True)
style(ax, "Combined overlay", "Mutually exclusive groups")

handles = [
    Line2D([], [], marker="o", ls="", ms=5, mfc=GREY, mec="none",
           label=f"Other cells (n={nO:,})"),
    Line2D([], [], marker="o", ls="", ms=5, mfc=RED, mec="none",
           label=f"EGFR-high (n={nEo:,})"),
    Line2D([], [], marker="o", ls="", ms=5, mfc=TEAL, mec="none",
           label=f"WBP5-high (n={nWo:,})"),
    Line2D([], [], marker="o", ls="", ms=5, mfc=NAVY, mec="none",
           label=f"WBP5-high / EGFR-high (n={nD:,})"),
]
fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.815, 0.30),
           frameon=True, edgecolor="#AAAAAA", fontsize=9)

fig.savefig("Fig4C_runB.png")
fig.savefig("Fig4C_runB.pdf")
print(f"[done] WBP5-high={nW}, EGFR-high={nE}, double={nD}; "
      f"combined {nO}+{nWo}+{nEo}+{nD}={nO+nWo+nEo+nD}")
