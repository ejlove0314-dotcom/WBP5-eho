#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. S3 - assembled from TIMER3 Gene_Correlation exports
=========================================================
TIMER3 returns a two-panel plot per query: one panel for the gene-gene
association and one for the gene-purity association. The x-axis titles are
assigned by facet position rather than by content, so they are transposed
whenever the queried gene sorts alphabetically before "Purity" (ITGA5, LAMC2,
KRT14, KRT17). The facet strip titles are correct.

This script crops the correct half of each export, discards the incorrect
axis title along with the strip band, and re-labels each panel. No data point,
LOESS fit, confidence band, axis tick or statistical annotation is altered.

Panel   Source file                      Half   New x-axis title
  A     genecorr_plot_ITGA5_HNSC.jpg     right  Tumor purity
  B     genecorr_plot_ITGA5_HNSC.jpg     left   ITGA5 expression (log2 TPM)
  C     genecorr_plot_LAMC2_HNSC.jpg     left   LAMC2 expression (log2 TPM)
  D     genecorr_plot_VIM_HNSC.jpg       right  VIM expression (log2 TPM)

Crop boundaries were determined by pixel profiling of the 4400 x 2400 exports
(panel borders, strip band, tick-label band, axis-title band).

Usage:  python make_figS3_from_timer3.py [input_dir] [output_dir]
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image

IN = sys.argv[1] if len(sys.argv) > 1 else "timer3"
OUT = sys.argv[2] if len(sys.argv) > 2 else "."

# pixel boundaries of the 4400 x 2400 TIMER3 export
R0, R1 = 261, 2185          # below the strip band, down to the tick labels
LEFT = (208, 2248)          # y tick labels + left panel
RIGHT = (2304, 4375)        # y tick labels + right panel

PANELS = [
    ("A", "Tumor purity",                   "genecorr_plot_ITGA5_HNSC.jpg", "R"),
    ("B", "ITGA5 expression (log$_2$ TPM)", "genecorr_plot_ITGA5_HNSC.jpg", "L"),
    ("C", "LAMC2 expression (log$_2$ TPM)", "genecorr_plot_LAMC2_HNSC.jpg", "L"),
    ("D", "VIM expression (log$_2$ TPM)",   "genecorr_plot_VIM_HNSC.jpg",   "R"),
]

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42

crops = []
for lab, xt, fn, side in PANELS:
    im = Image.open(os.path.join(IN, fn))
    c0, c1 = LEFT if side == "L" else RIGHT
    tmp = os.path.join(OUT, f"_panel_{lab}.png")
    im.crop((c0, R0, c1, R1)).save(tmp)
    crops.append((lab, xt, tmp))

fig, axes = plt.subplots(2, 2, figsize=(7.48, 6.6))   # 190 mm wide
for ax, (lab, xt, fn) in zip(axes.ravel(), crops):
    ax.imshow(mpimg.imread(fn))
    ax.set_axis_off()
    ax.text(-0.045, 1.02, f"({lab})", transform=ax.transAxes,
            fontsize=13, fontweight="bold", va="bottom", ha="left")
    ax.text(0.55, -0.055, xt, transform=ax.transAxes,
            fontsize=10.5, va="top", ha="center")

fig.text(0.012, 0.5, "TCEAL9 expression (log$_2$ TPM)", rotation=90,
         va="center", ha="left", fontsize=11)
fig.subplots_adjust(left=0.055, right=0.995, top=0.965, bottom=0.055,
                    wspace=0.06, hspace=0.16)
fig.savefig(os.path.join(OUT, "Fig_S3_TIMER3_purity.pdf"))
fig.savefig(os.path.join(OUT, "Fig_S3_TIMER3_purity.png"), dpi=300)
for _, _, fn in crops:
    os.remove(fn)
print("saved: Fig_S3_TIMER3_purity.pdf / .png")
