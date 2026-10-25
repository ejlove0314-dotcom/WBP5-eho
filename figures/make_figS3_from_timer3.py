#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. S3 - assembled from TIMER3 Gene_Correlation exports (revised 2026-10-25)
============================================================================
TIMER3 returns a two-panel plot per query: one panel for the gene-gene
association and one for the gene-purity association. The x-axis titles are
assigned by facet position rather than by content, so they are transposed
whenever the queried gene sorts alphabetically before "Purity" (ITGA5, LAMC2,
KRT14, KRT17). The facet strip titles are correct.

This script crops the correct half of each export, discards the incorrect
axis title along with the strip band, and re-labels each panel. No data point,
LOESS fit, confidence band, axis tick or statistical annotation is altered.
Crops are embedded at native resolution (no resampling).

2026-10-25 layout revision: the crop now starts above the top y-tick label so
the "9" is no longer clipped (strip remnant above the panel border is blanked);
each panel carries its own bold y-axis title,
x-axis titles are bold, and the panel letters sit in reserved space above
each panel so they no longer overlap the top y-tick label.

Panel   Source file                      Half   New x-axis title
  A     genecorr_plot_ITGA5_HNSC.jpg     right  Tumor purity
  B     genecorr_plot_ITGA5_HNSC.jpg     left   ITGA5 expression (log2 TPM)
  C     genecorr_plot_LAMC2_HNSC.jpg     left   LAMC2 expression (log2 TPM)
  D     genecorr_plot_VIM_HNSC.jpg       right  VIM expression (log2 TPM)

Usage:  python make_figS3_from_timer3.py [input_dir] [output_dir]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

IN = sys.argv[1] if len(sys.argv) > 1 else "timer3"
OUT = sys.argv[2] if len(sys.argv) > 2 else "."

# pixel boundaries of the 4400 x 2400 TIMER3 export
R0, R1 = 238, 2185          # from just above the top y-tick label ("9" starts at row 244)
STRIP_END = 250              # rows R0..249 right of the y-axis are facet-strip remnant -> blanked
AXIS = {"L": 330, "R": 2425} # y-axis line column of each half (pixel profiling)
LEFT = (208, 2248)          # y tick labels + left panel
RIGHT = (2304, 4375)        # y tick labels + right panel

YTITLE = "TCEAL9 expression level (log$_2$ TPM)"
PANELS = [
    ("A", "Tumor purity",                   "genecorr_plot_ITGA5_HNSC.jpg", "R"),
    ("B", "ITGA5 expression (log$_2$ TPM)", "genecorr_plot_ITGA5_HNSC.jpg", "L"),
    ("C", "LAMC2 expression (log$_2$ TPM)", "genecorr_plot_LAMC2_HNSC.jpg", "L"),
    ("D", "VIM expression (log$_2$ TPM)",   "genecorr_plot_VIM_HNSC.jpg",   "R"),
]

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Liberation Sans", "Arial", "DejaVu Sans"]
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "Liberation Sans"
matplotlib.rcParams["mathtext.bf"] = "Liberation Sans:bold"
matplotlib.rcParams["pdf.fonttype"] = 42

# ---- layout in inches (explicit, so labels never collide) ----------------
FIG_W = 7.48                                   # 190 mm
COL_L = [0.45, 4.22]                           # left edge of each image column
IMG_W = 3.20
ROW_TOP = [0.42, 4.12]                         # gap above each image row
# image height follows the crop aspect ratio
crop_h = R1 - R0
img_h = {s: IMG_W * crop_h / (c1 - c0) for s, (c0, c1) in (("L", LEFT), ("R", RIGHT))}
H_IMG = max(img_h.values())
XT_H = 0.30                                    # room for the x-axis title
FIG_H = ROW_TOP[1] + H_IMG + XT_H + 0.05

fig = plt.figure(figsize=(FIG_W, FIG_H))
checks = []
for k, (lab, xt, fn, side) in enumerate(PANELS):
    r, c = divmod(k, 2)
    c0, c1 = LEFT if side == "L" else RIGHT
    arr = np.array(Image.open(os.path.join(IN, fn)).convert("RGB").crop((c0, R0, c1, R1)))
    # blank the strip remnant above the panel border; the tick labels left of the axis are untouched
    arr[: STRIP_END - R0, AXIS[side] - c0 - 8:] = 255
    h = img_h[side]
    x0 = COL_L[c]
    y_top = ROW_TOP[r]
    ax = fig.add_axes([x0 / FIG_W, 1 - (y_top + h) / FIG_H, IMG_W / FIG_W, h / FIG_H])
    ax.imshow(arr, interpolation="none", aspect="auto")
    ax.set_axis_off()
    # panel letter: above the image, left of the y-axis title -> no overlap
    lt = fig.text((x0 - 0.40) / FIG_W, 1 - (y_top - 0.12) / FIG_H, f"({lab})",
                  fontsize=12, fontweight="bold", ha="left", va="bottom")
    # x-axis title under the image
    fig.text((x0 + IMG_W * 0.53) / FIG_W, 1 - (y_top + h + 0.04) / FIG_H, xt,
             fontsize=10, fontweight="bold", ha="center", va="top")
    # per-panel y-axis title
    yt = fig.text((x0 - 0.33) / FIG_W, 1 - (y_top + h / 2) / FIG_H, YTITLE,
                  rotation=90, fontsize=9, fontweight="bold", ha="left", va="center")
    checks.append((lab, lt, yt))

# guard: the panel letter must sit clearly above the y-axis title
fig.canvas.draw()
rend = fig.canvas.get_renderer()
for lab, lt, yt in checks:
    gap = (lt.get_window_extent(rend).y0 - yt.get_window_extent(rend).y1) / fig.dpi
    print(f"({lab}) vertical gap letter-to-y-title: {gap:.2f} in")
    assert gap > 0.08, f"panel {lab}: letter overlaps y-axis title"
fig.savefig(os.path.join(OUT, "Figure_S3.pdf"))
print("saved: Figure_S3.pdf  (%.2f x %.2f in)" % (FIG_W, FIG_H))
