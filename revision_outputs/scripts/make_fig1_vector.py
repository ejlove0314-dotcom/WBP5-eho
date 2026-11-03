#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_fig1_vector.py — Figure 1 (Computational workflow), vector PDF
===================================================================
Draws Fig. 1 directly with matplotlib and writes a vector PDF.

Why this script exists
----------------------
The journal states that it cannot accept "figures prepared in PowerPoint or
Word formats that have been converted to other acceptable formats such as
PostScript or PDF".  The previous generator (`make_fig1.py`) built a .pptx
with python-pptx and converted it with `soffice --headless --convert-to pdf`,
so the resulting file carried Creator "Impress" / Producer "LibreOffice" and
fell under that prohibition.  This script reproduces the same figure without
going through a presentation format at any stage.

Nothing about the figure's content changes.  The geometry (a 13.333 x 7.62 in
page), the colours, the box positions, the font sizes and every string are
carried over unchanged from make_fig1.py, including the corrections of
2026-11-03:

  * "WBP5 is expressed across the malignant and mesenchymal compartments ...
    (51.97% of fibroblasts and 68.89% of cancer cells positive)"
  * "full-length model; core 29-104"
  * the citation [14]
  * "Drugs-lib; 4,574 approved drugs" and "5 approved drugs + 5 de novo
    candidates"
  * the six tool boxes as two rows of three

Typography: the .pptx asked for Calibri and LibreOffice substituted Carlito,
its metric-compatible equivalent, so Carlito is requested directly here and
the type sets identically.

The only visible difference from the Impress render is the loss of the soft
drop shadows LibreOffice added to the two arrows and the two italic notes.
Those were a rendering artefact of the conversion, not a design element.

Output: Figure_1_20261103_vector.pdf (+ .png preview at 300 dpi)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon

matplotlib.rcParams["pdf.fonttype"] = 42        # embed TrueType, keep text live
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Carlito", "Calibri", "Liberation Sans"]

PAGE_W, PAGE_H = 13.333, 7.62                   # inches, as in the .pptx
BASELINE_FIX = 0.15                             # see text(); fraction of size

BLUE_FILL, BLUE_LINE, BLUE_TEXT = "#EAF3FB", "#2E75B6", "#1F4E79"
PUR_HEAD, PUR_LINE, PUR_FILL, PUR_TEXT = "#7030A0", "#7030A0", "#E9E1F4", "#4E2583"
PUR_SUBLN = "#B9A5D9"
GRN_HEAD, GRN_LINE, GRN_FILL, GRN_TEXT = "#2E7D32", "#2E7D32", "#E4F0DC", "#2F5B22"
GRN_SUBLN = "#A9D08E"
WHITE = "#FFFFFF"
ARROW = "#A6A6A6"

fig = plt.figure(figsize=(PAGE_W, PAGE_H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, PAGE_W)
ax.set_ylim(0, PAGE_H)
ax.axis("off")


def _y(y_top, h):
    """PowerPoint measures y downward from the top; matplotlib upward."""
    return PAGE_H - y_top - h


def box(x, y, w, h, fill, line, lw=1.0, adj=None):
    """A .pptx roundRect: its corner radius is `adj` x the shorter side."""
    r = 0.0 if adj is None else adj * min(w, h)
    if r > 0:
        p = FancyBboxPatch((x + r, _y(y, h) + r), w - 2 * r, h - 2 * r,
                           boxstyle=f"round,pad={r},rounding_size={r}",
                           linewidth=lw if line else 0,
                           edgecolor=line if line else "none",
                           facecolor=fill if fill else "none",
                           mutation_aspect=1)
    else:
        p = Rectangle((x, _y(y, h)), w, h,
                      linewidth=lw if line else 0,
                      edgecolor=line if line else "none",
                      facecolor=fill if fill else "none")
    ax.add_patch(p)
    return p


def text(x, y, w, h, lines, size=10, color="#000000", bold=False, italic=False,
         space=0):
    """Centre a stack of lines in the box, as PowerPoint's MIDDLE anchor does.

    `lines` is a string, or a list of strings, or a list of
    (text, bold, italic, size) tuples.  `space` is the paragraph spacing in pt.
    """
    if isinstance(lines, str):
        lines = [lines]
    norm = []
    for ln in lines:
        norm.append((ln, bold, italic, size) if isinstance(ln, str) else ln)

    heights = [1.2 * sz / 72.0 for (_, _, _, sz) in norm]
    total = sum(heights) + (len(norm) - 1) * space / 72.0
    cx = x + w / 2.0
    top = _y(y, h) + h / 2.0 + total / 2.0          # top edge of the stack

    for (t, b, it, sz), lh in zip(norm, heights):
        # matplotlib centres a run on its glyph bounding box; PowerPoint
        # centres it on the line box (ascent to descent).  Measured against
        # the LibreOffice render, the two differ by a constant fraction of
        # the point size, so the line is dropped by that amount.
        ax.text(cx, top - lh / 2.0 - BASELINE_FIX * sz / 72.0, t,
                ha="center", va="center", fontsize=sz,
                fontweight="bold" if b else "normal",
                fontstyle="italic" if it else "normal",
                color=color)
        top -= lh + space / 72.0


def down_arrow(xc, y0, y1, w=0.22):
    """The .pptx downArrow at its default adjustments: half-width stem, half-
    height head."""
    h = y1 - y0
    yt, yb = _y(y0, 0), _y(y1, 0)                   # top and bottom, mpl coords
    ym = yt - h / 2.0
    sw = w / 2.0
    ax.add_patch(Rectangle((xc - sw / 2.0, ym), sw, h / 2.0,
                           facecolor=ARROW, edgecolor="none"))
    ax.add_patch(Polygon([(xc - w / 2.0, ym), (xc + w / 2.0, ym), (xc, yb)],
                         closed=True, facecolor=ARROW, edgecolor="none"))


L, W = 0.35, 12.63
cw, gap = 5.97, 0.25
cx1, cx2 = L + 0.25, L + 0.25 + cw + gap

# -- 1. the premise box across the top -----------------------------------
box(L, 0.22, W, 0.52, BLUE_FILL, BLUE_LINE, 1.25, adj=0.22)
text(L, 0.22, W, 0.52,
     "WBP5 (TCEAL9): prognostic candidate carried forward from our prior "
     "bulk-tissue study [14] — not an unbiased discovery screen",
     size=12, color=BLUE_TEXT, bold=True)

# -- 2. Phase I ----------------------------------------------------------
p1_y, p1_h = 0.95, 3.56
box(L, p1_y, W, p1_h, None, PUR_LINE, 1.5, adj=0.05)

box(L + 0.08, p1_y + 0.07, W - 0.16, 0.42, PUR_HEAD, PUR_HEAD, 1.0, adj=0.30)
text(L + 0.08, p1_y + 0.07, W - 0.16, 0.42,
     "Phase I — Transcriptomic target nomination  (main text)",
     size=12.5, color=WHITE, bold=True)

r1y, r1h = 1.58, 0.74
box(cx1, r1y, cw, r1h, PUR_FILL, PUR_SUBLN, 1.0, adj=0.14)
text(cx1, r1y, cw, r1h,
     [("Discovery cohort: GSE103322 (HNSCC; 5,902 cells)", True, False, 10.5),
      ("Exploratory cohort: GSE198315 (HPV-negative OSCC)", True, False, 10.5)],
     color=PUR_TEXT, space=2)

box(cx2, r1y, cw, r1h, PUR_FILL, PUR_SUBLN, 1.0, adj=0.14)
text(cx2, r1y, cw, r1h,
     [("Within-patient co-expression", True, False, 10.5),
      ("+ compartment-resolved analysis", True, False, 10.5)],
     color=PUR_TEXT, space=2)

r2y, r2h = 2.42, 0.82
box(cx1, r2y, cw, r2h, PUR_FILL, PUR_SUBLN, 1.0, adj=0.13)
text(cx1, r2y, cw, r2h,
     [("Program-level p-EMT enrichment", True, False, 10.5),
      ("depth-adjusted; expression-matched null gene sets", False, False, 9.5)],
     color=PUR_TEXT, space=2)

box(cx2, r2y, cw, r2h, PUR_FILL, PUR_SUBLN, 1.0, adj=0.13)
text(cx2, r2y, cw, r2h,
     [("Robustness checks", True, False, 10.5),
      ("leave-one-patient-out; gene-family specificity; Xq22 module "
       "adjustment; negative and positive controls", False, False, 9.5)],
     color=PUR_TEXT, space=2)

tw = 6.55
box(L + (W - tw) / 2, 3.33, tw, 0.46, PUR_FILL, PUR_SUBLN, 1.0, adj=0.22)
text(L + (W - tw) / 2, 3.33, tw, 0.46,
     "TCGA HNSC bulk cohort (n = 520; purity-adjusted) — complementary",
     size=10.5, color=PUR_TEXT, bold=True)

text(L + 0.25, 3.86, W - 0.50, 0.56,
     ["Finding: WBP5 is expressed across the malignant and mesenchymal "
      "compartments and is largely absent from lymphoid cells "
      "(51.97% of fibroblasts and 68.89% of cancer cells positive);",
      "its association with the p-EMT program is detectable at the level "
      "of the program rather than at the level of individual markers"],
     size=10, color=PUR_TEXT, italic=True, space=1)

# -- 3. Phases II-III ----------------------------------------------------
p2_y, p2_h = 4.70, 2.82
box(L, p2_y, W, p2_h, None, GRN_LINE, 1.5, adj=0.07)

box(L + 0.08, p2_y + 0.07, W - 0.16, 0.42, GRN_HEAD, GRN_HEAD, 1.0, adj=0.30)
text(L + 0.08, p2_y + 0.07, W - 0.16, 0.42,
     "Phases II–III — Exploratory structural and chemical ligandability  "
     "(Supplementary Information)",
     size=12.5, color=WHITE, bold=True)

boxes = [
    [("AlphaFold2 model", True, False, 10.5),
     ("AF-Q9UHQ7-F1-model_v6", False, False, 9.0),
     ("full-length model; core 29–104", False, False, 9.0)],
    [("DoGSiteScorer", True, False, 10.5),
     ("candidate pocket P₀", False, False, 9.0),
     ("residues 38–58; Drug Score 0.84", False, False, 9.0)],
    [("MTiOpenScreen", True, False, 10.5),
     ("library-scale calibration", False, False, 9.0),
     ("Drugs-lib; 4,574 approved drugs", False, False, 9.0)],
    [("REINVENT4", True, False, 10.5),
     ("reinforcement-learning", False, False, 9.0),
     ("de novo generation", False, False, 9.0)],
    [("AutoDock Vina", True, False, 10.5),
     ("docking of 5 approved drugs", False, False, 9.0),
     ("+ 5 de novo candidates", False, False, 9.0)],
    [("SwissADME / ProTox", True, False, 10.5),
     ("drug-likeness and", False, False, 9.0),
     ("predicted toxicity", False, False, 9.0)],
]
cols, g = 3, 0.15
bw = (W - 0.50 - g * (cols - 1)) / cols
bh, vg, by0 = 0.80, 0.12, 5.16
for i, lines in enumerate(boxes):
    r, c = divmod(i, cols)
    bx, byy = L + 0.25 + c * (bw + g), by0 + r * (bh + vg)
    box(bx, byy, bw, bh, GRN_FILL, GRN_SUBLN, 1.0, adj=0.10)
    text(bx, byy, bw, bh, lines, color=GRN_TEXT, space=1)

text(L + 0.25, 6.98, W - 0.50, 0.50,
     ["Computationally prioritized hypotheses: these results do not validate "
      "the transcriptomic target nomination,",
      "and experimental validation is required"],
     size=10, color=GRN_TEXT, italic=True, space=1)

# -- 4. the arrows between the phases ------------------------------------
for y0, y1 in ((0.74, 0.95), (4.51, 4.70)):
    down_arrow(L + W / 2, y0, y1)

fig.savefig("Figure_1_20261103_vector.pdf")
fig.savefig("Figure_1_20261103_vector.png", dpi=300)
print("saved Figure_1_20261103_vector.pdf / .png")
