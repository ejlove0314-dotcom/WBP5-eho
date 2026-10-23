#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. S4 — TCEAL/WBP paralog specificity (300 dpi), 2026-10-24

Rebuilt from the deposited output because no generator for this figure existed
in the repository; the previously embedded image was a hand-made ~194 dpi
raster. Every plotted value is read from

    revision_outputs/credibility/gse103322/CR3_family_and_controls.csv

(A) matched-null percentile of each gene's own correlation with the p-EMT
    program, for the twelve evaluable TCEAL/BEX/WBP genes (column
    `pct_vs_reference`), sorted descending.
(B) detection rate in malignant cells for the five evaluable TCEALs
    (column `detection` x 100).

Colours and ordering reproduce the figure as published in the revision.
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

CSV = ("/home/claude/ejlove0314-dotcom/wbp5-eho/revision_outputs/"
       "credibility/gse103322/CR3_family_and_controls.csv")
OUT = os.path.dirname(os.path.abspath(__file__))

RED, DBLUE, LBLUE, GREY = "#D62728", "#1F77B4", "#9ECAE1", "#D2D2D2"

d = pd.read_csv(CSV)
d = d[d["evaluable"] == True].copy()
d["label"] = d["symbol"].str.replace(r"\s*\(TCEAL9\)", "", regex=True)

famA = d[d["role"].isin(["target", "family_core", "family_ext", "outgroup"])].copy()
famA = famA.sort_values("pct_vs_reference", ascending=False)


def colA(r):
    if r["role"] == "target":
        return RED
    if r["label"] == "TCEAL3":
        return DBLUE
    return GREY if r["role"] == "outgroup" else LBLUE


famA["c"] = famA.apply(colA, axis=1)

famB = d[((d["role"].isin(["family_core"])) & (d["label"].str.startswith("TCEAL")))
         | (d["role"] == "target")].copy()
famB["pct"] = famB["detection"] * 100
famB = famB.sort_values("pct", ascending=False)
famB["c"] = [RED if r == "target" else LBLUE for r in famB["role"]]

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42})

fig, (axA, axB) = plt.subplots(1, 2, figsize=(10.2, 2.50),
                               gridspec_kw={"width_ratios": [1.12, 1],
                                            "wspace": 0.42})

# ---------------- panel A ----------------
y = range(len(famA))[::-1]
axA.barh(list(y), famA["pct_vs_reference"], color=famA["c"], height=0.72)
axA.set_yticks(list(y))
axA.set_yticklabels(famA["label"], fontsize=8.5, fontweight="bold")
for yy, v in zip(y, famA["pct_vs_reference"]):
    axA.text(v + 1.5, yy, f"{v:.1f}", va="center", fontsize=8, zorder=5,
             bbox=dict(fc="white", ec="none", pad=0.6))
axA.axvline(95, color="k", ls="--", lw=1.0, zorder=0)
axA.text(95, len(famA) - 0.35, "95$^{th}$", ha="center", va="bottom", fontsize=6.5)
axA.set_xlim(0, 112)
axA.set_xticks([0, 20, 40, 60, 80, 100])
axA.set_xlabel("Matched-null percentile of each gene's correlation "
               "with the p-EMT program", fontsize=8.5, fontweight="bold")
axA.set_title("(A) p-EMT co-expression across TCEAL/WBP paralogs",
              fontsize=10.5, fontweight="bold", loc="left", pad=8)

handles = [Patch(fc=RED, label="WBP5 (TCEAL9)"),
           Patch(fc=DBLUE, label="TCEAL3 (TCEAL/BEX family; highest)"),
           Patch(fc=LBLUE, label="TCEAL/BEX family (other)"),
           Patch(fc=GREY, label="WBP outgroup")]
axA.legend(handles=handles, ncol=2, fontsize=7.5, frameon=False,
           loc="upper left", bbox_to_anchor=(-0.02, -0.19),
           handlelength=1.2, handleheight=0.9, columnspacing=1.4,
           labelspacing=0.45)

# ---------------- panel B ----------------
yb = range(len(famB))[::-1]
axB.barh(list(yb), famB["pct"], color=famB["c"], height=0.62)
axB.set_yticks(list(yb))
axB.set_yticklabels(famB["label"], fontsize=8.5, fontweight="bold")
for yy, v in zip(yb, famB["pct"]):
    axB.text(v + 1.2, yy, f"{v:.1f}%", va="center", fontsize=8)
axB.set_xlim(0, 80)
axB.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80])
axB.set_xlabel("Detection rate in malignant cells (%)",
               fontsize=8.5, fontweight="bold")
axB.set_title("(B) WBP5 has the highest detection among evaluable TCEALs",
              fontsize=10.5, fontweight="bold", loc="left", pad=8)

for ax in (axA, axB):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", length=0)

fig.savefig(f"{OUT}/Fig_S4_300dpi.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_S4.pdf", bbox_inches="tight")

print("panel A (sorted):")
for _, r in famA.iterrows():
    print(f"  {r['label']:<8} {r['pct_vs_reference']:6.1f}  {r['role']}")
print("panel B (sorted):")
for _, r in famB.iterrows():
    print(f"  {r['label']:<8} {r['pct']:6.1f}%")
print("saved", f"{OUT}/Fig_S4_300dpi.png")
