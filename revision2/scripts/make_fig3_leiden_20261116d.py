"""Fig. 3 (second revision) - GSE198315 with the revised annotation.
(A) UMAP colored by the Leiden / canonical-marker cell-type annotation, (B) anatomical region,
(C) WBP5 expression, (D) share of cells with detected WBP5 per annotated cell type.
Replaces the proxy-compartment Fig. 3 of the first revision, which moves to the Supplement.

Input : gse198315_checkpoint_annotated.h5ad (step 2 output of GSE198315_reannotation_20261015.py;
        80,000 cells, seed 42). UMAP = obsm['X_umap'] computed there on the Harmony-corrected
        neighborhood graph; expression = .raw, log1p(counts per 10,000).
Output: Figure_3.pdf / .eps / _preview.png and Figure_3_source.csv (panel D values).
Usage : python make_fig3_leiden_20261116d.py gse198315_checkpoint_annotated.h5ad figures_R2
Point layers are rasterized at 600 dpi; text and axes stay vector."""
import os, sys
import numpy as np, pandas as pd, anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
matplotlib.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.8, "xtick.major.width": 0.8,
    "ytick.major.width": 0.8, "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "axes.titlesize": 8.5})
MM = 1 / 25.4
SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
a = ad.read_h5ad(SRC)
assert a.n_obs == 80000, a.n_obs
assert "X_umap" in a.obsm, "no UMAP in this file"
Y = np.asarray(a.obsm["X_umap"])
raw = a.raw.to_adata()
gene = next(g for g in ("WBP5", "TCEAL9") if g in raw.var_names)
e = np.asarray(raw[:, gene].X.todense()).ravel() if hasattr(raw.X, "todense") else np.asarray(raw[:, gene].X).ravel()
ct = a.obs["cell_type"].astype(str).replace({"Ambiguous": "Unassigned"}).values
reg = a.obs["region"].astype(str).values
pat = a.obs["patient"].astype(str).values
# Okabe-Ito palette (color-vision-deficiency safe); checked with a CVD simulator: all pairs of the region
# palette separate under protan/deutan/tritan vision; for the nine cell types the closest pair is
# Myeloid/B-Plasma under deutan vision, so cell types are also labeled directly on the embedding.
CT = {"Epithelial": "#D55E00", "Fibroblast": "#0072B2", "Endothelial": "#56B4E9", "Muscle": "#000000",
      "T/NK": "#E69F00", "B/Plasma": "#CC79A7", "Myeloid": "#009E73", "Mast": "#F0E442", "Unassigned": "#CCCCCC"}
REG = {"NT": "#56B4E9", "TP": "#E69F00", "TC": "#0072B2", "mLN": "#CC79A7"}
assert set(ct) <= set(CT), set(ct) - set(CT)
assert set(reg) <= set(REG), set(reg) - set(REG)

fig, axs = plt.subplots(2, 2, figsize=(174 * MM, 132 * MM))
fig.subplots_adjust(hspace=0.45, wspace=0.85, left=0.07, right=0.88, top=0.93, bottom=0.09)
def letter(ax, s, x=-0.24, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
def fin(ax, t, L):
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_title(t, fontweight="bold")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    letter(ax, f"({L})")
def cat(ax, v, pal):
    rng = np.random.default_rng(42); o = rng.permutation(len(v))      # no category drawn on top of the others
    ax.scatter(Y[o, 0], Y[o, 1], s=0.25, c=[pal[k] for k in v[o]], lw=0, alpha=0.8, rasterized=True)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=pal[c], markersize=5, label=c)
                       for c in pal if c in set(v)], loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
              handletextpad=0.2, labelspacing=0.3, borderaxespad=0)
A, B, C, D = axs.ravel()
cat(A, ct, CT); fin(A, "Cell type (revised annotation)", "A")
# Direct labels. Positions are set by hand for this embedding (obsm['X_umap'] of the seed-42 file is fixed),
# because automatic placement put labels on minor islands or on top of one another.
LAB = {"Myeloid": (0.5, 16.0, "center"), "Endothelial": (10.0, 13.0, "center"), "Mast": (-8.3, 12.6, "center"),
       "T/NK": (-1.0, 5.0, "center"), "Unassigned": (11.0, 7.3, "center"), "Fibroblast": (15.3, 3.6, "center"),
       "Muscle": (9.7, 3.1, "right"), "Epithelial": (10.0, -3.0, "center"), "B/Plasma": (0.0, -9.0, "center")}
for k, (lx, ly, ha) in LAB.items():
    if (ct == k).any():
        A.text(lx, ly, k, fontsize=7.5, ha=ha, va="center", path_effects=[pe.withStroke(linewidth=2, foreground="white")])
cat(B, reg, REG); fin(B, "Region", "B")
o = np.argsort(e)
sc = C.scatter(Y[o, 0], Y[o, 1], s=0.25, c=e[o], cmap="viridis", lw=0, vmin=0,
               vmax=max(0.1, np.percentile(e, 99.5)), rasterized=True)
cb = fig.colorbar(sc, ax=C, fraction=0.05, pad=0.03); cb.outline.set_linewidth(0.6)
cb.set_label("WBP5, log1p(counts per 10,000)"); fin(C, "WBP5 expression", "C")

rows = []
for k in CT:
    m = ct == k
    if m.any():
        rows.append({"cell_type": k, "n_cells": int(m.sum()), "n_patients": int(pd.Series(pat[m]).nunique()),
                     "pct_WBP5_detected": 100 * float((e[m] > 0).mean()), "mean_WBP5": float(e[m].mean())})
src = pd.DataFrame(rows).sort_values("pct_WBP5_detected", ascending=False)
src.to_csv(f"{OUT}/Figure_3_source.csv", index=False)
y = np.arange(len(src))[::-1]                                      # highest detection at the top
D.barh(y, src.pct_WBP5_detected, 0.68, color=[CT[k] for k in src.cell_type], zorder=2)
for yi, p in zip(y, src.pct_WBP5_detected):
    D.text(p + 1.5, yi, f"{p:.1f}", ha="left", va="center", fontsize=7.5)
D.set_yticks(y); D.set_yticklabels(src.cell_type)
D.set_xlabel("Cells with detected WBP5 (%)"); D.set_xlim(0, 100)
D.set_title("WBP5 detection by cell type", fontweight="bold")
for s in ("top", "right"): D.spines[s].set_visible(False)
letter(D, "(D)", x=-0.42)
for ext in ("pdf", "eps"):
    fig.savefig(f"{OUT}/Figure_3.{ext}", dpi=600, bbox_inches="tight", pad_inches=0.03)
fig.savefig(f"{OUT}/Figure_3_preview.png", dpi=300, bbox_inches="tight", pad_inches=0.03, facecolor="white")
print("gene:", gene); print(src.to_string(index=False))
print("regions:", pd.Series(reg).value_counts().to_dict())
