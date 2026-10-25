"""Fig. 2 - GSE103322 UMAP (A cell type, B cancer vs non-cancer, C site, D WBP5).
Re-plot of fig_umap_gse103322.py as a vector figure for CSBJ submission.
Inputs (Drive): GSE103322/out/umap_coords.csv, obs.csv (5,902 cells, row-aligned)
                GSE103322_WBP5_cell_level_raw.csv (WBP5, log2(TPM/10+1); row-aligned with obs)
Changes vs the docx raster: colour-bar label placed on the bar and stated as
log2(TPM/10+1) (the old script wrote 'log1p counts'); legend 'Lymph' -> 'Lymph node'.
Data, colours, orientation and panel content are unchanged."""
import sys, numpy as np, pandas as pd
from matplotlib.lines import Line2D
sys.path.insert(0, sys.argv[1]); import _style as st; plt = st.plt
D, OUT = sys.argv[2], sys.argv[3]
umap = pd.read_csv(f"{D}/umap103322.csv"); obs = pd.read_csv(f"{D}/obs103322.csv")
raw = pd.read_csv(f"{D}/wbp5_103322_raw.csv")
assert len(umap) == len(obs) == len(raw) == 5902
assert (raw.cell_id.values == obs.cell.values).all()
Y = umap[["UMAP1", "UMAP2"]].values.copy(); expr = raw.WBP5.values.astype(float)
cm = obs.malignancy.values == "Cancer cells"
if Y[cm, 0].mean() < Y[~cm, 0].mean(): Y[:, 0] = -Y[:, 0]      # cancer cells on the right
COMP = {"Cancer cells": "#E8836C", "Fibroblast": "#6FD0C4", "Endothelial": "#2EA27F",
        "T cell": "#2C3E8A", "B cell": "#F2B79E", "Macrophage": "#9A9BC5",
        "Dendritic": "#A8D8C9", "Mast": "#D62728", "Myocyte": "#7B5E3A",
        "Fibroblast (excl.)": "#B4A7CC", "Unannotated": "#D9D9D9"}
MAL = {"Non-cancer": "#7FD7E0", "Cancer cells": "#E87A67"}
SITE = {"Primary": "#9AA2C0", "Lymph node": "#F2A68D"}
fig, axs = plt.subplots(2, 2, figsize=(174*st.MM, 128*st.MM))
fig.subplots_adjust(hspace=0.42, wspace=0.85, left=0.07, right=0.88, top=0.93, bottom=0.08)
def fin(ax, title, L):
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_title(title, fontweight="bold")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    st.letter(ax, f"({L})", x=-0.24, y=1.04)
def scat(ax, vals, pal):
    for c in pal:
        m = vals == c
        if m.any(): ax.scatter(Y[m, 0], Y[m, 1], s=1.2, c=pal[c], lw=0, alpha=0.9)
def leg(ax, pal, labels=None):
    h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=pal[c], markersize=5,
                label=(labels or {}).get(c, c)) for c in pal]
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False,
              handletextpad=0.2, labelspacing=0.3, borderaxespad=0)
a, b, c, d = axs.ravel()
scat(a, obs.compartment.values, COMP); fin(a, "Cell type", "A"); leg(a, COMP)
scat(b, obs.malignancy.values, MAL); fin(b, "Cancer vs non-cancer", "B")
leg(b, MAL, {"Cancer cells": "Cancer"})
scat(c, obs.site.values, SITE); fin(c, "Primary vs lymph node", "C"); leg(c, SITE)
o = np.argsort(expr)
sc = d.scatter(Y[o, 0], Y[o, 1], s=1.2, c=expr[o], cmap="viridis", lw=0, vmin=0,
               vmax=np.percentile(expr, 99.5))
cb = fig.colorbar(sc, ax=d, fraction=0.05, pad=0.03)
cb.set_label(r"WBP5, log$_2$(TPM/10+1)"); cb.outline.set_linewidth(0.6)
fin(d, "WBP5 expression", "D")
for ax in axs.ravel(): ax.set_rasterization_zorder(None)
st.save(fig, "Figure_2", OUT)
print("cells", len(Y), "WBP5>0", int((expr > 0).sum()), "cancer det%",
      round((expr[cm] > 0).mean()*100, 2))
