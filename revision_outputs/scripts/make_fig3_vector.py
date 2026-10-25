"""Fig. 3 - GSE198315 UMAP (80,000 cells): (A) proxy compartment, (B) region,
(C) WBP5, (D) EGFR. Port of fig_umap_gse198315.py (F:\\...\\figures_arial_unified\\_scripts).
Input: GSE198315/out/obs.csv + umap_coords.csv (run_gse198315.py output).
Expression = obs expr_<gene> = log1p of library-size-normalised values (target 10,000 per cell;
run_gse198315.py lines 210-211, 259-264). The old script's colour-bar label 'UMI count' was
wrong and is corrected. Point layers are rasterised at 600 dpi; text/axes stay vector."""
import sys, numpy as np, pandas as pd
from matplotlib.lines import Line2D
sys.path.insert(0, sys.argv[1]); import _style as st; plt = st.plt
D, OUT = sys.argv[2], sys.argv[3]
obs = pd.read_csv(f"{D}/obs.csv"); um = pd.read_csv(f"{D}/umap_coords.csv")
assert len(obs) == 80000 and (obs.barcode.values == um.barcode.values).all()
Y = um[["UMAP1", "UMAP2"]].values
COMP = {"Epithelial-like": "#E07B6B", "Immune-like": "#6FC7D1", "Stromal-like": "#36A48F",
        "Endothelial-like": "#3C4D9E", "Unclassified": "#C8B899"}
REG = {"NT": "#A3D9E8", "TP": "#3E998E", "TC": "#E86859", "mLN": "#8B5E46"}
fig, axs = plt.subplots(2, 2, figsize=(174*st.MM, 128*st.MM))
fig.subplots_adjust(hspace=0.42, wspace=0.85, left=0.07, right=0.88, top=0.93, bottom=0.08)
def fin(ax, t, L):
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_title(t, fontweight="bold")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    st.letter(ax, f"({L})", x=-0.24, y=1.04)
def cat(ax, v, pal):
    for c in pal:
        m = v == c
        ax.scatter(Y[m, 0], Y[m, 1], s=0.25, c=pal[c], lw=0, alpha=0.8, rasterized=True)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=pal[c], markersize=5, label=c)
                       for c in pal], loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False,
              handletextpad=0.2, labelspacing=0.3, borderaxespad=0)
a, b, c, d = axs.ravel()
cat(a, obs.compartment.values, COMP); fin(a, "Proxy compartment", "A")
cat(b, obs.region.values, REG); fin(b, "Region", "B")
for ax, g, L in ((c, "WBP5", "C"), (d, "EGFR", "D")):
    e = obs[f"expr_{g}"].values; o = np.argsort(e)
    sc = ax.scatter(Y[o, 0], Y[o, 1], s=0.25, c=e[o], cmap="viridis", lw=0, vmin=0,
                    vmax=max(0.1, np.percentile(e, 99.5)), rasterized=True)
    cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.03); cb.outline.set_linewidth(0.6)
    cb.set_label(f"{g}, log1p (counts per 10$^4$)"); fin(ax, f"{g} expression", L)
fig.savefig(f"{OUT}/Figure_3.pdf", dpi=600, bbox_inches="tight", pad_inches=0.03)
fig.savefig(f"{OUT}/Figure_3.eps", dpi=600, bbox_inches="tight", pad_inches=0.03)
fig.savefig(f"{OUT}/Figure_3_preview.png", dpi=300, bbox_inches="tight", pad_inches=0.03, facecolor="white")
print(obs.region.value_counts().to_dict(), obs.compartment.value_counts().to_dict())
