"""Fig. 5 - WBP5 co-expression with p-EMT/invasion and reference markers.
(A) GSE103322 cancer cells, n=2,215 -> values = manuscript Table 2 (Spearman rho vs WBP5).
(B) GSE198315 proxy tumour epithelial cells, n=53,192, depth-normalised
    (scale CPM1e4_log1p, cell_level) -> correlations_full.tsv (Drive, 2026-08-26).
Rebuilt because the docx Fig. 5 had corrupted y-axis tick labels (100/80/60...),
the typo 'GSE19831' and a mismatched bracket. Single bar colour, no class legend: the manuscript classes KRT17 differently in
Table 2 (pan-keratinocyte) and in the Fig. 5B/S1 legends (p-EMT/invasion marker), so the
figure does not encode a marker class."""
import sys
sys.path.insert(0, sys.argv[1]); import _style as st; plt = st.plt
from matplotlib.patches import Patch
OUT = sys.argv[2]
PEMT = {"LAMC2", "ITGA5", "KRT14", "VIM"}           # Table 2 class labels
A = [("LAMC2", .169), ("ITGA5", .169), ("KRT14", .125), ("VIM", .120),
     ("EGFR", .077), ("PTPRC", .015), ("KRT17", .002), ("EPCAM", -.081)]
B = [("LAMC2", .4375723567473858), ("ITGA5", .36663997035828005),
     ("KRT17", .2975980199128242), ("KRT14", .19346759986627327),
     ("VIM", -.05729310824915489)]
C1, C0 = "#4DB6D6", "#A9A9A9"
fig, (a, b) = plt.subplots(1, 2, figsize=(174*st.MM, 72*st.MM),
                           gridspec_kw=dict(width_ratios=[8, 5.6], wspace=0.34))
def panel(ax, data, title, sub, L, ylim):
    x = range(len(data))
    ax.bar(x, [v for _, v in data], 0.62, color=C1,
           zorder=2)
    ax.axhline(0, color="#666666", lw=0.8, ls=(0, (4, 3)), zorder=1)
    for i, (g, v) in enumerate(data):
        ax.text(i, v + (0.008 if v >= 0 else -0.008), f"{v:.3f}".replace("-", "−"),
                ha="center", va="bottom" if v >= 0 else "top", fontsize=6.5)
    ax.set_xticks(list(x)); ax.set_xticklabels([g for g, _ in data], rotation=45, ha="right",
                                               rotation_mode="anchor", fontstyle="italic")
    ax.set_ylabel("Spearman's ρ with WBP5"); ax.set_ylim(*ylim)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title(f"{title}\n", fontweight="bold", pad=2)
    ax.text(0.5, 1.0, sub, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=7, style="italic", color="#555555")
    st.letter(ax, f"({L})", x=-0.17, y=1.17)
panel(a, A, "GSE103322 cancer cells", "Discovery cohort, n = 2,215", "A", (-0.11, 0.20))
panel(b, B, "GSE198315 proxy tumor epithelium",
      "Exploratory cohort, n = 53,192 (depth-normalized)", "B", (-0.10, 0.49))
b.set_yticks([-0.1, 0.0, 0.1, 0.2, 0.3, 0.4])
fig.savefig(f"{OUT}/Figure_5.pdf", bbox_inches="tight", pad_inches=0.03)
