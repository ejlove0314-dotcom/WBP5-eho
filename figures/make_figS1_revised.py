"""
Fig. S1 - Multi-resolution WBP5-p-EMT correlation analysis
===========================================================
Panels A and B are unchanged (GSE198315, depth-corrected CPM1e4_log1p scale).
Panel C is regenerated: the TCGA HNSC bars now carry TIMER3 purity-adjusted
partial Spearman coefficients, replacing the earlier GEPIA3 zero-order values.

Sources
  A, B  : gse198315_recalc/correlations_full.tsv  (scale = CPM1e4_log1p)
  C     : GSE103322  -> WBP5_Phase1_recalc.ipynb
          GSE198315  -> correlations_full.tsv, proxy tumor epithelial cells
          TCGA HNSC  -> genecorr_table-*.csv (TIMER3, purity-adjusted)

Set FONT = "Arial" before final submission if Arial is installed.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

FONT = "Liberation Sans"          # metric-compatible Arial substitute
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = [FONT, "Arial", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["axes.linewidth"] = 0.9

MARKERS = ["LAMC2", "ITGA5", "KRT17", "KRT14", "VIM"]
CM = {"LAMC2": "#4FC3E8", "ITGA5": "#0E9E7E", "KRT17": "#8C8C8C",
      "KRT14": "#C9C9C9", "VIM": "#E8412E"}

# ---------------------------------------------------------------- panel A
STRATA = ["All\n(268,131)", "Proxy tumor\nepithelial\n(53,192)",
          "Non-proxy\ntumor\n(214,939)", "Primary proxy\nepithelial\n(48,108)",
          "LN proxy\nepithelial\n(5,084)"]
A = {  # correlations_full.tsv, CPM1e4_log1p, cell_level
    "LAMC2": [0.3000, 0.4376, 0.2016, 0.4238, 0.5382],
    "ITGA5": [0.4191, 0.3666, 0.4300, 0.3539, 0.4701],
    "KRT17": [0.1970, 0.2976, 0.1416, 0.2865, 0.3452],
    "KRT14": [0.1309, 0.1935, 0.1045, 0.1816, 0.2177],
    "VIM":   [0.2555, -0.0573, 0.3652, -0.0669, -0.1040],
}
# ---------------------------------------------------------------- panel B
B_RHO = [0.758, 0.891, -0.382, -0.608, 0.479]      # LN proxy, patient-mean
B_STAR = ["*", "**", "", "", ""]                    # p = .011 / .0005 / ns
# ---------------------------------------------------------------- panel C
COHORTS = ["GSE103322 (n=2,215)", "GSE198315 (n=53,192)",
           "TCGA HNSC bulk (n=520)"]
CC = ["#7B7FA8", "#4FC3E8", "#0E9E7E"]
C = {                       # marker: [GSE103322, GSE198315, TCGA purity-adj]
    "LAMC2": [0.169,  0.438,  0.288],
    "ITGA5": [0.169,  0.367,  0.304],
    "KRT17": [0.002,  0.298,  0.034],
    "KRT14": [0.125,  0.193, -0.045],
    "VIM":   [0.120, -0.057,  0.218],
}
C_NS_TCGA = {"KRT17", "KRT14"}      # q > 0.05 after purity adjustment

fig = plt.figure(figsize=(7.8, 9.4))
gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 1.0, 1.0],
                      hspace=0.62, left=0.115, right=0.975,
                      top=0.955, bottom=0.055)

def tidy(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.axhline(0, color="#999999", lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.tick_params(labelsize=8.5, length=3)

# ================================================================ A
axA = fig.add_subplot(gs[0])
w, x = 0.16, np.arange(len(STRATA))
for k, m in enumerate(MARKERS):
    axA.bar(x + (k - 2) * w, A[m], w, color=CM[m], zorder=2,
            edgecolor="none", label=m)
tidy(axA)
axA.set_xticks(x); axA.set_xticklabels(STRATA, fontsize=8)
axA.set_ylabel("Spearman's \u03c1", fontsize=9.5, fontweight="bold")
axA.set_ylim(-0.2, 0.62)
axA.set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6])
axA.legend(ncol=5, frameon=False, fontsize=8.5, loc="lower center",
           bbox_to_anchor=(0.5, 1.01), columnspacing=1.6, handlelength=1.1,
           handletextpad=0.5)
axA.text(-0.105, 1.14, "A", transform=axA.transAxes,
         fontsize=13, fontweight="bold", va="top")

# ================================================================ B
axB = fig.add_subplot(gs[1])
xb = np.arange(len(MARKERS))
cols = [CM[m] if s else "#A8A8A8" for m, s in zip(MARKERS, B_STAR)]
axB.bar(xb, B_RHO, 0.55, color=cols, zorder=2, edgecolor="none")
for i, (v, s) in enumerate(zip(B_RHO, B_STAR)):
    axB.text(i, v + (0.055 if v >= 0 else -0.055),
             f"\u03c1={v:+.3f}{s}".replace("+", ""),
             ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
tidy(axB)
axB.set_xticks(xb); axB.set_xticklabels(MARKERS, fontsize=9)
axB.set_ylabel("Patient-mean Spearman's \u03c1", fontsize=9.5, fontweight="bold")
axB.set_ylim(-0.95, 1.18)
axB.set_yticks([-0.5, 0.0, 0.5, 1.0])
axB.text(0.5, 1.055, "LN proxy tumor epithelial cells, n = 10 patients",
         transform=axB.transAxes, ha="center", fontsize=8.5,
         style="italic", color="#555555")
axB.text(-0.105, 1.14, "B", transform=axB.transAxes,
         fontsize=13, fontweight="bold", va="top")

# ================================================================ C
axC = fig.add_subplot(gs[2])
wc = 0.26
for j, coh in enumerate(COHORTS):
    vals = [C[m][j] for m in MARKERS]
    for i, (m, v) in enumerate(zip(MARKERS, vals)):
        ns = (j == 2 and m in C_NS_TCGA)
        axC.bar(i + (j - 1) * wc, v, wc, zorder=2,
                color="white" if ns else CC[j],
                edgecolor=CC[j] if ns else "none",
                hatch="////" if ns else None, linewidth=0.8)
        if ns:
            axC.text(i + (j - 1) * wc, v + (0.022 if v >= 0 else -0.022), "ns",
                     ha="center", va="bottom" if v >= 0 else "top",
                     fontsize=7, color="#555555")
# KRT17 in GSE103322 is 0.002 and invisible at this scale
axC.annotate("0.002", xy=(2 - wc, 0.004), xytext=(2 - wc, 0.085),
             ha="center", fontsize=7, color="#555555",
             arrowprops=dict(arrowstyle="-", lw=0.7, color="#999999"))
tidy(axC)
axC.set_xticks(np.arange(len(MARKERS)))
axC.set_xticklabels(MARKERS, fontsize=9)
axC.set_ylabel("Spearman's \u03c1", fontsize=9.5, fontweight="bold")
axC.set_ylim(-0.13, 0.52)
axC.set_yticks([-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
handles = [Patch(facecolor=c, edgecolor="none", label=l)
           for c, l in zip(CC, COHORTS)]
handles.append(Patch(facecolor="white", edgecolor=CC[2], hatch="////",
                     label="TCGA, q \u2265 0.05"))
axC.legend(handles=handles, frameon=False, fontsize=8,
           loc="upper right", bbox_to_anchor=(1.005, 1.03),
           handlelength=1.1, handletextpad=0.5, labelspacing=0.35)
axC.text(-0.105, 1.10, "C", transform=axC.transAxes,
         fontsize=13, fontweight="bold", va="top")

fig.savefig("/mnt/user-data/outputs/FigS1_revised.pdf")
fig.savefig("/mnt/user-data/outputs/FigS1_revised.png", dpi=600)
print("saved")
