"""
Fig. S2 — Per-residue AlphaFold2 pLDDT confidence for WBP5 (UniProt Q9UHQ7)

Regenerated 2026-10-24 from make_figS2.py (2026-08-11). Changes, and only these:
  1. The bracket label and the summary-box line that read "Structured core
     (res 35-100)" now read "Highest-confidence segment (res 35-100)", so that
     "structured core" refers only to residues 29-104 throughout the paper.
  2. Raster preview written at 300 dpi instead of 200 dpi.
No data, geometry, colour or numeric value is changed.

Values are the per-residue mean B-factor (pLDDT) of AF-Q9UHQ7-F1-model_v6.
P_0 pocket residues are those reported by DoGSiteScorer (job 0f462a25),
i.e. 38, 41, 42, 44, 45, 48, 49, 51, 52, 54, 55, 58 - NOT 64-74.
"""
import collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans",
                                          "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "Liberation Sans"
matplotlib.rcParams["mathtext.it"] = "Liberation Sans:italic"
matplotlib.rcParams["mathtext.bf"] = "Liberation Sans:bold"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
from matplotlib.patches import Rectangle

PDB = "/home/claude/ejlove0314-dotcom/wbp5-eho/WBP5_AlphaFold.pdb"
P0 = [38, 41, 42, 44, 45, 48, 49, 51, 52, 54, 55, 58]
CORE = range(35, 101)

LABEL = "Highest-confidence segment (res 35–100)"

b = collections.defaultdict(list)
for line in open(PDB):
    if line.startswith("ATOM"):
        b[int(line[22:26])].append(float(line[60:66]))
res = np.array(sorted(b))
pl = np.array([np.mean(b[r]) for r in res])

fig, ax = plt.subplots(figsize=(7.4, 4.1))

bands = [(0, 50, "#FF7D45", "Very low"), (50, 70, "#FFDB13", "Low"),
         (70, 90, "#65CBF3", "Confident"), (90, 100, "#0053D6", "Very confident")]
for lo, hi, c, lab in bands:
    ax.add_patch(Rectangle((0, lo), 105, hi - lo, color=c, alpha=0.16, zorder=0))
    ax.text(1.015, (lo + hi) / 2 / 100, lab, transform=ax.transAxes,
            va="center", ha="left", fontsize=7.5, color="#333333")

ax.axvspan(min(P0), max(P0), color="#B8860B", alpha=0.13, zorder=1)
ax.plot(res, pl, color="#1F3B99", lw=1.4, zorder=3)
ax.plot(P0, [np.mean(b[r]) for r in P0], "o", ms=5.2, color="#D2691E",
        mec="white", mew=0.6, zorder=4, label="P$_0$ pocket residues (n = 12)")

ax.annotate("", xy=(35, 6), xytext=(100, 6),
            arrowprops=dict(arrowstyle="|-|,widthA=0.4,widthB=0.4",
                            color="#1F3B99", lw=1.0))
ax.text(67.5, 10.5, LABEL,
        ha="center", fontsize=7.5, color="#1F3B99")

ax.text(0.02, -0.235,
        "P$_0$ pocket (res 38–58; 12 lining residues): "
        "mean pLDDT = 79.4 (range 60.9–91.0)\n"
        + LABEL + ": mean pLDDT = 78.4 (range 58.3–96.0)",
        transform=ax.transAxes, va="top", ha="left", fontsize=7.0,
        linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.4", fc="#FFFDF0",
                  ec="#B8860B", lw=0.9))

ax.set_xlim(0, 105); ax.set_ylim(0, 100)
ax.set_xlabel("Residue number", fontsize=9.5, fontweight="bold")
ax.set_ylabel("pLDDT", fontsize=9.5, fontweight="bold")
ax.set_xticks([1, 20, 38, 48, 58, 70, 85, 100])
ax.set_yticks([0, 25, 50, 70, 90, 100])
ax.tick_params(labelsize=8)
ax.legend(loc="lower left", bbox_to_anchor=(0.015, 0.03), frameon=False,
          handletextpad=0.4, prop={"size": 7.5})
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.subplots_adjust(left=0.085, right=0.845, top=0.82, bottom=0.24)
fig.savefig("FigS2.pdf")
fig.savefig("FigS2.svg")
fig.savefig("FigS2_300dpi.png", dpi=300)
print("saved FigS2.pdf / .svg / FigS2_300dpi.png")

# --- verification against the printed values ---
import statistics as st
core = [np.mean(b[r]) for r in CORE if r in b]
p0 = [np.mean(b[r]) for r in P0]
print(f"core 35-100: n={len(core)} mean={st.mean(core):.2f} "
      f"range {min(core):.1f}-{max(core):.1f}")
print(f"P0 twelve  : n={len(p0)} mean={st.mean(p0):.2f} "
      f"range {min(p0):.1f}-{max(p0):.1f}  >70: {sum(1 for x in p0 if x>70)}")
