"""Fig. 4 (second revision; Fig. 5 of the first revision) - WBP5 and the p-EMT signature in the two cohorts.
(A) GSE103322 cancer cells, n=2,215: cell-level Spearman rho, values of manuscript Table 2 (unchanged).
(B) GSE198315 tumor-region epithelial cells under the revised (Leiden, marker-annotated) annotation,
    n=5,359 cells from 10 patients: within-patient Spearman coefficient adjusted for sequencing depth,
    one point per patient, bar = median across patients.
    Source: revision_outputs/gse198315/D2_depth_adjusted_per_patient.csv (column rho_depth_partial).
(C) Gene-set enrichment of the 100-gene p-EMT signature among the transcriptome-wide correlates of WBP5:
    observed mean adjusted coefficient of the signature genes against the null (mean +/- 2 SD of 10,000
    matched random gene sets), under the expression-matched (primary) and the expression- and
    detection-rate-matched nulls, pooled profile, in GSE103322 and the three GSE198315 subsamples.
    Source: results_R2b/<run>/R2_2_null_and_weighting.csv written by R2_specificity_sensitivity_20261116b.py.
No value is recomputed here.
Usage: python make_fig4_enrichment_20261117.py D2_depth_adjusted_per_patient.csv results_R2b figures_R2"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.8, "xtick.major.width": 0.8,
    "ytick.major.width": 0.8, "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "axes.titlesize": 8.5})
MM = 1 / 25.4
SRC, R2B, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
A = [("LAMC2", .169), ("ITGA5", .169), ("KRT14", .125), ("VIM", .120),
     ("EGFR", .077), ("PTPRC", .015), ("KRT17", .002), ("EPCAM", -.081)]
ORDER = ["LAMC2", "ITGA5", "KRT17", "KRT14", "VIM"]                 # same marker order as the previous panel (B)
d = pd.read_csv(SRC)
d = d[d.feature.isin(ORDER)]
assert d.groupby("feature").patient.nunique().eq(10).all() and int(d[d.feature == "LAMC2"].n_cells.sum()) == 5359
C1, C0 = "#4DB6D6", "#444444"
fig = plt.figure(figsize=(174 * MM, 152 * MM))
gs = fig.add_gridspec(2, 2, width_ratios=[8, 5.6], wspace=0.36, hspace=0.78)
a, b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
c = fig.add_subplot(gs[1, 0])
def dress(ax, title, sub, L):
    ax.axhline(0, color="#666666", lw=0.8, ls=(0, (4, 3)), zorder=1)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_title(f"{title}\n", fontweight="bold", pad=2)
    ax.text(0.5, 1.0, sub, transform=ax.transAxes, ha="center", va="bottom", fontsize=7.5, style="italic", color="#555555")
    ax.text(-0.17, 1.17, f"({L})", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
x = range(len(A))
a.bar(x, [v for _, v in A], 0.62, color=C1, zorder=2)
for i, (g, v) in enumerate(A):
    a.text(i, v + (0.008 if v >= 0 else -0.008), f"{v:.3f}".replace("-", "−"), ha="center",
           va="bottom" if v >= 0 else "top", fontsize=7.5)
a.set_xticks(list(x)); a.set_xticklabels([g for g, _ in A], rotation=45, ha="right", rotation_mode="anchor", fontstyle="italic")
a.set_ylabel("Spearman's ρ with WBP5 (cell level)"); a.set_ylim(-0.11, 0.20)
dress(a, "GSE103322 cancer cells", "Discovery cohort, n = 2,215", "A")
rng = np.random.default_rng(42)
rows = []
YT = d.rho_depth_partial.max() + 0.035                              # one row of medians above all points
for i, g in enumerate(ORDER):
    v = d[d.feature == g].sort_values("patient").rho_depth_partial.values
    med = float(np.median(v))
    b.bar(i, med, 0.62, color=C1, zorder=2)
    b.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=9, color=C0, lw=0, alpha=0.85, zorder=3)
    b.text(i, YT, f"{med:.3f}".replace("-", "−"), ha="center", va="bottom", fontsize=7.5)
    rows.append({"marker": g, "median": med, "n_positive": int((v > 0).sum()), "n_patients": len(v), "min": v.min(), "max": v.max()})
b.set_xticks(range(len(ORDER))); b.set_xticklabels(ORDER, rotation=45, ha="right", rotation_mode="anchor", fontstyle="italic")
b.set_ylabel("Within-patient ρ with WBP5\n(depth-adjusted)")
lo, hi = d.rho_depth_partial.min(), d.rho_depth_partial.max()
b.set_ylim(min(-0.05, lo - 0.03), YT + 0.045)
b.text(len(ORDER) - 0.45, YT, "median", ha="left", va="bottom", fontsize=7.5, style="italic", color="#555555")
b.set_yticks([-0.2, -0.1, 0.0, 0.1, 0.2, 0.3])
dress(b, "GSE198315 tumor-region epithelium", "Exploratory cohort, 5,359 cells, 10 patients", "B")
# ---- panel C: gene-set enrichment, observed vs matched null ---------------------------------
RUNS = [("gse103322", "GSE103322"), ("seed42", "GSE198315\nseed 42"), ("seed2026", "GSE198315\nseed 2026"), ("seed7", "GSE198315\nseed 7")]
NULLS = [("expression-matched (primary)", "o", "expression-matched null"), ("expression+detection-matched", "s", "expression- and detection-matched null")]
crow = []
for i, (tag, lab) in enumerate(RUNS):
    t = pd.read_csv(os.path.join(R2B, tag, "R2_2_null_and_weighting.csv"))
    t = t[(t.profile == "pooled (primary)") & (t.signature == "PEMT100")]
    assert len(t) == 2 and t.mean_rho_set.nunique() == 1
    obs = float(t.mean_rho_set.iloc[0])
    for k, (nl, mk, nlab) in enumerate(NULLS):
        r = t[t.null == nl].iloc[0]
        xx = i + (-0.17 if k == 0 else 0.17)
        c.errorbar(xx, r.null_mean, yerr=2 * r.null_sd, fmt=mk, ms=4.2, mfc="white", mec=C0, ecolor=C0, elinewidth=0.9,
                   capsize=2.2, mew=0.9, zorder=2, label=nlab if i == 0 else None)
        crow.append({"run": tag, "null": nl, "observed_mean_rho": obs, "null_mean": r.null_mean, "null_sd": r.null_sd,
                     "z": r.z, "perm_p_upper": r.perm_p_upper, "n_genes": int(r.n_genes)})
    c.plot([i - 0.34, i + 0.34], [obs, obs], color=C1, lw=2.2, solid_capstyle="butt", zorder=3, label="signature (observed)" if i == 0 else None)
    zz = t.set_index("null").z
    c.text(i, obs + 0.004, f"z = {zz[NULLS[0][0]]:.1f} / {zz[NULLS[1][0]]:.1f}", ha="center", va="bottom", fontsize=7.5)
c.set_xticks(range(len(RUNS))); c.set_xticklabels([l for _, l in RUNS], fontsize=7.5)
c.set_xlim(-0.6, len(RUNS) - 0.4); c.set_ylim(0, 0.098); c.set_yticks([0, 0.02, 0.04, 0.06, 0.08])
c.set_ylabel("Mean adjusted ρ with WBP5\n(100 signature genes)")
for sp in ("top", "right"): c.spines[sp].set_visible(False)
c.set_title("p-EMT gene-set enrichment among the correlates of WBP5\n", fontweight="bold", pad=2)
c.text(0.5, 1.0, "100-gene signature; null shown as mean ± 2 SD of 10,000 matched gene sets", transform=c.transAxes, ha="center", va="bottom", fontsize=7.5, style="italic", color="#555555")
c.text(-0.17, 1.17, "(C)", transform=c.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
h, l = c.get_legend_handles_labels()
c.legend(h, l, loc="center left", bbox_to_anchor=(1.06, 0.62), frameon=False, fontsize=7.5, handlelength=1.6, borderaxespad=0, ncol=1)
c.text(1.06, 0.30, "z: expression-matched null /\nexpression- and detection-matched null", transform=c.transAxes, fontsize=7.5, va="top", ha="left", color="#555555", style="italic")
for ext in ("pdf", "eps"):
    fig.savefig(f"{OUT}/Figure_4.{ext}", bbox_inches="tight", pad_inches=0.03)
fig.savefig(f"{OUT}/Figure_4_preview.png", dpi=300, bbox_inches="tight", pad_inches=0.03, facecolor="white")
src = pd.DataFrame(rows); src.to_csv(f"{OUT}/Figure_4B_source.csv", index=False); print(src.to_string(index=False))

cs = pd.DataFrame(crow); cs.to_csv(f"{OUT}/Figure_4C_source.csv", index=False); print(cs.to_string(index=False))
