#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. S5 — 300 dpi regeneration (2026-10-24)

Based on fig_S5_regenerate_20261021.py. Two deliberate differences:

  1. Panel B histogram is read from the DEPOSITED null distribution
     (revision_outputs/reference/YAP_matched_null_distribution.csv) instead of
     being re-drawn, so the figure matches the deposited output exactly.
     The empirical p printed in the panel title is the PRIMARY-analysis value
     from P5_YAP_TAZ_permutation.csv (= 0.0002), as reported in the manuscript;
     the legend already states that the distribution shown is an independent
     re-draw.
  2. Panel labels match the SI v24 figure: no panel-A title, panel-B axis
     titles as approved, panel-C category labels as approved.

Everything else - data, geometry, colours, numbers - is unchanged.
"""
import os, re, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

MATRIX = "/mnt/user-data/uploads/WBP5 revision/GSE103322_HNSCC_all_data.txt.gz"
REPO   = "/home/claude/ejlove0314-dotcom/wbp5-eho/revision_outputs"
NULLCSV = f"{REPO}/reference/YAP_matched_null_distribution.csv"
P5CSV   = f"{REPO}/gse103322_puram/P5_YAP_TAZ_permutation.csv"
OUT     = os.path.dirname(os.path.abspath(__file__))

GENE_ALIASES = ["WBP5", "TCEAL9"]
YAP22 = ["MYOF","AMOTL2","LATS2","CTGF","CYR61","ANKRD1","ASAP1","AXL","F3","IGFBP3",
         "CRIM1","FJX1","FOXF2","GADD45A","CCDC80","NT5E","DOCK5","PTPN14","ARHGEF17",
         "NUAK2","TGFB2","RBMS3"]
PURAM30 = ["SERPINE1","TGFBI","MMP10","LAMC2","P4HA2","PDPN","ITGA5","LAMA3","CDH13","TNC",
           "MMP2","EMP3","INHBA","LAMB3","VIM","SEMA3C","PRKCDBP","ANXA5","DHRS7","ITGB1",
           "ACTN1","CXCR7","ITGB6","IGFBP7","THBS1","PTHLH","TNFRSF6B","PDLIM7","CAV1","DKK3"]


def parse_patient(name):
    name = str(name).strip().strip("'\"")
    if re.match(r"^HNSCC[_]?combo", name, re.I):
        return None
    m = re.match(r"^HN(?:SCC)?[_]?(\d+)", name, re.I)
    return f"HN{int(m.group(1))}" if m else None


def load_matrix_malignant(path):
    df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False, compression="gzip")
    df.index = [str(i).strip().strip("'\"").strip() for i in df.index]
    df.columns = [str(c).strip().strip("'\"").strip() for c in df.columns]
    nk = lambda s: " ".join(str(s).replace('"', '').replace("'", "").lower().split())
    crow = next((i for i in df.index if nk(i) == "classified as cancer cell"), None)
    if crow is None:
        raise SystemExit("'classified as cancer cell' row not found")
    is_cancer = pd.to_numeric(df.loc[crow], errors="coerce").fillna(0).astype(int)
    mkeys = {"processed by maxima enzyme", "maxima enriched", "lymph node",
             "classified as cancer cell", "classified as non-cancer cells",
             "non-cancer cell type"}
    expr = df.drop(index=[i for i in df.index if nk(i) in mkeys]).apply(pd.to_numeric, errors="coerce")
    pats = pd.Series([parse_patient(c) for c in expr.columns], index=expr.columns)
    keep = (is_cancer.reindex(expr.columns) == 1) & (pats.notna())
    mal = expr.columns[keep.values]
    sys.stderr.write(f"malignant cells {len(mal):,} / genes {expr.shape[0]:,}\n")
    return expr[mal]


def zmodule(expr_sub, genes):
    g = [x for x in genes if x in expr_sub.index]
    sub = expr_sub.loc[g].astype(float)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0, skipna=True).values, g


def partial_spearman(x, y, z):
    from scipy.stats import rankdata
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    res = lambda t, u: (t - np.c_[np.ones_like(u), u] @
                        np.linalg.lstsq(np.c_[np.ones_like(u), u], t, rcond=None)[0])
    return np.corrcoef(res(rx, rz), res(ry, rz))[0, 1]


CACHE = f"{OUT}/_figS5_cache.npz"
if os.path.exists(CACHE):
    z = np.load(CACHE)
    w, yz, pz, yap1 = z["w"], z["yz"], z["pz"], z["yap1"]
    nyap, npur = int(z["nyap"]), int(z["npur"])
else:
    expr = load_matrix_malignant(MATRIX)
    gene = next((g for g in GENE_ALIASES if g in expr.index), None)
    w = expr.loc[gene].astype(float).values
    yz, yused = zmodule(expr, YAP22)
    pz, pused = zmodule(expr, PURAM30)
    yap1 = expr.loc["YAP1"].astype(float).values
    nyap, npur = len(yused), len(pused)
    np.savez(CACHE, w=w, yz=yz, pz=pz, yap1=yap1, nyap=nyap, npur=npur)

rho_obs = stats.spearmanr(w, yz)[0]
rho_pemt = stats.spearmanr(w, pz)[0]
rho_yap1 = stats.spearmanr(w, yap1)[0]
partial_rho = partial_spearman(w, yz, pz)

null = pd.read_csv(NULLCSV)["null_rho"].values
p95 = np.percentile(np.abs(null), 95)
p5 = pd.read_csv(P5CSV)
emp_p = float(p5.loc[p5["matched_on_expression_bins"] == True, "perm_p_two_sided"].iloc[0])

print("=== recomputed ===")
print(f"YAP genes {nyap}/22, PURAM30 {npur}/30, cells {len(w):,}")
print(f"rho_obs    = {rho_obs:.4f}   [printed 0.173]")
print(f"rho_YAP1   = {rho_yap1:.4f}   [printed 0.019]")
print(f"rho_pEMT   = {rho_pemt:.4f}   [printed 0.250]")
print(f"partial    = {partial_rho:.4f}   [printed 0.040]")
print(f"null 95%|r|= {p95:.4f}   emp_p (P5) = {emp_p:.4f}")

plt.rcParams.update({"font.size": 11, "svg.fonttype": "none"})
fig = plt.figure(figsize=(15, 4.4), dpi=300)
gs = fig.add_gridspec(1, 3, wspace=0.42)

axA = fig.add_subplot(gs[0])
hb = axA.hexbin(w, yz, gridsize=32, cmap="viridis", mincnt=1)
axA.set_xlabel("WBP5  log2(TPM/10+1)")
axA.set_ylabel("YAP/TAZ target-signature score (z)")
cb = fig.colorbar(hb, ax=axA, fraction=0.046, pad=0.03)
cb.set_label("Cell count", fontsize=9)

axB = fig.add_subplot(gs[1])
axB.hist(null, bins=40, color="#bdbdbd")
axB.axvline(p95, color="k", ls="--", lw=1.1)
axB.axvline(rho_obs, color="crimson", lw=2)
axB.set_xlabel("Spearman $\\rho$ (WBP5 vs random 22-gene score)")
axB.set_ylabel("Number of gene sets (permutations)")
axB.set_title(f"Matched null (10,000);  empirical p = {emp_p:.4f}", fontsize=10)

axC = fig.add_subplot(gs[2])
labs = ["YAP-target", "YAP1\nmRNA", "p-EMT", "YAP-target\n/p-EMT (partial)"]
vals = [rho_obs, rho_yap1, rho_pemt, partial_rho]
cols = ["#2c7fb8", "#cccccc", "#31a354", "#756bb1"]
bars = axC.bar(labs, vals, color=cols)
axC.axhline(0, color="k", lw=.6)
axC.set_ylabel("Spearman $\\rho$ with WBP5")
axC.tick_params(axis="x", labelsize=8)
for b, v in zip(bars, vals):
    axC.text(b.get_x() + b.get_width() / 2, v + 0.004, f"{v:.3f}", ha="center", fontsize=8)

for ax, lab in zip([axA, axB, axC], ["(A)", "(B)", "(C)"]):
    ax.text(-0.16, 1.06, lab, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="top", ha="left")

fig.savefig(f"{OUT}/Fig_S5_300dpi.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig_S5.pdf", bbox_inches="tight")
print(f"saved -> {OUT}/Fig_S5_300dpi.png")
