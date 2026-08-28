#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Patient-resolved WBP5-marker correlation analysis (GSE103322)
==============================================================
Addresses the pseudoreplication limitation in Section 4.4: the cell-level
Spearman coefficients in Table 2 treat 2,215 malignant cells from 18 patients
as independent observations.

Two complementary analyses:

  (A) WITHIN-PATIENT   Spearman rho computed separately in each patient, then
                       summarised (median, IQR, number of patients positive,
                       two-sided sign test). This tests the manuscript's actual
                       claim - co-expression *within* malignant cells - and
                       shows whether the pooled coefficient is driven by a few
                       patients or is consistent across them.

  (B) PSEUDOBULK       Patient means, then Spearman across patients. This is
                       the between-patient axis, comparable in kind to the
                       TCGA bulk analysis.

Usage
  python patient_level_correlation.py [path/to/cell_level.csv]

Input: a cell-level table with, at minimum, columns WBP5, LAMC2, ITGA5,
KRT14, VIM, KRT17 and a malignant-cell indicator. A patient column is used if
present; otherwise the script attempts to parse the patient ID from the cell
barcode (GSE103322 barcodes begin with the patient label, e.g. "HN25_P1_...").

Dependencies: pandas, numpy, scipy
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest, false_discovery_control

CSV = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\Users\user\Desktop\wbp5\GSE103322_WBP5_EMT_markers_cell_level.csv"

MARKERS = ["LAMC2", "ITGA5", "KRT14", "VIM", "KRT17"]
CELL_LEVEL = {"LAMC2": 0.169, "ITGA5": 0.169, "KRT14": 0.125,
              "VIM": 0.120, "KRT17": 0.002}          # Table 2, n = 2,215
MIN_CELLS = 20            # patients with fewer malignant cells are dropped

# ---------------------------------------------------------------- load
d = pd.read_csv(CSV, encoding="utf-8-sig")
d.columns = [c.strip() for c in d.columns]
print(f"{len(d):,} rows")
print(f"columns: {list(d.columns)}\n")

for c in ["WBP5"] + MARKERS:
    if c not in d.columns:
        sys.exit(f"ERROR: column {c!r} not found.")
    d[c] = pd.to_numeric(d[c], errors="coerce")

# malignant-cell indicator
flag = next((c for c in ("is_cancer", "malignant", "is_malignant")
             if c in d.columns), None)
if flag is None:
    sys.exit("ERROR: no malignant-cell indicator column found.")

# patient identifier
pid = next((c for c in ("patient", "patient_id", "sample", "donor")
            if c in d.columns), None)
if pid is None:
    cellcol = next((c for c in ("cell_id", "cell", "barcode", d.columns[0])
                    if c in d.columns), None)
    print(f"No patient column; parsing patient ID from {cellcol!r}.")
    d["patient"] = d[cellcol].astype(str).str.split("_").str[0]
    pid = "patient"
    print("  parsed IDs:", sorted(d[pid].unique())[:20], "...\n")

m = d[d[flag] == 1].dropna(subset=["WBP5"] + MARKERS).copy()
sizes = m.groupby(pid).size().sort_values(ascending=False)
keep = sizes[sizes >= MIN_CELLS].index
print(f"malignant cells: {len(m):,} from {m[pid].nunique()} patients")
print(f"retained: {len(keep)} patients with >= {MIN_CELLS} cells "
      f"({sizes[keep].sum():,} cells)")
print("cells per retained patient:",
      ", ".join(f"{p}:{n}" for p, n in sizes[keep].items()), "\n")
m = m[m[pid].isin(keep)]

# ============================================================ (A) within
print("=" * 74)
print("(A) WITHIN-PATIENT Spearman rho")
print("=" * 74)
print(f"{'marker':<8}{'median':>9}{'IQR':>18}{'pos/n':>9}"
      f"{'sign p':>9}{'pooled':>9}")
rows_a, sign_p = [], []
for mk in MARKERS:
    rhos = []
    for p, g in m.groupby(pid):
        if g["WBP5"].nunique() > 1 and g[mk].nunique() > 1:
            r, _ = spearmanr(g["WBP5"], g[mk])
            if np.isfinite(r):
                rhos.append(r)
    rhos = np.array(rhos)
    npos = int((rhos > 0).sum())
    bt = binomtest(npos, len(rhos), 0.5, alternative="two-sided")
    q1, q3 = np.percentile(rhos, [25, 75])
    sign_p.append(bt.pvalue)
    rows_a.append((mk, np.median(rhos), q1, q3, npos, len(rhos), bt.pvalue))
    print(f"{mk:<8}{np.median(rhos):>+9.3f}"
          f"{f'{q1:+.3f} to {q3:+.3f}':>18}{f'{npos}/{len(rhos)}':>9}"
          f"{bt.pvalue:>9.4f}{CELL_LEVEL[mk]:>+9.3f}")
q_a = false_discovery_control(sign_p, method="bh")
print("\nBH-adjusted sign-test q: " +
      ", ".join(f"{mk}={q:.3f}" for mk, q in zip(MARKERS, q_a)))

# ============================================================ (B) pseudobulk
print("\n" + "=" * 74)
print("(B) PSEUDOBULK Spearman rho (patient means)")
print("=" * 74)
pb = m.groupby(pid)[["WBP5"] + MARKERS].mean()
print(f"n = {len(pb)} patients")
print(f"{'marker':<8}{'rho':>9}{'p':>10}{'q':>10}{'pooled':>9}{'same sign':>11}")
ps = [spearmanr(pb["WBP5"], pb[mk])[1] for mk in MARKERS]
q_b = false_discovery_control(ps, method="bh")
rows_b = []
for mk, p, q in zip(MARKERS, ps, q_b):
    r, _ = spearmanr(pb["WBP5"], pb[mk])
    same = "yes" if r * CELL_LEVEL[mk] > 0 else "NO"
    rows_b.append((mk, r, p, q))
    print(f"{mk:<8}{r:>+9.3f}{p:>10.4f}{q:>10.4f}"
          f"{CELL_LEVEL[mk]:>+9.3f}{same:>11}")

# ---------------------------------------------------------------- save
pd.DataFrame(rows_a, columns=["marker", "median_rho", "q1", "q3",
                              "n_positive", "n_patients", "sign_test_p"]
             ).assign(sign_test_q=q_a).round(4).to_csv(
    "WBP5_within_patient_rho.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(rows_b, columns=["marker", "rho", "p", "q"]).round(4).to_csv(
    "WBP5_pseudobulk_rho.csv", index=False, encoding="utf-8-sig")
pb.round(4).to_csv("WBP5_patient_means.csv", encoding="utf-8-sig")
print("\nsaved: WBP5_within_patient_rho.csv, WBP5_pseudobulk_rho.csv, "
      "WBP5_patient_means.csv")
