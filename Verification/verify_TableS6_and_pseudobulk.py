#!/usr/bin/env python3
"""
Local verification of the two quantities that could not be checked
in-session because the source file is too large to transfer.

Input : GSE103322_WBP5_EMT_markers_cell_level.csv   (~684 KB)
        columns include: cell_id, WBP5, LAMC2, ITGA5, VIM, KRT14, KRT17,
        is_cancer, site, non_cancer_type / compartment, patient

Outputs
  (A) Compartment-wise WBP5 positivity  -> checks Table S6
      manuscript: cancer 68.89, fibroblast 51.97, endothelial 51.15,
                  dendritic 3.92, T cell 1.21
  (B) Patient-level pseudobulk Spearman -> supports Section 4.4
      (pseudoreplication defence; not yet in the manuscript)

Usage:  python verify_TableS6_and_pseudobulk.py [path/to/csv]
"""
import sys
import pandas as pd
from scipy.stats import spearmanr

CSV = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\Users\user\Desktop\wbp5\GSE103322_WBP5_EMT_markers_cell_level.csv"

MARKERS = ["LAMC2", "ITGA5", "KRT14", "VIM", "KRT17"]
EXPECTED_PCT = {"cancer": 68.89, "fibroblast": 51.97, "endothelial": 51.15,
                "dendritic": 3.92, "T cell": 1.21}
CELL_LEVEL_RHO = {"LAMC2": 0.169, "ITGA5": 0.169, "KRT14": 0.125,
                  "VIM": 0.120, "KRT17": 0.002}

d = pd.read_csv(CSV, encoding="utf-8-sig")
d.columns = [c.strip() for c in d.columns]
print(f"{len(d):,} rows | columns: {list(d.columns)}\n")
for c in ["WBP5"] + MARKERS:
    if c in d.columns:
        d[c] = pd.to_numeric(d[c], errors="coerce")

# ── (A) Table S6 ──────────────────────────────────────────────────────
key = next((k for k in ("compartment", "non_cancer_type", "cell_type")
            if k in d.columns), None)
if key is None:
    print("No compartment column found; skipping part A.")
else:
    s6 = (d.assign(pos=d["WBP5"] > 0)
            .groupby(key)
            .agg(n_cells=("WBP5", "size"),
                 mean_expr=("WBP5", "mean"),
                 pct_pos=("pos", lambda x: 100 * x.mean()))
            .sort_values("pct_pos", ascending=False)
            .round(2))
    print("(A) Table S6 - compartment-wise WBP5")
    print(s6.to_string())
    print("\n    manuscript: " +
          ", ".join(f"{k} {v}" for k, v in EXPECTED_PCT.items()))

# ── (B) Patient-level pseudobulk ──────────────────────────────────────
if "patient" not in d.columns:
    print("\nNo patient column found; skipping part B.")
    sys.exit(0)
m = d[d.is_cancer == 1].dropna(subset=["WBP5"] + MARKERS)
sizes = m.groupby("patient").size()
pb = m.groupby("patient")[["WBP5"] + MARKERS].mean()
pb = pb[sizes >= 20]                      # drop patients with <20 malignant cells
print(f"\n(B) Patient-level pseudobulk Spearman (n={len(pb)} patients, "
      f">=20 malignant cells each)")
print(f"    {'marker':<8}{'rho':>9}{'p':>10}{'cell-level':>13}{'same sign':>11}")
for mk in MARKERS:
    r, p = spearmanr(pb["WBP5"], pb[mk])
    same = "yes" if r * CELL_LEVEL_RHO[mk] > 0 else "NO"
    print(f"    {mk:<8}{r:>9.3f}{p:>10.4f}{CELL_LEVEL_RHO[mk]:>13.3f}{same:>11}")
pb.round(4).to_csv("WBP5_patient_pseudobulk.csv", encoding="utf-8-sig")
print("\n    saved: WBP5_patient_pseudobulk.csv")
