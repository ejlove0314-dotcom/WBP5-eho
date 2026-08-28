#!/usr/bin/env python3
"""
WBP5 manuscript numerical verification — recomputation script
=============================================================
Reproduces every derived statistic in the manuscript that is not read
directly from a deposited results file. Run with no arguments.

    python verify_recompute.py

Dependencies: numpy, scipy
Verified 2026-08-27 against manuscript CSBJ_WBP5_v10_20260825.
"""
import numpy as np
from scipy.stats import chi2_contingency, false_discovery_control

FAIL = []

def check(label, got, expected, tol):
    ok = abs(got - expected) <= tol
    if not ok:
        FAIL.append(label)
    print(f"  [{'OK ' if ok else 'FAIL'}] {label:<46} {got:>12.4g}  (manuscript {expected})")

# ── 1. §3.1.2 chi-square tests ────────────────────────────────────────
# Source: gse198315_composition.tsv (counts for WBP5_high / EGFR_high /
# double_high). The three mutually exclusive subsets are obtained by
# subtracting the double-high cells from each single-marker group.
print("\n1. Section 3.1.2 chi-square tests")

# columns: Epithelial-like, Stromal-like, Immune-like, Endothelial-like, Unclassified
compartment = np.array([
    [450 - 36, 1553 - 46, 54 - 0, 270 - 1, 18 - 1],   # WBP5-only   n=2261
    [569 - 36,  269 - 46, 26 - 0,  11 - 1,  3 - 1],   # EGFR-only   n=794
    [      36,        46,      0,        1,      1],  # double-high n=84
])
# columns: NT, TP, TC, mLN
region = np.array([
    [512 - 15, 717 - 16, 1056 - 48, 60 - 5],
    [150 - 15, 241 - 16,  450 - 48, 37 - 5],
    [      15,       16,        48,      5],
])
assert list(compartment.sum(1)) == [2261, 794, 84]
assert list(region.sum(1))      == [2261, 794, 84]

c, p, d, e = chi2_contingency(compartment)
check("compartment chi2", c, 692.5, 0.05)
check("compartment df", d, 8, 0)
check("compartment expected-count<5 cells", (e < 5).sum(), 2, 0)
print(f"        p = {p:.3g}  (manuscript 3.0e-144)")

c, p, d, e = chi2_contingency(region)
check("region chi2", c, 26.8, 0.05)
check("region df", d, 6, 0)
check("region expected-count<5 cells", (e < 5).sum(), 1, 0)
print(f"        p = {p:.3g}  (manuscript 1.6e-04)")

# ── 2. Table S4 Benjamini-Hochberg FDR ────────────────────────────────
# Source: TIMER3 Gene_Correlation exports (genecorr_table-*.csv).
# BH is applied across the five markers within each cohort stratum and
# adjustment condition.
print("\n2. Table S4 Benjamini-Hochberg q values")
MARKERS = ["ITGA5", "LAMC2", "VIM", "KRT17", "KRT14"]
P = {
    ("adjusted", "HNSC"):  [5.82679620806112e-12, 7.84253907982802e-11,
                            1.07449897186814e-06, 0.448308579445154,
                            0.318925127177507],
    ("adjusted", "HPV-"):  [2.21761329320845e-10, 2.13431683335049e-10,
                            1.95157396512918e-06, 0.677766621854402,
                            0.143194210650613],
    ("adjusted", "HPV+"):  [0.0582271754691809, 0.321767777166349,
                            0.128339501474643, 0.949937403819684,
                            0.470750988199799],
    ("raw", "HNSC"):       [6.34287704644275e-12, 9.83221101992766e-11,
                            2.33280357631711e-05, 0.374373381040025,
                            0.162741194635822],
    ("raw", "HPV-"):       [8.30802237336247e-11, 1.09839888642425e-10,
                            1.27052015893235e-05, 0.646831097903963,
                            0.0615879743807213],
    ("raw", "HPV+"):       [0.132331981618752, 0.324359813177789,
                            0.689707382507937, 0.563406953739623,
                            0.33996090934726],
}
for key in sorted(P):
    q = false_discovery_control(P[key], method="bh")
    print(f"  {key[0]:<9} {key[1]:<5} " +
          "  ".join(f"{m}={v:.2g}" for m, v in zip(MARKERS, q)))

# ── 3. Table 4 / Section 3.5 / 3.7 docking margins ────────────────────
# Source: docking_replicates_P0_summary_v2.csv (mean of five fixed-seed
# replicates per ligand). All margins in the text are differences of
# these means.
print("\n3. Docking score margins (kcal/mol)")
MEAN = {"Pazopanib": -6.86, "Candidate_4": -6.26, "Candidate_2": -5.84,
        "Candidate_3": -5.64, "Verteporfin": -5.54, "Dasatinib": -5.48,
        "Candidate_1": -5.44, "Simvastatin": -5.42, "Dobutamine": -5.08,
        "Candidate_5": -5.08}
check("pazopanib - verteporfin",        MEAN["Verteporfin"] - MEAN["Pazopanib"],   1.32, 0.005)
check("verteporfin - dobutamine span",  MEAN["Dobutamine"] - MEAN["Verteporfin"],  0.46, 0.005)
check("Candidate_4 - pazopanib",        MEAN["Candidate_4"] - MEAN["Pazopanib"],   0.60, 0.005)
check("Candidate_2 - Candidate_4",      MEAN["Candidate_2"] - MEAN["Candidate_4"], 0.42, 0.005)
check("Candidate_3 - Candidate_2",      MEAN["Candidate_3"] - MEAN["Candidate_2"], 0.20, 0.005)

# ── 4. Table S3 depth-correction deltas ───────────────────────────────
# Source: gse198315_recalc/headline.tsv
print("\n4. Table S3 depth-correction deltas")
HEADLINE = {  # marker: (CPM1e4_log1p, as_deposited)
    "LAMC2": (0.4375723567473858, 0.5915163063215932),
    "ITGA5": (0.3666399703582800, 0.5165394615542203),
    "VIM":   (-0.0572931082491548, 0.0529238860581020),
    "KRT17": (0.2975980199128242, 0.4166487984041913),
    "KRT14": (0.1934675998662732, 0.2846975876897106),
}
EXPECTED_DELTA = {"LAMC2": -0.154, "ITGA5": -0.150, "VIM": -0.110,
                  "KRT17": -0.119, "KRT14": -0.091}
for m, (corr, raw) in HEADLINE.items():
    check(f"{m} delta", corr - raw, EXPECTED_DELTA[m], 0.0006)
lo, hi = min(raw - corr for corr, raw in HEADLINE.values()), \
         max(raw - corr for corr, raw in HEADLINE.values())
print(f"        inflation range {lo:.3f}-{hi:.3f}  (manuscript 0.09-0.15)")

print("\n" + ("ALL CHECKS PASSED" if not FAIL else f"FAILED: {FAIL}"))
