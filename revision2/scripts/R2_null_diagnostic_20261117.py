#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R2_null_diagnostic_20261117.py — why the expression-matched and the expression+detection-matched nulls differ
================================================================================================================
Descriptive check written on 17 November 2026, AFTER the pre-specified sensitivity analyses
(R2_specificity_sensitivity_20261116b.py) had been run and read. It tests nothing and changes no result;
it only describes the deposited profiles.

Input : CR1_profile_WBP5.csv.gz written by CRED_profiles_20261021.py (columns gene, rho_adj, detection,
        mean_expr, bin, in_PURAM100), one per run.
Output: R2_null_diagnostic.csv

For each run
  1. detection rate of the signature genes relative to non-signature genes in the same expression bin
     (percentile within the bin; fraction of signature genes above the bin median);
  2. Spearman correlation between detection rate and the adjusted coefficient with WBP5 among non-signature
     genes, computed within each expression bin and averaged with weights equal to the number of signature
     genes in the bin;
  3. the mean adjusted coefficient of non-signature genes drawn (a) from the same expression bin and
     (b) from the same expression bin and detection quartile, i.e. the analytic expectation of the two nulls.

Run:
  python R2_null_diagnostic_20261117.py --profile gse103322=results_cred\\gse103322\\CR1_profile_WBP5.csv.gz ^
     --profile seed42=results_cred\\seed42\\CR1_profile_WBP5.csv.gz ... --out R2_null_diagnostic.csv
"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats


def one(tag, path):
    p = pd.read_csv(path)
    S, B = p[p.in_PURAM100], p[~p.in_PURAM100]
    by_bin = {b: g for b, g in B.groupby("bin")}
    pct, exp_e, exp_ed = [], [], []
    for _, g in S.iterrows():
        b = by_bin[g.bin]
        pct.append(float((b.detection < g.detection).mean()))
        exp_e.append(b.rho_adj.mean())
        edges = b.detection.quantile([0.25, 0.5, 0.75]).values
        q = np.searchsorted(edges, b.detection.values)
        qg = min(int(np.searchsorted(edges, g.detection)), 3)
        exp_ed.append(b.rho_adj.values[q == qg].mean())
    n_sig = S.bin.value_counts()
    r = pd.Series({b: stats.spearmanr(g.detection, g.rho_adj)[0] for b, g in by_bin.items() if len(g) > 30})
    r_sig = r.reindex(n_sig.index)
    pct = np.array(pct)
    return {
        "run": tag, "n_genes": len(p), "n_signature": len(S),
        "signature_median_detection": S.detection.median(),
        "signature_median_within_bin_detection_percentile": 100 * np.median(pct),
        "signature_fraction_above_bin_median_detection": float((pct > 0.5).mean()),
        "wilcoxon_p_percentile_vs_50": stats.wilcoxon(pct - 0.5).pvalue,
        "within_bin_spearman_detection_vs_rho_adj_weighted_by_signature_bins":
            float(np.average(r_sig.values, weights=n_sig.values)),
        "within_bin_spearman_detection_vs_rho_adj_median_all_bins": r.median(),
        "signature_mean_rho_adj": S.rho_adj.mean(),
        "expected_null_mean_expression_matched": float(np.mean(exp_e)),
        "expected_null_mean_expression_and_detection_matched": float(np.mean(exp_ed)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", action="append", required=True, help="TAG=path, repeatable")
    ap.add_argument("--out", default="R2_null_diagnostic.csv")
    a = ap.parse_args()
    rows = [one(*x.split("=", 1)) for x in a.profile]
    d = pd.DataFrame(rows)
    d.to_csv(a.out, index=False)
    pd.set_option("display.width", 250)
    print(d.T.to_string(float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
