#!/usr/bin/env python3
"""
GSE198315 — region- and compartment-resolved statistics for WBP5/EGFR-high cells.

Reproduces every number reported in Results section 3.2 and Figure 4 of the
WBP5 (TCEAL9) manuscript, starting from the per-cell table written by
`run_gse198315.py` (obs.csv, 80,000 quality-filtered cells).

Definitions
-----------
WBP5-high / EGFR-high : top decile among cells with non-zero expression of that
                        gene, computed within the analysed subset.
double-high           : cells meeting both criteria.
WBP5-only / EGFR-only : mutually exclusive groups (high for one gene, not both).

Regional trend is reported as two contrasts rather than a single ordinal odds
ratio, because prevalence is not monotonic across NT -> TP -> TC -> mLN:
  (i)  ordinal logistic regression restricted to the primary site (NT, TP, TC)
  (ii) primary site versus lymph node metastasis (Fisher's exact test), with a
       patient-level Wilcoxon signed-rank test as a paired confirmation.

Composition differences among the three mutually exclusive groups (WBP5-only,
EGFR-only, double-high) are tested by chi-square.

Usage
-----
    python gse198315_region_stats.py --obs out/obs.csv --outdir out

Requires: numpy, pandas, scipy, statsmodels
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
from scipy.stats import chi2_contingency, fisher_exact, wilcoxon
from statsmodels.stats.proportion import proportion_confint

REGIONS = ["NT", "TP", "TC", "mLN"]
REGION_SCORE = {"NT": 1, "TP": 2, "TC": 3, "mLN": 4}
COMPARTMENTS = ["Epithelial-like", "Stromal-like", "Immune-like",
                "Endothelial-like", "Unclassified"]
GROUPS = ["WBP5_high", "EGFR_high", "double_high"]


def classify(df: pd.DataFrame, gene_a: str = "WBP5", gene_b: str = "EGFR",
             decile: float = 0.90) -> tuple[pd.DataFrame, float, float]:
    """Assign high / only / double-high labels by top-decile-of-positive rule."""
    a = df[f"expr_{gene_a}"]
    b = df[f"expr_{gene_b}"]
    thr_a = a[a > 0].quantile(decile)
    thr_b = b[b > 0].quantile(decile)
    df = df.copy()
    df[f"{gene_a}_high"] = a >= thr_a
    df[f"{gene_b}_high"] = b >= thr_b
    df["double_high"] = df[f"{gene_a}_high"] & df[f"{gene_b}_high"]
    df[f"{gene_a}_only"] = df[f"{gene_a}_high"] & ~df["double_high"]
    df[f"{gene_b}_only"] = df[f"{gene_b}_high"] & ~df["double_high"]
    df["score"] = df["region"].map(REGION_SCORE).astype(float)
    df["site"] = np.where(df["region"].eq("mLN"), "mLN", "Primary")
    return df, float(thr_a), float(thr_b)


def composition(df: pd.DataFrame, flag: str) -> list[dict]:
    """Region and compartment composition of one high-expression group."""
    sub = df[df[flag]]
    rows = []
    for key, levels in (("region", REGIONS), ("compartment", COMPARTMENTS)):
        for lv in levels:
            n = int((sub[key] == lv).sum())
            rows.append({"group": flag, "variable": key, "level": lv,
                         "n": n,
                         "pct_of_group": round(n / len(sub) * 100, 2) if len(sub) else np.nan})
    return rows


def prevalence(df: pd.DataFrame, flag: str) -> list[dict]:
    """Within-region prevalence with Wilson 95% confidence intervals."""
    rows = []
    for r in REGIONS:
        sub = df[df["region"] == r]
        k, n = int(sub[flag].sum()), len(sub)
        lo, hi = proportion_confint(k, n, method="wilson")
        rows.append({"group": flag, "region": r, "n_positive": k, "n_total": n,
                     "prevalence_pct": round(k / n * 100, 3),
                     "wilson_lo_pct": round(lo * 100, 3),
                     "wilson_hi_pct": round(hi * 100, 3)})
    return rows


def ordinal_trend(df: pd.DataFrame, flag: str, regions: list[str]) -> dict | None:
    """Logistic regression with region entered as an ordinal predictor."""
    sub = df[df["region"].isin(regions)]
    y = sub[flag].astype(int)
    if y.sum() < 10 or y.nunique() < 2:
        return None
    fit = sm.Logit(y, sm.add_constant(sub[["score"]])).fit(disp=0)
    beta, se = fit.params["score"], fit.bse["score"]
    return {"group": flag, "regions": "+".join(regions),
            "OR_per_step": round(float(np.exp(beta)), 4),
            "ci_lo": round(float(np.exp(beta - 1.96 * se)), 4),
            "ci_hi": round(float(np.exp(beta + 1.96 * se)), 4),
            "p_value": float(fit.pvalues["score"]),
            "n_positive": int(y.sum()), "n_total": int(len(sub))}


def primary_vs_mln(df: pd.DataFrame, flag: str) -> dict:
    """Fisher's exact test contrasting the primary site with lymph node."""
    tab = pd.crosstab(df["site"], df[flag])
    a, b = int(tab.loc["Primary", True]), int(tab.loc["Primary", False])
    c, d = int(tab.loc["mLN", True]), int(tab.loc["mLN", False])
    odds, p = fisher_exact([[a, b], [c, d]])
    return {"group": flag,
            "primary_pct": round(a / (a + b) * 100, 3),
            "mln_pct": round(c / (c + d) * 100, 3),
            "odds_ratio": round(float(odds), 3), "p_value": float(p),
            "n_primary_positive": a, "n_mln_positive": c}


def patient_level(df: pd.DataFrame, flag: str) -> dict:
    """Paired patient-level comparison of primary versus lymph node prevalence."""
    per = df.groupby(["patient", "region"])[flag].mean().unstack()
    per = per.reindex(columns=REGIONS)
    primary = per[["NT", "TP", "TC"]].mean(axis=1)
    mln = per["mLN"]
    ok = primary.notna() & mln.notna()
    out = {"group": flag, "n_patients": int(ok.sum()),
           "median_primary_pct": round(float(primary[ok].median() * 100), 3),
           "median_mln_pct": round(float(mln[ok].median() * 100), 3),
           "n_patients_primary_gt_mln": int((primary[ok] > mln[ok]).sum())}
    try:
        out["wilcoxon_p"] = float(wilcoxon(primary[ok], mln[ok]).pvalue)
    except ValueError:
        out["wilcoxon_p"] = None
    # within-primary consistency: does the tumour core exceed non-tumour tissue?
    both = per[["NT", "TC"]].notna().all(axis=1)
    out["n_patients_TC_gt_NT"] = int((per.loc[both, "TC"] > per.loc[both, "NT"]).sum())
    out["n_patients_NT_and_TC"] = int(both.sum())
    return out


def three_group_chisq(df: pd.DataFrame) -> list[dict]:
    """Chi-square tests across the three mutually exclusive high-expression groups."""
    grp = np.where(df["double_high"], "double_high",
          np.where(df["WBP5_only"], "WBP5_only",
          np.where(df["EGFR_only"], "EGFR_only", None)))
    sub = df.assign(group=grp).dropna(subset=["group"])
    rows = []
    for var in ("region", "compartment"):
        table = pd.crosstab(sub["group"], sub[var])
        chi2, p, dof, expected = chi2_contingency(table)
        rows.append({"variable": var, "chi2": round(float(chi2), 3), "dof": int(dof),
                     "p_value": float(p),
                     "n_cells_expected_below_5": int((expected < 5).sum()),
                     "n_table_cells": int(expected.size)})
    return rows


def enrichment(df: pd.DataFrame, flag: str) -> list[dict]:
    """Compartment enrichment of a group relative to the whole subset."""
    sub = df[df[flag]]
    rows = []
    for c in COMPARTMENTS:
        base = (df["compartment"] == c).mean()
        grp = (sub["compartment"] == c).mean()
        rows.append({"group": flag, "compartment": c,
                     "pct_all_cells": round(base * 100, 2),
                     "pct_in_group": round(grp * 100, 2),
                     "fold_enrichment": round(grp / base, 3) if base else np.nan})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--obs", required=True, type=Path,
                    help="obs.csv written by run_gse198315.py")
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.obs)
    df, thr_w, thr_e = classify(df)

    header = {
        "n_cells": int(len(df)),
        "wbp5_threshold_log1p": round(thr_w, 6),
        "egfr_threshold_log1p": round(thr_e, 6),
        "n_WBP5_high": int(df["WBP5_high"].sum()),
        "n_EGFR_high": int(df["EGFR_high"].sum()),
        "n_double_high": int(df["double_high"].sum()),
        "n_WBP5_only": int(df["WBP5_only"].sum()),
        "n_EGFR_only": int(df["EGFR_only"].sum()),
        "compartment_counts": {c: int((df["compartment"] == c).sum())
                               for c in COMPARTMENTS},
        "region_counts": {r: int((df["region"] == r).sum()) for r in REGIONS},
        "versions": {"numpy": np.__version__, "pandas": pd.__version__,
                     "scipy": scipy.__version__,
                     "statsmodels": statsmodels.__version__},
    }

    comp = pd.DataFrame([r for g in GROUPS for r in composition(df, g)])
    prev = pd.DataFrame([r for g in GROUPS for r in prevalence(df, g)])
    enr = pd.DataFrame([r for g in GROUPS + ["WBP5_only", "EGFR_only"]
                        for r in enrichment(df, g)])
    trend = pd.DataFrame([t for g in GROUPS
                          for t in [ordinal_trend(df, g, ["NT", "TP", "TC"])]
                          if t is not None])
    fisher = pd.DataFrame([primary_vs_mln(df, g) for g in GROUPS])
    patient = pd.DataFrame([patient_level(df, g) for g in GROUPS])
    chisq = pd.DataFrame(three_group_chisq(df))

    for name, table in [("composition", comp), ("prevalence_wilson", prev),
                        ("compartment_enrichment", enr),
                        ("trend_within_primary", trend),
                        ("primary_vs_mln_fisher", fisher),
                        ("patient_level", patient),
                        ("three_group_chisq", chisq)]:
        table.to_csv(args.outdir / f"gse198315_{name}.tsv", sep="\t", index=False)

    (args.outdir / "gse198315_region_stats_summary.json").write_text(
        json.dumps(header, indent=2))

    print(json.dumps(header, indent=2))
    for name, table in [("Within-region prevalence (Wilson 95% CI)", prev),
                        ("Ordinal trend within primary site (NT, TP, TC)", trend),
                        ("Primary versus lymph node (Fisher)", fisher),
                        ("Patient-level paired comparison", patient),
                        ("Compartment enrichment", enr),
                        ("Three-group chi-square (WBP5-only / EGFR-only / double-high)", chisq)]:
        print(f"\n=== {name} ===")
        print(table.to_string(index=False))
    print(f"\n[done] tables written to {args.outdir}")


if __name__ == "__main__":
    main()
