"""Patient-aware model of WBP5 on the five-marker p-EMT composite (GSE103322).

Reproduces the five values quoted in the response to Reviewer #1, point 1:
the ordinary least-squares coefficient, the mixed-model coefficient with its
95% confidence interval and p value, and the random-intercept variance.

Specification
-------------
Cells      GSE103322 malignant cells with an assignable patient identifier
           (n = 2,200; the fifteen cells on the pooled "combo" plate are
           excluded, as in Methods Section 2.2).
Composite  Mean of LAMC2, ITGA5, VIM, KRT14 and KRT17, each standardized to
           zero mean and unit variance first.  The standardization matters:
           on their original log2(TPM/10+1) scale the two keratins are an
           order of magnitude larger than VIM, so an unstandardized mean is
           close to a KRT14/KRT17 average and gives a different coefficient.
Model      WBP5 ~ composite, patient as a random intercept, fitted by
           restricted maximum likelihood (statsmodels MixedLM).

Usage
-----
    python mixedmodel_five_marker_composite.py \
        reproducibility/GSE103322_WBP5_EMT_markers_cell_level.csv

Requires pandas and statsmodels.  Verified with pandas 2.x and
statsmodels 0.15.0.
"""
import sys

import pandas as pd
import statsmodels.formula.api as smf

MARKERS = ["LAMC2", "ITGA5", "VIM", "KRT14", "KRT17"]
POOLED_PLATE_SUFFIX = "_comb"          # see note below

EXPECTED = {
    "ols": 0.597,
    "mixed": 0.505,
    "ci_lo": 0.320,
    "ci_hi": 0.689,
    "p": 8.2e-8,
    "group_var": 0.295,
}


def load(path):
    d = pd.read_csv(path)
    d = d[d["compartment"] == "Cancer cells"].copy()
    # The pooled "combo" plate cells cannot be assigned to an individual
    # patient and are excluded by Methods Section 2.2.  They are flagged in
    # the corrected deposit as pooled_plate; in the original file they are
    # identified by patient id alone, so the count is checked instead.
    if "pooled_plate" in d.columns:
        d = d[~d["pooled_plate"].astype(str).str.upper().eq("TRUE")]
    else:
        pooled = pd.read_csv(
            path.replace("_EMT_markers_cell_level", "_cell_level_raw"),
            usecols=["cell_id", "pooled_plate"],
        )
        flag = dict(zip(pooled["cell_id"], pooled["pooled_plate"].astype(str).str.upper()))
        d = d[~d["cell_id"].map(flag).eq("TRUE")]
    return d


def main(path):
    d = load(path)
    print(f"malignant cells with an assignable patient identifier: {len(d):,}")
    print(f"patients: {d['patient'].nunique()}")

    z = (d[MARKERS] - d[MARKERS].mean()) / d[MARKERS].std()
    df = pd.DataFrame({
        "WBP5": d["WBP5"].values,
        "composite": z.mean(axis=1).values,
        "patient": d["patient"].astype(str).values,
    }).dropna()

    ols = smf.ols("WBP5 ~ composite", data=df).fit()
    mix = smf.mixedlm("WBP5 ~ composite", data=df, groups=df["patient"]).fit(reml=True)

    b = mix.params["composite"]
    lo, hi = mix.conf_int().loc["composite"]
    got = {
        "ols": ols.params["composite"],
        "mixed": b,
        "ci_lo": lo,
        "ci_hi": hi,
        "p": mix.pvalues["composite"],
        "group_var": mix.cov_re.iloc[0, 0],
    }

    print(f"\n{'quantity':12s} {'reported':>12s} {'recomputed':>14s}")
    for k, want in EXPECTED.items():
        v = got[k]
        fmt = f"{v:14.4g}" if k == "p" else f"{v:14.4f}"
        wfmt = f"{want:12.4g}" if k == "p" else f"{want:12.3f}"
        print(f"{k:12s} {wfmt} {fmt}")
    print(f"\nmixed / OLS = {got['mixed'] / got['ols']:.4f} "
          f"({100 * (1 - got['mixed'] / got['ols']):.0f}% smaller)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1
         else "reproducibility/GSE103322_WBP5_EMT_markers_cell_level.csv")
