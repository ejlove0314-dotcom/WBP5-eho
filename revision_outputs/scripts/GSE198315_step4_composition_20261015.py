#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE198315 step4 — 전이 림프절 해석에서 '구성 교란'을 분리
==========================================================
생성: 2026-10-15
대응: Reviewer #5 Major comment 4 (mLN 고갈 재해석) + Reviewer #1 point 3 (구성 혼성)

[이 분석이 필요한 이유]
환자별로 전체 세포를 집계하면 10명 전원에서 mLN의 WBP5-high 비율이 최저로 나옵니다.
그러나 mLN 샘플은 림프구가 압도적이고(T/NK·B/Plasma) 림프구는 WBP5가 거의 0이므로,
이 차이의 상당 부분은 '세포 구성'이지 '세포 상태'가 아닙니다.
이는 Reviewer #1 point 3 이 지적한 것과 정확히 같은 구조의 오류입니다.

따라서 두 가지를 나란히 계산해 비교합니다.
  (A) 집계 비교  : 환자×부위의 전체 세포 기준  -> 교란된 값
  (B) 세포형 내부: 환자×세포형×부위 기준       -> 올바른 값
두 값의 차이 자체를 보고하면, 지적을 이해하고 전면 적용했음을 보일 수 있습니다.

실행:
    python GSE198315_step4_composition_20261015.py \
        --h5ad gse198315_checkpoint_full_annotated_cnv.h5ad \
        --outdir results_gse198315
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

# 차이가 전부 0인 세포형(예: 림프구)에서 wilcoxon 이 내는 무해한 경고를 가립니다.
warnings.filterwarnings("ignore", category=RuntimeWarning)

GENE_ALIASES = ["WBP5", "TCEAL9"]
PRIMARY = ["TC", "TP"]      # 원발 부위(종양 중심/변연)
MET = "mLN"                 # 전이 림프절
MIN_CELLS = 20              # 환자×세포형×부위 최소 세포 수
N_BOOT = 10000
SEED = 42


def log(m):
    sys.stderr.write(f"[{pd.Timestamp.now():%H:%M:%S}] {m}\n"); sys.stderr.flush()


def paired_test(df, label, rng):
    """
    df: columns = [patient, met, primary]  (각 값은 WBP5-high 비율 %)
    환자 단위 paired 검정 + 부트스트랩 95% CI.
    """
    d = (df["met"] - df["primary"]).to_numpy(float)
    n = len(d)
    if n < 3:
        return {"contrast": label, "n_patients": n, "note": "환자 3명 미만 — 검정 생략"}
    try:
        w_stat, w_p = stats.wilcoxon(d, alternative="two-sided")
    except ValueError:
        w_stat, w_p = np.nan, np.nan
    boot = np.array([np.median(rng.choice(d, size=n, replace=True)) for _ in range(N_BOOT)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return {
        "contrast": label,
        "n_patients": n,
        "median_pct_met": float(np.median(df["met"])),
        "median_pct_primary": float(np.median(df["primary"])),
        "median_diff_met_minus_primary": float(np.median(d)),
        "ci95_lo": float(lo), "ci95_hi": float(hi),
        "n_patients_met_lower": int((d < 0).sum()),
        "wilcoxon_p": float(w_p),
        "sign_test_p": float(stats.binomtest(int((d < 0).sum()), n, 0.5).pvalue),
    }


def main():
    ap = argparse.ArgumentParser(description="mLN 해석의 구성 교란 분리")
    ap.add_argument("--h5ad", required=True,
                    help="step3 산출물 (…_full_annotated_cnv.h5ad)")
    ap.add_argument("--outdir", default="./results_gse198315")
    ap.add_argument("--min-cells", type=int, default=MIN_CELLS)
    args = ap.parse_args()

    import scanpy as sc
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(SEED)

    log(f"적재: {args.h5ad}")
    a = sc.read_h5ad(args.h5ad)
    for col in ["cell_type", "patient", "region"]:
        if col not in a.obs.columns:
            raise SystemExit(f"오류: obs에 '{col}' 열이 없습니다. step2/3 산출물이 맞는지 확인하십시오.")

    gene = next((g for g in GENE_ALIASES if g in a.var_names), None)
    if gene is None:
        raise SystemExit(f"오류: {GENE_ALIASES} 를 찾지 못했습니다.")
    X = a[:, gene].X
    v = pd.Series(np.asarray(X.todense()).ravel() if hasattr(X, "todense")
                  else np.asarray(X).ravel(), index=a.obs_names)
    thr = v[v > 0].quantile(0.9)          # 원고 정의와 동일: 비영 세포 중 상위 10%
    high = (v > thr)
    log(f"  {a.shape}, WBP5-high 임계 {thr:.4f}, high 세포 {int(high.sum()):,}")

    obs = a.obs.copy()
    obs["WBP5"] = v.values
    obs["WBP5_high"] = high.values

    # ---- (0) 부위별 세포 구성 — 교란의 원인을 수치로 제시 -------------------
    comp = pd.crosstab(obs["region"], obs["cell_type"])
    comp_pct = (comp.T / comp.sum(axis=1)).T * 100
    comp.to_csv(f"{args.outdir}/I1_region_composition_counts.csv")
    comp_pct.round(2).to_csv(f"{args.outdir}/I1b_region_composition_pct.csv")
    lymph = [c for c in ["T/NK", "B/Plasma"] if c in comp_pct.columns]
    if lymph:
        log("  부위별 림프구 비율: " +
            ", ".join(f"{r} {comp_pct.loc[r, lymph].sum():.1f}%" for r in comp_pct.index))

    # ---- (A) 집계 비교 (교란된 값) ------------------------------------------
    agg = (obs.groupby(["patient", "region"], observed=True)
              .agg(n_cells=("WBP5_high", "size"),
                   pct_high=("WBP5_high", lambda s: s.mean() * 100),
                   mean_WBP5=("WBP5", "mean")).reset_index())
    agg.to_csv(f"{args.outdir}/I2_patient_region_aggregate.csv", index=False)

    rows_a = []
    for p, g in agg.groupby("patient", observed=True):
        m = g.loc[g.region == MET, "pct_high"]
        pr = g.loc[g.region.isin(PRIMARY)]
        if len(m) and len(pr):
            # 원발은 세포수 가중 평균
            w = pr["n_cells"].to_numpy(float)
            rows_a.append({"patient": p, "met": float(m.iloc[0]),
                           "primary": float((pr["pct_high"] * w).sum() / w.sum())})
    dfa = pd.DataFrame(rows_a)

    # ---- (B) 세포형 내부 비교 (올바른 값) -----------------------------------
    per = (obs.groupby(["cell_type", "patient", "region"], observed=True)
              .agg(n_cells=("WBP5_high", "size"),
                   n_high=("WBP5_high", "sum"),
                   pct_high=("WBP5_high", lambda s: s.mean() * 100),
                   mean_WBP5=("WBP5", "mean")).reset_index())
    per.to_csv(f"{args.outdir}/I3_patient_celltype_region.csv", index=False)

    results = []
    if len(dfa):
        r = paired_test(dfa, "ALL CELLS POOLED (composition-confounded)", rng)
        r["cell_type"] = "— aggregate —"
        results.append(r)

    for ct, g in per.groupby("cell_type", observed=True):
        rows = []
        for p, gp in g.groupby("patient", observed=True):
            m = gp[(gp.region == MET) & (gp.n_cells >= args.min_cells)]
            pr = gp[gp.region.isin(PRIMARY) & (gp.n_cells >= args.min_cells)]
            if len(m) and len(pr):
                w = pr["n_cells"].to_numpy(float)
                rows.append({"patient": p, "met": float(m["pct_high"].iloc[0]),
                             "primary": float((pr["pct_high"] * w).sum() / w.sum())})
        if len(rows) >= 3:
            r = paired_test(pd.DataFrame(rows), f"within {ct}", rng)
            r["cell_type"] = ct
            results.append(r)

    res = pd.DataFrame(results)
    cols = ["cell_type", "contrast", "n_patients", "median_pct_met", "median_pct_primary",
            "median_diff_met_minus_primary", "ci95_lo", "ci95_hi",
            "n_patients_met_lower", "wilcoxon_p", "sign_test_p"]
    res = res[[c for c in cols if c in res.columns]]
    res.to_csv(f"{args.outdir}/I4_mLN_paired_tests.csv", index=False)

    pd.set_option("display.width", 220)
    print("\n=== mLN vs primary (TC+TP), 환자 단위 paired 검정 ===")
    print(res.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("\n=== 부위별 세포 구성 (%) ===")
    print(comp_pct.round(1).to_string())

    with open(f"{args.outdir}/environment_step4.txt", "w") as fh:
        import scipy
        fh.write(f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n")
        fh.write(f"scipy {scipy.__version__}\nseed {SEED}\nn_boot {N_BOOT}\n")
        fh.write(f"min_cells_per_patient_celltype_region {args.min_cells}\n")
        fh.write(f"WBP5_high_threshold {thr}\ninput {os.path.abspath(args.h5ad)}\n")
    log("완료")


if __name__ == "__main__":
    main()
