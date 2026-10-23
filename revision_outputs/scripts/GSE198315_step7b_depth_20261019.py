#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE198315 step7b — 시퀀싱 깊이(depth) 보정 점검  (2026-10-19)
==============================================================
왜 필요한가
  step7 M3 에서 EPI_TUMOR 의 WBP5 전사체 상관 귀무분포 평균이 0.09–0.11 로 0 에서
  크게 벗어났습니다. 대부분의 유전자가 WBP5 와 양의 상관을 보인다는 뜻이고,
  이는 세포별 깊이(검출 유전자 수)가 공통 축으로 작용할 때의 전형적 패턴입니다.
  제출 원고 단계에서도 이 코호트에서 WBP5 ~ log depth ρ = 0.63 을 확인한 바 있습니다.
  따라서 M1/M2 의 양의 상관(LAMC2 0.17, ITGA5 0.18, Puram 0.20–0.24)이
  깊이 축을 넘어서는 신호인지 확인해야 응답서에 쓸 수 있습니다.

방법 (EPI_TUMOR = Epithelial ∩ TC/TP, 환자 내 순위 기반)
  깊이 지표: obs 의 total_counts / n_counts 가 있으면 log, 없으면 X 에서 검출 유전자 수.
  환자 내 Spearman 을 깊이에 대해 편상관으로 보정:
      ρ_partial = (ρ_gw − ρ_gd ρ_wd) / sqrt((1 − ρ_gd²)(1 − ρ_wd²))
  D1  WBP5 ~ 깊이 (환자별)
  D2  마커 5개 + Puram30/100 점수: 원래 ρ vs 깊이 보정 ρ, 환자별 sign test
  D3  깊이 보정 후 전사체 전체 귀무분포 (평균이 0 근처로 돌아오는지) + 마커의
      백분위 (전체 유전자 / 같은 발현구간 유전자 기준)
  D4  깊이 보정 상관에서 Puram 농축 재검정 + 상위 30 유전자

실행:
  python GSE198315_step7b_depth_20261019.py ^
     --h5ad gse198315_checkpoint_full_annotated_cnv.h5ad ^
     --pemt "p-EMT gene list.xlsx" --outdir results_gse198315
"""
import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import wbp5_rev_common as C
except ImportError:
    raise SystemExit("오류: wbp5_rev_common.py 를 이 스크립트와 같은 폴더에 두십시오.")

warnings.filterwarnings("ignore")
TUMOR_REGIONS = ["TC", "TP"]
MIN_CELLS = 20


def zrank_within(v, groups):
    """환자 내 순위를 환자 내 z-표준화."""
    out = np.zeros(len(v))
    for g in pd.unique(groups):
        s = groups == g
        r = stats.rankdata(v[s])
        sd = r.std()
        out[s] = (r - r.mean()) / sd if sd > 0 else 0.0
    return out


def zrank_within_rows(M, groups):
    R = np.zeros(M.shape, dtype=np.float32)
    for g in pd.unique(groups):
        s = np.where(groups == g)[0]
        r = stats.rankdata(M[:, s], axis=1).astype(np.float32)
        r -= r.mean(axis=1, keepdims=True)
        sd = r.std(axis=1, keepdims=True)
        sd[sd == 0] = 1.0
        R[:, s] = r / sd
    return R


def partial(r_gw, r_gd, r_wd):
    with np.errstate(invalid="ignore", divide="ignore"):
        return (r_gw - r_gd * r_wd) / np.sqrt((1 - r_gd ** 2) * (1 - r_wd ** 2))


def per_patient(x, w, d, groups):
    """환자별 원래 ρ 와 깊이 편상관 ρ."""
    rows = []
    for g in pd.unique(groups):
        s = groups == g
        if s.sum() < MIN_CELLS or np.std(x[s]) == 0:
            continue
        rgw = stats.spearmanr(x[s], w[s])[0]
        rgd = stats.spearmanr(x[s], d[s])[0]
        rwd = stats.spearmanr(w[s], d[s])[0]
        rows.append({"patient": g, "n_cells": int(s.sum()), "rho": rgw,
                     "rho_depth_partial": partial(rgw, rgd, rwd), "rho_with_depth": rgd})
    return pd.DataFrame(rows)


def summarize(per, col):
    r = per[col].dropna().values
    npos = int((r > 0).sum())
    return {f"{col}_median": float(np.median(r)), f"{col}_n_pos": npos, "n_patients": len(r),
            f"{col}_sign_p": stats.binomtest(npos, len(r), 0.5).pvalue}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--pemt", required=True)
    ap.add_argument("--outdir", default="results_gse198315")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--min-det", type=float, default=0.01)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(C.SEED)
    out = lambda f: os.path.join(args.outdir, f)  # noqa

    import anndata as ad
    C.log(f"적재: {args.h5ad}")
    a = ad.read_h5ad(args.h5ad)
    for c in ("cell_type", "region", "patient"):
        a.obs[c] = a.obs[c].astype(str)
    gene = C.resolve_one(C.WBP5_ALIASES, a.var_names)
    m = (a.obs.cell_type.values == "Epithelial") & np.isin(a.obs.region.values, TUMOR_REGIONS)
    e = a[m].copy()
    del a
    pats = e.obs.patient.values
    X = e.X.tocsc() if sp.issparse(e.X) else sp.csc_matrix(e.X)

    # ---- 깊이 지표 ----------------------------------------------------------
    src = None
    for c in ("total_counts", "n_counts", "nCount_RNA"):
        if c in e.obs:
            depth, src = np.log1p(e.obs[c].astype(float).values), f"log1p(obs.{c})"
            break
    if src is None:
        depth = np.asarray((X > 0).sum(axis=1)).ravel().astype(float)
        src = "검출 유전자 수 (X > 0)"
    C.log(f"  EPI_TUMOR {e.n_obs:,} 세포 / 깊이 지표: {src}  (obs 열: {list(e.obs.columns)})")

    w = np.asarray(X[:, e.var_names.get_loc(gene)].todense()).ravel()

    # ---- D1 ----------------------------------------------------------------
    d1 = []
    for g in pd.unique(pats):
        s = pats == g
        d1.append({"patient": g, "n_cells": int(s.sum()),
                   "rho_WBP5_depth": stats.spearmanr(w[s], depth[s])[0]})
    d1 = pd.DataFrame(d1)
    d1.to_csv(out("D1_WBP5_vs_depth.csv"), index=False)

    # ---- 서명 점수 (step7 과 동일한 z-평균) ------------------------------------
    plist = C.read_signature(args.pemt)
    used100, _ = C.resolve_genes(plist, e.var_names)
    used30, _ = C.resolve_genes(plist[:30], e.var_names)

    def block(gl):
        idx = [e.var_names.get_loc(g) for g in gl]
        return X[:, idx].toarray().T

    feats = {mk: np.asarray(X[:, e.var_names.get_loc(mk)].todense()).ravel()
             for mk in C.MARKER_PANEL if mk in e.var_names}
    feats["PURAM30_score"] = C.zmean_score(block(used30))
    feats["PURAM100_score"] = C.zmean_score(block(used100))

    # ---- D2 ----------------------------------------------------------------
    C.log("D2 마커·서명: 원래 ρ vs 깊이 보정 ρ (환자별)")
    rows, per_all = [], []
    for nm, x in feats.items():
        per = per_patient(x, w, depth, pats)
        per["feature"] = nm
        per_all.append(per)
        r = {"feature": nm}
        r.update(summarize(per, "rho"))
        r.update(summarize(per, "rho_depth_partial"))
        r["median_rho_feature_vs_depth"] = float(per.rho_with_depth.median())
        rows.append(r)
    d2 = pd.DataFrame(rows)
    d2.to_csv(out("D2_depth_adjusted_correlations.csv"), index=False)
    pd.concat(per_all).to_csv(out("D2_depth_adjusted_per_patient.csv"), index=False)

    # ---- D3 / D4 -----------------------------------------------------------
    C.log("D3 전사체 전체 (깊이 보정) — 수 분 소요")
    det = np.asarray((X > 0).mean(axis=0)).ravel()
    keep = np.where((det >= args.min_det) & (e.var_names.values != gene))[0]
    names = e.var_names.values[keep]
    M = X[:, keep].toarray().T.astype(np.float32)
    bins = C.expression_bins(M.mean(axis=1))
    R = zrank_within_rows(M, pats)
    del M
    tw, td = zrank_within(w, pats), zrank_within(depth, pats)
    n = len(tw)
    r_gw = (R.astype(np.float64) @ tw) / n
    r_gd = (R.astype(np.float64) @ td) / n
    r_wd = float(np.mean(tw * td))
    del R
    r_p = partial(r_gw, r_gd, r_wd)
    tab = pd.DataFrame({"gene": names, "rho_within": r_gw, "rho_depth_partial": r_p,
                        "rho_gene_vs_depth": r_gd, "detection": det[keep], "bin": bins,
                        "in_PURAM30": np.isin(names, used30), "in_PURAM100": np.isin(names, used100)})
    tab.sort_values("rho_depth_partial", ascending=False).to_csv(
        out("D3_transcriptome_depth_partial.csv.gz"), index=False)

    d3 = []
    for col in ("rho_within", "rho_depth_partial"):
        v = tab[col].values
        d3.append({"statistic": col, "median_all_genes": np.nanmedian(v),
                   "mean_all_genes": np.nanmean(v), "p95": np.nanpercentile(v, 95),
                   "p99": np.nanpercentile(v, 99)})
    d3 = pd.DataFrame(d3)
    mk_rows = []
    for mk in C.MARKER_PANEL:
        if mk not in set(names):
            continue
        rr = tab[tab.gene == mk].iloc[0]
        same = tab[tab.bin == rr.bin]
        r = {"marker": mk}
        for col in ("rho_within", "rho_depth_partial"):
            r[col] = rr[col]
            r[f"{col}_pct_all"] = 100 * np.nanmean(tab[col].values < rr[col])
            r[f"{col}_pct_same_bin"] = 100 * np.nanmean(same[col].values < rr[col])
        mk_rows.append(r)
    d3m = pd.DataFrame(mk_rows)
    d3.to_csv(out("D3_null_distribution.csv"), index=False)
    d3m.to_csv(out("D3_marker_percentiles.csv"), index=False)

    C.log("D4 Puram 농축 (깊이 보정 상관)")
    d4 = []
    for nm, gl in (("PURAM30", used30), ("PURAM100", used100)):
        present = [g for g in gl if g in set(names)]
        r, pg = C.set_enrichment(r_p, names, bins, present, rng, args.n_perm)
        r["signature"] = nm
        d4.append(r)
        pg.to_csv(out(f"D4_{nm}_gene_ranks_depth_partial.csv"), index=False)
    d4 = pd.DataFrame(d4)
    d4.to_csv(out("D4_Puram_enrichment_depth_partial.csv"), index=False)
    top = tab.sort_values("rho_depth_partial", ascending=False).gene.head(30).tolist()

    # ---- 출력 --------------------------------------------------------------
    pd.set_option("display.width", 220)
    f = lambda v: f"{v:.4g}"  # noqa
    print(f"\n=== D1 WBP5 ~ 깊이 [{src}] (환자별 Spearman) ===")
    print(f"  중앙값 {d1.rho_WBP5_depth.median():.3f}, 범위 {d1.rho_WBP5_depth.min():.3f} ~ {d1.rho_WBP5_depth.max():.3f}")
    print("\n=== D2 원래 ρ vs 깊이 보정 ρ (환자 내) ===")
    print(d2[["feature", "n_patients", "rho_median", "rho_n_pos", "rho_depth_partial_median",
              "rho_depth_partial_n_pos", "rho_depth_partial_sign_p", "median_rho_feature_vs_depth"]]
          .to_string(index=False, float_format=f))
    print("\n=== D3 전사체 전체 귀무분포: 보정 전 vs 후 ===")
    print(d3.to_string(index=False, float_format=f))
    print(d3m.to_string(index=False, float_format=f))
    print("\n=== D4 Puram 농축 (깊이 보정) ===")
    print(d4[["signature", "set_size_in_data", "mean_rho_set", "null_mean", "null_99pct", "perm_p_upper",
              "top100_overlap", "top100_expected", "top100_hypergeom_p", "top200_overlap",
              "top200_hypergeom_p"]].to_string(index=False, float_format=f))
    print(f"\n  깊이 보정 상위 30: {', '.join(top)}")
    C.log(f"완료 → {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
