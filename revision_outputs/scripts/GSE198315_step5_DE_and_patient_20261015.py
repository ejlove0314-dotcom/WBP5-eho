#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE198315 step5 — 유전자 선정 순위 + 환자 수준 검정 + inferCNV 타당성 점검
========================================================================
생성: 2026-10-15
대응:
  (A) inferCNV 타당성 점검      — 우리 스스로 확인해야 할 항목
  (B) Reviewer #4 point 1       — 악성 vs 정상 상피 차등발현 순위에서 WBP5 위치
  (C) Reviewer #5 Major 1       — 환자 수준 WBP5-마커 상관 (GSE103322 분석의 재현)
  (D) Reviewer #5 Major 1       — 환자 수준 악성 vs 비악성 농축 검정

[(A)를 넣은 이유 — 반드시 읽으십시오]
step3 결과에서 상피세포 8,496개 중 8,459개(99.6%)가 CNV 양성으로 판정됐습니다.
그런데 이 코호트에는 정상조직(NT) 유래 상피세포가 2,780개 포함되어 있습니다.
정상 상피까지 전부 CNV 양성으로 불렸다면, 그 신호는 실제 복제수 이상이 아니라
'상피세포와 면역세포의 발현 차이'를 반영한 것일 수 있습니다. inferCNV를 면역세포
reference로 돌릴 때 알려진 한계입니다.

(A)는 부위별 CNV score를 비교해 이를 직접 확인합니다.
  - NT 상피의 CNV score 가 TC/TP 상피보다 뚜렷이 낮으면 → 판정이 타당
  - 비슷하면 → 판정이 상피 정체성에 교란된 것이므로 답변서 문안을 수정해야 함
어느 쪽이든 결과대로 보고합니다. 리뷰어가 먼저 지적하기 전에 우리가 확인하는 편이
낫습니다.

실행:
    python GSE198315_step5_DE_and_patient_20261015.py \
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

warnings.filterwarnings("ignore")

GENE_ALIASES = ["WBP5", "TCEAL9"]
MARKERS = ["LAMC2", "ITGA5", "KRT14", "KRT17", "VIM"]
CONTEXT = ["EGFR", "EPCAM", "PTPRC"]
PRIMARY = ["TC", "TP"]
NORMAL = "NT"
MIN_CELLS = 20
SEED = 42


def log(m):
    sys.stderr.write(f"[{pd.Timestamp.now():%H:%M:%S}] {m}\n"); sys.stderr.flush()


def cliffs_delta(a, b):
    """Cliff's delta = P(a>b) - P(a<b), Mann-Whitney U 기반."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return np.nan
    r = stats.rankdata(np.concatenate([a, b]))
    u1 = r[:n1].sum() - n1 * (n1 + 1) / 2
    return (2 * u1) / (n1 * n2) - 1


def dense(adata, gene):
    X = adata[:, gene].X
    return np.asarray(X.todense()).ravel() if hasattr(X, "todense") else np.asarray(X).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="step3 산출물 (…_full_annotated_cnv.h5ad)")
    ap.add_argument("--outdir", default="./results_gse198315")
    ap.add_argument("--min-cells", type=int, default=MIN_CELLS)
    args = ap.parse_args()

    import scanpy as sc
    os.makedirs(args.outdir, exist_ok=True)

    log(f"적재: {args.h5ad}")
    a = sc.read_h5ad(args.h5ad)
    for c in ["cell_type", "patient", "region", "cnv_score", "cnv_malignant"]:
        if c not in a.obs.columns:
            raise SystemExit(f"오류: obs에 '{c}' 가 없습니다. step3 산출물이 맞는지 확인하십시오.")
    gene = next((g for g in GENE_ALIASES if g in a.var_names), None)
    if gene is None:
        raise SystemExit("오류: WBP5/TCEAL9 를 찾지 못했습니다.")
    log(f"  {a.shape}, 대상 유전자 {gene}")

    # =====================================================================
    # (A) inferCNV 타당성 점검
    # =====================================================================
    log("(A) inferCNV 타당성 점검 — 부위별 CNV score")
    tab = (a.obs.groupby(["cell_type", "region"], observed=True)
             .agg(n_cells=("cnv_score", "size"),
                  mean_cnv_score=("cnv_score", "mean"),
                  median_cnv_score=("cnv_score", "median"),
                  pct_called_malignant=("cnv_malignant", lambda s: s.mean() * 100),
                  n_patients=("patient", "nunique")).reset_index())
    tab.to_csv(f"{args.outdir}/K1_cnv_score_by_celltype_region.csv", index=False)

    ep = a.obs[a.obs.cell_type == "Epithelial"]
    verdict = "판정 불가 (정상조직 상피 부족)"
    if len(ep):
        nt = ep.loc[ep.region == NORMAL, "cnv_score"]
        tu = ep.loc[ep.region.isin(PRIMARY), "cnv_score"]
        print("\n=== (A) 상피세포 CNV score: 정상조직 vs 종양 ===")
        print(f"  NT   n={len(nt):>6}  mean {nt.mean():.4f}  median {nt.median():.4f}  "
              f"악성판정 {100*ep.loc[ep.region==NORMAL,'cnv_malignant'].mean():.1f}%")
        print(f"  TC+TP n={len(tu):>6}  mean {tu.mean():.4f}  median {tu.median():.4f}  "
              f"악성판정 {100*ep.loc[ep.region.isin(PRIMARY),'cnv_malignant'].mean():.1f}%")
        if len(nt) >= 20 and len(tu) >= 20:
            d = cliffs_delta(tu.values, nt.values)
            p = stats.mannwhitneyu(tu.values, nt.values, alternative="two-sided").pvalue
            print(f"  Cliff's delta (종양 - 정상) = {d:+.4f}, MWU p = {p:.3g}")
            if d > 0.3 and 100*ep.loc[ep.region==NORMAL,'cnv_malignant'].mean() < 60:
                verdict = "타당 — 정상 상피가 뚜렷이 낮고 악성 판정률도 낮음"
            elif d > 0.1:
                verdict = "부분적 — 차이는 있으나 정상 상피도 상당수 악성으로 불림. 답변서에 한계 명시 필요"
            else:
                verdict = "교란 의심 — CNV 신호가 상피 정체성을 반영할 가능성. 답변서 문안 수정 필요"
            print(f"  >>> 판정: {verdict}")

    # =====================================================================
    # (B) 악성 vs 정상 상피 차등발현 순위 (Reviewer #4 point 1)
    # =====================================================================
    log("(B) 차등발현 순위 — 종양 상피 vs 정상조직 상피")
    sub = a[(a.obs.cell_type == "Epithelial") &
            (a.obs.region.isin(PRIMARY + [NORMAL]))].copy()
    sub.obs["grp"] = np.where(sub.obs.region == NORMAL, "normal_epi", "tumor_epi")
    n_t = int((sub.obs.grp == "tumor_epi").sum()); n_n = int((sub.obs.grp == "normal_epi").sum())
    log(f"  종양 상피 {n_t:,} / 정상 상피 {n_n:,}")

    if n_n >= args.min_cells and n_t >= args.min_cells:
        sc.tl.rank_genes_groups(sub, "grp", groups=["tumor_epi"], reference="normal_epi",
                                method="wilcoxon", pts=True)
        r = sc.get.rank_genes_groups_df(sub, group="tumor_epi")
        r = r.rename(columns={"names": "gene", "logfoldchanges": "log2FC",
                              "pvals": "p", "pvals_adj": "p_adj"})
        r = r.sort_values("scores", ascending=False).reset_index(drop=True)
        r["rank_by_score"] = np.arange(1, len(r) + 1)
        r.to_csv(f"{args.outdir}/K2_DE_tumor_vs_normal_epithelium_full.csv", index=False)
        r.head(200).to_csv(f"{args.outdir}/K2b_DE_top200.csv", index=False)

        print("\n=== (B) 종양 상피 vs 정상 상피 차등발현 ===")
        print(f"  검정 유전자 {len(r):,}")
        hit = r[r.gene == gene]
        if len(hit):
            row = hit.iloc[0]
            pct = 100 * row.rank_by_score / len(r)
            print(f"  {gene}: 순위 {int(row.rank_by_score):,} / {len(r):,} (상위 {pct:.2f}%)")
            print(f"    log2FC {row.log2FC:+.4f}, p_adj {row.p_adj:.3g}")
        else:
            print(f"  경고: {gene} 가 결과에 없습니다")
        print("  상위 15개:")
        print(r.head(15)[["rank_by_score", "gene", "log2FC", "p_adj"]].to_string(
            index=False, float_format=lambda v: f"{v:.3f}"))
        print(f"  참고 — 마커 패널 순위:")
        for m in MARKERS + CONTEXT:
            h = r[r.gene == m]
            if len(h):
                print(f"    {m:<8} 순위 {int(h.iloc[0].rank_by_score):>6,}  "
                      f"log2FC {h.iloc[0].log2FC:+.3f}")
    else:
        log("  정상조직 상피가 부족해 (B) 생략")

    # =====================================================================
    # (C) 환자 내 WBP5-마커 상관 (Reviewer #5 Major 1)
    # =====================================================================
    log("(C) 환자 내 WBP5-마커 Spearman 상관 (악성세포)")
    mal_mask = a.obs.cnv_malignant.values.astype(bool)
    mal = a[mal_mask].copy()
    v = pd.Series(dense(mal, gene), index=mal.obs_names)
    pats = mal.obs.patient.values

    rows = []
    for m in MARKERS:
        if m not in mal.var_names:
            continue
        mv = pd.Series(dense(mal, m), index=mal.obs_names)
        for p in pd.unique(pats):
            sel = pats == p
            x, y = v.values[sel], mv.values[sel]
            if sel.sum() < args.min_cells or np.std(x) == 0 or np.std(y) == 0:
                continue
            rho, pv = stats.spearmanr(x, y)
            rows.append({"marker": m, "patient": p, "n_cells": int(sel.sum()),
                         "rho": rho, "p_within": pv})
    per = pd.DataFrame(rows)
    per.to_csv(f"{args.outdir}/K3_within_patient_rho_per_patient.csv", index=False)

    if len(per):
        from statsmodels.stats.multitest import multipletests
        summ = []
        for m, g in per.groupby("marker"):
            rh = g.rho.dropna().values
            npos = int((rh > 0).sum())
            summ.append({"marker": m, "median_rho": np.median(rh),
                         "q1": np.percentile(rh, 25), "q3": np.percentile(rh, 75),
                         "n_positive": npos, "n_patients": len(rh),
                         "sign_test_p": stats.binomtest(npos, len(rh), 0.5).pvalue})
        summ = pd.DataFrame(summ)
        summ["sign_test_q"] = multipletests(summ.sign_test_p, method="fdr_bh")[1]
        summ.to_csv(f"{args.outdir}/K4_within_patient_rho_summary.csv", index=False)
        print("\n=== (C) 환자 내 Spearman (CNV+ 악성세포) ===")
        print(summ.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    # =====================================================================
    # (D) 환자 수준 악성 vs 비악성 농축 (Reviewer #5 Major 1)
    # =====================================================================
    log("(D) 환자 수준 악성 vs 비악성 농축 검정")
    vv = pd.Series(dense(a, gene), index=a.obs_names)
    df = pd.DataFrame({"patient": a.obs.patient.values,
                       "malignant": a.obs.cnv_malignant.values.astype(bool),
                       "w": vv.values})
    rows = []
    for p, g in df.groupby("patient"):
        m, n = g[g.malignant], g[~g.malignant]
        if len(m) >= args.min_cells and len(n) >= args.min_cells:
            rows.append({"patient": p, "n_malignant": len(m), "n_nonmalignant": len(n),
                         "mean_malignant": m.w.mean(), "mean_nonmalignant": n.w.mean(),
                         "delta": m.w.mean() - n.w.mean(),
                         "cliffs_delta": cliffs_delta(m.w.values, n.w.values)})
    pl = pd.DataFrame(rows)
    pl.to_csv(f"{args.outdir}/K5_patient_level_enrichment.csv", index=False)
    if len(pl) >= 3:
        d = pl.delta.values
        w_s, w_p = stats.wilcoxon(d, alternative="two-sided")
        rng = np.random.default_rng(SEED)
        boot = np.array([np.mean(rng.choice(d, len(d), replace=True)) for _ in range(10000)])
        lo, hi = np.percentile(boot, [2.5, 97.5])
        print("\n=== (D) 환자 수준 악성 vs 비악성 ===")
        print(f"  환자 {len(pl)}명, 평균 차이 {d.mean():+.4f} "
              f"(95% CI {lo:+.4f} ~ {hi:+.4f})")
        print(f"  양수 환자 {int((d>0).sum())}/{len(d)}, Wilcoxon p = {w_p:.4g}")
        print(f"  환자별 Cliff's delta 중앙값 {pl.cliffs_delta.median():+.4f}")

    with open(f"{args.outdir}/environment_step5.txt", "w", encoding="utf-8") as fh:
        import scipy, statsmodels
        fh.write(f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n")
        fh.write(f"scipy {scipy.__version__}\nscanpy {sc.__version__}\n")
        fh.write(f"statsmodels {statsmodels.__version__}\nseed {SEED}\n")
        fh.write(f"min_cells {args.min_cells}\ninput {os.path.abspath(args.h5ad)}\n")
        fh.write(f"inferCNV_validity_verdict {verdict}\n")
    log("완료")


if __name__ == "__main__":
    main()
