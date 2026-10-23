#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE198315 step7 — inferCNV 없이 다시 정의한 세포 집합으로 재분석  (2026-10-19)
==========================================================================
배경: step6 판정 STILL CONFOUNDED (정상 상피 포함 reference 에서도 종양 vs 정상 상피
      CNV score Cliff's δ = −0.052). 따라서 inferCNV 는 악성 판정에 쓰지 않습니다.
      step5 (C) 의 10/10 결과는 immune-only inferCNV 의 CNV+ 10,203개
      (정상 상피 ~2,770 + 비상피 ~1,744 포함)로 계산된 것이라 그대로 보고할 수 없습니다.

세포 집합 정의
  EPI_TUMOR  (주분석) : cell_type == Epithelial  AND  region ∈ {TC, TP}
  EPI_ALL    (민감도) : cell_type == Epithelial  (모든 부위, NT 포함)
  OLD_CNV    (참고)   : step3 cnv_malignant == True  ← step5 (C) 가 쓴 집합, 비교용으로만

실행:
  python GSE198315_step7_20261019.py ^
     --h5ad gse198315_checkpoint_full_annotated_cnv.h5ad ^
     --pemt p-EMT_gene_list.xlsx ^
     --outdir results_gse198315

산출 (M = 이 스크립트 접두어)
  M0  집합 구성표 (OLD_CNV 가 실제로 무엇으로 이뤄졌는지 포함)
  M1  환자 내 WBP5–마커 Spearman (세 집합), sign test, BH
  M2  WBP5 ~ Puram 서명 (EPI_TUMOR) + ΔR²  ← GSE103322 에서 도출된 서명의 독립 코호트 검증
  M3  전사체 전체 WBP5 상관에서 Puram 유전자 농축 (EPI_TUMOR)
  M4  구획 대조: EPI_TUMOR vs 다른 세포형 (TC/TP 안에서), Cliff's δ + 환자 수준 paired
  M5  WBP5-high 임계 민감도 (상위 5/10/20% of nonzero, 상위 10% of all, 환자별 10%)
"""
import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import wbp5_rev_common as C
except ImportError:
    raise SystemExit("오류: wbp5_rev_common.py 를 이 스크립트와 같은 폴더에 두십시오.")

warnings.filterwarnings("ignore")
TUMOR_REGIONS = ["TC", "TP"]
MIN_CELLS = 20


def col(a, g):
    x = a[:, g].X
    return np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()


def block(a, genes):
    x = a[:, genes].X
    x = x.toarray() if sp.issparse(x) else np.asarray(x)
    return x.T.astype(np.float32)   # 유전자 x 세포


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--pemt", required=True)
    ap.add_argument("--outdir", default="results_gse198315")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--min-det", type=float, default=0.01,
                    help="M3 에서 쓸 유전자의 최소 검출률 (EPI_TUMOR 기준)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(C.SEED)
    out = lambda f: os.path.join(args.outdir, f)  # noqa

    import anndata as ad
    C.log(f"적재: {args.h5ad}")
    a = ad.read_h5ad(args.h5ad)
    for c in ("cell_type", "region", "patient"):
        if c not in a.obs:
            raise SystemExit(f"오류: obs 에 '{c}' 없음. 있는 열: {list(a.obs.columns)}")
        a.obs[c] = a.obs[c].astype(str)
    # 정규화 상태 확인 (step2 의 CPM 1e4 + log1p 여야 함)
    xs = a.X[:2000].toarray() if sp.issparse(a.X) else np.asarray(a.X[:2000])
    if xs.max() > 50 and np.allclose(xs, np.round(xs)):
        C.log("  X 가 raw count 로 보임 → CPM 1e4 + log1p 적용")
        import scanpy as sc
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    gene = C.resolve_one(C.WBP5_ALIASES, a.var_names)
    if gene is None:
        raise SystemExit("오류: WBP5/TCEAL9 없음")
    C.log(f"  세포 {a.n_obs:,} / 유전자 {a.n_vars:,} / 대상 유전자 {gene}")
    C.log(f"  region 값: {sorted(a.obs.region.unique())}")
    C.log(f"  cell_type 값: {sorted(a.obs.cell_type.unique())}")

    ct, rg = a.obs.cell_type.values, a.obs.region.values
    sets = {
        "EPI_TUMOR": (ct == "Epithelial") & np.isin(rg, TUMOR_REGIONS),
        "EPI_ALL": ct == "Epithelial",
    }
    if "cnv_malignant" in a.obs:
        sets["OLD_CNV"] = a.obs.cnv_malignant.astype(str).str.strip().str.lower().isin(["true", "1"]).values
    else:
        C.log("  cnv_malignant 없음 → OLD_CNV 비교 생략 (_cnv.h5ad 를 주십시오)")
    if sets["EPI_TUMOR"].sum() < 500:
        raise SystemExit(f"오류: EPI_TUMOR {sets['EPI_TUMOR'].sum()}개 — region/cell_type 표기 확인")

    # ---- M0 ----------------------------------------------------------------
    rows = []
    for nm, m in sets.items():
        sub = a.obs[m]
        comp = sub.groupby(["cell_type", "region"], observed=True).size().rename("n").reset_index()
        comp.insert(0, "set", nm)
        rows.append(comp)
        C.log(f"  {nm}: {m.sum():,} 세포, 환자 {sub.patient.nunique()}명")
    m0 = pd.concat(rows)
    m0.to_csv(out("M0_set_composition.csv"), index=False)
    w_all = col(a, gene)

    # ---- M1 ----------------------------------------------------------------
    C.log("M1 환자 내 WBP5–마커 상관")
    rows, per_rows = [], []
    for nm, m in sets.items():
        pats = a.obs.patient.values[m]
        w = w_all[m]
        summ = []
        for mk in C.MARKER_PANEL:
            if mk not in a.var_names:
                continue
            per, s = C.within_group_spearman(w, col(a, mk)[m], pats, MIN_CELLS)
            per["set"], per["marker"] = nm, mk
            per_rows.append(per)
            s.update({"set": nm, "marker": mk})
            summ.append(s)
        s = pd.DataFrame(summ)
        s["sign_test_q"] = multipletests(s.sign_test_p, method="fdr_bh")[1]
        rows.append(s)
    m1 = pd.concat(rows)
    m1.to_csv(out("M1_within_patient_markers.csv"), index=False)
    pd.concat(per_rows).to_csv(out("M1_within_patient_per_patient.csv"), index=False)

    # ---- 서명 ---------------------------------------------------------------
    plist = C.read_signature(args.pemt)
    used100, map100 = C.resolve_genes(plist, a.var_names)
    used30, _ = C.resolve_genes(plist[:30], a.var_names)
    map100.insert(0, "rank", range(1, len(map100) + 1))
    map100.to_csv(out("M2_signature_gene_mapping.csv"), index=False)
    miss = map100[map100.match == "missing"].list_symbol.tolist()
    C.log(f"  Puram 100개 중 {len(used100)}개 매핑, 누락 {miss}; "
          f"별칭 사용 {map100[map100.match.str.startswith('alias')].list_symbol.tolist()}")
    sig = {"PURAM30": used30, "PURAM100": used100,
           "PURAM30_minus_panel": [g for g in used30 if g not in C.MARKER_PANEL],
           "PURAM100_minus_panel": [g for g in used100 if g not in C.MARKER_PANEL]}
    for gl in sig.values():
        if gene in gl:
            raise SystemExit("오류: 서명에 WBP5 포함")

    # ---- M2 (EPI_TUMOR) ------------------------------------------------------
    C.log("M2 WBP5 ~ Puram 서명 (EPI_TUMOR)")
    m = sets["EPI_TUMOR"]
    e = a[m]
    pats = e.obs.patient.values
    w = w_all[m]
    scores = {}
    for nm, gl in sig.items():
        scores[nm] = C.zmean_score(block(e, gl))
    rows, per_rows = [], []
    for nm in ("PURAM30", "PURAM100"):
        r = {"signature": nm, "n_genes": len(sig[nm])}
        r.update(C.spearman_cluster_boot(w, scores[nm], pats, rng, args.n_boot))
        per, s = C.within_group_spearman(w, scores[nm], pats, MIN_CELLS)
        r.update({f"within_{k}": v for k, v in s.items()})
        r.update(C.mixed_model(w, scores[nm], pats))
        rows.append(r)
        per["signature"] = nm
        per_rows.append(per)
    m2 = pd.DataFrame(rows)
    preds = {mk: col(e, mk) for mk in C.MARKER_PANEL if mk in e.var_names}
    inc = []
    for nm in ("PURAM30_minus_panel", "PURAM100_minus_panel"):
        r = {"response_signature": nm, "n_genes": len(sig[nm])}
        r.update(C.incremental_information(scores[nm], preds, w, pats))
        inc.append(r)
    m2i = pd.DataFrame(inc)
    m2.to_csv(out("M2_WBP5_vs_Puram_signature.csv"), index=False)
    pd.concat(per_rows).to_csv(out("M2_within_patient_per_patient.csv"), index=False)
    m2i.to_csv(out("M2_incremental_information.csv"), index=False)
    pd.DataFrame({"cell": e.obs_names, "patient": pats, "region": e.obs.region.values,
                  "WBP5": w, **scores}).to_csv(out("M2_per_cell_scores.csv.gz"), index=False)

    # ---- M3 ----------------------------------------------------------------
    C.log("M3 전사체 전체 상관 + Puram 농축 (EPI_TUMOR, 수 분 소요)")
    X = e.X.tocsc() if sp.issparse(e.X) else sp.csc_matrix(e.X)
    det = np.asarray((X > 0).mean(axis=0)).ravel()
    keep = np.where((det >= args.min_det) & (e.var_names.values != gene))[0]
    names = e.var_names.values[keep]
    M = X[:, keep].toarray().T.astype(np.float32)
    bins = C.expression_bins(M.mean(axis=1))
    C.log(f"  검출률 ≥ {args.min_det:.0%} 유전자 {len(keep):,}개")
    res, gtabs = [], []
    for mode, grp in (("pooled", None), ("within_patient", pats)):
        rho = C.corr_with_all(w, M, grp)
        pd.DataFrame({"gene": names, "rho": rho, "detection": det[keep],
                      "in_PURAM100": np.isin(names, used100),
                      "in_PURAM30": np.isin(names, used30)}).sort_values(
            "rho", ascending=False).to_csv(out(f"M3_transcriptome_rho_{mode}.csv.gz"), index=False)
        for nm in ("PURAM30", "PURAM100"):
            present = [g for g in sig[nm] if g in set(names)]
            r, pg = C.set_enrichment(rho, names, bins, present, rng, args.n_perm)
            r.update({"mode": mode, "signature": nm, "n_genes_tested": int(np.isfinite(rho).sum()),
                      "all_genes_median": float(np.nanmedian(rho)),
                      "all_genes_99pct": float(np.nanpercentile(rho, 99))})
            res.append(r)
            pg["mode"], pg["signature"] = mode, nm
            gtabs.append(pg)
        top = names[np.argsort(-np.nan_to_num(rho, nan=-9))][:30]
        C.log(f"  [{mode}] 상위 30: {', '.join(top)}")
    del M
    m3 = pd.DataFrame(res)
    m3.to_csv(out("M3_Puram_enrichment.csv"), index=False)
    pd.concat(gtabs).to_csv(out("M3_Puram_gene_ranks.csv"), index=False)

    # ---- M4 ----------------------------------------------------------------
    C.log("M4 구획 대조 (TC/TP 안에서)")
    tumor_reg = np.isin(rg, TUMOR_REGIONS)
    epi = sets["EPI_TUMOR"]
    rows = []
    for t in sorted(set(ct[tumor_reg]) - {"Epithelial"}):
        mm = tumor_reg & (ct == t)
        d, p = C.cliffs_delta(w_all[epi], w_all[mm])
        rows.append({"contrast": f"EPI_TUMOR_vs_{t}", "n_epi": int(epi.sum()), "n_other": int(mm.sum()),
                     "mean_epi": w_all[epi].mean(), "mean_other": w_all[mm].mean(),
                     "det_epi": (w_all[epi] > 0).mean(), "det_other": (w_all[mm] > 0).mean(),
                     "cliffs_delta": d, "mwu_p": p})
    m4 = pd.DataFrame(rows)
    m4.to_csv(out("M4_compartment_contrasts_TCTP.csv"), index=False)
    # 환자 수준 paired
    pat_all = a.obs.patient.values
    rows, pp_rows = [], []
    for t in [x for x in ("Fibroblast", "Endothelial", "Ambiguous", "Myofibroblast") if x in set(ct)]:
        dd = []
        for p in pd.unique(pat_all):
            ia = epi & (pat_all == p)
            ib = tumor_reg & (ct == t) & (pat_all == p)
            if ia.sum() >= 10 and ib.sum() >= 10:
                dd.append(w_all[ia].mean() - w_all[ib].mean())
                pp_rows.append({"comparison": f"EPI_TUMOR_vs_{t}", "patient": p,
                                "n_epi": int(ia.sum()), "n_comp": int(ib.sum()),
                                "mean_epi": w_all[ia].mean(), "mean_comp": w_all[ib].mean()})
        dd = np.array(dd)
        r = {"comparison": f"EPI_TUMOR_vs_{t}", "n_patients": len(dd)}
        if len(dd) >= 3:
            boot = [rng.choice(dd, len(dd)).mean() for _ in range(10000)]
            r.update({"mean_diff": dd.mean(), "ci_lo": np.percentile(boot, 2.5),
                      "ci_hi": np.percentile(boot, 97.5), "n_positive": int((dd > 0).sum()),
                      "wilcoxon_p": stats.wilcoxon(dd).pvalue if len(dd) >= 5 else np.nan})
        rows.append(r)
    m4p = pd.DataFrame(rows)
    m4p.to_csv(out("M4_compartment_patient_level.csv"), index=False)
    pd.DataFrame(pp_rows).to_csv(out("M4_compartment_per_patient.csv"), index=False)

    # ---- M5 ----------------------------------------------------------------
    C.log("M5 WBP5-high 임계 민감도")
    nz = w_all > 0
    defs = {}
    for q in (5, 10, 20):
        thr = np.percentile(w_all[nz], 100 - q)
        defs[f"top{q}pct_of_nonzero_global"] = nz & (w_all >= thr)
    k = a.n_obs // 10
    thr_all = np.sort(w_all)[-k]
    defs["top10pct_of_all_global"] = (w_all >= thr_all) & nz
    if thr_all <= 0:
        C.log("  주의: 0 이 아닌 세포가 전체의 10% 미만 → 'top 10% of all' 은 nonzero 전체와 같음")
    hi = np.zeros(a.n_obs, bool)
    for p in pd.unique(pat_all):
        s = (pat_all == p) & nz
        if s.sum() >= 10:
            hi |= s & (w_all >= np.percentile(w_all[s], 90))
    defs["top10pct_of_nonzero_per_patient"] = hi
    rows = []
    for nm, h in defs.items():
        base = {"definition": nm, "n_high": int(h.sum()),
                "pct_of_all_cells": 100 * h.mean(), "pct_of_nonzero": 100 * h.sum() / nz.sum()}
        share = pd.Series(ct[h]).value_counts(normalize=True) * 100
        for t, v in share.items():
            base[f"share_of_high_%_{t}"] = v
        for t in sorted(set(ct)):
            mm = tumor_reg & (ct == t)
            if mm.sum():
                base[f"pct_high_in_TCTP_{t}"] = 100 * h[mm].mean()
        rows.append(base)
    m5 = pd.DataFrame(rows)
    m5.to_csv(out("M5_threshold_sensitivity.csv"), index=False)

    # ---- 요약 출력 ----------------------------------------------------------
    pd.set_option("display.width", 220)
    f = lambda v: f"{v:.4g}"  # noqa
    print("\n=== M0 OLD_CNV(step5 C 가 쓴 집합)의 실제 구성 ===")
    if "OLD_CNV" in sets:
        oc = a.obs[sets["OLD_CNV"]]
        tab = pd.crosstab(oc.cell_type, oc.region, margins=True)
        print(tab.to_string())
    print("\n=== M1 환자 내 WBP5–마커 Spearman ===")
    print(m1[["set", "marker", "median_rho", "n_positive", "n_patients", "sign_test_p", "sign_test_q",
              "lopo_median_min", "lopo_median_max"]].to_string(index=False, float_format=f))
    print("\n=== M2 WBP5 ~ Puram 서명 (EPI_TUMOR) — 독립 코호트 ===")
    print(m2[["signature", "n_genes", "rho_cell", "rho_ci_lo", "rho_ci_hi", "within_median_rho",
              "within_n_positive", "within_n_patients", "within_sign_test_p",
              "mixedlm_std_coef", "mixedlm_ci_lo", "mixedlm_ci_hi", "mixedlm_p"]]
          .to_string(index=False, float_format=f))
    for _, r in m2i.iterrows():
        print(f"  ΔR² [{r.response_signature}] {r.get('R2_markers', np.nan):.4f} → "
              f"{r.get('R2_markers_plus_WBP5', np.nan):.4f} (Δ {r.get('delta_R2_WBP5', np.nan):.4f}); "
              f"WBP5 부분계수(혼합) {r.get('WBP5_std_partial_coef_mixedlm', np.nan):.4f} "
              f"[{r.get('WBP5_partial_ci_lo_mixedlm', np.nan):.4f}, {r.get('WBP5_partial_ci_hi_mixedlm', np.nan):.4f}]")
    print("\n=== M3 Puram 유전자 농축 (EPI_TUMOR) ===")
    print(m3[["mode", "signature", "set_size_in_data", "mean_rho_set", "null_mean", "null_99pct",
              "perm_p_upper", "top100_overlap", "top100_expected", "top100_hypergeom_p"]]
          .to_string(index=False, float_format=f))
    print("\n=== M4 구획 대조 (TC/TP) ===")
    print(m4[["contrast", "n_epi", "n_other", "det_epi", "det_other", "cliffs_delta"]]
          .to_string(index=False, float_format=f))
    print(m4p.to_string(index=False, float_format=f))
    print("\n=== M5 임계 민감도 (WBP5-high 중 세포형 비율 %) ===")
    show = ["definition", "n_high", "pct_of_all_cells"] + [c for c in m5.columns if c.startswith("share_of_high_%_")]
    print(m5[show].to_string(index=False, float_format=lambda v: f"{v:.1f}"))

    import scipy, statsmodels  # noqa
    with open(out("environment_step7.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
                 f"scipy {scipy.__version__}\nstatsmodels {statsmodels.__version__}\n"
                 f"anndata {ad.__version__}\nseed {C.SEED}\nn_perm {args.n_perm}\n"
                 f"input {os.path.abspath(args.h5ad)}\npemt {os.path.abspath(args.pemt)}\n"
                 f"missing_pemt_symbols {miss}\n")
    C.log(f"완료 → {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
