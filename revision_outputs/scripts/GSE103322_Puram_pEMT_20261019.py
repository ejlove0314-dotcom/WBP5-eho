#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE103322 — Puram p-EMT 서명 기반 분석  (2026-10-19)
====================================================
대응: Reviewer 5 Major 1, 2, 9 / Reviewer 4 point 1 / Reviewer 1 point 4

필요 파일 (같은 폴더):
  wbp5_rev_common.py
  GSE103322_HNSCC_all_data.txt.gz
  p-EMT_gene_list.xlsx          (Puram 2017 Table S7 p-EMT 열, Rank/Gene, 100행)
  [선택] yap_taz_22genes.txt     (원고 Section 2.3 의 Wang et al. 22개, 한 줄 한 유전자)

실행:
  python GSE103322_Puram_pEMT_20261019.py ^
     --matrix GSE103322_HNSCC_all_data.txt.gz ^
     --pemt p-EMT_gene_list.xlsx ^
     --outdir results_gse103322_puram

분석 (악성세포 2,215개, 환자 10명; 환자 ID 규칙은 v2 스크립트와 동일)
  P1  WBP5 ~ p-EMT 서명 점수 : PURAM30(주분석, 원 논문 scoring 정의) / PURAM100(민감도)
      세포 수준 ρ + 환자 단위 부트스트랩 95% CI, 환자 내 ρ(중앙값·sign test·LOPO),
      혼합모형(환자 random intercept). 점수는 원고와 같은 z-평균(주) + 대조구간 보정(민감도).
  P2  증분 정보량: PURAM30에서 마커 패널(LAMC2, ITGA5, VIM)을 뺀 27개 점수를 반응변수로,
      마커 5개 ± WBP5 → ΔR², 표준화 부분계수(OLS·혼합모형), drop-one 비교, 상피-기질 PC1.
  P3  전사체 전체 WBP5 상관에서 Puram 유전자 농축 (발현구간 일치 무작위 집합 10,000개,
      상위 100/200/500 hypergeometric). pooled 및 환자 내 순위 두 방식.
  P4  구획별 환자 수준 비교: 악성 vs 섬유아세포 / 내피세포 (환자별 평균, paired).
  P5  [--yap 제공 시] YAP/TAZ: 무작위 집합 10,000개 (발현구간 일치 / 비일치), 양측 p,
      그리고 PURAM30 점수를 통제한 편상관.

주의 (답변서에 명시):
  - Puram 서명은 이 코호트(GSE103322)의 악성세포에서 도출됨 → P1·P3 는 in-sample.
    독립 검증은 GSE198315 (step7) 결과로 제시.
  - 데이터는 log2(TPM/10+1). count 기반 pseudobulk 아님.
"""
import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import wbp5_rev_common as C
except ImportError:
    raise SystemExit("오류: wbp5_rev_common.py 를 이 스크립트와 같은 폴더에 두십시오.")

warnings.filterwarnings("ignore")
MIN_CELLS = 20

# --- v2 스크립트(WBP5_patient_level_analysis_20261015_v2.py)와 동일한 적재 규칙 ---------
META_ROW_KEYS = {
    "processed by maxima enzyme": "maxima_enriched",
    "maxima enriched": "maxima_enriched",
    "lymph node": "lymph_node",
    "classified as cancer cell": "is_cancer",
    "classified as non-cancer cells": "is_noncancer",
    "non-cancer cell type": "noncancer_type",
}


def norm_key(s):
    return " ".join(str(s).replace('"', "").replace("'", "").lower().split())


def parse_patient(cell_name):
    """HN28 / HNSCC28 / HNSCC_17 → HN<번호>, combo 플레이트는 None(제외)."""
    name = str(cell_name).strip().strip("'\"")
    if re.match(r"^HNSCC[_]?combo", name, flags=re.I):
        return None
    m = re.match(r"^HN(?:SCC)?[_]?(\d+)", name, flags=re.I)
    return f"HN{int(m.group(1))}" if m else None


def load_gse103322(path):
    C.log(f"매트릭스 적재: {path}")
    df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    df.index = [str(i).strip().strip("'\"").strip() for i in df.index]
    df.columns = [str(c).strip().strip("'\"").strip() for c in df.columns]
    found = {}
    for idx in df.index[:20]:
        k = norm_key(idx)
        if k in META_ROW_KEYS:
            found[META_ROW_KEYS[k]] = idx
    if "is_cancer" not in found:
        raise SystemExit(f"오류: 'classified as cancer cell' 행 없음. 상단 행: {list(df.index[:8])}")
    meta = pd.DataFrame({n: df.loc[i] for n, i in found.items()})
    expr = df.drop(index=list(found.values())).apply(pd.to_numeric, errors="coerce")
    meta["patient"] = [parse_patient(c) for c in meta.index]
    drop = meta.index[meta["patient"].isna()]
    meta = meta[meta["patient"].notna()]
    expr = expr[meta.index]
    meta["is_cancer"] = pd.to_numeric(meta["is_cancer"], errors="coerce").fillna(0).astype(int)
    C.log(f"  유전자 {expr.shape[0]:,} / 세포 {expr.shape[1]:,} / 환자 {meta.patient.nunique()} "
          f"(combo 플레이트 {len(drop)}개 제외) / 악성 {int(meta.is_cancer.sum()):,}")
    return expr, meta


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--pemt", required=True, help="p-EMT_gene_list.xlsx (Rank, Gene)")
    ap.add_argument("--yap", default=None, help="원고의 YAP/TAZ 22개 유전자 목록 (선택)")
    ap.add_argument("--outdir", default="results_gse103322_puram")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-boot", type=int, default=2000)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(C.SEED)
    out = lambda f: os.path.join(args.outdir, f)  # noqa

    expr, meta = load_gse103322(args.matrix)
    gene = C.resolve_one(C.WBP5_ALIASES, expr.index)
    if gene is None:
        raise SystemExit("오류: WBP5/TCEAL9 없음")

    mal = meta.index[meta.is_cancer == 1]
    E = expr[mal]
    pats = meta.loc[mal, "patient"].values
    w = E.loc[gene].values.astype(float)

    # ---- 서명 정의 -----------------------------------------------------------
    plist = C.read_signature(args.pemt)
    if len(plist) != 100:
        C.log(f"  경고: 목록이 {len(plist)}개 (Table S7 은 100개)")
    used100, map100 = C.resolve_genes(plist, E.index)
    used30, _ = C.resolve_genes(plist[:30], E.index)
    map100.insert(0, "rank", range(1, len(map100) + 1))
    map100.to_csv(out("P0_signature_gene_mapping.csv"), index=False)
    miss = map100[map100.match == "missing"].list_symbol.tolist()
    C.log(f"  매핑: 100개 중 {len(used100)}개 사용, 누락 {miss}")
    for nm, gl in (("PURAM30", used30), ("PURAM100", used100)):
        if gene in gl or "TCEAL9" in gl or "WBP5" in gl:
            raise SystemExit(f"오류: {nm} 에 WBP5 포함 — 자기상관 문제")
    panel_in_30 = [m for m in C.MARKER_PANEL if m in used30]
    sig = {
        "PURAM30": used30,
        "PURAM100": used100,
        "PURAM30_minus_panel": [g for g in used30 if g not in C.MARKER_PANEL],
        "PURAM100_minus_panel": [g for g in used100 if g not in C.MARKER_PANEL],
    }
    C.log(f"  PURAM30 에 포함된 원고 마커: {panel_in_30} → *_minus_panel 서명에서 제외")

    # 점수
    scores = {}
    keep = E.index[E.std(axis=1) > 0]
    Eall = E.loc[keep].values.astype(np.float64)
    for nm, gl in sig.items():
        scores[nm] = C.zmean_score(E.loc[gl].values)
        scores[nm + "_ctrlbin"] = C.control_binned_score(Eall, keep, gl, rng)
    pd.DataFrame({"cell": mal, "patient": pats, "WBP5": w, **scores}).to_csv(
        out("P0_per_cell_scores.csv.gz"), index=False)

    # ---- P1 ----------------------------------------------------------------
    C.log("P1 WBP5 ~ p-EMT 서명")
    rows, per_all = [], []
    for nm in ("PURAM30", "PURAM100", "PURAM30_ctrlbin", "PURAM100_ctrlbin"):
        s = scores[nm]
        r = {"signature": nm, "n_genes": len(sig[nm.replace("_ctrlbin", "")])}
        r.update(C.spearman_cluster_boot(w, s, pats, rng, args.n_boot))
        per, summ = C.within_group_spearman(w, s, pats, MIN_CELLS)
        r.update({f"within_{k}": v for k, v in summ.items()})
        r.update(C.mixed_model(w, s, pats))
        rows.append(r)
        per["signature"] = nm
        per_all.append(per)
    p1 = pd.DataFrame(rows)
    p1.to_csv(out("P1_WBP5_vs_Puram_signature.csv"), index=False)
    pd.concat(per_all).to_csv(out("P1_within_patient_per_patient.csv"), index=False)

    # ---- P2 ----------------------------------------------------------------
    C.log("P2 증분 정보량 (ΔR²)")
    preds = {m: E.loc[m].values.astype(float) for m in C.MARKER_PANEL if m in E.index}
    rows = []
    for nm in ("PURAM30_minus_panel", "PURAM100_minus_panel"):
        r = {"response_signature": nm, "n_genes": len(sig[nm])}
        r.update(C.incremental_information(scores[nm], preds, w, pats))
        rows.append(r)
    # 상피-기질 축 (v2 와 동일 정의)
    try:
        from sklearn.decomposition import PCA
        ax = [g for g in ["EPCAM", "KRT14", "KRT17", "VIM", "PTPRC"] if g in E.index]
        pc1 = PCA(n_components=1, random_state=C.SEED).fit_transform(
            E.loc[ax].T.values.astype(float)).ravel()
        rr, pp = stats.spearmanr(w, pc1)
        rs, _ = stats.spearmanr(scores["PURAM30"], pc1)
        for r in rows:
            r.update({"rho_WBP5_vs_epi_stromal_PC1": rr, "p_WBP5_vs_PC1": pp,
                      "rho_PURAM30_vs_PC1": rs})
    except ImportError:
        pass
    p2 = pd.DataFrame(rows)
    p2.to_csv(out("P2_incremental_information.csv"), index=False)

    # ---- P3 ----------------------------------------------------------------
    C.log("P3 전사체 전체 상관 + Puram 농축 (수 분 소요)")
    Mt = E.loc[keep].drop(index=[gene], errors="ignore")
    names = Mt.index.values
    M = Mt.values.astype(np.float32)
    det = (M > 0).mean(axis=1)
    bins = C.expression_bins(M.mean(axis=1))
    res_rows, gene_tabs = [], []
    for mode, grp in (("pooled", None), ("within_patient", pats)):
        rho = C.corr_with_all(w, M, grp)
        pd.DataFrame({"gene": names, "rho": rho, "detection": det,
                      "in_PURAM100": np.isin(names, used100),
                      "in_PURAM30": np.isin(names, used30)}).sort_values(
            "rho", ascending=False).to_csv(out(f"P3_transcriptome_rho_{mode}.csv.gz"), index=False)
        for nm in ("PURAM30", "PURAM100"):
            r, pg = C.set_enrichment(rho, names, bins, sig[nm], rng, args.n_perm)
            r.update({"mode": mode, "signature": nm, "n_genes_tested": int(np.isfinite(rho).sum()),
                      "null_median_rho_all_genes": float(np.nanmedian(rho)),
                      "all_genes_95pct": float(np.nanpercentile(rho, 95)),
                      "all_genes_99pct": float(np.nanpercentile(rho, 99))})
            res_rows.append(r)
            pg["mode"], pg["signature"] = mode, nm
            gene_tabs.append(pg)
        top = names[np.argsort(-np.nan_to_num(rho, nan=-9))][:30]
        C.log(f"  [{mode}] 상위 30 상관 유전자: {', '.join(top)}")
    p3 = pd.DataFrame(res_rows)
    p3.to_csv(out("P3_Puram_enrichment.csv"), index=False)
    pd.concat(gene_tabs).to_csv(out("P3_Puram_gene_ranks.csv"), index=False)

    # ---- P4 ----------------------------------------------------------------
    C.log("P4 구획별 환자 수준 비교")
    nct = meta.get("noncancer_type")
    rows, rows_pp = [], []
    if nct is not None:
        nct = nct.astype(str).str.strip().str.lower()
        for label, key in (("Fibroblast", "fibroblast"), ("Endothelial", "endothelial")):
            comp = meta.index[(meta.is_cancer == 0) & nct.str.contains(key)]
            d, info = [], []
            for p in pd.unique(meta.patient):
                a = [c for c in mal if meta.at[c, "patient"] == p]
                b = [c for c in comp if meta.at[c, "patient"] == p]
                if len(a) >= 10 and len(b) >= 10:
                    va, vb = expr.loc[gene, a].astype(float), expr.loc[gene, b].astype(float)
                    d.append(va.mean() - vb.mean())
                    info.append({"comparison": f"Malignant_vs_{label}", "patient": p,
                                 "n_malignant": len(a), "n_comp": len(b),
                                 "mean_malignant": va.mean(), "mean_comp": vb.mean(),
                                 "det_malignant": (va > 0).mean(), "det_comp": (vb > 0).mean()})
            rows_pp.extend(info)
            d = np.array(d)
            r = {"comparison": f"Malignant_vs_{label}", "n_patients": len(d)}
            if len(d) >= 3:
                boot = [rng.choice(d, len(d)).mean() for _ in range(10000)]
                r.update({"mean_diff_log2": d.mean(), "ci_lo": np.percentile(boot, 2.5),
                          "ci_hi": np.percentile(boot, 97.5), "n_positive": int((d > 0).sum()),
                          "wilcoxon_p": stats.wilcoxon(d).pvalue if len(d) >= 5 else np.nan})
            rows.append(r)
    else:
        C.log("  noncancer_type 행 없음 → P4 생략")
    p4 = pd.DataFrame(rows)
    p4.to_csv(out("P4_compartment_patient_level.csv"), index=False)
    pd.DataFrame(rows_pp).to_csv(out("P4_compartment_per_patient.csv"), index=False)

    # ---- P5 ----------------------------------------------------------------
    p5 = None
    if args.yap and os.path.exists(args.yap):
        C.log("P5 YAP/TAZ 순열 (10,000회)")
        ylist = C.read_signature(args.yap)
        yused, ymap = C.resolve_genes(ylist, Mt.index)
        ymap.to_csv(out("P5_yap_gene_mapping.csv"), index=False)
        Z = C.zscore_rows(M)
        name_idx = {g: i for i, g in enumerate(names)}
        yidx = np.array([name_idx[g] for g in yused])
        rows = []
        for matched in (True, False):
            r = C.matched_null_setscore(w, Z, bins, yidx, rng, args.n_perm, matched)
            r["n_genes"] = len(yidx)
            rows.append(r)
        ys = np.nanmean(Z[yidx], axis=0)
        # 편상관 (순위 잔차): WBP5, YAP | PURAM30 , 그리고 | 원고 마커 패널 점수
        def partial(a, b, c):
            ra, rb, rc = (stats.rankdata(v) for v in (a, b, c))
            ea = ra - np.polyval(np.polyfit(rc, ra, 1), rc)
            eb = rb - np.polyval(np.polyfit(rc, rb, 1), rc)
            return stats.pearsonr(ea, eb)
        pr, pp = partial(w, ys, scores["PURAM30"])
        for r in rows:
            r.update({"partial_rho_given_PURAM30": pr, "partial_p_given_PURAM30": pp})
        p5 = pd.DataFrame(rows)
        p5.to_csv(out("P5_YAP_TAZ_permutation.csv"), index=False)
    else:
        C.log("P5 생략 (--yap 미제공)")

    # ---- 요약 출력 ----------------------------------------------------------
    pd.set_option("display.width", 200)
    f = lambda v: f"{v:.4g}"  # noqa
    print("\n=== P1 WBP5 ~ Puram p-EMT 서명 (악성세포) ===")
    print(p1[["signature", "n_genes", "rho_cell", "rho_ci_lo", "rho_ci_hi", "within_median_rho",
              "within_n_positive", "within_n_patients", "within_sign_test_p",
              "mixedlm_std_coef", "mixedlm_ci_lo", "mixedlm_ci_hi", "mixedlm_p"]]
          .to_string(index=False, float_format=f))
    print("\n=== P2 증분 정보량 (반응변수 = 마커 제외 Puram 점수) ===")
    for _, r in p2.iterrows():
        print(f"  [{r.response_signature}] R² 마커 {r.get('R2_markers', np.nan):.4f} → "
              f"+WBP5 {r.get('R2_markers_plus_WBP5', np.nan):.4f}  ΔR² {r.get('delta_R2_WBP5', np.nan):.4f}")
        print(f"     WBP5 표준화 부분계수 (혼합모형) {r.get('WBP5_std_partial_coef_mixedlm', np.nan):.4f} "
              f"[{r.get('WBP5_partial_ci_lo_mixedlm', np.nan):.4f}, {r.get('WBP5_partial_ci_hi_mixedlm', np.nan):.4f}] "
              f"p={r.get('WBP5_partial_p_mixedlm', np.nan):.3g}")
        print("     drop-one ΔR²: " + ", ".join(
            f"{k.replace('dropone_dR2_', '')} {v:.4f}" for k, v in r.items() if str(k).startswith("dropone_")))
    print("\n=== P3 Puram 유전자 농축 (WBP5 전사체 상관) ===")
    print(p3[["mode", "signature", "set_size_in_data", "mean_rho_set", "null_mean", "null_99pct",
              "perm_p_upper", "top100_overlap", "top100_expected", "top100_hypergeom_p",
              "top200_overlap", "top200_hypergeom_p"]].to_string(index=False, float_format=f))
    if len(p4):
        print("\n=== P4 구획별 환자 수준 (악성 - 비교군, log2(TPM/10+1)) ===")
        print(p4.to_string(index=False, float_format=f))
    if p5 is not None:
        print("\n=== P5 YAP/TAZ ===")
        print(p5.to_string(index=False, float_format=f))

    import scipy, statsmodels  # noqa
    with open(out("environment.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
                 f"scipy {scipy.__version__}\nstatsmodels {statsmodels.__version__}\n"
                 f"seed {C.SEED}\nn_perm {args.n_perm}\nn_boot {args.n_boot}\n"
                 f"pemt_file {os.path.abspath(args.pemt)}\nyap_file {args.yap}\n"
                 f"missing_pemt_symbols {miss}\n")
    C.log(f"완료 → {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
