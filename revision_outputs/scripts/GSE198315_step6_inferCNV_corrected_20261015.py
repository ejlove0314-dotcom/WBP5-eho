#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE198315 step6 — inferCNV 재실행 (정상 상피를 reference 에 포함)
================================================================
생성: 2026-10-15
대응: Reviewer #1 point 3 / Reviewer #5 Major comment 3

[왜 다시 돌리는가]
step3 은 면역·골수 계열만 reference 로 써서 inferCNV 를 돌렸습니다. 그 결과
상피세포의 99.6% 가 CNV 양성으로 나왔지만, step5 의 점검에서 정상조직(NT) 유래
상피세포도 99.7% 가 같은 기준을 넘었습니다(Cliff's delta = -0.058).

즉 그 점수는 '복제수 이상'이 아니라 '상피세포라는 계통 정체성'을 잡고 있었습니다.
상피세포는 면역세포와 유전체 전반에 걸쳐 발현 양상이 달라서, 면역세포만을
reference 로 삼으면 그 차이가 CNV 유사 신호로 나타납니다. 알려진 한계입니다.

이 스크립트는 **정상조직 유래 상피세포를 reference 에 추가**합니다. 그러면
상피 계통 성분이 상쇄되고, 종양 상피에서 남는 신호만이 실제 복제수 변화 후보가
됩니다. 악성 판정 임계도 면역세포가 아니라 **정상 상피의 분포**를 기준으로 잡습니다.

[주의 — 이 설계의 한계도 함께 보고해야 합니다]
HNSCC 의 'normal tissue' 는 보통 종양 인접 조직이라 field cancerization 으로
이미 변화가 있을 수 있고, 미량의 종양세포가 섞여 있을 수도 있습니다. 따라서
'정상 상피 대비' 판정은 '진짜 정상 대비' 판정보다 보수적입니다. 이 점을 답변서
Methods 에 명시하십시오.

실행:
    python GSE198315_step6_inferCNV_corrected_20261015.py \
        --h5ad gse198315_checkpoint_full_annotated.h5ad \
        --gtf gencode.v50.basic.annotation.gtf.gz \
        --outdir results_gse198315
"""

import argparse
import gzip
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

GENE_ALIASES = ["WBP5", "TCEAL9"]
MARKERS = ["LAMC2", "ITGA5", "KRT14", "KRT17", "VIM"]
IMMUNE = ["T/NK", "B/Plasma", "Myeloid", "Mast"]
NORMAL_REGION = "NT"
TUMOR_REGIONS = ["TC", "TP"]
SEED = 42


def log(m):
    sys.stderr.write(f"[{pd.Timestamp.now():%H:%M:%S}] {m}\n"); sys.stderr.flush()


def parse_gtf_gene_positions(gtf_path):
    """GTF 에서 gene 단위 좌표를 읽어 DataFrame(index=gene_name) 반환 (의존성 없음)."""
    op = gzip.open if str(gtf_path).endswith(".gz") else open
    recs = {}
    with op(gtf_path, "rt") as fh:
        for ln in fh:
            if not ln or ln[0] == "#":
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            i = p[8].find('gene_name "')
            if i < 0:
                continue
            name = p[8][i + 11:p[8].find('"', i + 11)]
            try:
                s, e = int(p[3]), int(p[4])
            except ValueError:
                continue
            prev = recs.get(name)
            if prev is None or (e - s) > (prev[2] - prev[1]):
                recs[name] = (p[0], s, e)
    if not recs:
        raise SystemExit("오류: GTF 에서 gene_name 을 가진 'gene' 행을 찾지 못했습니다.")
    return pd.DataFrame.from_dict(recs, orient="index",
                                  columns=["chromosome", "start", "end"])


def cliffs_delta(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return np.nan
    r = stats.rankdata(np.concatenate([a, b]))
    u1 = r[:n1].sum() - n1 * (n1 + 1) / 2
    return (2 * u1) / (n1 * n2) - 1


def dense(ad, gene):
    X = ad[:, gene].X
    return np.asarray(X.todense()).ravel() if hasattr(X, "todense") else np.asarray(X).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True,
                    help="step2 산출물 (…_full_annotated.h5ad). _cnv 가 아닌 쪽입니다")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--outdir", default="./results_gse198315")
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--cnv-sd", type=float, default=2.0,
                    help="정상 상피 분포의 평균 + N*SD 를 악성 임계로 사용")
    ap.add_argument("--min-mapped", type=int, default=5000)
    ap.add_argument("--min-ref-epi", type=int, default=200,
                    help="reference 로 쓸 정상조직 상피의 최소 세포 수")
    args = ap.parse_args()

    import scanpy as sc
    import infercnvpy as cnv
    os.makedirs(args.outdir, exist_ok=True)

    log(f"적재: {args.h5ad}")
    a = sc.read_h5ad(args.h5ad)
    for c in ["cell_type", "patient", "region"]:
        if c not in a.obs.columns:
            raise SystemExit(f"오류: obs 에 '{c}' 가 없습니다. step2 산출물이 맞는지 확인하십시오.")

    # dtype 정규화 (anndata 0.13+ 의 arrow 문자열이 infercnvpy 내부에서 실패)
    a.var.index = pd.Index(np.asarray(a.var_names, dtype=object), dtype=object)
    a.obs.index = pd.Index(np.asarray(a.obs_names, dtype=object), dtype=object)
    for c in a.var.columns:
        if str(a.var[c].dtype) in ("str", "string", "string[pyarrow]"):
            a.var[c] = a.var[c].astype(object)
    for c in a.obs.columns:
        if str(a.obs[c].dtype) in ("str", "string", "string[pyarrow]"):
            a.obs[c] = a.obs[c].astype(object)

    # ---- reference 그룹 정의: 면역 + 정상조직 상피 -------------------------
    ct = a.obs["cell_type"].astype(str)
    rg = a.obs["region"].astype(str)
    grp = ct.copy()
    is_norm_epi = (ct == "Epithelial") & (rg == NORMAL_REGION)
    grp[is_norm_epi] = "Epithelial_normal"
    a.obs["cnv_ref_group"] = grp.values

    n_norm = int(is_norm_epi.sum())
    log(f"  정상조직 상피 {n_norm:,} 개를 reference 에 추가")
    if n_norm < args.min_ref_epi:
        raise SystemExit(
            f"오류: 정상조직 상피가 {n_norm}개뿐입니다. 이 설계로는 보정이 불가능합니다.")

    ref = [c for c in IMMUNE if c in set(a.obs["cnv_ref_group"])] + ["Epithelial_normal"]
    log(f"  reference 구획: {ref}")

    # ---- 유전자 좌표 -------------------------------------------------------
    log(f"GTF 좌표 읽는 중: {args.gtf}")
    pos = parse_gtf_gene_positions(args.gtf)
    a.var["chromosome"] = a.var_names.map(pos["chromosome"]).values
    a.var["start"] = pd.to_numeric(a.var_names.map(pos["start"]), errors="coerce").values
    a.var["end"] = pd.to_numeric(a.var_names.map(pos["end"]), errors="coerce").values
    ok = a.var["chromosome"].notna()
    frac = float(ok.sum()) / max(a.n_vars, 1)
    log(f"  좌표 부여 {int(ok.sum()):,} / {a.n_vars:,} ({frac:.1%})")
    if frac < 0.30 or ok.sum() < args.min_mapped:
        raise SystemExit("오류: 좌표 매핑이 부족합니다.")

    # ---- inferCNV ----------------------------------------------------------
    log("inferCNV 실행 중 (정상 상피 포함 reference)")
    cnv.tl.infercnv(a, reference_key="cnv_ref_group", reference_cat=ref,
                    window_size=args.window)
    cnv.tl.pca(a); cnv.pp.neighbors(a); cnv.tl.leiden(a); cnv.tl.cnv_score(a)

    s = a.obs["cnv_score"]
    norm_epi = s[is_norm_epi.values]
    thr = norm_epi.mean() + args.cnv_sd * norm_epi.std()
    a.obs["cnv_malignant"] = (s > thr).values
    log(f"  임계 {thr:.5f} (정상 상피 평균 + {args.cnv_sd}SD)")

    # ---- 타당성 점검 (step5 (A) 와 동일 설계) -------------------------------
    ep = a.obs[ct.values == "Epithelial"]
    nt = ep.loc[ep.region == NORMAL_REGION, "cnv_score"]
    tu = ep.loc[ep.region.isin(TUMOR_REGIONS), "cnv_score"]
    print("\n=== 타당성 점검: 상피 CNV score, 정상조직 vs 종양 ===")
    print(f"  NT    n={len(nt):>6}  mean {nt.mean():.5f}  median {nt.median():.5f}  "
          f"악성판정 {100*ep.loc[ep.region==NORMAL_REGION,'cnv_malignant'].mean():.1f}%")
    print(f"  TC+TP n={len(tu):>6}  mean {tu.mean():.5f}  median {tu.median():.5f}  "
          f"악성판정 {100*ep.loc[ep.region.isin(TUMOR_REGIONS),'cnv_malignant'].mean():.1f}%")
    d = cliffs_delta(tu.values, nt.values)
    p = stats.mannwhitneyu(tu.values, nt.values, alternative="two-sided").pvalue
    print(f"  Cliff's delta (종양 - 정상) = {d:+.4f}, MWU p = {p:.3g}")
    pct_nt_mal = 100 * ep.loc[ep.region == NORMAL_REGION, "cnv_malignant"].mean()
    if d > 0.3 and pct_nt_mal < 30:
        verdict = "OK: normal epithelium now separates from tumour epithelium"
    elif d > 0.15:
        verdict = "PARTIAL: separation improved but incomplete; state as a limitation"
    else:
        verdict = "STILL CONFOUNDED: do not use inferCNV as an independent malignancy call"
    print(f"  >>> 판정: {verdict}")

    # ---- 산출물 ------------------------------------------------------------
    pd.crosstab(a.obs["cell_type"], a.obs["cnv_malignant"]).to_csv(
        f"{args.outdir}/L1_celltype_vs_inferCNV_corrected.csv")
    (a.obs.groupby(["cell_type", "region"], observed=True)
       .agg(n_cells=("cnv_score", "size"),
            mean_cnv_score=("cnv_score", "mean"),
            pct_malignant=("cnv_malignant", lambda x: x.mean()*100),
            n_patients=("patient", "nunique")).reset_index()
       .to_csv(f"{args.outdir}/L2_cnv_by_celltype_region_corrected.csv", index=False))
    a.obs[["cell_type", "region", "patient", "cnv_score", "cnv_malignant"]].to_csv(
        f"{args.outdir}/L3_percell_cnv_corrected.csv")

    # ---- 보정된 악성 판정으로 WBP5 구획 대조 재계산 -------------------------
    gene = next((g for g in GENE_ALIASES if g in a.var_names), None)
    if gene:
        v = pd.Series(dense(a, gene), index=a.obs_names)
        mal = v[a.obs["cnv_malignant"].values]
        rows = []
        for c in sorted(set(ct)):
            other = v[(ct.values == c) & (~a.obs["cnv_malignant"].values)]
            if len(other) >= 20 and len(mal) >= 20:
                rows.append({"contrast": f"CNV+ malignant vs {c} (CNV-)",
                             "n_malignant": len(mal), "n_other": len(other),
                             "mean_malignant": mal.mean(), "mean_other": other.mean(),
                             "cliffs_delta": cliffs_delta(mal.values, other.values),
                             "MWU_p": stats.mannwhitneyu(mal.values, other.values,
                                                         alternative="two-sided").pvalue})
        res = pd.DataFrame(rows)
        res.to_csv(f"{args.outdir}/L4_compartment_contrasts_corrected.csv", index=False)
        print("\n=== 보정된 악성 판정 기준 구획 대조 ===")
        print(res.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

        # 환자 내 상관도 보정된 악성세포로 재계산
        mal_mask = a.obs["cnv_malignant"].values.astype(bool)
        pats = a.obs["patient"].values[mal_mask]
        vm = v.values[mal_mask]
        rows = []
        for m in MARKERS:
            if m not in a.var_names:
                continue
            mv = dense(a, m)[mal_mask]
            for pt in pd.unique(pats):
                sel = pats == pt
                if sel.sum() < 20 or np.std(vm[sel]) == 0 or np.std(mv[sel]) == 0:
                    continue
                rho, _ = stats.spearmanr(vm[sel], mv[sel])
                rows.append({"marker": m, "patient": pt, "n_cells": int(sel.sum()), "rho": rho})
        per = pd.DataFrame(rows)
        if len(per):
            from statsmodels.stats.multitest import multipletests
            summ = []
            for m, g in per.groupby("marker"):
                rh = g.rho.dropna().values
                npos = int((rh > 0).sum())
                summ.append({"marker": m, "median_rho": np.median(rh),
                             "n_positive": npos, "n_patients": len(rh),
                             "sign_test_p": stats.binomtest(npos, len(rh), 0.5).pvalue})
            summ = pd.DataFrame(summ)
            summ["sign_test_q"] = multipletests(summ.sign_test_p, method="fdr_bh")[1]
            summ.to_csv(f"{args.outdir}/L5_within_patient_rho_corrected.csv", index=False)
            per.to_csv(f"{args.outdir}/L5b_within_patient_rho_per_patient_corrected.csv",
                       index=False)
            print("\n=== 보정된 악성세포 기준 환자 내 Spearman ===")
            print(summ.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    a.write_h5ad(args.h5ad.replace(".h5ad", "_cnv_corrected.h5ad"), compression="gzip")

    with open(f"{args.outdir}/environment_step6.txt", "w", encoding="utf-8") as fh:
        import scipy, statsmodels
        fh.write(f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n")
        fh.write(f"scipy {scipy.__version__}\nscanpy {sc.__version__}\n")
        fh.write(f"infercnvpy {cnv.__version__}\nstatsmodels {statsmodels.__version__}\n")
        fh.write(f"seed {SEED}\nwindow {args.window}\ncnv_sd {args.cnv_sd}\n")
        fh.write(f"reference_groups {ref}\nn_normal_epithelium_in_reference {n_norm}\n")
        fh.write(f"threshold {thr}\nvalidity_verdict {verdict}\n")
        fh.write(f"gtf {os.path.abspath(args.gtf)}\ninput {os.path.abspath(args.h5ad)}\n")
    log("완료")


if __name__ == "__main__":
    main()
