#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CRED_crosscohort_20261021.py — 신뢰도 분석 1번: 교차 코호트 유전자 이름 중복
=========================================================================
사전 판정 규칙: claude/CSBJ_신뢰도분석4종_사전규칙_20261021.md (1번)

입력: CRED_profiles_20261021.py 가 만든 CR1_profile_WBP5.csv.gz
      GSE103322 1개 + GSE198315 서브샘플 3개(seed 42/2026/7)

주 결과
  GSE103322 보정 프로파일 상위 N 그리고 GSE198315 세 seed 모두의 상위 N 에 드는 Puram 100 유전자 이름.
  N = 200 (주), 100·500 (민감도). 상위 N 은 각 코호트 자체 유전자 집합에서 정한다(D4·P3 와 같은 기준).
유의성 (두 코호트 공통 유전자 범위 U 안에서)
  (a) Puram 유전자 중복: U 안의 Puram 유전자 K개 중, 코호트 1 목록에 든 a개와 코호트 2 목록에 든 b개가
      겹치는 수의 초기하 검정 → "같은 p-EMT 유전자가 반복 회수되는가"
  (b) 전체 유전자 상위 N 중복의 초기하 검정, 공통 유전자 전체 프로파일의 코호트 간 Spearman

실행:
  python CRED_crosscohort_20261021.py --g103322 results_cred\\gse103322 ^
     --g198315 results_cred\\seed42 results_cred\\seed2026 results_cred\\seed7 ^
     --outdir results_cred\\cross
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wbp5_rev_common as C  # noqa: E402

TOP_NS = (200, 100, 500)
# 2017 기호 → 현재 기호 (공통 모듈 ALIASES + 이 분석에서 필요한 것)
CANON = dict(C.ALIASES)
CANON.update({"NGFRAP1": "BEX3", "C15orf24": "EMC7"})


def load(d, col):
    p = pd.read_csv(os.path.join(d, "CR1_profile_WBP5.csv.gz"))
    p["canon"] = [CANON.get(g, g) for g in p.gene]
    if p.canon.duplicated().any():
        p = p.sort_values(col, ascending=False).drop_duplicates("canon")
    p = p.sort_values(col, ascending=False).reset_index(drop=True)
    p["rank"] = np.arange(1, len(p) + 1)
    return p.set_index("canon")


def hyper_overlap(N, A, B):
    """N 개 중 크기 |A|, |B| 인 두 무작위 부분집합의 중복이 관측값 이상일 확률."""
    x = len(A & B)
    return x, len(A) * len(B) / N if N else np.nan, stats.hypergeom.sf(x - 1, N, len(A), len(B))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g103322", required=True)
    ap.add_argument("--g198315", required=True, nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--col", default="rho_adj", help="rho_adj(주) 또는 rho_within_unadj")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda f: os.path.join(args.outdir, f)  # noqa

    A = load(args.g103322, args.col)
    Bs = {os.path.basename(os.path.normpath(d)): load(d, args.col) for d in args.g198315}
    U = set(A.index)
    for b in Bs.values():
        U &= set(b.index)
    pur_A = set(A.index[A.in_PURAM100])
    pur_B = set.intersection(*[set(b.index[b.in_PURAM100]) for b in Bs.values()])
    K = sorted(pur_A & pur_B & U)
    C.log(f"공통 유전자 {len(U):,} / 공통 Puram 100 {len(K)}개 / GSE198315 서브샘플 {list(Bs)}")

    summ, named = [], []
    for N in TOP_NS:
        topA = set(A.index[A["rank"] <= N]) & U
        topB_each = {k: set(b.index[b["rank"] <= N]) & U for k, b in Bs.items()}
        topB_all = set.intersection(*topB_each.values())
        topB_2of3 = {g for g in U if sum(g in s for s in topB_each.values()) >= 2}
        for lab, topB in (("all_seeds", topB_all), ("ge2_seeds", topB_2of3)):
            a, b = topA & set(K), topB & set(K)
            x, ex, p = hyper_overlap(len(K), a, b)
            xg, exg, pg = hyper_overlap(len(U), topA, topB)
            summ.append({"top_N": N, "GSE198315_rule": lab, "n_shared_universe": len(U), "K_puram_shared": len(K),
                         "puram_in_GSE103322_top": len(a), "puram_in_GSE198315_top": len(b),
                         "puram_overlap": x, "puram_overlap_expected": ex, "puram_overlap_hypergeom_p": p,
                         "allgene_GSE103322_top": len(topA), "allgene_GSE198315_top": len(topB),
                         "allgene_overlap": xg, "allgene_overlap_expected": exg, "allgene_overlap_hypergeom_p": pg,
                         "shared_puram_genes": ";".join(sorted(a & b)),
                         "shared_top_genes": ";".join(sorted(topA & topB))})
            if N == 200 and lab == "all_seeds":
                for g in sorted(a & b, key=lambda g: A.at[g, "rank"]):
                    row = {"gene_current_symbol": g, "GSE103322_symbol": A.at[g, "gene"],
                           "GSE103322_rank": int(A.at[g, "rank"]), "GSE103322_rho": A.at[g, args.col]}
                    for k, bb in Bs.items():
                        row[f"{k}_rank"] = int(bb.at[g, "rank"])
                        row[f"{k}_rho"] = bb.at[g, args.col]
                    named.append(row)
    summ = pd.DataFrame(summ)
    named = pd.DataFrame(named)
    summ.to_csv(out("CR1x_crosscohort_overlap.csv"), index=False)
    named.to_csv(out("CR1x_shared_puram_genes_top200.csv"), index=False)

    conc = []
    Ul = sorted(U)
    for k, b in Bs.items():
        r, p = stats.spearmanr(A.loc[Ul, args.col], b.loc[Ul, args.col])
        rp, pp = stats.spearmanr(A.loc[K, args.col], b.loc[K, args.col]) if len(K) > 2 else (np.nan, np.nan)
        conc.append({"GSE198315": k, "spearman_all_shared_genes": r, "p_all": p, "n_all": len(Ul),
                     "spearman_puram_genes": rp, "p_puram": pp, "n_puram": len(K)})
    seeds = list(Bs)
    for i in range(len(seeds)):
        for j in range(i + 1, len(seeds)):
            r, _ = stats.spearmanr(Bs[seeds[i]].loc[Ul, args.col], Bs[seeds[j]].loc[Ul, args.col])
            conc.append({"GSE198315": f"{seeds[i]} vs {seeds[j]}", "spearman_all_shared_genes": r, "n_all": len(Ul)})
    conc = pd.DataFrame(conc)
    conc.to_csv(out("CR1x_profile_concordance.csv"), index=False)

    pd.set_option("display.width", 240)
    f = lambda v: f"{v:.4g}"  # noqa
    print(f"\n######## CR1x 교차 코호트 ({args.col}) ########")
    print(summ.drop(columns=["shared_puram_genes", "shared_top_genes"]).to_string(index=False, float_format=f))
    print("\n=== 주 결과: 상위 200, GSE198315 세 seed 모두 — 공통 Puram 유전자 ===")
    print(named.to_string(index=False, float_format=f) if len(named) else "  없음")
    main_row = summ[(summ.top_N == 200) & (summ.GSE198315_rule == "all_seeds")].iloc[0]
    verdict = "반복 회수됨 (p < 0.05)" if main_row.puram_overlap_hypergeom_p < 0.05 else "반복성 주장 불가 (이름만 보고)"
    print(f"  사전 규칙 판정: {verdict}")
    print("\n=== 프로파일 일치도 (공통 유전자) ===")
    print(conc.to_string(index=False, float_format=f))
    C.log(f"완료 → {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
