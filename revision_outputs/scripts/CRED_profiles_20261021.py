#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CRED_profiles_20261021.py — 신뢰도 분석 4종 중 2·3·4번 (+ 1번 입력 파일 생성)
==========================================================================
사전 판정 규칙: 프로젝트 문서 claude/CSBJ_신뢰도분석4종_사전규칙_20261021.md
(실행 전 기록. 이 스크립트는 그 규칙을 그대로 구현하며 결과를 보고 기준을 바꾸지 않는다.)

기존 분석과의 연속성
  GSE198315: step7b(D3/D4)와 같은 세포(Epithelial ∩ TC/TP), 같은 유전자 집합(검출률 ≥ 1%),
             같은 깊이 지표(log1p total_counts), 같은 환자 내 순위 z, 같은 편상관 식(분모 n).
  GSE103322: P3와 같은 세포(악성, combo 제외), 같은 유전자 집합(분산 > 0),
             같은 상관 정의(C.corr_with_all, 분모 = 제곱합). 기술 축 = log1p(검출 유전자 수).
  → CR0 에서 기존 결과 파일(--check-dir)을 다시 계산해 일치하는지 먼저 확인한다.
    불일치하면 "재현 실패"를 출력하고, 이후 결과는 응답서에 쓰지 않는다.

분석
  CR0  재현 점검 (D3·D4 또는 P3)
  CR1  WBP5 프로파일 저장 (교차 코호트 분석의 입력)
  CR2  유전자 단위 농축 leave-one-patient-out (환자 ≥ 20세포)
  CR3  TCEAL/BEX/WBP 계열 특이성 — 전체 유전자를 각각 표적으로 삼아 같은 지표를 계산하고,
       WBP5와 발현이 가장 가까운 기준 유전자 1,000개 대비 백분위로 비교
  CR3b Xq22 모듈 보정 후 WBP5의 Puram 농축 (기술 축 + 모듈 점수 동시 보정)
  CR4  음성 대조(하우스키핑 11개)·양성 대조(LAMC2, ITGA5, Puram 100 구성원)

실행 (Windows cmd, WBP5 revision 폴더):
  python CRED_profiles_20261021.py --cohort gse198315 ^
     --h5ad gse198315_checkpoint_full_annotated_cnv.h5ad --pemt "p-EMT gene list.xlsx" ^
     --tag seed42 --check-dir results_gse198315 --outdir results_cred\\seed42
  python CRED_profiles_20261021.py --cohort gse103322 ^
     --matrix GSE103322_HNSCC_all_data.txt.gz --pemt "p-EMT gene list.xlsx" ^
     --tag gse103322 --check-dir results_gse103322_puram --outdir results_cred\\gse103322

필요: wbp5_rev_common.py, (gse103322) GSE103322_Puram_pEMT_20261019.py — 같은 폴더.
"""
import argparse
import hashlib
import importlib
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import wbp5_rev_common as C
except ImportError:
    raise SystemExit("오류: wbp5_rev_common.py 를 이 스크립트와 같은 폴더에 두십시오.")

warnings.filterwarnings("ignore")

TUMOR_REGIONS = ["TC", "TP"]
MIN_CELLS = 20          # LOPO 대상 환자 기준 (step7b 와 동일)
EVAL_DET = 0.10         # 계열·대조 유전자 평가 가능 기준 (사전 규칙)
TOP_N = (100, 200, 500)
CHUNK = 256

# ---- 사전 지정 유전자 (사전 규칙 문서와 동일) -------------------------------------
FAMILY_CORE = [f"TCEAL{i}" for i in range(1, 9)]                 # TCEAL1–8
FAMILY_EXT = ["BEX1", "BEX2", "BEX3", "BEX4", "BEX5"]           # Xq22 인접 (BEX3 = NGFRAP1)
OUTGROUP = ["WBP1", "WBP2", "WBP4", "WBP11"]
NEG_CTRL = ["C1orf43", "CHMP2A", "EMC7", "GPI", "PSMB2", "PSMB4",
            "RAB7A", "REEP5", "SNRPD3", "VCP", "VPS29"]          # Eisenberg & Levanon 2013
POS_CTRL = ["LAMC2", "ITGA5"]
# 이 스크립트에서만 쓰는 기호 대응 (공통 모듈의 ALIASES 는 건드리지 않음)
EXTRA_ALIASES = {"BEX3": ["NGFRAP1"], "EMC7": ["C15orf24"], "C1orf43": ["C1ORF43"]}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def resolve_symbol(sym, available):
    """사전 지정 기호를 데이터의 실제 행 이름으로 대응. 없으면 None."""
    avail = set(available)
    for cand in [sym] + EXTRA_ALIASES.get(sym, []):
        if cand in avail:
            return cand
    if sym in C.ALIASES and C.ALIASES[sym] in avail:
        return C.ALIASES[sym]
    if sym in C.ALIASES_REV and C.ALIASES_REV[sym] in avail:
        return C.ALIASES_REV[sym]
    return None


# ---------------------------------------------------------------------------
# 순위 z (step7b 와 동일 코드)
# ---------------------------------------------------------------------------

def zrank_within(v, groups):
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


def matvec64(R, v, chunk=4096):
    """R(float32) @ v 를 float64 로, 메모리 사본 없이."""
    out = np.empty(R.shape[0], dtype=np.float64)
    for i in range(0, R.shape[0], chunk):
        out[i:i + chunk] = R[i:i + chunk].astype(np.float64) @ v
    return out


# ---------------------------------------------------------------------------
# 데이터 적재
# ---------------------------------------------------------------------------

def load_198315(args):
    import anndata as ad
    import scipy.sparse as sp
    C.log(f"적재: {args.h5ad}")
    a = ad.read_h5ad(args.h5ad)
    for c in ("cell_type", "region", "patient"):
        a.obs[c] = a.obs[c].astype(str)
    xs = a.X[:2000]
    xs = xs.toarray() if sp.issparse(xs) else np.asarray(xs)
    looks_raw = bool(xs.max() > 50 and np.allclose(xs, np.round(xs)))
    C.log(f"  X 상태 점검: max {xs.max():.2f}, 정수형 {np.allclose(xs, np.round(xs))} → "
          f"{'원시 count 로 보임 (step7b 와 동일하게 변환하지 않음)' if looks_raw else '정규화 값'}")
    gene = C.resolve_one(C.WBP5_ALIASES, a.var_names)
    m = (a.obs.cell_type.values == "Epithelial") & np.isin(a.obs.region.values, TUMOR_REGIONS)
    e = a[m].copy()
    del a
    pats = e.obs.patient.values
    X = e.X.tocsc() if sp.issparse(e.X) else sp.csc_matrix(e.X)
    src = None
    for c in ("total_counts", "n_counts", "nCount_RNA"):
        if c in e.obs:
            tech, src = np.log1p(e.obs[c].astype(float).values), f"log1p(obs.{c})"
            break
    if src is None:
        tech = np.asarray((X > 0).sum(axis=1)).ravel().astype(float)
        src = "검출 유전자 수 (X > 0)"
    w = np.asarray(X[:, e.var_names.get_loc(gene)].todense()).ravel()
    det_all = np.asarray((X > 0).mean(axis=0)).ravel()
    keep = np.where((det_all >= args.min_det) & (e.var_names.values != gene))[0]
    names = e.var_names.values[keep].astype(str)
    M = X[:, keep].toarray().T.astype(np.float32)
    var_names = e.var_names.values.astype(str)
    det_map = dict(zip(var_names, det_all))
    C.log(f"  EPI_TUMOR {e.n_obs:,} 세포 / 환자 {len(pd.unique(pats))} / 유전자 집합 {len(names):,} "
          f"(검출률 ≥ {args.min_det:.0%}) / 기술 축 {src}")
    return dict(w=w, pats=pats, tech=tech, tech_src=src, names=names, M=M, gene=gene,
                var_names=var_names, det_map=det_map, mode="n", looks_raw=looks_raw)


def load_103322(args):
    here = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(os.path.join(here, "GSE103322_Puram_pEMT_20261019.py")):
        raise SystemExit("오류: GSE103322_Puram_pEMT_20261019.py 를 같은 폴더에 두십시오 (적재 함수 재사용).")
    P = importlib.import_module("GSE103322_Puram_pEMT_20261019")
    expr, meta = P.load_gse103322(args.matrix)
    gene = C.resolve_one(C.WBP5_ALIASES, expr.index)
    mal = meta.index[meta.is_cancer == 1]
    E = expr[mal]
    pats = meta.loc[mal, "patient"].values
    w = E.loc[gene].values.astype(float)
    ndet = (E.values > 0).sum(axis=0).astype(float)
    tech, src = np.log1p(ndet), "log1p(검출 유전자 수, 전 유전자)"
    keep = E.index[E.std(axis=1) > 0]                      # P3 와 동일
    Mt = E.loc[keep].drop(index=[gene], errors="ignore")
    names = Mt.index.values.astype(str)
    M = Mt.values.astype(np.float32)
    det_all = (E.values > 0).mean(axis=1)
    var_names = E.index.values.astype(str)
    det_map = dict(zip(var_names, det_all))
    C.log(f"  악성세포 {E.shape[1]:,} / 환자 {len(pd.unique(pats))} / 유전자 집합 {len(names):,} (분산 > 0) "
          f"/ 기술 축 {src}")
    return dict(w=w, pats=pats, tech=tech, tech_src=src, names=names, M=M, gene=gene,
                var_names=var_names, det_map=det_map, mode="ss", looks_raw=False)


# ---------------------------------------------------------------------------
# 환자별 분해 — LOPO 와 전체 값을 같은 식으로 계산하기 위함
# ---------------------------------------------------------------------------

class Decomp:
    """R(유전자 × 세포, 환자 내 순위 z)을 환자별로 쪼개 둔다."""

    def __init__(self, R, M, pats, tw, td, mode):
        self.mode = mode
        self.pids = list(pd.unique(pats))
        self.ncell = np.array([(pats == p).sum() for p in self.pids])
        G, P = R.shape[0], len(self.pids)
        self.NUMw = np.zeros((G, P))
        self.NUMd = np.zeros((G, P))
        self.SSR = np.zeros((G, P))
        self.SUM = np.zeros((G, P))         # 발현 합 (LOPO 발현구간 재계산용)
        self.ssw = np.zeros(P)
        self.ssd = np.zeros(P)
        self.wd = np.zeros(P)
        for j, p in enumerate(self.pids):
            s = np.where(pats == p)[0]
            Rs = R[:, s].astype(np.float64)
            self.NUMw[:, j] = Rs @ tw[s]
            self.NUMd[:, j] = Rs @ td[s]
            self.SSR[:, j] = (Rs ** 2).sum(axis=1)
            self.SUM[:, j] = M[:, s].astype(np.float64).sum(axis=1)
            self.ssw[j] = (tw[s] ** 2).sum()
            self.ssd[j] = (td[s] ** 2).sum()
            self.wd[j] = (tw[s] * td[s]).sum()

    def profile(self, keep_mask):
        k = np.asarray(keep_mask, bool)
        nw, nd = self.NUMw[:, k].sum(1), self.NUMd[:, k].sum(1)
        wd = self.wd[k].sum()
        if self.mode == "n":
            n = self.ncell[k].sum()
            r_gw, r_gd, r_wd = nw / n, nd / n, wd / n
        else:
            ss = self.SSR[:, k].sum(1)
            sw, sd = self.ssw[k].sum(), self.ssd[k].sum()
            with np.errstate(invalid="ignore", divide="ignore"):
                r_gw = nw / np.sqrt(ss * sw)
                r_gd = nd / np.sqrt(ss * sd)
            r_wd = wd / np.sqrt(sw * sd)
        mean_expr = self.SUM[:, k].sum(1) / self.ncell[k].sum()
        return r_gw, r_gd, r_wd, partial(r_gw, r_gd, r_wd), mean_expr


# ---------------------------------------------------------------------------
# 전체 유전자를 표적으로 한 Puram 농축 지표 (해석적 발현구간 맞춤 귀무)
# ---------------------------------------------------------------------------

def target_scan(R, mode, td_vec, r_gd, bins, is_puram, targets, chunk=CHUNK, topn=200):
    """
    targets: 유전자 집합 안의 행 번호 배열. 각 표적 j 에 대해
      P_gj = 표적 j 의 기술 축 편상관 프로파일 (자기 자신 제외)
      E1   = Puram 유전자의 평균 P (자기 자신 제외)
      mu, sd = 같은 발현구간의 비-Puram 유전자에서 무작위로 뽑을 때의 평균·표준편차 (해석적)
      z    = (E1 - mu) / sd,  top = 상위 topn 중 Puram 수
    """
    G, n = R.shape
    ssR = (R.astype(np.float64) ** 2).sum(axis=1)
    nb = int(bins.max()) + 1
    pool = [np.where((bins == b) & ~is_puram)[0] for b in range(nb)]
    pur_idx = np.where(is_puram)[0]
    pur_bin_counts = np.bincount(bins[pur_idx], minlength=nb).astype(float)
    r_td_all = r_gd                                     # 표적 j 의 (j, 기술 축) 상관 = r_gd[j]
    res = {k: np.full(len(targets), np.nan) for k in ("E1", "mu", "sd", "z", "top")}
    for c0 in range(0, len(targets), chunk):
        J = targets[c0:c0 + chunk]
        k = len(J)
        Cm = (R @ R[J].T).astype(np.float64)            # G × k
        if mode == "n":
            Cm /= n
        else:
            Cm /= np.sqrt(np.outer(ssR, ssR[J]))
        rj = r_td_all[J]
        with np.errstate(invalid="ignore", divide="ignore"):
            Pm = (Cm - np.outer(r_gd, rj)) / np.sqrt(np.outer(1 - r_gd ** 2, 1 - rj ** 2))
        del Cm
        Pm = np.nan_to_num(Pm, nan=0.0, posinf=0.0, neginf=0.0)   # 환자 내 분산 0 인 행 (드묾)
        cols = np.arange(k)
        Pm[J, cols] = 0.0                               # 자기 자신: 합계에서 빠지도록 0, 아래에서 개수 보정
        self_pur = is_puram[J].astype(float)
        K = len(pur_idx) - self_pur
        E1 = Pm[pur_idx].sum(0) / K
        S1 = np.vstack([Pm[pool[b]].sum(0) for b in range(nb)])        # nb × k
        S2 = np.vstack([(Pm[pool[b]] ** 2).sum(0) for b in range(nb)])
        cnt = np.array([len(pool[b]) for b in range(nb)], float)[:, None] * np.ones((1, k))
        own_bin = bins[J]
        nonpur = ~is_puram[J]
        cnt[own_bin[nonpur], cols[nonpur]] -= 1         # 표적이 비-Puram 이면 자기 구간 풀에서 제외
        with np.errstate(invalid="ignore", divide="ignore"):
            mu_b = S1 / cnt
            var_b = S2 / cnt - mu_b ** 2
        pc = pur_bin_counts[:, None] * np.ones((1, k))
        pc[own_bin[~nonpur], cols[~nonpur]] -= 1        # 표적이 Puram 이면 자기 자신을 서명에서 제외
        mu = np.nansum(pc * mu_b, 0) / K
        sd = np.sqrt(np.nansum(pc * var_b, 0)) / K
        Pm[J, cols] = -np.inf
        top_idx = np.argpartition(-Pm, topn, axis=0)[:topn]
        top = is_puram[top_idx].sum(0)
        res["E1"][c0:c0 + k], res["mu"][c0:c0 + k], res["sd"][c0:c0 + k] = E1, mu, sd
        res["z"][c0:c0 + k] = (E1 - mu) / sd
        res["top"][c0:c0 + k] = top
        if (c0 // chunk) % 10 == 0:
            C.log(f"    표적 {c0 + k:,}/{len(targets):,}")
    return res


def single_target_stats(P, bins, is_puram, topn=200):
    """WBP5 처럼 유전자 집합 밖의 표적 하나 (프로파일 P 가 이미 계산됨)."""
    nb = int(bins.max()) + 1
    pur_idx = np.where(is_puram)[0]
    K = len(pur_idx)
    E1 = float(np.nanmean(P[pur_idx]))
    mu, var = 0.0, 0.0
    for b in range(nb):
        pool = np.where((bins == b) & ~is_puram)[0]
        nbp = int(is_puram[bins == b].sum())
        if nbp == 0:
            continue
        v = P[pool]
        mu += nbp * np.nanmean(v)
        var += nbp * np.nanvar(v)
    mu, sd = mu / K, np.sqrt(var) / K
    order = np.argsort(-np.nan_to_num(P, nan=-np.inf))
    return {"E1": E1, "mu": mu, "sd": sd, "z": (E1 - mu) / sd, "top": int(is_puram[order[:topn]].sum())}


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=["gse198315", "gse103322"])
    ap.add_argument("--h5ad")
    ap.add_argument("--matrix")
    ap.add_argument("--pemt", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--check-dir", default=None, help="기존 결과 폴더 (D3/D4 또는 P3 파일)")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-ref", type=int, default=1000)
    ap.add_argument("--min-det", type=float, default=0.01)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda f: os.path.join(args.outdir, f)  # noqa

    me = os.path.abspath(__file__)
    with open(out("CR_run_info.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"script {os.path.basename(me)} sha256 {sha256(me)}\n"
                 f"common sha256 {sha256(C.__file__)}\n"
                 f"started {pd.Timestamp.now()}\ncohort {args.cohort}\ntag {args.tag}\n"
                 f"input {args.h5ad or args.matrix}\npemt {os.path.abspath(args.pemt)}\n"
                 f"n_perm {args.n_perm}\nn_ref {args.n_ref}\nmin_det {args.min_det}\nseed {C.SEED}\n"
                 "preregistration claude/CSBJ_신뢰도분석4종_사전규칙_20261021.md\n")

    D = load_198315(args) if args.cohort == "gse198315" else load_103322(args)
    w, pats, tech, names, M, mode = D["w"], D["pats"], D["tech"], D["names"], D["M"], D["mode"]
    name_idx = {g: i for i, g in enumerate(names)}

    plist = C.read_signature(args.pemt)
    used100, _ = C.resolve_genes(plist, D["var_names"])
    used30, _ = C.resolve_genes(plist[:30], D["var_names"])
    pres100 = [g for g in used100 if g in name_idx]
    pres30 = [g for g in used30 if g in name_idx]
    is_p100 = np.isin(names, pres100)

    # ---- 순위 행렬 · 전체 프로파일 ------------------------------------------------
    C.log("순위 행렬 계산")
    mean_expr = M.mean(axis=1)
    bins = C.expression_bins(mean_expr)
    det = (M > 0).mean(axis=1)
    R = zrank_within_rows(M, pats)
    tw, td = zrank_within(w, pats), zrank_within(tech, pats)
    n = len(tw)
    if mode == "n":
        r_gw = matvec64(R, tw) / n
        r_gd = matvec64(R, td) / n
        r_wd = float(np.mean(tw * td))
        r_unadj = r_gw
    else:
        ssR = (R.astype(np.float64) ** 2).sum(axis=1)
        r_gw = matvec64(R, tw) / np.sqrt(ssR * (tw ** 2).sum())
        r_gd = matvec64(R, td) / np.sqrt(ssR * (td ** 2).sum())
        r_wd = float((tw * td).sum() / np.sqrt((tw ** 2).sum() * (td ** 2).sum()))
        r_unadj = C.corr_with_all(w, M, pats)            # P3 와 동일 경로 (재현 점검용)
    r_adj = partial(r_gw, r_gd, r_wd)
    C.log(f"  WBP5 ~ 기술 축 (환자 내 순위 상관, 전체) {r_wd:.3f}")

    # ---- CR0 재현 점검 ----------------------------------------------------------
    cr0 = []
    if args.check_dir:
        cd = args.check_dir
        if args.cohort == "gse198315":
            f3 = os.path.join(cd, "D3_transcriptome_depth_partial.csv.gz")
            if os.path.exists(f3):
                old = pd.read_csv(f3).set_index("gene")
                common = [g for g in names if g in old.index]
                a1 = pd.Series(r_unadj, index=names)[common]
                a2 = pd.Series(r_adj, index=names)[common]
                cr0.append({"check": "D3 rho_within", "n_genes_old": len(old), "n_genes_new": len(names),
                            "n_common": len(common),
                            "max_abs_diff": float(np.nanmax(np.abs(a1.values - old.loc[common, "rho_within"].values)))})
                cr0.append({"check": "D3 rho_depth_partial", "n_genes_old": len(old), "n_genes_new": len(names),
                            "n_common": len(common),
                            "max_abs_diff": float(np.nanmax(np.abs(a2.values - old.loc[common, "rho_depth_partial"].values)))})
            f4 = os.path.join(cd, "D4_Puram_enrichment_depth_partial.csv")
            if os.path.exists(f4):
                old4 = pd.read_csv(f4).set_index("signature")
                rng = np.random.default_rng(C.SEED)            # step7b 와 같은 rng 순서
                for nm, gl in (("PURAM30", pres30), ("PURAM100", pres100)):
                    r, _ = C.set_enrichment(r_adj, names, bins, gl, rng, args.n_perm)
                    if nm in old4.index:
                        for col in ("perm_p_upper", "top200_overlap", "mean_rho_set"):
                            cr0.append({"check": f"D4 {nm} {col}", "old": old4.at[nm, col], "new": r[col],
                                        "max_abs_diff": abs(float(old4.at[nm, col]) - float(r[col]))})
        else:
            f3 = os.path.join(cd, "P3_transcriptome_rho_within_patient.csv.gz")
            if os.path.exists(f3):
                old = pd.read_csv(f3).set_index("gene")
                common = [g for g in names if g in old.index]
                a1 = pd.Series(r_unadj, index=names)[common]
                cr0.append({"check": "P3 within_patient rho", "n_genes_old": len(old), "n_genes_new": len(names),
                            "n_common": len(common),
                            "max_abs_diff": float(np.nanmax(np.abs(a1.values - old.loc[common, "rho"].values)))})
            f4 = os.path.join(cd, "P3_Puram_enrichment.csv")
            if os.path.exists(f4):
                old4 = pd.read_csv(f4)
                old4 = old4[old4["mode"] == "within_patient"].set_index("signature")
                rng = np.random.default_rng(C.SEED)
                for nm, gl in (("PURAM30", pres30), ("PURAM100", pres100)):
                    r, _ = C.set_enrichment(r_unadj, names, bins, gl, rng, 100)   # 결정적 항목만 비교
                    if nm in old4.index:
                        for col in ("top200_overlap", "top200_hypergeom_p", "mean_rho_set"):
                            cr0.append({"check": f"P3 {nm} {col}", "old": old4.at[nm, col], "new": r[col],
                                        "max_abs_diff": abs(float(old4.at[nm, col]) - float(r[col]))})
    cr0 = pd.DataFrame(cr0)
    tol_ok = bool(len(cr0)) and bool((cr0.max_abs_diff.fillna(0) < 1e-6).all())
    cr0.to_csv(out("CR0_reproduction_check.csv"), index=False)

    # ---- CR1 프로파일 ----------------------------------------------------------
    prof = pd.DataFrame({"gene": names, "rho_within_unadj": r_unadj, "rho_adj": r_adj,
                         "rho_gene_vs_tech": r_gd, "detection": det, "mean_expr": mean_expr, "bin": bins,
                         "in_PURAM30": np.isin(names, pres30), "in_PURAM100": is_p100})
    prof.sort_values("rho_adj", ascending=False).to_csv(out("CR1_profile_WBP5.csv.gz"), index=False)
    full = {}
    for lab, vec in (("adj", r_adj), ("unadj", r_unadj)):
        rng = np.random.default_rng(C.SEED)
        for nm, gl in (("PURAM30", pres30), ("PURAM100", pres100)):
            r, _ = C.set_enrichment(vec, names, bins, gl, rng, args.n_perm)
            full[(lab, nm)] = r

    # ---- CR2 LOPO --------------------------------------------------------------
    C.log("CR2 leave-one-patient-out")
    dec = Decomp(R, M, pats, tw, td, mode)
    rows = []
    for j, p in enumerate(dec.pids):
        if dec.ncell[j] < MIN_CELLS:
            continue
        keepm = np.ones(len(dec.pids), bool)
        keepm[j] = False
        g_w, g_d, g_wd, g_p, mexp = dec.profile(keepm)
        b2 = C.expression_bins(mexp)
        for lab, vec in (("adj", g_p), ("unadj", g_w)):
            rng = np.random.default_rng(C.SEED)
            for nm, gl in (("PURAM30", pres30), ("PURAM100", pres100)):
                r, _ = C.set_enrichment(vec, names, b2, gl, rng, args.n_perm)
                rows.append({"left_out": p, "n_cells_left_out": int(dec.ncell[j]), "profile": lab,
                             "signature": nm, "rho_wbp5_tech": g_wd, **{k: r[k] for k in (
                                 "mean_rho_set", "null_mean", "perm_p_upper", "top100_overlap",
                                 "top200_overlap", "top200_expected", "top200_hypergeom_p", "top500_overlap")}})
    cr2 = pd.DataFrame(rows)
    cr2.to_csv(out("CR2_LOPO.csv"), index=False)
    # 전체(제외 없음)를 분해식으로 다시 계산해 위 전체 값과 같은지 확인
    _, _, _, chk_p, _ = dec.profile(np.ones(len(dec.pids), bool))
    decomp_diff = float(np.nanmax(np.abs(chk_p - r_adj)))

    # ---- CR3 / CR4 전체 유전자 표적 스캔 ---------------------------------------------
    C.log(f"CR3 전체 유전자 표적 스캔 ({len(names):,}개) — 가장 오래 걸리는 단계")
    targets = np.arange(len(names))
    scan = target_scan(R, mode, td, r_gd, bins, is_p100, targets)
    tab = pd.DataFrame({"gene": names, "detection": det, "mean_expr": mean_expr, "bin": bins,
                        "in_PURAM100": is_p100, **{f"enr_{k}": v for k, v in scan.items()}})
    wst = single_target_stats(r_adj, bins, is_p100)
    w_mean, w_det = float(np.mean(w)), float(np.mean(w > 0))

    # 사전 지정 유전자 대응
    roles = {}
    for lab, gl in (("family_core", FAMILY_CORE), ("family_ext", FAMILY_EXT), ("outgroup", OUTGROUP),
                    ("neg_ctrl", NEG_CTRL), ("pos_ctrl", POS_CTRL)):
        for s in gl:
            roles[s] = (lab, resolve_symbol(s, D["var_names"]))
    special = {d for _, d in roles.values() if d is not None}
    # 기준 유전자: WBP5 와 평균 발현이 가장 가까운 n_ref 개 (Puram 100, 계열, 외군, 대조 제외)
    elig = tab[~tab.in_PURAM100 & ~tab.gene.isin(special)].copy()
    elig["dist"] = (elig.mean_expr - w_mean).abs()
    ref = elig.nsmallest(args.n_ref, "dist")
    ref.to_csv(out("CR3_reference_genes.csv"), index=False)
    ref_z = ref.enr_z.values
    allz = tab.loc[~tab.in_PURAM100, "enr_z"].values

    def pct(z):
        return (100 * np.mean(ref_z < z), 100 * np.mean(allz < z)) if np.isfinite(z) else (np.nan, np.nan)

    fam_rows = [{"symbol": "WBP5 (TCEAL9)", "role": "target", "data_symbol": D["gene"], "in_universe": False,
                 "detection": w_det, "mean_expr": w_mean, "evaluable": w_det >= EVAL_DET,
                 **{f"enr_{k}": v for k, v in wst.items()},
                 "pct_vs_reference": pct(wst["z"])[0], "pct_vs_all_nonPuram": pct(wst["z"])[1]}]
    for s, (lab, ds) in roles.items():
        r = {"symbol": s, "role": lab, "data_symbol": ds}
        if ds is None:
            r.update({"in_universe": False, "detection": np.nan, "evaluable": False, "note": "데이터에 없음"})
        elif ds not in name_idx:
            r.update({"in_universe": False, "detection": D["det_map"].get(ds, np.nan), "evaluable": False,
                      "note": "유전자 집합 기준 미달 (검출률/분산)"})
        else:
            t = tab.iloc[name_idx[ds]]
            r.update({"in_universe": True, "detection": t.detection, "mean_expr": t.mean_expr,
                      "evaluable": bool(t.detection >= EVAL_DET), "in_PURAM100": bool(t.in_PURAM100),
                      **{c: t[c] for c in tab.columns if c.startswith("enr_")},
                      "pct_vs_reference": pct(t.enr_z)[0], "pct_vs_all_nonPuram": pct(t.enr_z)[1]})
        fam_rows.append(r)
    fam = pd.DataFrame(fam_rows)
    fam.to_csv(out("CR3_family_and_controls.csv"), index=False)
    pur_tab = tab[tab.in_PURAM100].copy()
    pur_tab["pct_vs_reference"] = [pct(z)[0] for z in pur_tab.enr_z]
    tab["pct_vs_reference"] = [pct(z)[0] for z in tab.enr_z]
    tab.to_csv(out("CR3_all_targets.csv.gz"), index=False)

    # ---- CR3b Xq22 모듈 보정 -----------------------------------------------------
    C.log("CR3b Xq22 모듈 보정")
    mod_members = [ds for s, (lab, ds) in roles.items()
                   if lab in ("family_core", "family_ext") and ds in name_idx
                   and det[name_idx[ds]] >= EVAL_DET]
    cr3b = []
    if mod_members:
        mod = np.mean(np.vstack([R[name_idx[g]].astype(np.float64) for g in mod_members]), axis=0)
        tm = zrank_within(mod, pats)
        ssR = (R.astype(np.float64) ** 2).sum(axis=1)
        for lab, covs in (("tech_only", [td]), ("tech_plus_Xq22_module", [td, tm])):
            Q, _ = np.linalg.qr(np.column_stack(covs))
            w_res = tw - Q @ (Q.T @ tw)
            RQ = np.column_stack([matvec64(R, Q[:, i]) for i in range(Q.shape[1])])
            ss_res = ssR - (RQ ** 2).sum(1)
            with np.errstate(invalid="ignore", divide="ignore"):
                rr = matvec64(R, w_res) / np.sqrt(ss_res * (w_res ** 2).sum())
            rng = np.random.default_rng(C.SEED)
            for nm, gl in (("PURAM30", pres30), ("PURAM100", pres100)):
                r, _ = C.set_enrichment(rr, names, bins, gl, rng, args.n_perm)
                cr3b.append({"adjustment": lab, "signature": nm, "module_members": ";".join(mod_members),
                             "rho_wbp5_module": float(np.mean(tw * tm) / np.sqrt(np.mean(tw ** 2) * np.mean(tm ** 2))),
                             **{k: r[k] for k in ("mean_rho_set", "null_mean", "perm_p_upper", "top200_overlap",
                                                  "top200_expected", "top200_hypergeom_p")}})
    cr3b = pd.DataFrame(cr3b)
    cr3b.to_csv(out("CR3b_Xq22_module_adjusted.csv"), index=False)

    # ---- 출력 ---------------------------------------------------------------------
    pd.set_option("display.width", 240)
    f = lambda v: f"{v:.4g}"  # noqa
    print(f"\n######## CRED [{args.tag}] {args.cohort} ########")
    print(f"세포 {n:,} / 환자 {len(dec.pids)} / 유전자 {len(names):,} / 기술 축 {D['tech_src']}"
          f" / WBP5~기술 축 {r_wd:.3f} / 분해식 검산 차이 {decomp_diff:.2e}")
    print("\n=== CR0 재현 점검 ===")
    if len(cr0):
        print(cr0.to_string(index=False, float_format=f))
        print("  판정:", "일치 (허용오차 1e-6)" if tol_ok else "**불일치 — 이후 결과를 응답서에 쓰지 마십시오**")
    else:
        print("  --check-dir 미지정 또는 파일 없음 → 재현 점검 안 함")
    print("\n=== CR1 전체 농축 (이번 실행) ===")
    for (lab, nm), r in full.items():
        print(f"  [{lab:5s}] {nm:8s} 평균ρ {r['mean_rho_set']:.4f} (귀무 {r['null_mean']:.4f}) "
              f"p {r['perm_p_upper']:.2g} / 상위200 {r['top200_overlap']} (기대 {r['top200_expected']:.2f}, "
              f"p {r['top200_hypergeom_p']:.2g})")
    print("\n=== CR2 LOPO (PURAM100) ===")
    for lab in ("adj", "unadj"):
        s = cr2[(cr2.profile == lab) & (cr2.signature == "PURAM100")]
        if len(s):
            worst = s.loc[s.perm_p_upper.idxmax()]
            print(f"  [{lab}] {len(s)}회: p 최대 {s.perm_p_upper.max():.2g} (제외 환자 {worst.left_out}), "
                  f"상위200 {int(s.top200_overlap.min())}–{int(s.top200_overlap.max())}, "
                  f"p<0.05 {int((s.perm_p_upper < 0.05).sum())}/{len(s)}, p<0.01 {int((s.perm_p_upper < 0.01).sum())}/{len(s)}")
    print(cr2[cr2.signature == "PURAM100"][["left_out", "n_cells_left_out", "profile", "mean_rho_set", "null_mean",
                                            "perm_p_upper", "top200_overlap", "top200_hypergeom_p"]]
          .to_string(index=False, float_format=f))
    print("\n=== CR3 계열 특이성 · CR4 대조 (PURAM100, 해석적 z / 기준 1,000개 대비 백분위) ===")
    print(f"  WBP5 해석적 귀무 평균 {wst['mu']:.4f} vs 순열 귀무 평균 {full[('adj', 'PURAM100')]['null_mean']:.4f} (검산)")
    cols = ["symbol", "role", "data_symbol", "detection", "evaluable", "enr_E1", "enr_z", "enr_top",
            "pct_vs_reference", "pct_vs_all_nonPuram"]
    print(fam[[c for c in cols if c in fam.columns]].to_string(index=False, float_format=f))
    print(f"  기준 유전자 z: 중앙값 {np.median(ref_z):.2f}, 95백분위 {np.percentile(ref_z, 95):.2f}, "
          f"z>1.96 비율 {100 * np.mean(ref_z > 1.96):.1f}%")
    print(f"  Puram 100 구성원(자기 제외) 기준 대비 백분위 중앙값 {np.nanmedian(pur_tab.pct_vs_reference):.1f}")
    ev_core = fam[(fam.role == "family_core") & (fam.evaluable == True)]  # noqa
    wz = wst["z"]
    print(f"  평가 가능 핵심 TCEAL {len(ev_core)}개: WBP5 보다 z 가 높은 구성원 "
          f"{', '.join(ev_core.loc[ev_core.enr_z >= wz, 'symbol']) or '없음'}; "
          f"95백분위 이상 {int((ev_core.pct_vs_reference >= 95).sum())}/{len(ev_core)}")
    neg = fam[(fam.role == "neg_ctrl") & (fam.evaluable == True)]  # noqa
    if len(neg):
        print(f"  음성 대조 {len(neg)}개: 백분위 중앙값 {neg.pct_vs_reference.median():.1f}, "
              f"95 이상 {int((neg.pct_vs_reference >= 95).sum())}개")
    print("\n=== CR3b Xq22 모듈 보정 ===")
    if len(cr3b):
        print(cr3b.drop(columns=["module_members"]).to_string(index=False, float_format=f))
        print(f"  모듈 구성: {', '.join(mod_members)}")
    else:
        print("  평가 가능한 계열 구성원 없음 → 생략")

    import scipy
    with open(out("CR_run_info.txt"), "a", encoding="utf-8") as fh:
        fh.write(f"finished {pd.Timestamp.now()}\npython {sys.version}\nnumpy {np.__version__}\n"
                 f"pandas {pd.__version__}\nscipy {scipy.__version__}\n"
                 f"reproduction_ok {tol_ok}\ndecomp_diff {decomp_diff}\n"
                 f"x_looks_raw {D['looks_raw']}\ntech {D['tech_src']}\n")
    C.log(f"완료 → {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
