# -*- coding: utf-8 -*-
"""
wbp5_rev_common.py — CSBJ-D-26-00923 리비전 공통 함수 (2026-10-19)

GSE103322_Puram_pEMT_20261019.py 와 GSE198315_step7_20261019.py 가 함께 씁니다.
세 파일을 같은 폴더(WBP5 revision)에 두십시오.

설계 원칙
  - 유전자 목록은 파일에서만 읽습니다(하드코딩 금지).
  - 2017년 기호와 현재 HGNC 기호 차이는 ALIASES 로 양방향 해결하고 로그에 남깁니다.
  - 모든 무작위 절차는 seed 고정.
"""
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats

SEED = 42
warnings.filterwarnings("ignore")

# Puram 2017 (hg19 시절 기호) <-> 현재 HGNC 승인 기호
ALIASES = {
    "PRKCDBP": "CAVIN3",
    "CXCR7": "ACKR3",
    "DFNA5": "GSDME",
    "LEPREL1": "P3H2",
}
ALIASES_REV = {v: k for k, v in ALIASES.items()}

WBP5_ALIASES = ["WBP5", "TCEAL9"]
MARKER_PANEL = ["LAMC2", "ITGA5", "KRT14", "VIM", "KRT17"]   # 원고 Table 2, 변경 금지


def log(msg):
    sys.stderr.write(f"[{pd.Timestamp.now():%H:%M:%S}] {msg}\n")
    sys.stderr.flush()


# ---------------------------------------------------------------------------
# 유전자 목록
# ---------------------------------------------------------------------------

def read_signature(path):
    """p-EMT 목록 읽기. xlsx/xls (Rank, Gene 열) 또는 txt (한 줄 한 유전자, 순서=순위)."""
    if str(path).lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(path)
        cols = {c.lower().strip(): c for c in df.columns}
        if "gene" not in cols:
            raise SystemExit(f"오류: {path} 에 'Gene' 열이 없습니다. 열: {list(df.columns)}")
        genes = df[cols["gene"]].astype(str).str.strip()
        if "rank" in cols:
            ranks = pd.to_numeric(df[cols["rank"]], errors="coerce")
            df = pd.DataFrame({"rank": ranks, "gene": genes}).dropna()
            df = df.sort_values("rank")
            genes = df["gene"]
        genes = [g for g in genes.tolist() if g and g.lower() != "nan"]
    else:
        with open(path, encoding="utf-8") as fh:
            genes = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    genes = [g.upper() for g in genes]
    dup = sorted({g for g in genes if genes.count(g) > 1})
    if dup:
        raise SystemExit(f"오류: 목록에 중복 유전자 {dup}")
    log(f"p-EMT 목록 {len(genes)}개 읽음 (1위 {genes[0]}, 마지막 {genes[-1]})")
    return genes


def resolve_genes(genes, available):
    """목록 유전자를 데이터의 실제 행 이름으로 대응. 반환: (used 목록, 대응표 DataFrame)."""
    avail = set(available)
    rows, used = [], []
    for g in genes:
        hit, how = None, "missing"
        if g in avail:
            hit, how = g, "exact"
        elif g in ALIASES and ALIASES[g] in avail:
            hit, how = ALIASES[g], "alias(old->new)"
        elif g in ALIASES_REV and ALIASES_REV[g] in avail:
            hit, how = ALIASES_REV[g], "alias(new->old)"
        rows.append({"list_symbol": g, "data_symbol": hit, "match": how})
        if hit is not None and hit not in used:
            used.append(hit)
    return used, pd.DataFrame(rows)


def resolve_one(candidates, available):
    avail = set(available)
    for g in candidates:
        if g in avail:
            return g
    return None


# ---------------------------------------------------------------------------
# 점수 계산
# ---------------------------------------------------------------------------

def zscore_rows(M):
    """행(유전자)별 z-표준화. 표준편차 0 인 행은 NaN."""
    M = np.asarray(M, dtype=np.float64)
    mu = M.mean(axis=1, keepdims=True)
    sd = M.std(axis=1, keepdims=True)
    sd[sd == 0] = np.nan
    return (M - mu) / sd


def zmean_score(M_sig):
    """원고 Section 2.3 과 동일: 유전자별 z 의 평균 (세포 단위)."""
    return np.nanmean(zscore_rows(M_sig), axis=0)


def control_binned_score(M_all, gene_names, sig_genes, rng, n_bins=25, n_ctrl=100):
    """
    Seurat AddModuleScore / Tirosh 2016 방식 (원고에서 민감도 분석으로 쓴 방식과 동일 설정).
    평균 발현으로 유전자를 n_bins 구간에 나누고, 서명 유전자마다 같은 구간의
    비서명 유전자 n_ctrl 개를 대조로 뽑아 (서명 평균 - 대조 평균)을 점수로 함.
    """
    names = np.asarray(gene_names)
    means = M_all.mean(axis=1)
    order = np.argsort(np.argsort(means))
    bins = np.floor(order / len(means) * n_bins).astype(int)
    idx_of = {g: i for i, g in enumerate(names)}
    sig_idx = [idx_of[g] for g in sig_genes if g in idx_of]
    sig_set = set(sig_idx)
    ctrl = []
    for i in sig_idx:
        pool = np.where((bins == bins[i]))[0]
        pool = np.array([p for p in pool if p not in sig_set])
        k = min(n_ctrl, len(pool))
        ctrl.extend(rng.choice(pool, size=k, replace=False).tolist())
    return M_all[sig_idx].mean(axis=0) - M_all[ctrl].mean(axis=0)


# ---------------------------------------------------------------------------
# 상관 / 검정
# ---------------------------------------------------------------------------

def cliffs_delta(x, y):
    """Cliff's delta = 2U/(n1 n2) - 1  (x 가 크면 양수)."""
    x, y = np.asarray(x), np.asarray(y)
    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 2 * u / (len(x) * len(y)) - 1, p


def spearman_cluster_boot(x, y, groups, rng, n_boot=2000):
    """세포 수준 Spearman + 환자 단위 재표본 부트스트랩 95% CI."""
    x, y, groups = np.asarray(x, float), np.asarray(y, float), np.asarray(groups)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, groups = x[ok], y[ok], groups[ok]
    rho, p = stats.spearmanr(x, y)
    pats = np.unique(groups)
    idx_by = {g: np.where(groups == g)[0] for g in pats}
    boot = []
    for _ in range(n_boot):
        pick = rng.choice(pats, size=len(pats), replace=True)
        ii = np.concatenate([idx_by[g] for g in pick])
        if np.std(x[ii]) > 0 and np.std(y[ii]) > 0:
            boot.append(stats.spearmanr(x[ii], y[ii])[0])
    lo, hi = (np.percentile(boot, [2.5, 97.5]) if boot else (np.nan, np.nan))
    return {"rho_cell": rho, "p_cell": p, "rho_ci_lo": lo, "rho_ci_hi": hi,
            "n_cells": int(ok.sum()), "n_patients": len(pats), "n_boot": len(boot)}


def within_group_spearman(x, y, groups, min_n=20):
    """환자별 Spearman → 중앙값, IQR, 양수 환자 수, sign test, LOPO 범위."""
    x, y, groups = np.asarray(x, float), np.asarray(y, float), np.asarray(groups)
    rows = []
    for g in pd.unique(groups):
        s = groups == g
        xs, ys = x[s], y[s]
        ok = np.isfinite(xs) & np.isfinite(ys)
        if ok.sum() < min_n or np.std(xs[ok]) == 0 or np.std(ys[ok]) == 0:
            continue
        r, p = stats.spearmanr(xs[ok], ys[ok])
        rows.append({"patient": g, "n_cells": int(ok.sum()), "rho": r, "p_within": p})
    per = pd.DataFrame(rows)
    if per.empty:
        return per, {"n_patients": 0}
    r = per["rho"].values
    npos = int((r > 0).sum())
    lopo = [np.median(np.delete(r, i)) for i in range(len(r))] if len(r) >= 4 else [np.nan]
    summ = {
        "median_rho": float(np.median(r)),
        "q1": float(np.percentile(r, 25)), "q3": float(np.percentile(r, 75)),
        "n_positive": npos, "n_patients": len(r),
        "sign_test_p": stats.binomtest(npos, len(r), 0.5).pvalue,
        "lopo_median_min": float(np.nanmin(lopo)), "lopo_median_max": float(np.nanmax(lopo)),
    }
    return per, summ


def mixed_model(yv, xv, groups):
    """y ~ x + (1|patient), 두 변수 모두 z-표준화 → 계수는 표준화 기울기."""
    import statsmodels.formula.api as smf
    d = pd.DataFrame({"y": yv, "x": xv, "g": groups}).dropna()
    for c in ("y", "x"):
        d[c] = (d[c] - d[c].mean()) / d[c].std()
    out = {}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = smf.mixedlm("y ~ x", d, groups=d["g"]).fit(reml=True)
        out["mixedlm_converged"] = bool(m.converged)
        ci = m.conf_int().loc["x"]
        out.update({"mixedlm_std_coef": float(m.params["x"]),
                    "mixedlm_ci_lo": float(ci[0]), "mixedlm_ci_hi": float(ci[1]),
                    "mixedlm_p": float(m.pvalues["x"]),
                    "mixedlm_group_var": float(m.cov_re.iloc[0, 0])})
        o = smf.ols("y ~ x", d).fit()
        out.update({"naive_ols_std_coef": float(o.params["x"]), "naive_ols_p": float(o.pvalues["x"])})
    except Exception as e:  # noqa
        out["mixedlm_error"] = str(e)
    return out


def incremental_information(resp, predictors, wbp5, groups):
    """
    resp ~ markers (+ WBP5).  모든 변수 z-표준화.
    반환: R²(마커), R²(마커+WBP5), ΔR², WBP5 표준화 부분계수(OLS, 혼합모형),
          그리고 비교 기준으로 각 예측변수를 하나씩 뺐을 때의 ΔR² (drop-one).
    """
    import statsmodels.formula.api as smf
    d = pd.DataFrame(predictors)
    d["resp"] = resp
    d["wbp5"] = wbp5
    d["g"] = groups
    d = d.dropna()
    for c in d.columns:
        if c != "g":
            sd = d[c].std()
            d[c] = (d[c] - d[c].mean()) / sd if sd > 0 else 0.0
    preds = [c for c in predictors.keys() if d[c].std() > 0]
    base = "resp ~ " + " + ".join(preds)
    full = base + " + wbp5"
    out = {"predictors": ";".join(preds), "n_cells": len(d), "n_patients": d["g"].nunique()}
    try:
        mb, mf = smf.ols(base, d).fit(), smf.ols(full, d).fit()
        out.update({"R2_markers": mb.rsquared, "R2_markers_plus_WBP5": mf.rsquared,
                    "delta_R2_WBP5": mf.rsquared - mb.rsquared,
                    "WBP5_std_partial_coef_ols": mf.params["wbp5"],
                    "WBP5_partial_p_ols": mf.pvalues["wbp5"]})
        # 비교: 각 변수를 뺐을 때 R² 감소분
        allv = preds + ["wbp5"]
        for v in allv:
            rest = [u for u in allv if u != v]
            r2 = smf.ols("resp ~ " + " + ".join(rest), d).fit().rsquared
            out[f"dropone_dR2_{v}"] = mf.rsquared - r2
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mm = smf.mixedlm(full, d, groups=d["g"]).fit(reml=True)
        out["mixedlm_converged"] = bool(mm.converged)
        ci = mm.conf_int().loc["wbp5"]
        out.update({"WBP5_std_partial_coef_mixedlm": float(mm.params["wbp5"]),
                    "WBP5_partial_ci_lo_mixedlm": float(ci[0]),
                    "WBP5_partial_ci_hi_mixedlm": float(ci[1]),
                    "WBP5_partial_p_mixedlm": float(mm.pvalues["wbp5"])})
    except Exception as e:  # noqa
        out["increment_error"] = str(e)
    return out


# ---------------------------------------------------------------------------
# 전사체 전체 상관 + 집합 농축
# ---------------------------------------------------------------------------

def _rank_rows(M):
    return stats.rankdata(M, axis=1)


def corr_with_all(target, M, groups=None):
    """
    target(세포 벡터) 과 M(유전자 x 세포)의 모든 행 사이 Spearman.
    groups 를 주면 환자 내 순위(환자별 z-표준화한 순위)로 계산 → 환자 간 차이 제거.
    반환: np.ndarray (유전자 수)
    """
    target = np.asarray(target, float)
    if groups is None:
        R = _rank_rows(M)
        t = stats.rankdata(target)
        R = R - R.mean(axis=1, keepdims=True)
        t = t - t.mean()
        num = R @ t
        den = np.sqrt((R ** 2).sum(axis=1) * (t ** 2).sum())
    else:
        groups = np.asarray(groups)
        R = np.empty_like(M, dtype=np.float32)
        t = np.empty(len(target), dtype=np.float64)
        for g in pd.unique(groups):
            s = np.where(groups == g)[0]
            r = _rank_rows(M[:, s]).astype(np.float32)
            r -= r.mean(axis=1, keepdims=True)
            sd = r.std(axis=1, keepdims=True)
            sd[sd == 0] = np.nan
            R[:, s] = r / sd
            tt = stats.rankdata(target[s])
            t[s] = (tt - tt.mean()) / (tt.std() if tt.std() > 0 else np.nan)
        R = np.nan_to_num(R)
        t = np.nan_to_num(t)
        num = R.astype(np.float64) @ t
        den = np.sqrt((R.astype(np.float64) ** 2).sum(axis=1) * (t ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def expression_bins(mean_expr, n_bins=25):
    order = np.argsort(np.argsort(mean_expr))
    return np.floor(order / len(mean_expr) * n_bins).astype(int)


def set_enrichment(rho, gene_names, bins, set_genes, rng, n_perm=10000, top_ks=(100, 200, 500)):
    """
    집합 유전자의 평균 rho 를 같은 발현 구간에서 뽑은 무작위 집합 n_perm 개와 비교 (단측, 상향).
    + 상위 k 개 상관 유전자와의 겹침 hypergeometric.
    """
    names = np.asarray(gene_names)
    idx_of = {g: i for i, g in enumerate(names)}
    sidx = np.array([idx_of[g] for g in set_genes if g in idx_of])
    ok = np.isfinite(rho)
    obs = float(np.nanmean(rho[sidx]))
    pools = {}
    for i in sidx:
        b = bins[i]
        if b not in pools:
            pools[b] = np.setdiff1d(np.where((bins == b) & ok)[0], sidx)
    picks = np.column_stack([rng.choice(pools[bins[i]], size=n_perm) for i in sidx])
    null = np.nanmean(rho[picks], axis=1)
    p = (1 + np.sum(null >= obs)) / (1 + n_perm)
    out = {"set_size_in_data": len(sidx), "mean_rho_set": obs,
           "null_mean": float(null.mean()), "null_95pct": float(np.percentile(null, 95)),
           "null_99pct": float(np.percentile(null, 99)), "perm_p_upper": p,
           "perm_resolution": 1 / (1 + n_perm), "n_perm": n_perm}
    # 순위 기반
    valid = np.where(ok)[0]
    order = valid[np.argsort(-rho[valid])]
    rank_of = {g: r + 1 for r, g in enumerate(order)}
    N, K = len(valid), len(set(sidx) & set(valid))
    for k in top_ks:
        top = set(order[:k])
        x = len(top & set(sidx))
        out[f"top{k}_overlap"] = x
        out[f"top{k}_expected"] = k * K / N
        out[f"top{k}_hypergeom_p"] = stats.hypergeom.sf(x - 1, N, K, k)
    per_gene = pd.DataFrame({
        "gene": names[sidx], "rho": rho[sidx],
        "rank_among_all": [rank_of.get(i, np.nan) for i in sidx],
        "percentile": [100 * (1 - (rank_of.get(i, np.nan) - 1) / N) for i in sidx],
    }).sort_values("rank_among_all")
    return out, per_gene


def matched_null_setscore(wbp5, Z, bins, set_idx, rng, n_perm=10000, matched=True):
    """
    YAP/TAZ 등: ρ(WBP5, 집합 z-평균 점수) 의 귀무분포.
    matched=True → 각 집합 유전자와 같은 발현 구간에서 대체 유전자 추출.
    Z: 유전자 x 세포 z-행렬 (NaN 행 제외된 것).  양측(|ρ|) p.
    """
    t = stats.rankdata(wbp5)
    t = (t - t.mean()) / t.std()
    def rho_of(idx):
        s = np.nanmean(Z[idx], axis=0)
        r = stats.rankdata(s)
        r = (r - r.mean()) / r.std()
        return float(np.mean(r * t))
    obs = rho_of(set_idx)
    all_idx = np.arange(Z.shape[0])
    others = np.setdiff1d(all_idx, set_idx)
    pools = {b: np.intersect1d(np.where(bins == b)[0], others) for b in np.unique(bins[set_idx])}
    if matched:
        picks = np.column_stack([rng.choice(pools[bins[i]], size=n_perm) for i in set_idx])
    else:
        picks = np.vstack([rng.choice(others, size=len(set_idx), replace=False)
                           for _ in range(n_perm)])
    null = np.array([rho_of(row) for row in picks])
    p = (1 + np.sum(np.abs(null) >= abs(obs))) / (1 + n_perm)
    return {"rho_obs": obs, "null_mean": float(null.mean()),
            "null_95pct_abs": float(np.percentile(np.abs(null), 95)),
            "perm_p_two_sided": p, "resolution": 1 / (1 + n_perm), "n_perm": n_perm,
            "matched_on_expression_bins": matched}
