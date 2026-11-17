#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R2_specificity_sensitivity_20261116c.py   (revision c, 16 Nov 2026)
======================================
Revision c: adds --extra LABEL=FILE (repeatable). With --extra the script runs the specificity blocks
(R2-5 to R2-8) for the extra gene sets only, as a sensitivity analysis added after the pre-specified
analysis had been read (plan addendum: reproducibility/pre_analysis_plan_R2_addendum_20261116.md);
blocks R2-1 to R2-4 are skipped and the pre-specified outputs are not recomputed.
Revision b: the sensitivity profiles (down-sampling, conditional analysis) now use the same
normalization as the primary profile in both cohorts; in revision a they used the sums-of-squares
normalization throughout, which differs from the GSE198315 primary profile for sparsely detected genes.
Second-round revision analyses (CSBJ-D-26-00923R1), editor/reviewer points 1 and 3.
Decision rules were written down before the first run:
  reproducibility/pre_analysis_plan_R2_20261116.md

Continuity with the first-round analysis
  Same cells, genes, technical covariate, within-patient z-ranks and partial-correlation
  formula as CRED_profiles_20261021.py (its loaders and Decomp class are imported, not
  re-implemented). Block R2-0 recomputes the primary profile and checks it against the
  CR1 profile file of the first round when --check-dir is given.

Blocks
  R2-0  reproduction of the primary profile
  R2-1  contribution of each patient to the pooled estimate (cells, share of weight)
  R2-2  p-EMT enrichment under (a) expression-matched null [primary, as before],
        (b) expression- AND detection-rate-matched null; each for the pooled profile and
        for the patient-balanced profile (equal-weight Fisher-z mean of per-patient
        partial coefficients, patients with >= 20 cells)
  R2-3  patient-level test: per-patient enrichment (set mean minus matched expectation),
        sign test and Wilcoxon signed-rank test across patients
  R2-4  capped down-sampling: every patient capped at the median cell count, 100 repeats
  R2-5  overlap among the p-EMT signature and the related programs
  R2-6  enrichment of each program, full and mutually exclusive ("disjoint") gene sets
  R2-7  conditional enrichment: p-EMT-only genes after adjusting for each competitor's
        module score, and the reverse
  R2-8  module-level partial correlations (WBP5 ~ p-EMT | technical + competitor; reverse)

Run (Windows cmd, "WBP5 revision" folder; same inputs as CRED_profiles):
  python R2_specificity_sensitivity_20261116.py --cohort gse103322 ^
     --matrix GSE103322_HNSCC_all_data.txt.gz --pemt "p-EMT gene list.xlsx" ^
     --tag gse103322 --outdir results_R2\\gse103322 --derive-caf ^
     --check-dir results_cred\\gse103322
  python R2_specificity_sensitivity_20261116.py --cohort gse198315 ^
     --h5ad gse198315_checkpoint_full_annotated_cnv.h5ad --pemt "p-EMT gene list.xlsx" ^
     --tag seed42 --outdir results_R2\\seed42 ^
     --caf results_R2\\gse103322\\R2_CAF_derived_GSE103322.txt --check-dir results_cred\\seed42
  (repeat the second command for the seed2026 and seed7 h5ad files)

Gene-set libraries: if --gmt is not given, the two Enrichr libraries MSigDB_Hallmark_2020
and GO_Biological_Process_2023 are downloaded once into --outdir\\..\\genesets and reused.
Needs: wbp5_rev_common.py, CRED_profiles_20261021.py, GSE103322_Puram_pEMT_20261019.py.
"""
import argparse
import hashlib
import os
import sys
import urllib.request
import warnings

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wbp5_rev_common as C                    # noqa: E402
import CRED_profiles_20261021 as CR            # noqa: E402

warnings.filterwarnings("ignore")

MIN_CELLS = 20
N_EXPR_BINS = 25
N_DET_BINS = 4
MIN_POOL = 10
TOPK = 200
ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={}"
LIBS = ("MSigDB_Hallmark_2020", "GO_Biological_Process_2023")
# label -> (library, term matched case-insensitively as a substring)
DEFAULT_SETS = {
    "HALLMARK_EMT": ("MSigDB_Hallmark_2020", "epithelial mesenchymal transition"),
    "HALLMARK_TGFB": ("MSigDB_Hallmark_2020", "tgf-beta signaling"),
    "GO_ECM_ORG": ("GO_Biological_Process_2023", "(go:0030198)"),
}


# ---------------------------------------------------------------------------
# gene sets
# ---------------------------------------------------------------------------

def read_gmt(path):
    out = {}
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            genes = [g.split(",")[0].strip().upper() for g in f[1:] if g.strip()]
            out[f[0].strip()] = [g for g in dict.fromkeys(genes) if g and not g.startswith("HTTP")]
    return out


def get_library(lib, cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    p = os.path.join(cache_dir, lib + ".gmt")
    if not os.path.exists(p):
        C.log(f"downloading {lib} from Enrichr")
        try:
            urllib.request.urlretrieve(ENRICHR.format(lib), p)
        except Exception as e:                                        # noqa: BLE001
            raise SystemExit(f"ERROR: could not download {lib} ({e}). Download it in a browser from\n"
                             f"  {ENRICHR.format(lib)}\nand save it as {p}")
    return p


def pick_term(lib_sets, needle):
    hits = [k for k in lib_sets if needle.lower() in k.lower()]
    if len(hits) != 1:
        raise SystemExit(f"ERROR: term '{needle}' matched {len(hits)} entries: {hits[:5]}")
    return hits[0], lib_sets[hits[0]]


def derive_caf(args, out):
    """Fibroblast-identity signature from the author-annotated GSE103322 cells:
    genes detected in >= 25% of fibroblasts, ranked by Cliff's delta (fibroblast vs malignant);
    top 100. Derived on cell-type labels only, never on WBP5."""
    P = CR.importlib.import_module("GSE103322_Puram_pEMT_20261019")
    expr, meta = P.load_gse103322(args.matrix)
    if "noncancer_type" not in meta:
        raise SystemExit("ERROR: 'non-cancer cell type' row not found; cannot derive the fibroblast signature")
    fib = meta.index[(meta.is_cancer == 0) & meta.noncancer_type.astype(str).str.lower().str.contains("fibro")]
    mal = meta.index[meta.is_cancer == 1]
    A = np.nan_to_num(expr[fib].values.astype(np.float32))
    B = np.nan_to_num(expr[mal].values.astype(np.float32))
    det = (A > 0).mean(axis=1)
    r = stats.rankdata(np.hstack([A, B]), axis=1)
    na, nb = A.shape[1], B.shape[1]
    u = r[:, :na].sum(axis=1) - na * (na + 1) / 2
    delta = 2 * u / (na * nb) - 1
    wb = set(C.WBP5_ALIASES)
    tab = pd.DataFrame({"gene": expr.index.astype(str), "cliffs_delta_fib_vs_malignant": delta, "det_fibroblast": det})
    tab = tab[(tab.det_fibroblast >= 0.25) & ~tab.gene.isin(wb)].sort_values("cliffs_delta_fib_vs_malignant", ascending=False)
    top = tab.head(100)
    top.to_csv(out("R2_CAF_derived_GSE103322_table.csv"), index=False)
    with open(out("R2_CAF_derived_GSE103322.txt"), "w", encoding="utf-8") as fh:
        fh.write("# fibroblast-identity signature derived in GSE103322 (see script docstring)\n")
        fh.write("\n".join(top.gene.str.upper()) + "\n")
    C.log(f"fibroblast signature: {na:,} fibroblasts vs {nb:,} malignant cells; top 100 written "
          f"(delta {top.cliffs_delta_fib_vs_malignant.min():.3f} to {top.cliffs_delta_fib_vs_malignant.max():.3f})")
    return top.gene.str.upper().tolist()


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def make_strata(mean_expr, det):
    """(a) 25 expression bins as in the first round; (b) each split into detection-rate quartiles."""
    eb = C.expression_bins(mean_expr, N_EXPR_BINS)
    s2 = np.zeros(len(eb), dtype=int)
    for b in np.unique(eb):
        i = np.where(eb == b)[0]
        q = np.floor(np.argsort(np.argsort(det[i], kind="stable"), kind="stable") / len(i) * N_DET_BINS).astype(int)
        s2[i] = b * N_DET_BINS + q
    return eb, s2


def enrich(rho, sidx, strata, fallback, rng, n_perm, exclude=None):
    """One-sided upper-tail test of the mean coefficient of a gene set against sets drawn
    gene-for-gene from the same stratum. `exclude`: genes never used as null draws."""
    sidx = np.asarray([i for i in sidx if np.isfinite(rho[i])], dtype=int)
    if len(sidx) < 5:
        return None
    ok = np.isfinite(rho)
    ban = np.zeros(len(rho), bool)
    ban[sidx] = True
    if exclude is not None:
        ban[np.asarray(list(exclude), dtype=int)] = True
    pools, n_fb = {}, 0
    cols = []
    for i in sidx:
        key = ("s", strata[i])
        if key not in pools:
            p = np.where((strata == strata[i]) & ok & ~ban)[0]
            if len(p) < MIN_POOL:
                p = np.where((fallback == fallback[i]) & ok & ~ban)[0]
                pools[key] = (p, True)
            else:
                pools[key] = (p, False)
        p, fb = pools[key]
        n_fb += fb
        cols.append(rng.choice(p, size=n_perm))
    null = rho[np.column_stack(cols)].mean(axis=1)
    obs = float(rho[sidx].mean())
    valid = np.where(ok)[0]
    order = valid[np.argsort(-rho[valid])]
    top = set(order[:TOPK].tolist())
    x, N, K = len(top & set(sidx.tolist())), len(valid), len(sidx)
    return {"n_genes": K, "mean_rho_set": obs, "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "z": float((obs - null.mean()) / null.std()), "perm_p_upper": float((1 + np.sum(null >= obs)) / (1 + n_perm)),
            f"top{TOPK}_overlap": x, f"top{TOPK}_expected": TOPK * K / N,
            f"top{TOPK}_hypergeom_p": float(stats.hypergeom.sf(x - 1, N, K, TOPK)), "n_fallback_to_expr_bin": int(n_fb)}


def stratum_expectation(rho, sidx, strata, fallback):
    """Analytic matched expectation of the set mean (no sampling)."""
    ok = np.isfinite(rho)
    ban = np.zeros(len(rho), bool)
    ban[sidx] = True
    mu = []
    for i in sidx:
        p = np.where((strata == strata[i]) & ok & ~ban)[0]
        if len(p) < MIN_POOL:
            p = np.where((fallback == fallback[i]) & ok & ~ban)[0]
        mu.append(rho[p].mean())
    return float(np.mean(mu))


def per_patient_partials(dec, mode):
    """G x P matrix of per-patient technical-covariate-adjusted coefficients."""
    out = np.full(dec.NUMw.shape, np.nan)
    for j in range(len(dec.pids)):
        if mode == "n":
            n = dec.ncell[j]
            r_gw, r_gd, r_wd = dec.NUMw[:, j] / n, dec.NUMd[:, j] / n, dec.wd[j] / n
        else:
            with np.errstate(invalid="ignore", divide="ignore"):
                r_gw = dec.NUMw[:, j] / np.sqrt(dec.SSR[:, j] * dec.ssw[j])
                r_gd = dec.NUMd[:, j] / np.sqrt(dec.SSR[:, j] * dec.ssd[j])
                r_wd = dec.wd[j] / np.sqrt(dec.ssw[j] * dec.ssd[j])
        p = CR.partial(r_gw, r_gd, r_wd)
        p[dec.SSR[:, j] <= 0] = np.nan                 # gene constant within this patient
        out[:, j] = p
    return out


def fisher_mean(PP, cols):
    z = np.arctanh(np.clip(PP[:, cols], -0.999, 0.999))
    nfin = np.isfinite(z).sum(axis=1)
    m = np.tanh(np.nanmean(z, axis=1))
    m[nfin < int(np.ceil(len(cols) / 2))] = np.nan       # evaluable in at least half of the patients
    return m


def residual_profile(R, ssR, tw, covs, mode="ss"):
    """Pooled within-patient coefficient between WBP5 and every gene, adjusted for the covariates
    (already within-patient z-ranks). Partial correlation computed from the correlation matrix with
    the same normalization as the primary profile of the first round: mode "ss" divides by sums of
    squares (GSE103322), mode "n" by the number of cells (GSE198315). With one covariate this is
    exactly the first-order formula of CRED_profiles_20261021.partial()."""
    X = np.column_stack(covs)
    n = len(tw)
    G = np.column_stack([CR.matvec64(R, tw)] + [CR.matvec64(R, X[:, i]) for i in range(X.shape[1])])
    if mode == "n":
        G = G / n
        r_wc = X.T @ tw / n
        S = X.T @ X / n
    else:
        sx = np.sqrt((X ** 2).sum(axis=0))
        sw = np.sqrt((tw ** 2).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            G = G / np.sqrt(ssR)[:, None] / np.concatenate([[sw], sx])[None, :]
        r_wc = X.T @ tw / (sw * sx)
        S = X.T @ X / np.outer(sx, sx)
    Si = np.linalg.inv(S)
    r_gw, r_gc = G[:, 0], G[:, 1:]
    with np.errstate(invalid="ignore", divide="ignore"):
        rr = (r_gw - r_gc @ Si @ r_wc) / np.sqrt((1 - np.einsum("ij,jk,ik->i", r_gc, Si, r_gc)) * (1 - r_wc @ Si @ r_wc))
    rr[ssR <= 0] = np.nan
    return rr


def partial_vec(a, b, covs):
    Q, _ = np.linalg.qr(np.column_stack(covs))
    ra, rb = a - Q @ (Q.T @ a), b - Q @ (Q.T @ b)
    return float((ra * rb).sum() / np.sqrt((ra ** 2).sum() * (rb ** 2).sum()))


def module_score(R, idx, pats):
    return CR.zrank_within(R[idx].astype(np.float64).mean(axis=0), pats)


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=["gse198315", "gse103322"])
    ap.add_argument("--h5ad")
    ap.add_argument("--matrix")
    ap.add_argument("--pemt", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--check-dir", default=None)
    ap.add_argument("--gmt", nargs="*", default=None, help="local GMT files replacing the Enrichr download")
    ap.add_argument("--caf", default=None, help="text file, one gene per line (fibroblast/CAF signature)")
    ap.add_argument("--derive-caf", action="store_true", help="GSE103322 only: derive the fibroblast signature")
    ap.add_argument("--extra", action="append", default=[], metavar="LABEL=FILE",
                    help="additional gene set (one gene per line); switches to the extra-set sensitivity mode")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-rep", type=int, default=100)
    ap.add_argument("--n-perm-rep", type=int, default=2000)
    ap.add_argument("--min-det", type=float, default=0.01)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda f: os.path.join(args.outdir, f)  # noqa: E731

    with open(out("R2x_run_info.txt" if args.extra else "R2_run_info.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"started {pd.Timestamp.now()}\nargs {vars(args)}\n")
        for f in (__file__, CR.__file__, C.__file__):
            fh.write(f"sha256 {os.path.basename(f)} {hashlib.sha256(open(f, 'rb').read()).hexdigest()}\n")

    # ---- gene sets ---------------------------------------------------------------
    raw_sets, src = {}, {}
    if args.gmt:
        libs = {os.path.splitext(os.path.basename(g))[0]: read_gmt(g) for g in args.gmt}
    else:
        cache = os.path.join(os.path.dirname(os.path.abspath(args.outdir)), "genesets")
        libs = {lib: read_gmt(get_library(lib, cache)) for lib in LIBS}
    for lab, (lib, needle) in DEFAULT_SETS.items():
        pool = libs.get(lib) or {k: v for d in libs.values() for k, v in d.items()}
        term, genes = pick_term(pool, needle)
        raw_sets[lab], src[lab] = genes, f"{lib}: {term}"
    if args.derive_caf:
        if args.cohort != "gse103322":
            raise SystemExit("ERROR: --derive-caf is for --cohort gse103322")
        raw_sets["FIBROBLAST_GSE103322"], src["FIBROBLAST_GSE103322"] = derive_caf(args, out), "derived (R2_CAF_derived_GSE103322_table.csv)"
    elif args.caf:
        raw_sets["FIBROBLAST_GSE103322"] = [ln.strip().upper() for ln in open(args.caf, encoding="utf-8")
                                            if ln.strip() and not ln.startswith("#")]
        src["FIBROBLAST_GSE103322"] = args.caf
    else:
        C.log("NOTE: no fibroblast/CAF signature supplied (--caf or --derive-caf)")
    plist = C.read_signature(args.pemt)
    extra_labels = []
    for item in args.extra:
        lab, path = item.split("=", 1)
        raw_sets[lab] = [ln.strip().upper() for ln in open(path, encoding="utf-8") if ln.strip() and not ln.startswith("#")]
        src[lab] = open(path, encoding="utf-8").readline().lstrip("# ").strip() or path
        extra_labels.append(lab)

    # ---- data and primary profile ---------------------------------------------------
    D = CR.load_198315(args) if args.cohort == "gse198315" else CR.load_103322(args)
    M, names, pats, mode = D["M"], D["names"], np.asarray(D["pats"]), D["mode"]
    wb = set(C.WBP5_ALIASES)
    idx_of = {g: i for i, g in enumerate(names)}

    def resolve(genes):
        used, _ = C.resolve_genes([g for g in genes if g not in wb], names)
        return np.array(sorted({idx_of[g] for g in used}), dtype=int)

    S = {"PEMT100": resolve(plist), "PEMT30": resolve(plist[:30])}
    for lab, genes in raw_sets.items():
        S[lab] = resolve(genes)
    prespecified = [k for k in S if not k.startswith("PEMT") and k not in extra_labels]
    comp = extra_labels if extra_labels else prespecified

    mean_expr, det = M.mean(axis=1), (M > 0).mean(axis=1)
    eb, s2 = make_strata(mean_expr, det)
    R = CR.zrank_within_rows(M, pats)
    ssR = (R.astype(np.float64) ** 2).sum(axis=1)
    tw, td = CR.zrank_within(D["w"], pats), CR.zrank_within(D["tech"], pats)
    dec = CR.Decomp(R, M, pats, tw, td, mode)
    allp = np.ones(len(dec.pids), bool)
    r_gw, r_gd, r_wd, r_adj, _ = dec.profile(allp)
    elig = np.where(dec.ncell >= MIN_CELLS)[0]
    PP = per_patient_partials(dec, mode)
    r_bal = fisher_mean(PP, elig)
    r_res = residual_profile(R, ssR, tw, [td], mode)

    # ---- R2-0 reproduction ------------------------------------------------------------
    rep = [{"check": "primary profile recomputed by the general partial-correlation routine",
            "max_abs_diff": float(np.nanmax(np.abs(r_adj - r_res)))}]
    if args.check_dir:
        for fn in os.listdir(args.check_dir):
            if fn.startswith("CR1") and "profile" in fn.lower() and fn.endswith((".csv", ".csv.gz")):
                old = pd.read_csv(os.path.join(args.check_dir, fn), index_col=0)
                col = [c for c in old.columns if "adj" in c.lower() and "unadj" not in c.lower()]
                if col:
                    common = [g for g in names if g in old.index]
                    a = pd.Series(r_adj, index=names).loc[common].values
                    rep.append({"check": f"primary profile vs {fn}[{col[0]}] ({len(common):,} genes)",
                                "max_abs_diff": float(np.nanmax(np.abs(a - old.loc[common, col[0]].values)))})
    pd.DataFrame(rep).to_csv(out("R2_0_reproduction.csv"), index=False)

    profiles = {"pooled (primary)": r_adj, "patient-balanced (equal-weight Fisher z)": r_bal}
    nulls = {"expression-matched (primary)": (eb, eb), "expression+detection-matched": (s2, eb)}
    if not extra_labels:
        # ---- R2-1 patient weights ----------------------------------------------------------
        wt = dec.ncell if mode == "n" else dec.ssw
        w1 = pd.DataFrame({"patient": dec.pids, "n_cells": dec.ncell, "share_of_cells": dec.ncell / dec.ncell.sum(),
                           "share_of_pooled_weight": wt / wt.sum(), "eligible_ge20": dec.ncell >= MIN_CELLS,
                           "weight_in_balanced_profile": np.where(dec.ncell >= MIN_CELLS, 1 / len(elig), 0.0)})
        for nm in ("PEMT100", "PEMT30"):
            w1[f"patient_mean_rho_{nm}"] = np.where(dec.ncell >= MIN_CELLS, np.nanmean(PP[S[nm]], axis=0), np.nan)
        w1.to_csv(out("R2_1_patient_weights.csv"), index=False)

        # ---- R2-2 null model x weighting ------------------------------------------------------
        rows = []
        for pl, rho in profiles.items():
            for nl, (st, fb) in nulls.items():
                for nm in ("PEMT100", "PEMT30"):
                    r = enrich(rho, S[nm], st, fb, np.random.default_rng(C.SEED), args.n_perm)
                    rows.append({"profile": pl, "null": nl, "signature": nm, **r})
        pd.DataFrame(rows).to_csv(out("R2_2_null_and_weighting.csv"), index=False)

        # ---- R2-3 patient-level test ---------------------------------------------------------------
        rows = []
        for nm in ("PEMT100", "PEMT30"):
            for nl, (st, fb) in nulls.items():
                d = []
                for j in elig:
                    rho = PP[:, j]
                    sidx = np.array([i for i in S[nm] if np.isfinite(rho[i])], dtype=int)
                    d.append(float(rho[sidx].mean() - stratum_expectation(rho, sidx, st, fb)))
                d = np.array(d)
                npos = int((d > 0).sum())
                rows.append({"signature": nm, "null": nl, "n_patients": len(d), "n_positive": npos,
                             "median_excess": float(np.median(d)), "min_excess": float(d.min()), "max_excess": float(d.max()),
                             "sign_test_p": float(stats.binomtest(npos, len(d), 0.5).pvalue),
                             "wilcoxon_p": float(stats.wilcoxon(d).pvalue),
                             "per_patient_excess": ";".join(f"{dec.pids[j]}:{v:.4f}" for j, v in zip(elig, d))})
        pd.DataFrame(rows).to_csv(out("R2_3_patient_level.csv"), index=False)

        # ---- R2-4 capped down-sampling -----------------------------------------------------------------
        cap = int(np.median(dec.ncell[elig]))
        rng_ds = np.random.default_rng(C.SEED)
        rows = []
        for rpt in range(args.n_rep):
            keep = np.concatenate([(lambda s: s if len(s) <= cap else np.sort(rng_ds.choice(s, cap, replace=False)))(np.where(pats == p)[0])
                                   for p in dec.pids])
            Ms, ps = M[:, keep], pats[keep]
            Rs = CR.zrank_within_rows(Ms, ps)
            tws, tds = CR.zrank_within(D["w"][keep], ps), CR.zrank_within(D["tech"][keep], ps)
            rho = residual_profile(Rs, (Rs.astype(np.float64) ** 2).sum(axis=1), tws, [tds], mode)
            ebs, s2s = make_strata(Ms.mean(axis=1), (Ms > 0).mean(axis=1))
            for nm in ("PEMT100", "PEMT30"):
                for nl, (st, fb) in {"expression-matched (primary)": (ebs, ebs), "expression+detection-matched": (s2s, ebs)}.items():
                    r = enrich(rho, S[nm], st, fb, np.random.default_rng(C.SEED + rpt), args.n_perm_rep)
                    rows.append({"repeat": rpt, "cap": cap, "n_cells": len(keep), "signature": nm, "null": nl, **r})
            if rpt % 10 == 0:
                C.log(f"  down-sampling {rpt + 1}/{args.n_rep}")
        ds = pd.DataFrame(rows)
        ds.to_csv(out("R2_4_downsample_all.csv"), index=False)
        dss = ds.groupby(["signature", "null"]).agg(
            n_repeats=("repeat", "size"), cap=("cap", "first"), n_cells=("n_cells", "first"),
            median_mean_rho=("mean_rho_set", "median"), median_null_mean=("null_mean", "median"),
            median_z=("z", "median"), min_z=("z", "min"), median_p=("perm_p_upper", "median"), max_p=("perm_p_upper", "max"),
            frac_p_lt_0_05=("perm_p_upper", lambda v: float((v < 0.05).mean())),
            frac_p_lt_0_01=("perm_p_upper", lambda v: float((v < 0.01).mean()))).reset_index()
        dss.to_csv(out("R2_4_downsample_summary.csv"), index=False)

    # ---- R2-5 overlap among programs ---------------------------------------------------------------------
    keys = ["PEMT100", "PEMT30"] + comp
    tagx = "x" if extra_labels else ""                       # extra-mode outputs get their own file names
    rows = []
    for a in keys:
        for b in keys:
            A, B = set(S[a].tolist()), set(S[b].tolist())
            rows.append({"set_a": a, "set_b": b, "n_a": len(A), "n_b": len(B), "shared": len(A & B),
                         "jaccard": len(A & B) / len(A | B) if A | B else np.nan,
                         "shared_genes": ";".join(sorted(names[i] for i in A & B)) if a < b else ""})
    pd.DataFrame(rows).to_csv(out("R2_5_program_overlap.csv".replace("R2_", "R2" + tagx + "_")), index=False)
    pd.DataFrame([{"set": k, "source": src.get(k, args.pemt), "n_listed": len(raw_sets.get(k, plist if k == "PEMT100" else plist[:30])),
                   "n_evaluable": len(S[k])} for k in keys]).to_csv(out("R2_5_program_sources.csv".replace("R2_", "R2" + tagx + "_")), index=False)

    # ---- R2-6 enrichment of each program, full and disjoint -----------------------------------------------------
    union_comp = set().union(*[set(S[k].tolist()) for k in comp]) if comp else set()
    union_all = set().union(*[set(S[k].tolist()) for k in prespecified + extra_labels]) if (prespecified or extra_labels) else set()
    pe = set(S["PEMT100"].tolist())
    variants = {"PEMT100 (full)": S["PEMT100"], "PEMT30 (full)": S["PEMT30"],
                ("PEMT100 only (in no extra gene set)" if extra_labels else "PEMT100 only (in no related program)"): np.array(sorted(pe - union_comp), dtype=int)}
    for k in comp:
        variants[f"{k} (full)"] = S[k]
        variants[f"{k} only (not in PEMT100)"] = np.array(sorted(set(S[k].tolist()) - pe), dtype=int)
        variants[f"PEMT100 and {k} (shared)"] = np.array(sorted(set(S[k].tolist()) & pe), dtype=int)
    all_prog = np.array(sorted(pe | union_all), dtype=int)       # never used as null draws
    rows = []
    for pl, rho in profiles.items():
        for nl, (st, fb) in nulls.items():
            for vn, sidx in variants.items():
                r = enrich(rho, sidx, st, fb, np.random.default_rng(C.SEED), args.n_perm, exclude=all_prog)
                rows.append({"profile": pl, "null": nl, "gene_set": vn,
                             **(r or {"n_genes": int(np.isfinite(rho[sidx]).sum()), "note": "fewer than 5 evaluable genes"})})
    pd.DataFrame(rows).to_csv(out("R2_6_program_comparison.csv".replace("R2_", "R2" + tagx + "_")), index=False)

    # ---- R2-7 conditional enrichment / R2-8 module partials ----------------------------------------------------------
    rows7, rows8 = [], []
    for k in comp:
        ck = set(S[k].tolist())
        p_only, k_only = np.array(sorted(pe - ck), dtype=int), np.array(sorted(ck - pe), dtype=int)
        if len(p_only) < 5 or len(k_only) < 5:
            continue
        m_p, m_k = module_score(R, p_only, pats), module_score(R, k_only, pats)
        for lab, covm, test in ((f"p-EMT-only genes | technical + {k}-only module", m_k, p_only),
                                (f"{k}-only genes | technical + p-EMT-only module", m_p, k_only)):
            rr = residual_profile(R, ssR, tw, [td, covm], mode)
            for nl, (st, fb) in nulls.items():
                r0 = enrich(r_adj, test, st, fb, np.random.default_rng(C.SEED), args.n_perm, exclude=all_prog)
                r1 = enrich(rr, test, st, fb, np.random.default_rng(C.SEED), args.n_perm, exclude=all_prog)
                rows7.append({"competitor": k, "test": lab, "null": nl, "n_genes": r1["n_genes"],
                              "mean_rho_before": r0["mean_rho_set"], "z_before": r0["z"], "p_before": r0["perm_p_upper"],
                              "mean_rho_after": r1["mean_rho_set"], "null_mean_after": r1["null_mean"], "z_after": r1["z"],
                              "p_after": r1["perm_p_upper"], f"top{TOPK}_after": r1[f"top{TOPK}_overlap"],
                              f"top{TOPK}_expected": r1[f"top{TOPK}_expected"], f"top{TOPK}_hypergeom_p_after": r1[f"top{TOPK}_hypergeom_p"]})
        pos_p = pos_k = 0
        for j in elig:
            s = np.where(pats == dec.pids[j])[0]
            pos_p += partial_vec(tw[s], m_p[s], [td[s], m_k[s]]) > 0
            pos_k += partial_vec(tw[s], m_k[s], [td[s], m_p[s]]) > 0
        rows8.append({"competitor": k, "n_pemt_only": len(p_only), "n_competitor_only": len(k_only),
                      "rho_modules": partial_vec(m_p, m_k, [td]),
                      "rho_wbp5_pemt_only | tech": partial_vec(tw, m_p, [td]),
                      "rho_wbp5_competitor_only | tech": partial_vec(tw, m_k, [td]),
                      "rho_wbp5_pemt_only | tech+competitor": partial_vec(tw, m_p, [td, m_k]),
                      "rho_wbp5_competitor_only | tech+pemt": partial_vec(tw, m_k, [td, m_p]),
                      "patients_positive_pemt | tech+competitor": f"{int(pos_p)}/{len(elig)}",
                      "patients_positive_competitor | tech+pemt": f"{int(pos_k)}/{len(elig)}"})
    pd.DataFrame(rows7).to_csv(out("R2_7_conditional_enrichment.csv".replace("R2_", "R2" + tagx + "_")), index=False)
    pd.DataFrame(rows8).to_csv(out("R2_8_module_partials.csv".replace("R2_", "R2" + tagx + "_")), index=False)

    # ---- console summary ---------------------------------------------------------------------------------------
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    f = lambda v: f"{v:.4g}"  # noqa: E731
    print(f"\n######## R2 [{args.tag}] {args.cohort}: {len(pats):,} cells / {len(dec.pids)} patients "
          f"({len(elig)} with >= {MIN_CELLS} cells) / {len(names):,} genes / technical covariate: {D['tech_src']} ########")
    shown = [("R2-0 reproduction", "R2_0_reproduction.csv")]
    if not extra_labels:
        shown += [("R2-1 patient weights", "R2_1_patient_weights.csv"), ("R2-2 null model x weighting", "R2_2_null_and_weighting.csv"),
                  ("R2-3 patient-level", "R2_3_patient_level.csv"), ("R2-4 capped down-sampling", "R2_4_downsample_summary.csv")]
    shown += [(f"R2{tagx}-5 program sources", f"R2{tagx}_5_program_sources.csv"), (f"R2{tagx}-6 program comparison", f"R2{tagx}_6_program_comparison.csv"),
              (f"R2{tagx}-7 conditional enrichment", f"R2{tagx}_7_conditional_enrichment.csv"), (f"R2{tagx}-8 module partials", f"R2{tagx}_8_module_partials.csv")]
    for nm, fn in shown:
        t = pd.read_csv(out(fn))
        t = t.drop(columns=[c for c in ("per_patient_excess", "shared_genes", "source") if c in t.columns])
        print(f"\n=== {nm} ===\n{t.to_string(index=False, float_format=f)}")
    import scipy
    with open(out("R2x_run_info.txt" if args.extra else "R2_run_info.txt"), "a", encoding="utf-8") as fh:
        fh.write(f"finished {pd.Timestamp.now()}\npython {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n"
                 f"scipy {scipy.__version__}\ntech {D['tech_src']}\nmode {mode}\n"
                 + "".join(f"set {k}: {src.get(k, args.pemt)} / evaluable {len(S[k])}\n" for k in keys))
    C.log(f"done -> {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
