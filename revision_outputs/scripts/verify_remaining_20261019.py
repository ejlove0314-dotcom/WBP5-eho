#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
응답서 미검증 수치 마감 스크립트  (2026-10-19)
==============================================
목적
  J4_PAE_summary.csv 로 구조 수치 대부분이 검증되었으나, 다음 다섯 항목이 남았습니다.
    (1) 구간 간 PAE 최솟값 12.00          (J4 에는 최댓값만 기록됨)
    (2) metapredict 구간 42-57 / 72-96 평균 0.283 / 0.115
    (3) metapredict 무질서(>0.5) 잔기 비율 100% / 17% / 19% / 0% / 0%
    (4) GSE103322 선형 스케일 평균비 0.85
    (5) GSE103322 섬유아세포 선형 평균 23.5 / 중앙값 0.52

  동시에 J4 에 이미 있는 값(블록 평균·중앙값·최댓값, 분율, pLDDT, metapredict 평균)도
  원본 파일에서 다시 계산해 대조합니다. J4 를 '믿고 쓰는' 것이 아니라 '확인하고 쓰기'
  위함입니다.

계산 규약 (J4 역산으로 확정)
  구간 내/간 PAE 통계는 대각선을 포함한 전체 블록 PAE[A][:, B] 에 대해 계산합니다.
    42-57  : 16 x 16 = 256 쌍  -> 3.0078125 x 256 = 770.0   (정수, 규약 일치)
    72-96  : 25 x 25 = 625 쌍  -> 3.6384    x 625 = 2274.0  (정수, 규약 일치)
    N(1-28) vs core(29-104) : 28 x 76 = 2128 쌍 -> 28.5197 x 2128 = 60690.0 (정수)
  임계 초과 분율은 한쪽 방향(A->B)과 양방향(A->B + B->A) 두 가지로 모두 보고합니다.

실행 (파트 A: 구조. PAE json 만 있으면 됩니다)
  python verify_remaining_20261019.py ^
     --pae AF-Q9UHQ7-F1-predicted_aligned_error_v6.json ^
     --pdb AF-Q9UHQ7-F1-model_v6.pdb ^
     --outdir results_verify

실행 (파트 B: GSE103322 선형 스케일. 매트릭스가 필요합니다)
  python verify_remaining_20261019.py ^
     --matrix GSE103322_HNSCC_all_data.txt ^
     --outdir results_verify

  둘을 한 번에 돌려도 됩니다. 파트 A 만 먼저 돌리는 편을 권합니다 (수 초).
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd

# --- 아미노산 3문자 -> 1문자 -------------------------------------------------
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

# 응답서에서 쓰는 구간 (1-based, 양끝 포함)
SEGMENTS = {
    "Nterm_1_28": (1, 28),
    "core_29_104": (29, 104),
    "pocket_38_58": (38, 58),
    "seg_42_57": (42, 57),
    "seg_72_96": (72, 96),
}

# 응답서에 적힌 값 (대조용). None 이면 아직 응답서에 없는 값.
EXPECTED = {
    "PAE_block_mean|seg_42_57|seg_42_57": 3.0078125,
    "PAE_block_mean|seg_72_96|seg_72_96": 3.6384,
    "PAE_block_mean|seg_42_57|seg_72_96": 23.17875,
    "PAE_block_mean|Nterm_1_28|core_29_104": 28.519736842105264,
    "PAE_block_median|seg_42_57|seg_42_57": 2.0,
    "PAE_block_median|seg_72_96|seg_72_96": 2.0,
    "PAE_block_median|seg_42_57|seg_72_96": 24.0,
    "PAE_block_median|Nterm_1_28|core_29_104": 29.0,
    "PAE_block_max|seg_42_57|seg_42_57": 12.0,
    "PAE_block_max|seg_72_96|seg_72_96": 18.0,
    "PAE_block_max|seg_42_57|seg_72_96": 29.0,
    "PAE_block_max|Nterm_1_28|core_29_104": 31.0,
    "PAE_block_min|seg_42_57|seg_72_96": 12.0,          # (1) 미검증
    "max_predicted_aligned_error": 31.75,
    "PAE_pair|45|69": 18.5,
    "PAE_pair|45|72": 21.0,
    "PAE_pair|45|48": 1.5,
    "PAE_pair|69|72": 4.0,
    "PAE_contactblock_41_48_vs_69_72_mean": 21.97,
    "frac_gt20_twoway|seg_42_57|seg_72_96": 0.771,
    "frac_gt15_twoway|seg_42_57|seg_72_96": 0.948,
    "pLDDT_mean|pocket_38_58": 79.09,
    "pLDDT_min|pocket_38_58": 60.88,
    "pLDDT_max|pocket_38_58": 91.00,
    "pLDDT_mean|Nterm_1_28": 44.10,
    "pLDDT_mean|core_29_104": 75.45,
    "pLDDT_mean|full_1_104": 67.01,
    "metapredict_mean|Nterm_1_28": 0.970,
    "metapredict_mean|core_29_104": 0.263,
    "metapredict_mean|pocket_38_58": 0.330,
    "metapredict_mean|seg_42_57": 0.283,                # (2) 미검증
    "metapredict_mean|seg_72_96": 0.115,                # (2) 미검증
    "metapredict_fracgt0.5|Nterm_1_28": 1.00,           # (3) 미검증
    "metapredict_fracgt0.5|core_29_104": 0.17,          # (3) 미검증
    "metapredict_fracgt0.5|pocket_38_58": 0.19,         # (3) 미검증
    "metapredict_fracgt0.5|seg_42_57": 0.00,            # (3) 미검증
    "metapredict_fracgt0.5|seg_72_96": 0.00,            # (3) 미검증
    "metapredict_IDR_boundary": 41,
    "GSE103322_log2_mean_ratio_malig_over_nonmalig": 1.80,
    "GSE103322_linear_mean_ratio_malig_over_nonmalig": 0.85,   # (4) 미검증
    "GSE103322_linear_mean|Fibroblast": 23.5,                  # (5) 미검증
    "GSE103322_linear_median|Fibroblast": 0.52,                # (5) 미검증
}

rows = []


def rec(key, value, note=""):
    """결과 한 줄을 기록하고 응답서 값과 대조한다."""
    exp = EXPECTED.get(key)
    if exp is None or value is None:
        status = "" if exp is None else "계산실패"
        delta = None
    else:
        delta = float(value) - float(exp)
        # 응답서 표기 자릿수를 고려한 허용오차
        tol = 0.5 * 10 ** (-_decimals(exp))
        status = "일치" if abs(delta) <= max(tol, 1e-9) else "불일치"
    rows.append({"quantity": key, "recomputed": value, "in_letter": exp,
                 "diff": delta, "status": status, "note": note})
    tag = f"[{status}]" if status else "[신규]"
    e = "" if exp is None else f"  (응답서 {exp})"
    print(f"  {tag:<8} {key:<52} {value}{e}  {note}")


def _decimals(x):
    s = f"{float(x):.10f}".rstrip("0")
    return len(s.split(".")[1]) if "." in s else 0


def log(m):
    print(f"[{pd.Timestamp.now():%H:%M:%S}] {m}", flush=True)


# =============================================================================
# 파트 A-1 : PAE
# =============================================================================
def load_pae(path):
    """AFDB v6(행렬형)과 구형(residue_1/residue_2/distance) 양쪽을 읽는다."""
    with open(path, "r") as fh:
        obj = json.load(fh)
    if isinstance(obj, list):
        obj = obj[0]
    cap = obj.get("max_predicted_aligned_error")
    if "predicted_aligned_error" in obj:
        M = np.asarray(obj["predicted_aligned_error"], dtype=float)
    elif "pae" in obj:
        M = np.asarray(obj["pae"], dtype=float)
    elif "distance" in obj:
        r1 = np.asarray(obj["residue_1"], dtype=int)
        r2 = np.asarray(obj["residue_2"], dtype=int)
        d = np.asarray(obj["distance"], dtype=float)
        n = int(max(r1.max(), r2.max()))
        M = np.zeros((n, n))
        M[r1 - 1, r2 - 1] = d
    else:
        raise SystemExit(f"오류: PAE 키를 찾을 수 없습니다. 키 목록: {list(obj)[:10]}")
    return M, cap


def block(M, a, b):
    """1-based 폐구간 a=(s,e), b=(s,e) 의 PAE[A][:, B] 블록 (대각선 포함)."""
    return M[a[0] - 1:a[1], b[0] - 1:b[1]]


def part_a_pae(pae_path):
    log(f"PAE 적재: {pae_path}")
    M, cap = load_pae(pae_path)
    n = M.shape[0]
    print(f"  PAE 행렬 {M.shape}, 비대칭도 max|M - M.T| = {np.abs(M - M.T).max():.4f}")
    if n != 104:
        print(f"  주의: 잔기 수가 104 가 아닙니다 ({n}). 구간 정의를 확인하십시오.")
    rec("max_predicted_aligned_error", float(cap) if cap is not None else None)

    pairs = [("seg_42_57", "seg_42_57"), ("seg_72_96", "seg_72_96"),
             ("seg_42_57", "seg_72_96"), ("Nterm_1_28", "core_29_104")]
    for A, B in pairs:
        blk = block(M, SEGMENTS[A], SEGMENTS[B])
        rec(f"PAE_block_mean|{A}|{B}", round(float(blk.mean()), 10),
            f"{blk.shape[0]}x{blk.shape[1]}={blk.size}쌍")
        rec(f"PAE_block_median|{A}|{B}", float(np.median(blk)))
        rec(f"PAE_block_min|{A}|{B}", float(blk.min()))
        rec(f"PAE_block_max|{A}|{B}", float(blk.max()))

    # 임계 초과 분율 — 한 방향 / 양방향 모두
    A, B = "seg_42_57", "seg_72_96"
    ab = block(M, SEGMENTS[A], SEGMENTS[B]).ravel()
    ba = block(M, SEGMENTS[B], SEGMENTS[A]).ravel()
    two = np.concatenate([ab, ba])
    for thr in (20, 15):
        rec(f"frac_gt{thr}_oneway|{A}|{B}", round(float((ab > thr).mean()), 6),
            f"{int((ab > thr).sum())}/{ab.size}")
        rec(f"frac_gt{thr}_twoway|{A}|{B}", round(float((two > thr).mean()), 6),
            f"{int((two > thr).sum())}/{two.size}")

    # 개별 잔기쌍 — 비대칭이므로 양방향 모두 기록
    for i, j in [(45, 69), (45, 72), (45, 48), (69, 72)]:
        v_ij, v_ji = float(M[i - 1, j - 1]), float(M[j - 1, i - 1])
        rec(f"PAE_pair|{i}|{j}", v_ij, f"역방향 M[{j},{i}] = {v_ji}")

    cb = block(M, (41, 48), (69, 69))
    cb2 = block(M, (41, 48), (72, 72))
    both = np.concatenate([cb.ravel(), cb2.ravel()])
    rec("PAE_contactblock_41_48_vs_69_72_mean", round(float(both.mean()), 6),
        f"{both.size}쌍 (41-48 x {{69,72}})")
    return M


# =============================================================================
# 파트 A-2 : pLDDT + 서열 (PDB)
# =============================================================================
def part_a_pdb(pdb_path):
    log(f"PDB 적재: {pdb_path}")
    plddt, seq = {}, {}
    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                r = int(line[22:26])
                plddt[r] = float(line[60:66])
                seq[r] = AA3.get(line[17:20].strip(), "X")
    idx = np.array(sorted(plddt))
    vals = np.array([plddt[i] for i in idx])
    sequence = "".join(seq[i] for i in idx)
    print(f"  잔기 {len(idx)}개, 서열 길이 {len(sequence)}")

    def sub(a, b):
        m = (idx >= a) & (idx <= b)
        return vals[m]

    for nm, (a, b) in SEGMENTS.items():
        v = sub(a, b)
        rec(f"pLDDT_mean|{nm}", round(float(v.mean()), 2))
        if nm == "pocket_38_58":
            rec(f"pLDDT_min|{nm}", round(float(v.min()), 2))
            rec(f"pLDDT_max|{nm}", round(float(v.max()), 2))
    rec("pLDDT_mean|full_1_104", round(float(vals.mean()), 2))
    return sequence


# =============================================================================
# 파트 A-3 : metapredict
# =============================================================================
def part_a_metapredict(sequence, outdir):
    try:
        import metapredict as meta
    except ImportError:
        print("  [건너뜀] metapredict 미설치. 설치 후 다시 실행하십시오:")
        print("           pip install metapredict")
        for nm in SEGMENTS:
            rec(f"metapredict_mean|{nm}", None, "metapredict 미설치")
            rec(f"metapredict_fracgt0.5|{nm}", None, "metapredict 미설치")
        rec("metapredict_IDR_boundary", None, "metapredict 미설치")
        return
    ver = getattr(meta, "__version__", "unknown")
    log(f"metapredict {ver} 실행 (서열 {len(sequence)} aa)")
    d = np.asarray(meta.predict_disorder(sequence), dtype=float)
    if len(d) != len(sequence):
        print(f"  주의: 점수 길이 {len(d)} != 서열 길이 {len(sequence)}")
    pos = np.arange(1, len(d) + 1)
    pd.DataFrame({"residue": pos, "aa": list(sequence)[:len(d)],
                  "disorder": d}).to_csv(
        os.path.join(outdir, "V_metapredict_percell.csv"), index=False)

    for nm, (a, b) in SEGMENTS.items():
        m = (pos >= a) & (pos <= b)
        rec(f"metapredict_mean|{nm}", round(float(d[m].mean()), 3))
        rec(f"metapredict_fracgt0.5|{nm}", round(float((d[m] > 0.5).mean()), 4),
            f"{int((d[m] > 0.5).sum())}/{int(m.sum())} 잔기")

    # IDR 경계: 앞에서부터 0.5 이상이 연속으로 유지되다가 처음 끊기는 지점
    above = d > 0.5
    boundary = None
    for i in range(len(above)):
        if not above[i]:
            boundary = i  # 1-based 로는 i (마지막 무질서 잔기), 경계 표기는 i+1 직전
            break
    rec("metapredict_IDR_boundary", boundary,
        "앞쪽 연속 무질서 구간이 끝나는 1-based 잔기 번호")

    # metapredict 자체 IDR 호출도 함께 기록 (참고용)
    try:
        idrs = meta.predict_disorder_domains(sequence).disordered_domain_boundaries
        rec("metapredict_called_IDR_domains", str(idrs),
            "metapredict predict_disorder_domains (0-based, 반열린구간)")
    except Exception as ex:  # noqa
        print(f"  (predict_disorder_domains 생략: {ex})")


# =============================================================================
# 파트 B : GSE103322 선형 스케일
# =============================================================================
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
    name = str(cell_name).strip().strip("'\"")
    if re.match(r"^HNSCC[_]?combo", name, flags=re.I):
        return None
    m = re.match(r"^HN(?:SCC)?[_]?(\d+)", name, flags=re.I)
    return f"HN{int(m.group(1))}" if m else None


def part_b(matrix_path, outdir):
    log(f"GSE103322 적재: {matrix_path}  (수 분 걸립니다)")
    df = pd.read_csv(matrix_path, sep="\t", index_col=0, low_memory=False)
    df.index = [str(i).strip().strip("'\"").strip() for i in df.index]
    df.columns = [str(c).strip().strip("'\"").strip() for c in df.columns]
    found = {}
    for idx in df.index[:20]:
        k = norm_key(idx)
        if k in META_ROW_KEYS:
            found[META_ROW_KEYS[k]] = idx
    if "is_cancer" not in found:
        raise SystemExit(f"오류: 'classified as cancer cell' 행 없음. 상단: {list(df.index[:8])}")
    meta = pd.DataFrame({n: df.loc[i] for n, i in found.items()})
    meta["patient"] = [parse_patient(c) for c in meta.index]
    keep = meta.index[meta["patient"].notna()]
    n_combo = len(meta) - len(keep)
    meta = meta.loc[keep]
    meta["is_cancer"] = pd.to_numeric(meta["is_cancer"], errors="coerce").fillna(0).astype(int)

    # WBP5 / TCEAL9 행만 숫자로 변환 (전체 변환은 불필요하게 무겁습니다)
    expr_index = [g for g in df.index if g not in found.values()]
    target = None
    for alias in ("WBP5", "TCEAL9"):
        if alias in expr_index:
            target = alias
            break
    if target is None:
        raise SystemExit("오류: WBP5/TCEAL9 행을 찾을 수 없습니다.")
    v_log = pd.to_numeric(df.loc[target, keep], errors="coerce").astype(float)
    del df
    log(f"  세포 {len(v_log):,} (combo 플레이트 {n_combo}개 제외), 유전자 행 '{target}'")

    # 원 데이터는 log2(TPM/10 + 1) → 선형 TPM = (2^x - 1) * 10
    v_lin = (np.power(2.0, v_log.values) - 1.0) * 10.0
    is_c = meta["is_cancer"].values == 1
    ntype = meta.get("noncancer_type", pd.Series("", index=meta.index)).astype(str).values
    is_nc = ~is_c

    def block_stats(label, mask):
        if mask.sum() == 0:
            return
        lg, ln = v_log.values[mask], v_lin[mask]
        rec(f"GSE103322_n|{label}", int(mask.sum()))
        rec(f"GSE103322_log2_mean|{label}", round(float(lg.mean()), 4))
        rec(f"GSE103322_log2_median|{label}", round(float(np.median(lg)), 4))
        rec(f"GSE103322_linear_mean|{label}", round(float(ln.mean()), 4))
        rec(f"GSE103322_linear_median|{label}", round(float(np.median(ln)), 4))
        rec(f"GSE103322_detrate|{label}", round(float((lg > 0).mean()), 4))

    block_stats("Malignant", is_c)
    block_stats("NonMalignant_pooled", is_nc)
    for t in sorted(set(ntype[is_nc])):
        if t and t.lower() not in ("nan", "0", "none", ""):
            block_stats(t, is_nc & (ntype == t))

    lg_c, lg_n = v_log.values[is_c].mean(), v_log.values[is_nc].mean()
    ln_c, ln_n = v_lin[is_c].mean(), v_lin[is_nc].mean()
    rec("GSE103322_log2_mean_ratio_malig_over_nonmalig", round(float(lg_c / lg_n), 4),
        "로그 스케일 평균의 비 — 원고에 쓰인(부적절한) 수치")
    rec("GSE103322_linear_mean_ratio_malig_over_nonmalig", round(float(ln_c / ln_n), 4),
        "선형 TPM 평균의 비 — 응답서에 쓰는 수치")

    pd.DataFrame({"cell": v_log.index, "patient": meta["patient"].values,
                  "is_cancer": is_c.astype(int), "noncancer_type": ntype,
                  "log2_tpm10p1": v_log.values, "linear_tpm": v_lin}).to_csv(
        os.path.join(outdir, "V_GSE103322_WBP5_percell.csv.gz"), index=False)


# =============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pae", default=None, help="AF-Q9UHQ7-F1-predicted_aligned_error_v6.json")
    ap.add_argument("--pdb", default=None, help="AF-Q9UHQ7-F1-model_v6.pdb")
    ap.add_argument("--matrix", default=None, help="GSE103322_HNSCC_all_data.txt")
    ap.add_argument("--outdir", default="results_verify")
    args = ap.parse_args()
    if not any([args.pae, args.pdb, args.matrix]):
        ap.error("--pae / --pdb / --matrix 중 적어도 하나는 필요합니다.")
    os.makedirs(args.outdir, exist_ok=True)

    print("\n=== 응답서 미검증 수치 재계산 ===")
    print("    [일치]=응답서와 같음  [불일치]=다름(확인 필요)  [신규]=응답서에 없는 참고값\n")

    if args.pae:
        print("--- PAE ---")
        part_a_pae(args.pae)
    sequence = None
    if args.pdb:
        print("\n--- pLDDT / 서열 ---")
        sequence = part_a_pdb(args.pdb)
    if sequence:
        print("\n--- metapredict ---")
        part_a_metapredict(sequence, args.outdir)
    elif args.pae:
        print("\n--- metapredict ---\n  [건너뜀] --pdb 가 없어 서열을 얻을 수 없습니다.")
    if args.matrix:
        print("\n--- GSE103322 선형 스케일 ---")
        part_b(args.matrix, args.outdir)

    res = pd.DataFrame(rows)
    path = os.path.join(args.outdir, "V_verification_20261019.csv")
    res.to_csv(path, index=False)

    n_ok = int((res.status == "일치").sum())
    n_bad = int((res.status == "불일치").sum())
    n_fail = int((res.status == "계산실패").sum())
    n_new = int((res.status == "").sum())
    print(f"\n=== 요약 ===")
    print(f"  일치 {n_ok} / 불일치 {n_bad} / 계산실패 {n_fail} / 신규(참고) {n_new}")
    if n_bad:
        print("\n  [불일치 항목 — 응답서에서 수정해야 합니다]")
        for _, r in res[res.status == "불일치"].iterrows():
            print(f"    {r.quantity:<52} 재계산 {r.recomputed}  vs  응답서 {r.in_letter}"
                  f"  (차 {r['diff']:+.6g})")
    print(f"\n  저장: {os.path.abspath(path)}")


if __name__ == "__main__":
    main()
