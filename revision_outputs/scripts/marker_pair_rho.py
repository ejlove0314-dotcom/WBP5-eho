#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
marker_pair_rho.py — GSE103322 악성세포에서 마커 쌍 Spearman 상관 계산
=====================================================================

목적
----
원고 §3.2의 다음 참조 문장을 검증한다.

    "for reference, canonical p-EMT marker pairs (e.g., LAMC2-ITGA5)
     in the same cohort yielded rho-values in the 0.30-0.45 range."

LAMC2-ITGA5 한 쌍만이 아니라 패널 전체의 쌍별 상관을 계산해, 문장에
적을 실제 범위를 확정한다.

실행
----
    cd "C:\\Users\\user\\Desktop\\WBP5 revision"
    python marker_pair_rho.py

필요 패키지: pandas, numpy, scipy
입력: GSE103322_HNSCC_all_data.txt.gz  (regen_tableS6.py 실행 시 이미 받음)

출력
----
1) marker_pair_rho_<날짜>.csv  — 모든 쌍의 rho / p / n
2) 화면에 LAMC2-ITGA5 값, p-EMT 마커 쌍의 실제 범위, 문장 제안

주의
----
- 대상은 악성세포 2,215개. Table 2 및 §3.2와 동일한 모집단이다.
- 행렬 값은 이미 log2(TPM/10+1) 스케일. 추가 변환 없음.
- Spearman은 단조 변환에 불변이므로 스케일 자체는 결과에 영향을 주지 않는다.
"""

import os
import sys
import gzip
import datetime

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

FNAME = "GSE103322_HNSCC_all_data.txt.gz"
STAMP = datetime.date.today().strftime("%Y%m%d")

# Table 2의 마커 패널 + WBP5
PANEL = ["WBP5", "LAMC2", "ITGA5", "KRT14", "KRT17", "VIM", "EGFR", "EPCAM", "PTPRC"]

# 원고가 "canonical p-EMT marker pairs"로 지칭하는 대상
PEMT = ["LAMC2", "ITGA5", "KRT14", "KRT17", "VIM"]


def read_annotations():
    """상단 주석 행을 읽어 세포 바코드와 악성 플래그를 돌려준다."""
    with gzip.open(FNAME, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = [c.strip().strip('"') for c in header[1:]]
        ann = {}
        for _ in range(12):
            line = fh.readline()
            if not line:
                break
            parts = line.rstrip("\n").split("\t")
            key = parts[0].strip().strip('"')
            vals = [v.strip().strip('"') for v in parts[1:]]
            try:
                if max(float(v) for v in vals[:200]) > 1.0001:
                    break
            except ValueError:
                pass
            ann[key] = vals
    for k, v in ann.items():
        low = " ".join(k.lower().split())
        if "cancer" in low and "cell" in low and "non-cancer" not in low:
            return cells, np.array([str(x).strip() in ("1", "1.0") for x in v])
    sys.exit("ERROR: 악성 플래그 주석 행을 찾지 못했습니다. 주석 행: %s"
             % list(ann))


def read_genes(n_cells, wanted):
    """패널 유전자 행만 뽑아 DataFrame으로 돌려준다."""
    want = {w.upper(): w for w in wanted}
    out = {}
    with gzip.open(FNAME, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        for line in fh:
            tab = line.find("\t")
            name = line[:tab].strip().strip("'\"").upper()
            if name in want and want[name] not in out:
                vals = np.array([float(v) for v in
                                 line.rstrip("\n").split("\t")[1:]], dtype=float)
                if len(vals) != n_cells:
                    sys.exit("ERROR: %s 행 길이 불일치 (%d vs %d)"
                             % (name, len(vals), n_cells))
                out[want[name]] = vals
    missing = [w for w in wanted if w not in out]
    if missing:
        print("  경고: 행렬에서 찾지 못한 유전자 -> %s" % ", ".join(missing))
    return pd.DataFrame(out)


def main():
    if not os.path.exists(FNAME):
        sys.exit("ERROR: %s 가 현재 폴더에 없습니다.\n"
                 "       regen_tableS6.py 를 돌린 폴더에서 실행해 주세요." % FNAME)

    print("[1/3] 세포 주석 읽기")
    cells, is_malignant = read_annotations()
    print("      전체 세포 %d, 악성 %d (기대 2215)  %s"
          % (len(cells), is_malignant.sum(),
             "OK" if is_malignant.sum() == 2215 else "*** 불일치 ***"))

    print("[2/3] 마커 발현 행 추출")
    df = read_genes(len(cells), PANEL)
    mal = df.loc[is_malignant].reset_index(drop=True)
    print("      악성세포 행렬: %d x %d" % mal.shape)

    print("[3/3] 쌍별 Spearman 상관")
    genes = [g for g in PANEL if g in mal.columns]
    rows = []
    for i, a in enumerate(genes):
        for b in genes[i + 1:]:
            rho, p = spearmanr(mal[a], mal[b])
            rows.append({"gene_A": a, "gene_B": b,
                         "rho": round(float(rho), 4),
                         "p": float(p), "n": len(mal)})
    res = pd.DataFrame(rows).sort_values("rho", ascending=False)
    f = "marker_pair_rho_%s.csv" % STAMP
    res.to_csv(f, index=False, encoding="utf-8-sig")

    # ── 검증 대상 ────────────────────────────────────────────────────
    print("\n── 원고 문장 검증 ──")
    li = res[((res.gene_A == "LAMC2") & (res.gene_B == "ITGA5")) |
             ((res.gene_A == "ITGA5") & (res.gene_B == "LAMC2"))]
    if len(li):
        r = float(li.iloc[0].rho)
        print("  LAMC2-ITGA5 rho = %.4f" % r)
        print("  원고 서술 '0.30-0.45 range' 안에 들어가는가 : %s"
              % ("예" if 0.30 <= r <= 0.45 else "아니오  *** 수정 필요 ***"))

    sub = res[res.gene_A.isin(PEMT) & res.gene_B.isin(PEMT)]
    if len(sub):
        lo, hi = sub.rho.min(), sub.rho.max()
        print("\n  p-EMT 마커 쌍 %d개의 실제 범위 : %.3f ~ %.3f" % (len(sub), lo, hi))
        print("  (중앙값 %.3f)" % sub.rho.median())
        print("\n  쌍별 상세:")
        for _, r in sub.iterrows():
            print("    %-6s - %-6s  rho = %7.4f" % (r.gene_A, r.gene_B, r.rho))
        print("\n  문장 제안: \"...in the same cohort yielded \u03c1-values in the "
              "%.2f\u2013%.2f range.\"" % (np.floor(lo * 20) / 20,
                                          np.ceil(hi * 20) / 20))

    print("\n── WBP5 와 각 마커 (Table 2 대조용) ──")
    w = res[(res.gene_A == "WBP5") | (res.gene_B == "WBP5")]
    for _, r in w.iterrows():
        other = r.gene_B if r.gene_A == "WBP5" else r.gene_A
        print("    WBP5 - %-6s  rho = %7.4f" % (other, r.rho))

    print("\n저장: %s" % f)


if __name__ == "__main__":
    main()
