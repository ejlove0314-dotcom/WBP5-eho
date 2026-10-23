#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
regen_tableS6.py — Table S6 재생성 (GSE103322 원본 주석 기준)
=============================================================

목적
----
SI Table S6의 `n cells`, `Mean log2(TPM/10+1)`, `PCT+ (%)` 세 열을
GSE103322 원본 발현 행렬에서 직접 재계산한다.

현재 Table S6은 세 열이 서로 다른 분석에서 나온 것으로 의심되는 상태다
(PCT+는 11행 중 9행이 맞는데 n과 mean은 대부분 어긋남). 이 스크립트는
세 열을 하나의 동일한 계산에서 산출하여 그 문제를 제거한다.

실행
----
    cd "C:\\Users\\user\\Desktop\\WBP5 revision"
    python regen_tableS6.py

필요 패키지: pandas, numpy  (둘 다 이미 설치돼 있을 것)

입력
----
GSE103322_HNSCC_all_data.txt.gz
  - 같은 폴더에 있으면 그대로 쓴다.
  - 없으면 GEO에서 자동으로 내려받는다(약 90 MB).

출력
----
1) TableS6_regenerated_<날짜>.csv   ← 새 Table S6 (본문에 넣을 값)
2) TableS6_comparison_<날짜>.csv    ← 현재 Table S6 선언값과의 대조표
3) 화면에 요약 출력 + 자가검증 결과

주의
----
행렬 값은 이미 log2(TPM/10+1) 스케일이다. 추가 변환을 하지 않는다.
"""

import os
import sys
import gzip
import datetime
import urllib.request

import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────
# 설정
# ─────────────────────────────────────────────────────────────────────
FNAME = "GSE103322_HNSCC_all_data.txt.gz"
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/"
       "GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz")

GENE = "WBP5"          # 대상 유전자. 행렬에 없으면 TCEAL9로 자동 재시도
STAMP = datetime.date.today().strftime("%Y%m%d")

# 현재 SI Table S6에 적혀 있는 값 (대조용). 바꾸지 말 것.
CURRENT = {
    "Cancer cells":       (2215, 2.38, 68.89),
    "Fibroblast":         (1440, 1.61, 51.97),
    "Endothelial":        (280,  1.48, 51.15),
    "Fibroblast (excl.)": (370,  1.12, 38.92),
    "Unannotated":        (135,  0.95, 34.07),
    "Macrophage":         (220,  0.61, 23.64),
    "Myocyte":            (25,   0.41, 16.00),
    "Mast":               (55,   0.38, 14.55),
    "B cell":             (145,  0.29, 11.72),
    "Dendritic":          (190,  0.15,  3.92),
    "T cell":             (827,  0.04,  1.21),
}


# ─────────────────────────────────────────────────────────────────────
# 1. 입력 파일 확보
# ─────────────────────────────────────────────────────────────────────
def ensure_matrix():
    if os.path.exists(FNAME):
        print("[1/5] 입력 파일 확인: %s (%.1f MB)"
              % (FNAME, os.path.getsize(FNAME) / 1e6))
        return
    print("[1/5] 파일이 없어 GEO에서 내려받습니다 (약 90 MB, 수 분 소요)")
    urllib.request.urlretrieve(URL, FNAME)
    print("      완료: %.1f MB" % (os.path.getsize(FNAME) / 1e6))


# ─────────────────────────────────────────────────────────────────────
# 2. 헤더(세포 주석) 파싱
# ─────────────────────────────────────────────────────────────────────
# GSE103322의 상단 5행은 유전자가 아니라 세포별 주석이다.
#   row0 : 세포 바코드 (열 이름)
#   이후 : 'Maxima enzyme', 'classified  as cancer cell',
#          'classified as non-cancer cells', 'non-cancer cell type' 등
# 주석 행은 첫 칸에 그 이름이 들어 있고, 값이 숫자가 아니거나 0/1이다.
def read_annotations():
    print("[2/5] 세포 주석 행 파싱")
    ann = {}
    with gzip.open(FNAME, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = [c.strip().strip('"') for c in header[1:]]
        for _ in range(12):                       # 상단 12행 안에 주석이 다 있다
            line = fh.readline()
            if not line:
                break
            parts = line.rstrip("\n").split("\t")
            key = parts[0].strip().strip('"')
            vals = [v.strip().strip('"') for v in parts[1:]]
            # 값이 전부 숫자이고 0/1을 벗어나면 유전자 행으로 간주하고 중단
            try:
                nums = [float(v) for v in vals[:200]]
                if max(nums) > 1.0001:
                    break
            except ValueError:
                pass
            ann[key] = vals
    print("      세포 수: %d" % len(cells))
    print("      주석 행: %s" % ", ".join(repr(k) for k in ann))
    return cells, ann


def pick(ann, *keywords):
    """주석 행 이름에 keywords가 모두 들어간 행을 찾는다 (공백/대소문자 무시)."""
    for k, v in ann.items():
        low = " ".join(k.lower().split())
        if all(w in low for w in keywords):
            return k, v
    return None, None


# ─────────────────────────────────────────────────────────────────────
# 3. 구획(compartment) 라벨 생성
# ─────────────────────────────────────────────────────────────────────
def build_compartments(cells, ann):
    print("[3/5] 구획 라벨 생성")
    k_can, cancer = pick(ann, "cancer", "cell")          # 'classified as cancer cell'
    k_typ, ctype = pick(ann, "non-cancer", "type")       # 'non-cancer cell type'
    if cancer is None or ctype is None:
        sys.exit("ERROR: 주석 행을 찾지 못했습니다. 위에 출력된 주석 행 이름을 "
                 "확인해 pick() 키워드를 맞춰 주세요.")
    print("      악성 플래그 행 : %r" % k_can)
    print("      세포유형 행   : %r" % k_typ)

    # 원본 라벨 → Table S6 표기로 정규화
    NORM = {
        "fibroblast":  "Fibroblast",
        "-fibroblast": "Fibroblast (excl.)",   # 원 분석에서 제외된 섬유아세포
        "endothelial": "Endothelial",
        "macrophage":  "Macrophage",
        "mast":        "Mast",
        "b cell":      "B cell",
        "dendritic":   "Dendritic",
        "t cell":      "T cell",
        "myocyte":     "Myocyte",
    }

    labels, raw_seen = [], {}
    for i in range(len(cells)):
        is_cancer = str(cancer[i]).strip() in ("1", "1.0")
        if is_cancer:
            labels.append("Cancer cells")
            continue
        raw = str(ctype[i]).strip()
        raw_seen[raw] = raw_seen.get(raw, 0) + 1
        key = raw.lower()
        labels.append(NORM.get(key, "Unannotated" if raw in ("", "0", "nan") else raw))

    print("      원본 non-cancer cell type 값의 분포:")
    for k, v in sorted(raw_seen.items(), key=lambda x: -x[1]):
        print("        %-24r %5d" % (k, v))
    return pd.Series(labels, index=cells, name="compartment")


# ─────────────────────────────────────────────────────────────────────
# 4. 대상 유전자 행 추출
# ─────────────────────────────────────────────────────────────────────
def read_gene_row(n_cells, gene=GENE):
    print("[4/5] %s 발현 행 추출" % gene)
    targets = [gene, gene.strip("'\""), "TCEAL9"]
    with gzip.open(FNAME, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()                       # 헤더
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            name = parts[0].strip().strip("'\"")
            if name in targets:
                vals = np.array([float(v) for v in parts[1:]], dtype=float)
                if len(vals) != n_cells:
                    sys.exit("ERROR: 유전자 행 길이(%d)가 세포 수(%d)와 다릅니다."
                             % (len(vals), n_cells))
                print("      찾음: %r  (비영 세포 %d / %d)"
                      % (name, int((vals > 0).sum()), n_cells))
                return vals
    sys.exit("ERROR: %s 행을 행렬에서 찾지 못했습니다." % targets)


# ─────────────────────────────────────────────────────────────────────
# 5. 집계 · 대조 · 저장
# ─────────────────────────────────────────────────────────────────────
def main():
    ensure_matrix()
    cells, ann = read_annotations()
    comp = build_compartments(cells, ann)
    expr = pd.Series(read_gene_row(len(cells)), index=cells, name="wbp5")

    print("[5/5] 구획별 집계")
    df = pd.DataFrame({"compartment": comp.values, "wbp5": expr.values})
    out = (df.groupby("compartment")["wbp5"]
             .agg(n="size",
                  mean="mean",
                  pct=lambda s: 100.0 * (s > 0).sum() / len(s))
             .reset_index()
             .sort_values("pct", ascending=False)
             .reset_index(drop=True))
    out["mean"] = out["mean"].round(3)
    out["pct"] = out["pct"].round(2)

    # ── 자가검증 ────────────────────────────────────────────────────
    print("\n── 자가검증 ──")
    total = int(out["n"].sum())
    print("  총 세포 수            : %d   (기대 5902)  %s"
          % (total, "OK" if total == 5902 else "*** 불일치 ***"))
    can = out.loc[out.compartment == "Cancer cells", "n"]
    can = int(can.iloc[0]) if len(can) else -1
    print("  Cancer cells n        : %d   (기대 2215)  %s"
          % (can, "OK" if can == 2215 else "*** 불일치 ***"))
    cm = out.loc[out.compartment == "Cancer cells", "mean"]
    cp = out.loc[out.compartment == "Cancer cells", "pct"]
    if len(cm):
        print("  Cancer cells mean/PCT+: %.3f / %.2f  (기대 2.38 / 68.89)"
              % (cm.iloc[0], cp.iloc[0]))

    # ── 현재 Table S6과 대조 ────────────────────────────────────────
    rows = []
    for _, r in out.iterrows():
        cur = CURRENT.get(r.compartment, (None, None, None))
        rows.append({
            "compartment": r.compartment,
            "n_new": r.n,        "n_TableS6": cur[0],
            "mean_new": r["mean"], "mean_TableS6": cur[1],
            "pct_new": r.pct,    "pct_TableS6": cur[2],
            "n_match": (cur[0] == r.n) if cur[0] is not None else None,
            "mean_match": (abs(cur[1] - r["mean"]) < 0.005) if cur[1] is not None else None,
            "pct_match": (abs(cur[2] - r.pct) < 0.05) if cur[2] is not None else None,
        })
    cmp_df = pd.DataFrame(rows)

    f1 = "TableS6_regenerated_%s.csv" % STAMP
    f2 = "TableS6_comparison_%s.csv" % STAMP
    out.to_csv(f1, index=False, encoding="utf-8-sig")
    cmp_df.to_csv(f2, index=False, encoding="utf-8-sig")

    print("\n── 새 Table S6 ──")
    print(out.to_string(index=False))
    print("\n── 현재 Table S6과의 대조 ──")
    print(cmp_df[["compartment", "n_new", "n_TableS6", "mean_new",
                  "mean_TableS6", "pct_new", "pct_TableS6"]].to_string(index=False))
    print("\n저장: %s\n      %s" % (f1, f2))


if __name__ == "__main__":
    main()
