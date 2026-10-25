#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figS7_PAE.py — Fig. S7 (PAE 매트릭스) 생성
================================================

SI Fig. S7의 자리표시자를 채울 그림을 만든다.
캡션이 약속하는 내용을 그대로 그린다: 전체 PAE 매트릭스 + 두 고신뢰 세그먼트
(42-57, 72-96)의 대각 블록 + 그 사이 비대각 블록 + P0 포켓 잔기 표시.

준비
----
AlphaFold DB에서 PAE JSON 파일을 내려받아 이 스크립트와 같은 폴더에 둔다.

  1) https://alphafold.ebi.ac.uk/entry/Q9UHQ7 접속
  2) Downloads 항목에서 "Predicted aligned error (JSON)" 클릭
  3) 받은 파일을 이 스크립트와 같은 폴더에 저장

  직접 주소로 받으려면 (브라우저 주소창에 붙여넣기):
    https://alphafold.ebi.ac.uk/files/AF-Q9UHQ7-F1-predicted_aligned_error_v6.json
  v6가 없다면 끝의 v6를 v4로 바꿔 시도한다.

실행
----
    python make_figS7_PAE.py

    # 파일명이 다르면 인자로 지정
    python make_figS7_PAE.py 내려받은파일.json

필요 패키지: numpy, matplotlib  (pip install numpy matplotlib)

출력
----
  Fig_S7_PAE_matrix.pdf   벡터 (논문 제출용)
  Fig_S7_PAE_matrix.png   300 dpi (확인용)
  화면에 검증 수치 출력 -> 응답서의 23.18 / 3.0 / 3.6 Å 과 대조할 것
"""

import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"]="sans-serif"
matplotlib.rcParams["font.sans-serif"]=["Liberation Sans","Arial","DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"]=42
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ─────────────────────────────────────────────────────────────────────
# 원고/응답서가 사용하는 구간 정의 (1-기반, 양끝 포함)
# ─────────────────────────────────────────────────────────────────────
SEG1 = (42, 57)          # 고신뢰 세그먼트 1
SEG2 = (72, 96)          # 고신뢰 세그먼트 2
P0 = (38, 58)            # P0 영역
POCKET = [38, 41, 42, 44, 45, 48, 49, 51, 52, 54, 55, 58]   # 포켓 구성 12잔기


def load_pae(path):
    """AlphaFold DB의 PAE JSON을 2차원 배열로 읽는다. 두 포맷을 모두 처리."""
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)
    obj = raw[0] if isinstance(raw, list) else raw

    # 신형: {"predicted_aligned_error": [[...], ...]}
    for key in ("predicted_aligned_error", "pae"):
        if key in obj:
            return np.asarray(obj[key], dtype=float)

    # 구형: {"residue1": [...], "residue2": [...], "distance": [...]}
    if "distance" in obj and "residue1" in obj:
        r1 = np.asarray(obj["residue1"], dtype=int)
        r2 = np.asarray(obj["residue2"], dtype=int)
        d = np.asarray(obj["distance"], dtype=float)
        n = int(max(r1.max(), r2.max()))
        m = np.zeros((n, n), dtype=float)
        m[r1 - 1, r2 - 1] = d
        return m

    sys.exit("ERROR: PAE 배열을 찾지 못했습니다. JSON 키: %s" % list(obj)[:8])


def block(mat, a, b):
    """1-기반 구간 a, b 의 부분행렬 (대각 포함)."""
    return mat[a[0] - 1:a[1], b[0] - 1:b[1]]


def main():
    # ---- 입력 파일 찾기 ------------------------------------------------
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        cand = sorted(glob.glob("*predicted_aligned_error*.json")) or \
               sorted(glob.glob("*pae*.json")) or sorted(glob.glob("*.json"))
        if not cand:
            sys.exit("ERROR: PAE JSON 파일이 이 폴더에 없습니다.\n"
                     "       스크립트 상단의 안내대로 내려받아 주세요.")
        path = cand[0]
    print("[1/3] 읽는 중: %s" % path)

    pae = load_pae(path)
    n = pae.shape[0]
    print("      매트릭스 %d x %d  (기대 104 x 104)" % pae.shape)
    if n != 104:
        print("      경고: 잔기 수가 104가 아닙니다. 구간 표시가 어긋날 수 있습니다.")

    # ---- 검증 수치 ------------------------------------------------------
    # PAE 행렬은 비대칭이다. 검증 파일(J4_PAE_summary.csv)은 두 방향을 합쳐
    # 계산했으므로 여기서도 양쪽 블록을 이어붙여 쓴다.
    s1 = block(pae, SEG1, SEG1)
    s2 = block(pae, SEG2, SEG2)
    offA = block(pae, SEG1, SEG2).ravel()      # 행 42-57, 열 72-96
    offB = block(pae, SEG2, SEG1).ravel()      # 행 72-96, 열 42-57
    off = np.concatenate([offA, offB])

    def pair(i, j):
        """잔기 i, j 의 PAE (양방향 평균). 1-기반."""
        return 0.5 * (pae[i - 1, j - 1] + pae[j - 1, i - 1])

    print("\n[2/3] 검증 수치 (J4_PAE_summary.csv 와 대조)")
    print("      세그먼트 42-57 내부  mean %8.5f   (검증 3.0078125)" % s1.mean())
    print("      세그먼트 72-96 내부  mean %8.5f   (검증 3.6384)" % s2.mean())
    print("      (참고) 방향별 평균 %.4f / %.4f -> 합산 평균 사용"
          % (offA.mean(), offB.mean()))
    print("      두 세그먼트 사이     mean %8.5f   (검증 23.17875)" % off.mean())
    print("      두 세그먼트 사이   median %8.2f   (검증 24.0)" % np.median(off))
    print("      두 세그먼트 사이     max  %8.2f   (검증 29.0)" % off.max())
    print("      20 A 초과 비율       %8.3f   (검증 0.771)" % (off > 20).mean())
    print("      15 A 초과 비율       %8.3f   (검증 0.948)" % (off > 15).mean())
    print("      Phe45 - Arg72        %8.2f   (검증 21.0)" % pair(45, 72))
    print("      Phe45 - Asp69        %8.2f   (검증 18.5)" % pair(45, 69))
    print("      pLDDT 포켓 평균은 Table S18 참조")

    # ---- 그리기 ---------------------------------------------------------
    print("\n[3/3] 그림 생성")
    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    im = ax.imshow(pae, cmap="Greens_r", vmin=0, vmax=31.75,
                   extent=[0.5, n + 0.5, n + 0.5, 0.5], interpolation="nearest")

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Predicted aligned error (Å)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    def box(a, b, color, lw, ls="-"):
        ax.add_patch(Rectangle((b[0] - 0.5, a[0] - 0.5),
                               b[1] - b[0] + 1, a[1] - a[0] + 1,
                               fill=False, edgecolor=color, lw=lw, ls=ls))

    # 두 세그먼트의 대각 블록
    box(SEG1, SEG1, "#1f4e79", 1.4)
    box(SEG2, SEG2, "#1f4e79", 1.4)
    # 그 사이 비대각 블록 (양쪽 모두)
    box(SEG1, SEG2, "#c00000", 1.6)
    box(SEG2, SEG1, "#c00000", 1.6)
    # P0 영역
    box(P0, P0, "#ed7d31", 1.2, ls="--")

    # 포켓 잔기 위치 표시
    for r in POCKET:
        ax.plot(r, n + 2.5, marker="v", color="#ed7d31", ms=3.2,
                clip_on=False, zorder=5)

    ax.set_xlabel("Residue (scored)", fontsize=11.5, fontweight="bold")
    ax.set_ylabel("Residue (aligned)", fontsize=11.5, fontweight="bold")
    ax.set_xticks([1, 20, 40, 60, 80, 104])
    ax.set_yticks([1, 20, 40, 60, 80, 104])
    ax.tick_params(labelsize=8)

    ax.annotate("inter-segment\nmean %.1f \u00c5" % off.mean(),
                xy=((SEG2[0] + SEG2[1]) / 2, (SEG1[0] + SEG1[1]) / 2),
                xytext=(0.62, 0.30), textcoords="axes fraction",
                fontsize=8, color="#c00000", ha="left", va="top",
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec="#c00000", lw=0.6, alpha=0.9),
                arrowprops=dict(arrowstyle="->", color="#c00000", lw=0.8))

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = "Fig_S7_PAE_matrix.%s" % ext
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print("      저장: %s" % out)

    print("\n위 검증 수치가 기대값과 일치하면 그림을 SI Fig. S7 자리에 넣으세요.")


if __name__ == "__main__":
    main()
