#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
find_params.py — GSE198315 재분석 파이프라인 파라미터 찾기
===========================================================

Methods 2.2에 기술할 전처리/배치보정/클러스터링 설정을 로컬 스크립트에서
직접 뽑아낸다. 기억이나 추정에 의존하지 않기 위한 도구다.

실행
----
    cd "C:\\Users\\user\\Desktop\\WBP5 revision"
    python find_params.py

    # 다른 폴더까지 함께 뒤지려면 경로를 인자로 준다
    python find_params.py "C:\\Users\\user\\Desktop\\wbp5" "C:\\Users\\user\\Downloads"

하는 일
-------
.py / .ipynb / .txt / .md 파일을 재귀적으로 훑어, scanpy·harmony·leiden
관련 호출과 수치 설정이 들어간 줄을 파일별로 모아 출력하고
params_found_<날짜>.txt 로 저장한다.

출력 파일을 그대로 올려주면 Methods 문단으로 옮길 수 있다.
"""

import io
import os
import re
import sys
import json
import datetime

STAMP = datetime.date.today().strftime("%Y%m%d")
OUT = "params_found_%s.txt" % STAMP

EXTS = (".py", ".ipynb", ".txt", ".md", ".log")

# 찾을 패턴: (라벨, 정규식)
PATTERNS = [
    ("품질필터",   r"min_genes|min_cells|min_counts|pct_counts_mt|percent\.mt|n_genes_by_counts|filter_cells|filter_genes"),
    ("정규화",     r"normalize_total|target_sum|log1p|CPM|1e4|10000\b"),
    ("HVG",        r"highly_variable|n_top_genes|flavor\s*=|seurat_v3"),
    ("스케일/PCA", r"scale\(|max_value|pca\(|n_comps|n_pcs|svd_solver"),
    ("Harmony",    r"harmony|batch_key|max_iter|theta\s*=|lamb\s*="),
    ("이웃/클러스터", r"neighbors\(|n_neighbors|use_rep|leiden|louvain|resolution\s*=|n_clusters"),
    ("시드",       r"random_state|random\.seed|np\.random|seed\s*=|RANDOM_SEED"),
    ("서브샘플",   r"subsample|sample\(|n_obs|80000|80_000|5359|5,359"),
    ("마커/주석",  r"EPCAM|KRT14|KRT17|LAMC2|PTPRC|COL1A1|DCN|PECAM1|CD3D|LYZ|marker|annot"),
    ("구획명",     r"Epithelial|Fibroblast|Endothelial|Immune|Stromal|tumour-region|tumor-region"),
]
COMPILED = [(lab, re.compile(rx, re.I)) for lab, rx in PATTERNS]

SKIP_DIRS = {".git", "__pycache__", ".ipynb_checkpoints", "node_modules",
             "deprecated", ".venv", "venv"}


def iter_lines(path):
    """일반 텍스트는 그대로, .ipynb 는 코드 셀만 펼쳐서 (줄번호, 내용) 생성."""
    if path.lower().endswith(".ipynb"):
        try:
            nb = json.load(io.open(path, encoding="utf-8", errors="replace"))
        except Exception:
            return
        n = 0
        for cell in nb.get("cells", []):
            if cell.get("cell_type") != "code":
                continue
            for src in cell.get("source", []):
                n += 1
                yield n, src.rstrip("\n")
    else:
        try:
            with io.open(path, encoding="utf-8", errors="replace") as fh:
                for n, line in enumerate(fh, 1):
                    yield n, line.rstrip("\n")
        except Exception:
            return


def scan(root):
    hits = {}
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in filenames:
            if not fn.lower().endswith(EXTS):
                continue
            full = os.path.join(dirpath, fn)
            if os.path.getsize(full) > 20 * 1024 * 1024:   # 20 MB 초과는 건너뜀
                continue
            found = []
            for n, line in iter_lines(full):
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                labs = [lab for lab, rx in COMPILED if rx.search(s)]
                if labs:
                    found.append((n, labs, s[:200]))
            if found:
                hits[full] = found
    return hits


def main():
    roots = sys.argv[1:] or ["."]
    buf = []

    def emit(s=""):
        print(s)
        buf.append(s)

    emit("GSE198315 파이프라인 파라미터 탐색  (%s)" % STAMP)
    emit("=" * 72)

    total = 0
    for root in roots:
        if not os.path.isdir(root):
            emit("\n[건너뜀] 폴더가 없습니다: %s" % root)
            continue
        emit("\n### 탐색 경로: %s" % os.path.abspath(root))
        hits = scan(root)
        if not hits:
            emit("    (해당 파일 없음)")
            continue
        # 관련 줄이 많은 파일부터
        for full, found in sorted(hits.items(), key=lambda kv: -len(kv[1])):
            emit("\n--- %s   (%d줄)" % (os.path.relpath(full, root), len(found)))
            for n, labs, s in found:
                emit("  %5d  [%s]  %s" % (n, "/".join(labs), s))
                total += 1

    emit("\n" + "=" * 72)
    emit("총 %d줄" % total)
    emit("""
다음으로 할 일
  1) 위 출력에서 GSE198315 재분석 스크립트가 어느 파일인지 확인
  2) %s 파일을 그대로 올리기
  3) 해당 스크립트 파일 자체도 함께 올리면 가장 정확함

아무것도 안 나왔다면
  - 다른 폴더를 인자로 지정:  python find_params.py "C:\\Users\\user\\Desktop\\wbp5"
  - 또는 GitHub 저장소 phase1_singlecell/ 폴더를 확인
""" % OUT)

    with io.open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(buf))
    print("저장: %s" % OUT)


if __name__ == "__main__":
    main()
