#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE198315 재주석 파이프라인 — Reviewer #1 point 3 / Reviewer #5 Major comment 3 대응
====================================================================================
proxy 마커 기반 분류를 Leiden 클러스터링 + 마커 주석 + inferCNV 악성 판정으로 교체합니다.

생성: 2026-10-15
입력: GSE198315 10x triplet (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
출력: WBP5/REVISION/results/ 에 업로드할 CSV 일체

[왜 CellTypist/Azimuth를 쓰지 않는가]
리뷰어는 "Leiden clustering with marker-based annotation OR reference-based label
transfer (CellTypist, Azimuth)"라고 OR로 요구했습니다. CellTypist의 공개 모델은
면역세포 중심이고 HNSCC 상피/기질에 대해 검증된 레퍼런스가 없습니다. 부적합한
레퍼런스로 label transfer를 하면 마커 기반 주석보다 부정확해지므로 앞쪽을
선택합니다. 이 판단 근거는 답변서에 명시합니다.

[실행 순서]
  # 0) 환경
  pip install scanpy infercnvpy leidenalg igraph

  # 1) GENCODE GTF 내려받기 (inferCNV의 유전자 좌표용, 약 50 MB)
  #    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/
  #        gencode.v44.basic.annotation.gtf.gz
  #    ※ GSE198315이 정렬된 참조 유전체 버전과 맞추십시오 (보통 GRCh38)

  # 2) 1단계 — 적재 + 서브샘플 + 체크포인트 저장 (메모리 가장 많이 씀)
  python GSE198315_reannotation_20261015.py step1 \
      --mtx GSE198315_matrix.mtx.gz \
      --barcodes GSE198315_barcodes.tsv.gz \
      --features GSE198315_features.tsv.gz \
      --n-cells 80000 --seed 42 \
      --out gse198315_checkpoint.h5ad

  # 3) 2단계 — 클러스터링 + 마커 주석
  python GSE198315_reannotation_20261015.py step2 \
      --h5ad gse198315_checkpoint.h5ad --outdir ./results_gse198315

  # 4) 3단계 — inferCNV 악성 판정 + 최종 대조표
  python GSE198315_reannotation_20261015.py step3 \
      --h5ad gse198315_checkpoint.h5ad --gtf gencode.v44.basic.annotation.gtf.gz \
      --outdir ./results_gse198315

단계를 나눈 이유: step1이 메모리 피크이고, 여기서 끊기면 처음부터 다시 해야 하므로
h5ad로 저장해 step2/step3을 몇 번이든 다시 돌릴 수 있게 했습니다.

[메모리 가이드]
  268,131 세포 전체 → 32 GB 이상 권장
  80,000 세포 서브샘플 → 16 GB 정도면 충분
  --n-cells 0 을 주면 전체를 사용합니다.
"""

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

RANDOM_SEED = 42

# ---------------------------------------------------------------------------
# 마커 패널 — 클러스터 주석용 (표준 canonical marker)
# 원고의 106-gene proxy 패널과 별개입니다. 이 패널은 "클러스터가 무엇인가"를
# 판정하는 데만 쓰고, 세포 단위 점수 매기기에는 쓰지 않습니다.
# ---------------------------------------------------------------------------
MARKERS = {
    "Epithelial":  ["EPCAM", "KRT5", "KRT14", "KRT17", "KRT19", "SFN", "CDH1"],
    "Fibroblast":  ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "THY1", "FAP"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "RAMP2"],
    "T/NK":        ["CD3D", "CD3E", "CD2", "IL7R", "NKG7", "GNLY"],
    "B/Plasma":    ["CD79A", "MS4A1", "IGHG1", "MZB1", "JCHAIN"],
    "Myeloid":     ["LYZ", "CD68", "CD14", "AIF1", "ITGAX", "C1QA"],
    "Mast":        ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Muscle":      ["ACTA2", "MYH11", "DES", "TAGLN"],
}
IMMUNE_LABELS = ["T/NK", "B/Plasma", "Myeloid", "Mast"]   # inferCNV reference

# 원고에서 쓰는 p-EMT / 침습 마커 (대조 분석용)
PEMT_MARKERS = ["LAMC2", "ITGA5", "KRT14", "KRT17", "VIM"]
GENE = "WBP5"
GENE_ALIASES = ["WBP5", "TCEAL9"]


def log(msg):
    sys.stderr.write(f"[{pd.Timestamp.now():%H:%M:%S}] {msg}\n")
    sys.stderr.flush()


def resolve_gene(adata, aliases):
    for g in aliases:
        if g in adata.var_names:
            return g
    return None


# ===========================================================================
# STEP 1 — 적재 / 서브샘플 / 체크포인트
# ===========================================================================

def read_features(path):
    """CellRanger features.tsv(.gz) 를 읽어 DataFrame 반환. 2열/3열 모두 처리."""
    op = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            rows.append(parts)
    ncol = max(len(r) for r in rows)
    if ncol >= 3:
        df = pd.DataFrame(rows, columns=["gene_id", "gene_symbol", "feature_type"][:ncol])
    elif ncol == 2:
        df = pd.DataFrame(rows, columns=["gene_id", "gene_symbol"])
    else:
        df = pd.DataFrame(rows, columns=["gene_symbol"])
        df["gene_id"] = df["gene_symbol"]
    return df


def _mtx_header(path):
    """MatrixMarket 헤더를 읽어 (n_rows, n_cols, nnz, 건너뛸 줄 수) 반환."""
    op = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with op(path, "rt") as fh:
        for ln in fh:
            n += 1
            if ln.startswith("%"):
                continue
            p = ln.split()
            return int(p[0]), int(p[1]), int(p[2]), n
    raise SystemExit("오류: MTX 헤더를 읽지 못했습니다.")


def _mtx_chunks(path, skiprows, chunksize):
    """MTX 본문을 청크 단위로 흘려보냅니다. (row0, col0, val) 0-기반."""
    for ch in pd.read_csv(path, sep=r"\s+", header=None, names=["r", "c", "v"],
                          skiprows=skiprows, chunksize=chunksize,
                          dtype={"r": np.int32, "c": np.int32, "v": np.float32},
                          comment="%"):
        yield (ch["r"].to_numpy() - 1, ch["c"].to_numpy() - 1, ch["v"].to_numpy())


def step1(args):
    """
    적재 + QC + 서브샘플 + 체크포인트.

    [메모리 설계 — 중요]
    scipy.io.mmread 로 전체 행렬을 한 번에 읽으면 COO(약 10 GB) + CSR 변환(+5 GB)
    + 전치(+5 GB) 로 피크가 20 GB에 달합니다. GSE198315은 nnz 가 6.2억이라
    16 GB 장비에서는 실패합니다.

    그래서 파일을 두 번 흘려 읽습니다.
      1차: 세포별 통계(엔트리 수, 총 count, MT count)만 누적 → QC 통과 세포 선정
           → 서브샘플 → 남길 엔트리 수를 정확히 계산
      2차: 선정된 세포의 엔트리만 미리 할당한 배열에 채워 넣음
    피크 메모리는 최종 CSR 크기(80,000 세포 기준 약 2 GB) 수준으로 떨어집니다.
    대신 파일을 두 번 읽으므로 시간이 조금 더 걸립니다.
    """
    import anndata as ad
    from scipy.sparse import coo_matrix

    log(f"features 읽는 중: {args.features}")
    feat = read_features(args.features)
    log(f"  features {len(feat):,}행 / 열 {list(feat.columns)}")
    log(f"  예시: {feat.head(2).to_dict('records')}")

    log(f"barcodes 읽는 중: {args.barcodes}")
    bc = pd.read_csv(args.barcodes, sep="\t", header=None)
    barcodes = bc.iloc[:, 0].astype(str).tolist()
    log(f"  barcodes {len(barcodes):,}")

    nr, nc, nnz, hdr = _mtx_header(args.mtx)
    log(f"  mtx {nr} x {nc}, nnz {nnz:,} (헤더 {hdr}줄)")
    if nr == len(feat) and nc == len(barcodes):
        genes_are_rows = True
    elif nr == len(barcodes) and nc == len(feat):
        genes_are_rows = False
    else:
        raise SystemExit(
            f"오류: mtx ({nr},{nc}) 가 features({len(feat)})/barcodes({len(barcodes)})와 맞지 않습니다.")
    n_genes_tot, n_cells_tot = (nr, nc) if genes_are_rows else (nc, nr)

    # pandas 문자열 메서드로 처리합니다. np.char 는 object dtype 배열에서
    # 'string operation on non-string array' 로 실패합니다.
    mt_gene = (feat["gene_symbol"].astype(str).str.upper()
               .str.startswith("MT-").to_numpy())

    # ---- 1차 통과: 세포별 통계만 누적 -------------------------------------
    log("1차 통과 — 세포별 QC 통계 누적 중 (파일 전체를 한 번 읽습니다)")
    entries = np.zeros(n_cells_tot, dtype=np.int64)
    tot_counts = np.zeros(n_cells_tot, dtype=np.float64)
    mt_counts = np.zeros(n_cells_tot, dtype=np.float64)
    seen = 0
    for r, c, v in _mtx_chunks(args.mtx, hdr, args.chunksize):
        gi, ci = (r, c) if genes_are_rows else (c, r)
        entries += np.bincount(ci, minlength=n_cells_tot)
        tot_counts += np.bincount(ci, weights=v, minlength=n_cells_tot)
        m = mt_gene[gi]
        if m.any():
            mt_counts += np.bincount(ci[m], weights=v[m], minlength=n_cells_tot)
        seen += len(v)
        if seen % (args.chunksize * 5) < args.chunksize:
            log(f"  {seen:,} / {nnz:,} ({seen/nnz:.0%})")
    log(f"  1차 완료: {seen:,} 엔트리")

    with np.errstate(invalid="ignore", divide="ignore"):
        pct_mt = np.where(tot_counts > 0, mt_counts / tot_counts * 100.0, 0.0)

    # ---- 세포 선정: QC -> (바코드 목록 | 서브샘플) -------------------------
    keep = (entries >= args.min_genes) & (pct_mt <= args.max_mt)
    log(f"  QC 통과 {int(keep.sum()):,} / {n_cells_tot:,} "
        f"(min_genes={args.min_genes}, max_mt={args.max_mt}%)")

    bc_arr = np.asarray(barcodes, dtype=object)
    cand = np.flatnonzero(keep)
    if args.barcode_list and os.path.exists(args.barcode_list):
        want = {l.strip() for l in open(args.barcode_list) if l.strip()}
        cand = np.array([i for i in cand if bc_arr[i] in want], dtype=np.int64)
        log(f"  바코드 목록 적용: {len(cand):,} 세포")
    elif args.n_cells and args.n_cells > 0 and len(cand) > args.n_cells:
        rng = np.random.default_rng(args.seed)
        cand = np.sort(rng.choice(cand, size=args.n_cells, replace=False))
        log(f"  서브샘플 seed={args.seed}: {len(cand):,} 세포")
    if len(cand) == 0:
        raise SystemExit("오류: 남은 세포가 없습니다. QC 기준을 확인하십시오.")

    new_index = np.full(n_cells_tot, -1, dtype=np.int64)
    new_index[cand] = np.arange(len(cand))
    keep_nnz = int(entries[cand].sum())
    log(f"  남길 엔트리 {keep_nnz:,} (전체의 {keep_nnz/nnz:.1%}) "
        f"→ 예상 메모리 약 {keep_nnz*12/1e9:.1f} GB")

    # ---- 2차 통과: 선정 세포의 엔트리만 미리 할당한 배열에 채움 ------------
    log("2차 통과 — 선정된 세포의 엔트리만 적재 중")
    out_r = np.empty(keep_nnz, dtype=np.int32)   # 새 세포 인덱스
    out_c = np.empty(keep_nnz, dtype=np.int32)   # 유전자 인덱스
    out_v = np.empty(keep_nnz, dtype=np.float32)
    pos, seen = 0, 0
    for r, c, v in _mtx_chunks(args.mtx, hdr, args.chunksize):
        gi, ci = (r, c) if genes_are_rows else (c, r)
        ni = new_index[ci]
        m = ni >= 0
        k = int(m.sum())
        if k:
            out_r[pos:pos+k] = ni[m]
            out_c[pos:pos+k] = gi[m]
            out_v[pos:pos+k] = v[m]
            pos += k
        seen += len(v)
        if seen % (args.chunksize * 5) < args.chunksize:
            log(f"  {seen:,} / {nnz:,} ({seen/nnz:.0%})")
    if pos != keep_nnz:
        log(f"  경고: 채워진 엔트리 {pos:,} != 예상 {keep_nnz:,} — 배열을 잘라냅니다")
        out_r, out_c, out_v = out_r[:pos], out_c[:pos], out_v[:pos]
    log(f"  2차 완료: {pos:,} 엔트리")

    X = coo_matrix((out_v, (out_r, out_c)),
                   shape=(len(cand), n_genes_tot), dtype=np.float32).tocsr()
    del out_r, out_c, out_v

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=pd.Index([str(bc_arr[i]) for i in cand], name="barcode")),
        var=feat.set_index("gene_symbol"))
    adata.var_names_make_unique()
    adata.obs["n_genes_by_counts"] = entries[cand]
    adata.obs["total_counts"] = tot_counts[cand]
    adata.obs["pct_counts_mt"] = pct_mt[cand]
    log(f"  AnnData {adata.shape} (cells x genes)")

    # 샘플 토큰 추출: '...-1' 접미사 또는 '샘플_바코드' 접두사
    def sample_of(b):
        """
        바코드에서 샘플 라벨 추출. 두 관례를 모두 처리합니다.
          (a) 접두사형: 'P04-mLN_AAACCTGAGA...'  -> 'P04-mLN'   (GSE198315)
          (b) 접미사형: 'AAACCTGAGA...-1'        -> 'S1'        (CellRanger aggr)
        접두사를 먼저 봅니다. 'P01-TC_AAAC...-1' 처럼 둘이 섞여 있을 때
        접미사를 먼저 보면 환자/부위 정보를 잃기 때문입니다.
        """
        b = str(b).strip()
        if "_" in b:
            pre = b.split("_")[0]
            # 접두사가 순수 염기서열이면 샘플 라벨이 아니라 바코드입니다
            if pre and set(pre.upper()) - set("ACGTN"):
                return pre
        if "-" in b:
            tail = b.rsplit("-", 1)[1]
            if tail.isdigit():
                return f"S{tail}"
        return "ALL"
    adata.obs["sample"] = [sample_of(b) for b in adata.obs_names]
    vc = adata.obs["sample"].value_counts()
    log(f"  감지된 샘플 토큰 {len(vc)}개: {dict(list(vc.items())[:12])}")
    if len(vc) == 1:
        log("  경고: 바코드에서 샘플을 구분할 수 없습니다. --sample-meta 로 매핑을 주십시오.")

    if args.sample_meta and os.path.exists(args.sample_meta):
        meta = pd.read_csv(args.sample_meta)
        log(f"  샘플 메타 적용: {list(meta.columns)}")
        meta = meta.set_index(meta.columns[0])
        for col in meta.columns:
            adata.obs[col] = adata.obs["sample"].map(meta[col]).astype(str)

    adata.var["mt"] = (pd.Series(adata.var_names.astype(str))
                       .str.upper().str.startswith("MT-").to_numpy())

    # 사용한 바코드를 기록 (재현성 — Reviewer #5 Major comment 10)
    pd.Series(adata.obs_names).to_csv(
        os.path.splitext(args.out)[0] + "_barcodes_used.csv", index=False, header=["barcode"])

    adata.layers["counts"] = adata.X.copy()      # inferCNV/재현용 원본 count 보존
    adata.write_h5ad(args.out, compression="gzip")
    log(f"저장 완료: {args.out}  {adata.shape}")


# ===========================================================================
# STEP 2 — 정규화 / Leiden / 마커 주석
# ===========================================================================

def step2(args):
    import scanpy as sc
    os.makedirs(args.outdir, exist_ok=True)
    sc.settings.n_jobs = args.threads

    log(f"적재: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    log(f"  {adata.shape}")

    # --- 샘플 토큰을 환자/부위로 분해
    # GSE198315의 바코드는 'P04-mLN_...' 형태로 환자와 조직 부위를 함께 담고 있습니다.
    #   P01..P09+ = 환자, TC/TP/NT/mLN = 조직 부위
    # 환자 단위 분석(Reviewer #5 Major 1)과 전이 림프절 재해석(Major 4)에 필요합니다.
    s = adata.obs["sample"].astype(str)
    if s.str.contains("-").any():
        adata.obs["patient"] = s.str.split("-").str[0]
        adata.obs["region"] = s.str.split("-").str[1:].str.join("-")
        log(f"  환자 {adata.obs['patient'].nunique()}명 / "
            f"부위 {sorted(adata.obs['region'].unique())}")
        log(f"  환자별 세포수: {adata.obs['patient'].value_counts().to_dict()}")
    else:
        adata.obs["patient"] = s
        adata.obs["region"] = "NA"
        log("  경고: 샘플 토큰에서 환자/부위를 분리하지 못했습니다")

    # Run B와 동일한 정규화 (CPM 1e4 + log1p) — 원고 §2.3과 일치시킵니다
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, batch_key=
                                "sample" if adata.obs["sample"].nunique() > 1 else None)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=args.n_pcs, svd_solver="arpack", random_state=RANDOM_SEED)

    # 배치(샘플) 보정 — harmony가 있으면 사용, 없으면 생략하고 기록
    rep = "X_pca"
    if adata.obs["sample"].nunique() > 1:
        try:
            import scanpy.external as sce
            sce.pp.harmony_integrate(adata, "sample", random_state=RANDOM_SEED)
            rep = "X_pca_harmony"
            log("  harmony 배치보정 적용")
        except Exception as e:
            log(f"  harmony 미적용 ({type(e).__name__}) — X_pca 사용. 답변서에 명시할 것")

    sc.pp.neighbors(adata, n_neighbors=15, use_rep=rep, random_state=RANDOM_SEED)
    sc.tl.leiden(adata, resolution=args.resolution, key_added="leiden",
                 random_state=RANDOM_SEED, flavor="igraph", n_iterations=2, directed=False)
    log(f"  Leiden 클러스터 {adata.obs['leiden'].nunique()}개 (resolution={args.resolution})")
    sc.tl.umap(adata, random_state=RANDOM_SEED)

    # --- 클러스터별 마커 점수 → 최고 점수 세포형으로 주석
    full = adata.raw.to_adata()
    full.obs["leiden"] = adata.obs["leiden"].values
    score_tbl = {}
    for ct, genes in MARKERS.items():
        present = [g for g in genes if g in full.var_names]
        if not present:
            log(f"  경고: {ct} 마커가 매트릭스에 없음 — 건너뜀")
            continue
        # h5ad(HDF5) 키에는 '/' 를 쓸 수 없으므로 열 이름만 치환합니다.
        # 표시용 라벨("T/NK")은 그대로 두고 저장 키만 "score_T_NK" 가 됩니다.
        key = "score_" + ct.replace("/", "_")
        import scanpy as _sc
        _sc.tl.score_genes(full, present, score_name=key, random_state=RANDOM_SEED)
        score_tbl[ct] = full.obs.groupby("leiden", observed=True)[key].mean()
    S = pd.DataFrame(score_tbl)
    S.to_csv(f"{args.outdir}/G1_cluster_marker_scores.csv")

    ann = S.idxmax(axis=1)
    margin = S.max(axis=1) - S.apply(lambda r: r.nlargest(2).iloc[-1], axis=1)
    # 1·2위 점수 차가 작은 클러스터는 Ambiguous 로 두고 수동 확인 대상으로 표시
    ann[margin < args.min_margin] = "Ambiguous"
    adata.obs["cell_type"] = adata.obs["leiden"].map(ann).astype(str)
    full.obs["cell_type"] = adata.obs["cell_type"].values

    # 클러스터가 세포 상태가 아니라 환자를 반영하는지 점검합니다.
    # 배치 보정 여부와 무관하게, 각 클러스터에 몇 명의 환자가 기여했는지와
    # 최다 기여 환자의 점유율을 보고하면 "환자 교란" 지적에 직접 답할 수 있습니다.
    pat_by_clu = pd.crosstab(adata.obs["leiden"], adata.obs["patient"])
    pat_by_clu.to_csv(f"{args.outdir}/G2b_cluster_by_patient.csv")
    n_pat = pat_by_clu.gt(0).sum(axis=1)
    top_frac = pat_by_clu.max(axis=1) / pat_by_clu.sum(axis=1)

    summary = pd.DataFrame({
        "leiden": S.index, "assigned": ann.values, "score_margin": margin.values,
        "n_cells": adata.obs["leiden"].value_counts().reindex(S.index).values,
        "n_patients": n_pat.reindex(S.index).values,
        "max_single_patient_frac": top_frac.reindex(S.index).values})
    summary.to_csv(f"{args.outdir}/G2_cluster_annotation.csv", index=False)
    dom = summary[summary.max_single_patient_frac > 0.5]
    if len(dom):
        log(f"  ※ 단일 환자가 50% 넘게 차지한 클러스터 {len(dom)}개: "
            f"{dom.leiden.tolist()} — 답변서에서 다룰 것")
    else:
        log("  모든 클러스터가 여러 환자에 걸쳐 있습니다 (단일 환자 점유율 ≤50%)")
    log(f"  주석 결과:\n{adata.obs['cell_type'].value_counts().to_string()}")
    if (ann == "Ambiguous").any():
        log(f"  ※ Ambiguous 클러스터 {(ann=='Ambiguous').sum()}개 — G1 점수표를 직접 확인하십시오")

    # --- proxy 분류와의 대조 (있는 경우)
    for col in ["proxy_class", "proxy_classification"]:
        if col in adata.obs.columns:
            ct = pd.crosstab(adata.obs[col], adata.obs["cell_type"])
            ct.to_csv(f"{args.outdir}/G3_proxy_vs_new_annotation.csv")
            log(f"  proxy 대조표 저장")
            break

    # --- WBP5 분포를 새 주석 기준으로 재집계
    g = resolve_gene(full, GENE_ALIASES)
    if g:
        v = pd.Series(np.asarray(full[:, g].X.todense()).ravel(), index=full.obs_names)
        thr = v[v > 0].quantile(0.9)          # 원고 정의: 비영 세포 중 상위 10%
        full.obs["WBP5_high"] = (v > thr).values
        tab = pd.DataFrame({
            "n_cells": full.obs.groupby("cell_type", observed=True).size(),
            "mean_WBP5": v.groupby(full.obs["cell_type"]).mean(),
            "pct_positive": v.groupby(full.obs["cell_type"]).apply(lambda s: (s > 0).mean()*100),
            "n_WBP5_high": full.obs.groupby("cell_type", observed=True)["WBP5_high"].sum(),
        })
        tab["pct_of_all_WBP5_high"] = tab.n_WBP5_high / tab.n_WBP5_high.sum() * 100
        # 각 세포형에 몇 명의 환자가 기여했는지 — 세포 수만 보고하지 말라는 요구(#5 Major 1)
        tab["n_patients"] = full.obs.groupby("cell_type", observed=True)["patient"].nunique()
        tab.to_csv(f"{args.outdir}/G4_WBP5_by_new_celltype.csv")
        log(f"  WBP5 재집계 저장 (threshold={thr:.4f})")

        # --- 세포형 × 조직부위 교차 (전이 림프절 재해석, #5 Major 4)
        g5 = pd.DataFrame({
            "n_cells": full.obs.groupby(["cell_type", "region"], observed=True).size(),
            "n_WBP5_high": full.obs.groupby(["cell_type", "region"], observed=True)["WBP5_high"].sum(),
            "n_patients": full.obs.groupby(["cell_type", "region"], observed=True)["patient"].nunique(),
            "mean_WBP5": v.groupby([full.obs["cell_type"], full.obs["region"]]).mean(),
        })
        g5["pct_WBP5_high"] = g5.n_WBP5_high / g5.n_cells * 100
        g5.to_csv(f"{args.outdir}/G5_WBP5_by_celltype_region.csv")

        # --- 환자별 WBP5-high 비율 (환자 단위 추론의 입력)
        g6 = pd.DataFrame({
            "n_cells": full.obs.groupby(["patient", "region"], observed=True).size(),
            "n_WBP5_high": full.obs.groupby(["patient", "region"], observed=True)["WBP5_high"].sum(),
            "mean_WBP5": v.groupby([full.obs["patient"], full.obs["region"]]).mean(),
        })
        g6["pct_WBP5_high"] = g6.n_WBP5_high / g6.n_cells * 100
        g6.to_csv(f"{args.outdir}/G6_WBP5_by_patient_region.csv")
        log("  세포형×부위, 환자×부위 교차표 저장")
    else:
        log("  경고: WBP5/TCEAL9 를 찾지 못했습니다")

    adata.write_h5ad(args.h5ad.replace(".h5ad", "_annotated.h5ad"), compression="gzip")
    full.write_h5ad(args.h5ad.replace(".h5ad", "_full_annotated.h5ad"), compression="gzip")
    log("STEP 2 완료")


# ===========================================================================
# STEP 3 — inferCNV 악성 판정
# ===========================================================================

def parse_gtf_gene_positions(gtf_path):
    """
    GTF에서 gene 단위 좌표를 읽어 DataFrame(index=gene_name,
    columns=[chromosome, start, end]) 반환.

    외부 패키지를 쓰지 않습니다. 표준 GTF 포맷만 가정합니다:
      seqname  source  feature  start  end  score  strand  frame  attributes
    feature == 'gene' 인 행에서 attributes 의 gene_name 을 키로 씁니다.
    같은 심볼이 여러 번 나오면 가장 긴 구간을 채택합니다(패러로그/PAR 영역 대비).
    """
    op = gzip.open if str(gtf_path).endswith(".gz") else open
    recs = {}
    with op(gtf_path, "rt") as fh:
        for ln in fh:
            if not ln or ln[0] == "#":
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            attrs = p[8]
            i = attrs.find('gene_name "')
            if i < 0:
                continue
            j = attrs.find('"', i + 11)
            name = attrs[i + 11:j]
            try:
                s, e = int(p[3]), int(p[4])
            except ValueError:
                continue
            prev = recs.get(name)
            if prev is None or (e - s) > (prev[2] - prev[1]):
                recs[name] = (p[0], s, e)
    if not recs:
        raise SystemExit(
            "오류: GTF에서 gene_name 속성을 가진 'gene' 행을 찾지 못했습니다. "
            "GENCODE/Ensembl 표준 GTF인지 확인하십시오.")
    return pd.DataFrame.from_dict(
        recs, orient="index", columns=["chromosome", "start", "end"])

def step3(args):
    import scanpy as sc
    import infercnvpy as cnv
    os.makedirs(args.outdir, exist_ok=True)

    path = args.h5ad.replace(".h5ad", "_full_annotated.h5ad")
    if not os.path.exists(path):
        raise SystemExit(f"오류: {path} 가 없습니다. step2 를 먼저 실행하십시오.")
    log(f"적재: {path}")
    adata = sc.read_h5ad(path)

    # --- dtype 정규화 (필수)
    # anndata 0.13+ 는 문자열을 arrow 기반 dtype('str')으로 읽습니다.
    # infercnvpy 내부가 var 인덱스에 .take() 를 수행하는데, arrow 배열에서
    # 'ArrowInvalid: only handle 1-dimensional arrays' 로 죽습니다.
    # object dtype 으로 되돌리면 해결됩니다. 버전 올라가면 불필요해질 수 있지만
    # 두어도 해롭지 않습니다.
    adata.var.index = pd.Index(np.asarray(adata.var_names, dtype=object), dtype=object)
    adata.obs.index = pd.Index(np.asarray(adata.obs_names, dtype=object), dtype=object)
    for c in adata.var.columns:
        if str(adata.var[c].dtype) in ("str", "string", "string[pyarrow]"):
            adata.var[c] = adata.var[c].astype(object)
    for c in adata.obs.columns:
        if str(adata.obs[c].dtype) in ("str", "string", "string[pyarrow]"):
            adata.obs[c] = adata.obs[c].astype(object)

    # 유전자 좌표 부여 (inferCNV 필수)
    #
    # infercnvpy.io.genomic_position_from_gtf 는 쓰지 않습니다.
    # 해당 함수는 polars 기반 gtfparse 를 기대하는데, pandas 를 반환하는 버전이
    # 설치되면 'DataFrame has no attribute to_pandas' 로 죽습니다(버전 충돌).
    # 아래 파서는 표준 GTF 9번째 열만 읽으므로 의존성이 없고 버전에 영향받지 않습니다.
    log(f"GTF에서 유전자 좌표 읽는 중: {args.gtf}")
    pos = parse_gtf_gene_positions(args.gtf)
    log(f"  GTF에서 유전자 {len(pos):,}개 좌표 확보")
    adata.var["chromosome"] = adata.var_names.map(pos["chromosome"]).values
    adata.var["start"] = pd.to_numeric(adata.var_names.map(pos["start"]), errors="coerce").values
    adata.var["end"] = pd.to_numeric(adata.var_names.map(pos["end"]), errors="coerce").values
    ok = adata.var["chromosome"].notna()
    frac = float(ok.sum()) / max(adata.n_vars, 1)
    log(f"  좌표 부여된 유전자 {int(ok.sum()):,} / {adata.n_vars:,} ({frac:.1%})")
    # 실제 데이터에서 매핑률이 낮으면 GTF 어셈블리/심볼 표기가 어긋난 것입니다.
    # 비율과 절대수를 함께 보되, 소규모 테스트를 막지 않도록 --min-mapped 로 조정 가능.
    if frac < 0.30 or ok.sum() < args.min_mapped:
        raise SystemExit(
            f"오류: 좌표 매핑이 부족합니다 ({ok.sum():,}개, {frac:.1%}). "
            "GTF 어셈블리 버전과 유전자 심볼 표기를 확인하거나 --min-mapped 로 기준을 조정하십시오.")

    ref = [c for c in IMMUNE_LABELS if c in set(adata.obs["cell_type"])]
    if not ref:
        raise SystemExit("오류: reference 로 쓸 면역 클러스터가 없습니다.")
    log(f"  reference 구획: {ref}")

    cnv.tl.infercnv(adata, reference_key="cell_type", reference_cat=ref,
                    window_size=args.window)
    cnv.tl.pca(adata); cnv.pp.neighbors(adata); cnv.tl.leiden(adata)
    cnv.tl.cnv_score(adata)

    # 악성 판정: 면역 reference 의 CNV score 분포 대비 임계
    s = adata.obs["cnv_score"]
    ref_mask = adata.obs["cell_type"].isin(ref)
    thr = s[ref_mask].mean() + args.cnv_sd * s[ref_mask].std()
    adata.obs["cnv_malignant"] = (s > thr)
    log(f"  CNV score 임계 {thr:.4f} (면역 평균 + {args.cnv_sd}SD)")
    log(f"  악성 판정 {int(adata.obs['cnv_malignant'].sum()):,} / {adata.n_obs:,}")

    pd.crosstab(adata.obs["cell_type"], adata.obs["cnv_malignant"]).to_csv(
        f"{args.outdir}/H1_celltype_vs_inferCNV.csv")
    adata.obs[["cell_type", "cnv_score", "cnv_malignant", "sample"]].to_csv(
        f"{args.outdir}/H2_percell_cnv.csv")

    # --- 최종 대조: WBP5-high 세포의 정체
    g = resolve_gene(adata, GENE_ALIASES)
    if g:
        v = pd.Series(np.asarray(adata[:, g].X.todense()).ravel(), index=adata.obs_names)
        thr_w = v[v > 0].quantile(0.9)
        high = v > thr_w
        grp = adata.obs["cell_type"].astype(str) + np.where(
            adata.obs["cnv_malignant"], " (CNV+)", " (CNV-)")
        out = pd.DataFrame({
            "n_cells": grp.groupby(grp).size(),
            "n_WBP5_high": pd.Series(high.values, index=grp.index).groupby(grp).sum(),
            "mean_WBP5": v.groupby(grp).mean(),
            "pct_positive": v.groupby(grp).apply(lambda s: (s > 0).mean()*100),
        })
        out["pct_of_all_WBP5_high"] = out.n_WBP5_high / out.n_WBP5_high.sum() * 100
        out.to_csv(f"{args.outdir}/H3_WBP5_high_identity.csv")
        log(f"  WBP5-high 정체 분해 저장")

        # --- 구획별 Cliff's delta (GSE103322 분석과 동일 설계)
        from scipy import stats as st
        mal = v[adata.obs["cnv_malignant"].values]
        rows = []
        for ct in sorted(set(adata.obs["cell_type"])):
            if ct in ref:
                pass
            other = v[(adata.obs["cell_type"] == ct).values & (~adata.obs["cnv_malignant"].values)]
            if len(other) < 20 or len(mal) < 20:
                continue
            a, b = mal.values, other.values
            n1, n2 = len(a), len(b)
            r = st.rankdata(np.concatenate([a, b]))
            u1 = r[:n1].sum() - n1*(n1+1)/2
            d = (2*u1)/(n1*n2) - 1
            rows.append({"contrast": f"CNV+ malignant vs {ct} (CNV-)", "n_other": n2,
                         "mean_malignant": a.mean(), "mean_other": b.mean(),
                         "cliffs_delta": d,
                         "MWU_p": st.mannwhitneyu(a, b, alternative="two-sided").pvalue})
        pd.DataFrame(rows).to_csv(f"{args.outdir}/H4_compartment_contrasts.csv", index=False)
        log("  구획별 Cliff's delta 저장")

    adata.write_h5ad(path.replace(".h5ad", "_cnv.h5ad"), compression="gzip")

    with open(f"{args.outdir}/environment.txt", "w") as fh:
        import scipy, statsmodels
        fh.write(f"python {sys.version}\nnumpy {np.__version__}\npandas {pd.__version__}\n")
        fh.write(f"scanpy {sc.__version__}\ninfercnvpy {cnv.__version__}\n")
        fh.write(f"seed {RANDOM_SEED}\nresolution {args.resolution if hasattr(args,'resolution') else 'n/a'}\n")
        fh.write(f"gtf {args.gtf}\ncnv_window {args.window}\ncnv_threshold_sd {args.cnv_sd}\n")
        fh.write(f"reference_compartments {ref}\n")
    log("STEP 3 완료")


# ===========================================================================

def main():
    ap = argparse.ArgumentParser(description="GSE198315 재주석 (Reviewer #1 p3 / #5 Major 3)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("step1", help="적재 + QC + 서브샘플 + h5ad 저장")
    p1.add_argument("--mtx", required=True); p1.add_argument("--barcodes", required=True)
    p1.add_argument("--features", required=True); p1.add_argument("--out", default="gse198315_checkpoint.h5ad")
    p1.add_argument("--n-cells", type=int, default=80000, help="0 이면 전체 사용")
    p1.add_argument("--seed", type=int, default=RANDOM_SEED)
    p1.add_argument("--barcode-list", default=None, help="원고의 80k subset 바코드 목록(있으면 우선)")
    p1.add_argument("--sample-meta", default=None, help="sample,patient,region 매핑 CSV")
    p1.add_argument("--min-genes", type=int, default=200)
    p1.add_argument("--max-mt", type=float, default=20.0)
    p1.add_argument("--chunksize", type=int, default=20_000_000,
                    help="MTX를 한 번에 읽어들일 엔트리 수. 메모리가 빠듯하면 줄이십시오(예: 5000000)")
    p1.set_defaults(func=step1)

    p2 = sub.add_parser("step2", help="정규화 + Leiden + 마커 주석")
    p2.add_argument("--h5ad", required=True); p2.add_argument("--outdir", default="./results_gse198315")
    p2.add_argument("--n-hvg", type=int, default=2000); p2.add_argument("--n-pcs", type=int, default=50)
    p2.add_argument("--resolution", type=float, default=1.0)
    p2.add_argument("--min-margin", type=float, default=0.05,
                    help="1·2위 마커점수 차가 이보다 작으면 Ambiguous 처리")
    p2.add_argument("--threads", type=int, default=4)
    p2.set_defaults(func=step2)

    p3 = sub.add_parser("step3", help="inferCNV 악성 판정 + 최종 대조표")
    p3.add_argument("--h5ad", required=True); p3.add_argument("--gtf", required=True)
    p3.add_argument("--outdir", default="./results_gse198315")
    p3.add_argument("--window", type=int, default=100)
    p3.add_argument("--cnv-sd", type=float, default=2.0)
    p3.add_argument("--min-mapped", type=int, default=5000,
                    help="GTF 좌표가 부여되어야 할 최소 유전자 수 (실데이터 기본 5000)")
    p3.set_defaults(func=step3)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
