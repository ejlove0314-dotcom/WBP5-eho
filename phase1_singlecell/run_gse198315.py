#!/usr/bin/env python3
"""
GSE198315 (OSCC — Oral Squamous Cell Carcinoma) single-cell UMAP reconstruction.
Marker-panel subsampled analysis.

Usage
-----
    python run_gse198315.py --indir . --outdir out

Expected input (in --indir):
    GSE198315_matrix.mtx.gz       (~6.8 GB, sparse 10x matrix)
    GSE198315_barcodes.tsv.gz     (~1.2 MB, cell barcodes with P##-REGION_ prefix)
    GSE198315_features.tsv.gz     (~239 KB, gene list — 10x format)

Pipeline
--------
1. Load 10x triplet via scanpy.read_10x_mtx
2. Parse patient & region from barcode prefix (P##-{NT|TP|TC|mLN}_...)
3. QC filter (min_genes, min_cells, max_mt%)
4. Optional cell subsampling (--max_cells; default 80,000 for speed)
5. Marker-panel subsampling: keep only a curated list of 4-compartment markers
6. Compute compartment score per cell → assign proxy compartment
7. PCA + UMAP on marker panel (faster, cleaner than whole transcriptome)
8. 4-panel figure (compartments / regions / WBP5 / EGFR)

Requires
--------
    scanpy, anndata, numpy, pandas, matplotlib  (already in `scrna` conda env)
"""
from __future__ import annotations
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# -------------------------------------------------------------------
# Marker panel for 4 proxy compartments (curated, OSCC-relevant)
# -------------------------------------------------------------------
MARKERS = {
    "Epithelial-like": [
        "EPCAM", "KRT5", "KRT6A", "KRT6B", "KRT14", "KRT15", "KRT16", "KRT17",
        "KRT1", "KRT10", "KRT4", "KRT13", "KRT19", "SFN", "KRT8", "KRT18",
        "CDH1", "CLDN4", "CLDN7", "TP63", "SOX2", "PERP", "S100A2", "S100A8", "S100A9",
    ],
    "Immune-like": [
        "PTPRC",  # CD45 — pan-immune
        "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "TRAC", "TRBC1", "TRBC2",  # T cells
        "CD79A", "CD79B", "MS4A1", "IGHM", "IGHG1", "IGKC", "IGLC2", "JCHAIN",  # B/plasma
        "CD68", "CD163", "CD14", "C1QA", "C1QB", "C1QC", "LYZ", "AIF1",  # Myeloid/macs
        "FCGR3A", "FCGR3B", "S100A12",  # Monocytes/neutrophils
        "TPSAB1", "TPSB2", "CPA3", "KIT",  # Mast
        "CLEC9A", "CLEC10A", "CD1C", "CLEC4C",  # DC
        "NKG7", "GNLY", "KLRD1", "GZMB", "GZMA", "PRF1",  # NK / cytotoxic
    ],
    "Stromal-like": [
        "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3",
        "DCN", "LUM", "PDPN", "PDGFRA", "PDGFRB", "THY1",
        "ACTA2", "MYH11", "TAGLN",  # smooth muscle / myofibroblasts
        "FAP", "POSTN", "CXCL12",  # CAF
        "MCAM", "NOTCH3", "RGS5",  # pericyte
    ],
    "Endothelial-like": [
        "PECAM1", "VWF", "CDH5", "CLDN5", "ENG", "KDR", "FLT1",
        "CLEC14A", "ESAM", "ERG", "EMCN", "TIE1", "TEK", "PLVAP",
        "PROX1", "LYVE1", "CCL21",  # lymphatic endothelial
    ],
}

COMP_COLORS = {
    "Epithelial-like":  "#E07B6B",   # salmon
    "Immune-like":      "#6FC7D1",   # cyan
    "Stromal-like":     "#36A48F",   # teal-green
    "Endothelial-like": "#3C4D9E",   # deep blue
    "Unclassified":     "#C8B899",   # beige
}
COMP_ORDER = ["Epithelial-like", "Immune-like", "Stromal-like",
              "Endothelial-like", "Unclassified"]

REGION_COLORS = {
    "NT":  "#A3D9E8",   # light cyan
    "TP":  "#3E998E",   # teal
    "TC":  "#E86859",   # red
    "mLN": "#8B5E46",   # brown
}
REGION_ORDER = ["NT", "TP", "TC", "mLN"]


# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, type=Path,
                    help="Folder containing matrix.mtx.gz + barcodes.tsv.gz + features.tsv.gz")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--max_cells", type=int, default=80000,
                    help="Random subsample to this many cells after QC (default 80000)")
    ap.add_argument("--min_genes", type=int, default=500,
                    help="Min genes per cell for QC filter")
    ap.add_argument("--max_mt_pct", type=float, default=20.0,
                    help="Max mitochondrial percentage per cell")
    ap.add_argument("--feature_gene", type=str, default="WBP5",
                    help="Gene for Panel C feature plot")
    ap.add_argument("--feature_gene2", type=str, default="EGFR",
                    help="Gene for Panel D feature plot")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rng = np.random.default_rng(args.seed)

    # -------- 1. Load 10x triplet --------
    import scanpy as sc
    import anndata as ad

    # Rename files temporarily if needed: scanpy.read_10x_mtx expects
    #   matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz (without GSE prefix)
    # We'll use scanpy's read_mtx + manual attach for reliability.
    from scipy.io import mmread
    import gzip

    print(f"[load] reading matrix from {args.indir}")
    mtx_path = args.indir / "GSE198315_matrix.mtx.gz"
    bc_path  = args.indir / "GSE198315_barcodes.tsv.gz"
    ft_path  = args.indir / "GSE198315_features.tsv.gz"
    for p in [mtx_path, bc_path, ft_path]:
        if not p.exists():
            raise SystemExit(f"ERROR: missing file {p}")

    print("[load] reading sparse matrix (this may take 2-5 minutes)...")
    with gzip.open(mtx_path, "rb") as f:
        X = mmread(f).tocsr().astype(np.float32)
    print(f"[load] matrix shape: {X.shape}, nnz={X.nnz:,}, {time.time()-t0:.1f}s")

    print("[load] reading barcodes & features...")
    barcodes = pd.read_csv(bc_path, sep="\t", header=None).iloc[:, 0].astype(str).values
    features = pd.read_csv(ft_path, sep="\t", header=None)
    # features columns: gene_id, gene_symbol, feature_type
    if features.shape[1] >= 2:
        gene_names = features.iloc[:, 1].astype(str).values
    else:
        gene_names = features.iloc[:, 0].astype(str).values
    print(f"[load] n_cells={len(barcodes):,}  n_genes={len(gene_names):,}")

    # Matrix orientation: 10x convention is genes x cells. Transpose to cells x genes.
    if X.shape[0] == len(gene_names) and X.shape[1] == len(barcodes):
        X = X.T.tocsr()
    elif X.shape[0] == len(barcodes) and X.shape[1] == len(gene_names):
        pass
    else:
        raise SystemExit(
            f"Matrix shape {X.shape} doesn't match barcodes ({len(barcodes)}) "
            f"or features ({len(gene_names)})"
        )

    # -------- 2. Parse barcode metadata (patient + region) --------
    # barcode format: P##-{REGION}_{bases} e.g., P01-NT_AAACCCAAG...
    print("[parse] extracting patient & region from barcodes...")
    prefixes = pd.Series(barcodes).str.split("_").str[0]  # P01-NT
    patients = prefixes.str.split("-").str[0].values      # P01
    regions  = prefixes.str.split("-").str[1].values      # NT/TP/TC/mLN

    obs = pd.DataFrame({
        "barcode":  barcodes,
        "patient":  patients,
        "region":   regions,
    }).set_index("barcode")
    print(f"[parse] region counts:\n{obs['region'].value_counts().to_string()}")
    print(f"[parse] patient counts:\n{obs['patient'].value_counts().head(15).to_string()}")

    # -------- 3. Build AnnData & QC filter --------
    var = pd.DataFrame(index=gene_names)
    var.index.name = "gene"
    # Drop duplicate gene symbols (keep first)
    dup_mask = var.index.duplicated(keep="first")
    if dup_mask.any():
        print(f"[qc] dropping {dup_mask.sum()} duplicate gene symbols")
        X = X[:, ~dup_mask]
        var = var.loc[~dup_mask]

    adata = ad.AnnData(X=X, obs=obs, var=var)
    del X  # free memory
    print(f"[qc] adata: {adata.shape}")

    # QC metrics (scanpy handles sparse efficiently)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                                log1p=False, inplace=True)
    print(f"[qc] pre-filter: {adata.n_obs:,} cells")
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=10)
    adata = adata[adata.obs["pct_counts_mt"] < args.max_mt_pct].copy()
    print(f"[qc] post-filter: {adata.n_obs:,} cells, {adata.n_vars:,} genes "
          f"({time.time()-t0:.1f}s)")

    # -------- 4. Optional cell subsampling --------
    if args.max_cells and adata.n_obs > args.max_cells:
        print(f"[subsample] {adata.n_obs:,} -> {args.max_cells:,} cells (random)")
        idx = rng.choice(adata.n_obs, size=args.max_cells, replace=False)
        adata = adata[np.sort(idx)].copy()

    # -------- 5. Normalize & log1p (for marker scoring + feature plots) --------
    print("[norm] normalize_total + log1p")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # adata.X is now log1p(CPM-like) values

    # -------- 6. Compartment scoring (marker-based) --------
    print("[score] computing compartment scores")
    genes_present = set(adata.var_names)
    score_cols = []
    for comp, markers in MARKERS.items():
        use = [g for g in markers if g in genes_present]
        print(f"  {comp:20s}: {len(use)}/{len(markers)} markers present")
        sc.tl.score_genes(adata, gene_list=use, score_name=f"score_{comp}",
                          random_state=args.seed)
        score_cols.append(f"score_{comp}")

    # Assign compartment = argmax score, with "Unclassified" if best score is low
    scores = adata.obs[score_cols].values
    best_idx = scores.argmax(axis=1)
    best_val = scores[np.arange(len(scores)), best_idx]
    comp_names = [c.replace("score_", "") for c in score_cols]
    comp_assigned = np.array([comp_names[i] for i in best_idx])
    # Threshold: if best score is below a small positive value → Unclassified
    comp_assigned[best_val < 0.05] = "Unclassified"
    adata.obs["compartment"] = pd.Categorical(comp_assigned, categories=COMP_ORDER)

    print("[score] compartment assignment:")
    print(adata.obs["compartment"].value_counts().to_string())

    # -------- 7. Build marker-panel subsampled view for UMAP --------
    all_markers = sorted(set(m for ms in MARKERS.values() for m in ms))
    panel = [g for g in all_markers if g in genes_present]
    print(f"[panel] marker-panel genes: {len(panel)}")
    ad_panel = adata[:, panel].copy()

    sc.pp.scale(ad_panel, max_value=10)
    n_pcs = min(50, len(panel) - 1)
    sc.tl.pca(ad_panel, n_comps=n_pcs, random_state=args.seed)
    print("[umap] neighbors + UMAP...")
    sc.pp.neighbors(ad_panel, n_neighbors=20, use_rep="X_pca",
                    random_state=args.seed)
    sc.tl.umap(ad_panel, min_dist=0.3, spread=1.0, random_state=args.seed)

    adata.obsm["X_umap"] = ad_panel.obsm["X_umap"]
    Y = adata.obsm["X_umap"]
    print(f"[umap] done, Y shape {Y.shape}, {time.time()-t0:.1f}s")

    # Save outputs
    pd.DataFrame(Y, columns=["UMAP1", "UMAP2"],
                 index=adata.obs.index).to_csv(args.outdir / "umap_coords.csv")
    # Extract feature gene expression and save to obs so replotting doesn't need matrix
    for gene in [args.feature_gene, args.feature_gene2]:
        if gene in adata.var_names:
            expr = adata[:, gene].X
            expr = np.asarray(expr.todense()).ravel() if hasattr(expr, "todense") else np.asarray(expr).ravel()
            adata.obs[f"expr_{gene}"] = expr
    adata.obs.to_csv(args.outdir / "obs.csv")

    # -------- 8. 4-panel figure --------
    print("[plot] building 4-panel figure")
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.42,
                          left=0.06, right=0.97, top=0.96, bottom=0.06)
    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[1, 0])
    ax_D = fig.add_subplot(gs[1, 1])

    def finalize(ax, title, panel_letter):
        ax.set_xlabel("UMAP1", fontsize=12)
        ax.set_ylabel("UMAP2", fontsize=12)
        ax.set_title(title, fontsize=14)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.12, 1.04, panel_letter, transform=ax.transAxes,
                fontsize=20, fontweight="bold", va="top")

    def scatter_cat(ax, values, palette, order, s=3):
        for c in order:
            m = (values == c)
            if not m.any():
                continue
            ax.scatter(Y[m, 0], Y[m, 1], s=s, c=palette.get(c, "#cccccc"),
                       linewidth=0, alpha=0.75)

    # --- A: Proxy compartments ---
    scatter_cat(ax_A, adata.obs["compartment"].values, COMP_COLORS, COMP_ORDER)
    finalize(ax_A, "Proxy compartments", "A")
    ax_A.legend(handles=[
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=COMP_COLORS[c], markersize=9, label=c)
        for c in COMP_ORDER
    ], loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True,
       fontsize=9, borderpad=0.5, labelspacing=0.45)

    # --- B: Regions (simple inline legend, frameless, GSE103322 Panel B style) ---
    scatter_cat(ax_B, adata.obs["region"].values, REGION_COLORS, REGION_ORDER)
    finalize(ax_B, "Regions", "B")
    ax_B.legend(handles=[
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=REGION_COLORS[c], markersize=9, label=c)
        for c in REGION_ORDER
    ], loc="upper right", bbox_to_anchor=(1.02, 1.0),
       frameon=False, fontsize=10)

    # --- C/D: feature genes ---
    for ax, gene, letter in [(ax_C, args.feature_gene, "C"),
                              (ax_D, args.feature_gene2, "D")]:
        if gene in adata.var_names:
            expr = adata[:, gene].X
            expr = np.asarray(expr.todense()).ravel() if hasattr(expr, "todense") else np.asarray(expr).ravel()
            order_idx = np.argsort(expr)
            sc_pl = ax.scatter(Y[order_idx, 0], Y[order_idx, 1], s=3,
                               c=expr[order_idx], cmap="viridis", linewidth=0,
                               alpha=0.9, vmin=0,
                               vmax=max(0.1, np.percentile(expr, 99.5)))
            cbar = plt.colorbar(sc_pl, ax=ax, fraction=0.045, pad=0.02)
            cbar.set_label(f"{gene} (log1p counts)", fontsize=9)
        else:
            ax.text(0.5, 0.5, f"Gene '{gene}' not found",
                    ha="center", va="center", transform=ax.transAxes)
        finalize(ax, f"{gene} expression", letter)

    # No suptitle / footer — match clean GSE103322 style
    out_png = args.outdir / "GSE198315_umap_4panel.png"
    out_pdf = args.outdir / "GSE198315_umap_4panel.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"[done] saved:\n  {out_png}\n  {out_pdf}")
    print(f"[done] total time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
