# Phase I — Single-cell analysis (GSE198315)

Scripts used for the single-cell transcriptomic component of the WBP5 (TCEAL9)
study, covering the multiregional OSCC cohort GSE198315.

## Input data

Download from GEO accession **GSE198315**:

| Script | Required input |
|---|---|
| `run_gse198315.py` | `GSE198315_matrix.mtx.gz`, `GSE198315_barcodes.tsv.gz`, `GSE198315_features.tsv.gz` (10x triplet, placed in one folder) |
| `make_gse198315_wbp5_marker_pipeline.py` | `GSE198315_OSCC_UMI_count_matrix.txt.gz` (gene x cell, tab-separated) |

If the supplementary files are distributed uncompressed, gzip them first, or
adjust the filenames at the top of the corresponding script.

Cell barcodes are expected in the format `P##-REGION_BASES`
(e.g. `P01-NT_AAACCCAAGGGAGAAT`), from which patient identity and anatomical
region (NT, TP, TC, mLN) are parsed.

## Scripts

### `run_gse198315.py`

Quality control, proxy compartment assignment, UMAP embedding and the
four-panel landscape figure.

```
python run_gse198315.py --indir /path/to/GSE198315 --outdir out
```

Key defaults (all overridable): `--min_genes 500`, `--max_mt_pct 20`,
`--max_cells 80000`, `--seed 42`.

Pipeline: QC filtering (>= 500 genes per cell; mitochondrial fraction < 20%;
genes detected in >= 10 cells) -> random subsampling to 80,000 cells ->
`normalize_total(target_sum=1e4)` + `log1p` -> control-gene-corrected module
scoring over a curated 106-gene, four-programme marker panel -> compartment
assignment by argmax (unclassified if the maximum score does not exceed 0.05)
-> scaling, PCA, neighbourhood graph and UMAP computed on the marker-panel
genes only.

Outputs: `obs.csv` (per-cell metadata, module scores, compartment labels and
feature-gene expression), `umap_coords.csv`, and
`GSE198315_umap_4panel.{png,pdf}`.

### `make_gse198315_wbp5_marker_pipeline.py`

WBP5 co-expression with p-EMT and invasion markers, at cell level and at
patient-mean level, across five cell-population strata.

```
python make_gse198315_wbp5_marker_pipeline.py \
    --data /path/to/GSE198315_OSCC_UMI_count_matrix.txt.gz \
    --outdir out
```

Outputs: `gse198315_wbp5_marker_correlations.tsv`,
`gse198315_wbp5_marker_correlations_summary.tsv`,
`gse198315_wbp5_marker_scatter_filtered_stats.tsv`,
`gse198315_wbp5_combined_summary.tsv`,
`gse198315_wbp5_marker_pipeline_report.json`, and the corresponding figures.

## Two proxy schemes

The two scripts use different, deliberately separate marker-based proxies,
because neither dataset ships with author-provided malignant-cell annotation:

- **Compartment proxy** (`run_gse198315.py`): four-programme module scoring over
  a 106-gene panel, assigned by argmax. Used for the compartment-resolved
  analyses.
- **Proxy tumor epithelial** (`make_gse198315_wbp5_marker_pipeline.py`):
  non-zero expression of any of EPCAM, KRT14, KRT17 or LAMC2, restricted to
  tumour-associated regions (TP, TC, mLN). Used for the marker co-expression
  analyses.

Both are approximations and do not establish compartment identity.

## Environment

Python 3.10+ with `scanpy`, `anndata`, `numpy`, `pandas`, `scipy`,
`matplotlib`, `seaborn` and `scienceplots`. Random seeds are fixed at 42
throughout.
