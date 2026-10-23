# Revision analysis outputs — WBP5 (TCEAL9) in HNSCC

Analysis outputs, barcode lists, scripts and structural descriptors supporting the
revised manuscript (CSBJ-D-26-00923). All files here were produced during the
revision; nothing in the pre-existing top-level folders of this repository has
been modified.

## Contents

| Folder | Contents |
|---|---|
| `scripts/` | Analysis scripts for every result reported in the revision. `wbp5_rev_common.py` holds the shared statistical routines; the `GSE198315_step*` scripts run in numerical order. |
| `gse103322_puram/` | Discovery-cohort (GSE103322) signature scoring, within-patient correlations, transcriptome-wide correlation profiles, Puram gene-set enrichment, and the YAP/TAZ permutation (`P5_YAP_TAZ_permutation.csv`). |
| `gse198315/` | Secondary-cohort (GSE198315) re-annotation, compartment composition, inferCNV, differential expression and depth-adjusted analyses. Primary subsample, random seed 42. |
| `gse198315_seed2026/`, `gse198315_seed7/` | The same pipeline repeated under two further subsampling seeds. |
| `credibility/` | Leave-one-patient-out resampling, gene-family and Xq22 module-adjustment tests, cross-cohort profile concordance, and run metadata (`CR_run_info.txt`). |
| `verify/` | Independent recomputation of reported values against the raw tool outputs. |
| `barcodes/` | Cell barcodes retained in each 80,000-cell subsample, one file per seed. |
| `reference/` | Input gene lists, resolved parameters, AlphaFold predicted aligned error, and regenerated comparison tables. |
| `structure/` | Per-residue relative solvent accessibility of the AlphaFold model, with the script that computes it. |

## Key files

- `gse198315/K2_DE_tumor_vs_normal_epithelium_full.csv` — the complete ranked
  differential-expression table for tumour-region versus normal-tissue epithelial
  cells, all 33,694 genes, with effect sizes, adjusted p values and ranks.
- `gse103322_puram/P5_YAP_TAZ_permutation.csv` — the 10,000-permutation
  WBP5–YAP/TAZ test against expression- and detection-matched null gene sets.
- `gse198315/D2_depth_adjusted_correlations.csv` — within-patient marker
  correlations before and after adjustment for per-cell sequencing depth.
- `credibility/*/CR3_family_and_controls.csv` — TCEAL/WBP family specificity with
  negative and positive control genes.
- `structure/WBP5_RSA_per_residue.csv` — per-residue SASA and RSA for
  AF-Q9UHQ7-F1-model_v6, with the P₀ pocket span (residues 38–58) flagged.

## Data not deposited here

The following are redistributed by their original sources and are not duplicated
in this repository:

- GSE103322 expression matrix — NCBI GEO, accession GSE103322.
- GSE198315 count matrix, barcodes and features — NCBI GEO, accession GSE198315.
- GENCODE v50 basic annotation — GENCODE release 50.
- The ChEMBL-25 pre-trained prior used for de novo generation is distributed with
  REINVENT4 and is subject to its licence; it is not redistributed here.

Intermediate AnnData checkpoints (`.h5ad`) are not deposited because of their size.
They are reproducible from the deposited scripts together with the barcode lists in
`barcodes/`, which fix the exact cells used in each subsample.

## Reproducing the analyses

Package versions for each stage are recorded in the `environment*.txt` files inside
the corresponding output folders. Random seeds are fixed in the scripts; the primary
subsample uses seed 42, with seeds 2026 and 7 used for the reproducibility checks.
