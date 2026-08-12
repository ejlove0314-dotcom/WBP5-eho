# Reproducibility — manuscript tables

Scripts and verified output files that reproduce the numerical values reported in
Tables 1-4 of the manuscript.

## Contents

| File | Description |
|------|-------------|
| `WBP5_Phase1_recalc.ipynb` | Colab notebook. Reproduces Table 1, Table 2, and the P0 contact analysis of Table 3. Despite the file name, it covers both the single-cell and structural analyses. |
| `WBP5_Table1_results_v2.csv` | Table 1 - WBP5 expression, cancer vs. non-cancer and primary vs. lymph node. |
| `WBP5_Table2_results_v2.csv` | Table 2 - Spearman correlations between WBP5 and eight markers. |
| `WBP5_Table2_meanPCT_v2.csv` | Table 2 - group means and detection rates for WBP5-positive and WBP5-negative malignant cells. |
| `WBP5_P0_contacts_seed42.csv` | Table 3 - residues within 4 A of each ligand in the top-scoring pose (seed 42). |
| `docking_replicates_P0_summary_v2.csv` | Table 4 - AutoDock Vina binding energies across five seeds, with the corrected Candidate_1 value. |
| `GSE103322_WBP5_cell_level_raw.csv` | Input. Per-cell WBP5 expression with cell-type, site and patient annotation. |
| `GSE103322_WBP5_EMT_markers_cell_level.csv` | Input. As above, with ITGA5, LAMC2, VIM, KRT14 and KRT17 added. |

EGFR, PTPRC and EPCAM are not present in the extracted input files; the notebook
retrieves them directly from the original GSE103322 expression matrix
(`GSE103322_HNSCC_all_data.txt.gz`, available from GEO).

## Notes

**Docking pocket.** All docking reported in the manuscript was performed at pocket P0
(grid center -23.835, -0.451, 23.910; box 25 x 25 x 25 A) in five replicates per ligand
(seeds 42, 1042, 2042, 3042, 4042). See `docking_P0/`.

**Candidate_1.** Candidate_1 was re-docked after its SMILES was corrected to the
pyridin-4-yl variant (MW 347.44). `docking_replicates_P0_summary_v2.csv` contains the
corrected value (-5.44 +/- 0.06 kcal/mol).

## Data availability

Source data are available from GEO under accessions GSE103322 and GSE198315.
