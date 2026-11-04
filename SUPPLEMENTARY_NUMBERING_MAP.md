# Supplementary numbering: submitted version versus revised version

Supplementary tables and figures were renumbered during revision so that they follow the order
in which they are first cited in the main text. Eighteen of the twenty tables and four of the
seven figures changed number.

**Anything in this repository that names a supplementary table or figure refers to the revised
numbering, except where a file explicitly describes the submitted version.** Script and data
filenames were deliberately *not* renamed: `CR_run_info.txt` records the SHA-256 of the scripts
it points to, so renaming a script would break a verification a reader can run. The filename
table below gives the correspondence instead.

## Tables

| Submitted | Revised | Content |
|---|---|---|
| Table S1 | **Table S1** | Marker-gene panel for proxy compartment assignment |
| Table S2 | **Table S4** | Summary of transcriptomic datasets |
| Table S3 | **Table S5** | Per-cell sequencing-depth correction |
| Table S4 | **Table S6** | Bulk-cohort co-expression, TCGA HNSC |
| Table S5 | **Table S8** | SMILES, SwissADME and ProTox profiles |
| Table S6 | **Table S12** | Compartment-wise distribution of WBP5 |
| Table S7 | **Table S3** | GSE198315 high-expression stratification |
| Table S8 | **Table S13** | Approved-drug docking at P₀ |
| — | Table S2 | Leiden cluster annotation *(added at revision)* |
| — | Table S7 | Per-residue solvent accessibility and pLDDT *(added at revision)* |
| — | Table S9 | mLN versus primary by compartment *(added at revision)* |
| — | Table S10 | Leave-one-patient-out stability *(added at revision)* |
| — | Table S11 | Transcriptome-wide null percentiles *(added at revision)* |
| — | Table S14 | P₀ structural features *(added at revision)* |
| — | Table S15 | De novo candidate profiles *(added at revision)* |
| — | Table S16 | TCEAL family specificity *(added at revision)* |
| — | Table S17 | Xq22 module adjustment *(added at revision)* |
| — | Table S18 | Per-region intrinsic disorder *(added at revision)* |
| — | Table S19 | PAE within and between segments *(added at revision)* |
| — | Table S20 | Pairwise PAE *(added at revision)* |

## Figures

| Submitted | Revised | Content |
|---|---|---|
| Fig. S1 | **Fig. S1** | Multi-resolution WBP5–p-EMT correlation analysis |
| Fig. S2 | **Fig. S3** | AlphaFold2 per-residue confidence |
| Fig. S3 | **Fig. S2** | Purity-adjusted co-expression, TCGA HNSC |
| Fig. S4 | **Fig. S5** | TCEAL/WBP family paralogs |
| Fig. S5 | **Fig. S4** | YAP/TAZ target signature |
| — | Fig. S6 | Structural modeling and de novo design *(relocated from main-text Fig. 6)* |
| — | Fig. S7 | Predicted aligned error matrix *(added at revision)* |

## Filenames that carry a number

These names were kept as they are. The table says which table or figure each one now corresponds to.

| File | Now corresponds to |
|---|---|
| `tables/TableS3_gse198315_depth_correction.csv` | Table S5 |
| `tables/TableS4_timer3_hnsc_correlations.csv` | Table S6 |
| `verification/verify_TableS6_and_pseudobulk.py` | Table S12 |
| `revision_outputs/scripts/regen_tableS6.py` | Table S12 |
| `revision_outputs/reference/TableS6_regenerated_20261022.csv` | Table S12 |
| `revision_outputs/reference/TableS6_comparison_20261022.csv` | Table S12 |
| `figures/make_figS1_revised.py` | Fig. S1 |
| `figures/make_figS3_from_timer3.py` | Fig. S2 |
| `revision_outputs/scripts/make_figS2_v2.py` | Fig. S3 |
| `revision_outputs/scripts/make_figS4.py` | Fig. S5 |
| `deprecated/make_figS4.py` | Fig. S5 (superseded script) |
| `revision_outputs/scripts/figS5_300dpi.py` | Fig. S4 |
| `revision_outputs/scripts/fig_S5_regenerate_20261021.py` | Fig. S4 |
| `revision_outputs/scripts/make_figS6.py` | Fig. S6 |
| `revision_outputs/scripts/make_figS7_PAE.py` | Fig. S7 |

One filename is **not** part of this mapping:
`revision_outputs/reference/pEMT_gene_list_Puram_TableS7.xlsx`. The "TableS7" there is
Table S7 of Puram et al., the published source of the p-EMT signature. It has nothing to do
with the supplementary tables of this manuscript.

## Convention going forward

Each supplementary table legend names the deposited file or files from which its values are
derived, so that a reader can go from a printed number to the data product without consulting
this map. Table S11 is the reference example: its legend names
`revision_outputs/gse103322_puram/P3_transcriptome_rho_pooled.csv.gz` and
`P3_Puram_enrichment.csv`, and every value in the table can be recomputed from those two files.
