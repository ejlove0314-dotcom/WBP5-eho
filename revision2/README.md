# Second-round revision outputs — WBP5 (TCEAL9) in HNSCC

Scripts, decision rules, results and figure files added for the second revision of
manuscript CSBJ-D-26-00923 (editor and reviewer comments 1, 3 and 4). Everything here is
new; no file deposited with the first revision has been modified.

The two decision documents are in the top-level `reproducibility/` folder:

- `reproducibility/pre_analysis_plan_R2_20261116.md` — rules for reading the specificity
  and sensitivity analyses, written before the script was run on either cohort.
- `reproducibility/pre_analysis_plan_R2_addendum_20261116.md` — addendum for the published
  CAF signatures. It was written after the pre-specified analysis had been read and before
  the added analysis was run, so that analysis is a post hoc sensitivity analysis and is
  reported as one.

## Contents of this folder

| Folder | Contents |
|---|---|
| `scripts/` | `R2_specificity_sensitivity_20261116b.py` (pre-specified analysis), `R2_specificity_sensitivity_20261116c.py` (the same script with the `--extra` option used for the CAF sensitivity analysis), `R2_null_diagnostic_20261117.py`, the figure scripts, and `make_full_from_annotated.py`. |
| `specificity/<run>/` | Output of revision b, one folder per run: `gse103322` (discovery cohort) and `seed42`, `seed2026`, `seed7` (the three GSE198315 subsamples). |
| `caf_sensitivity/<run>/` | Output of revision c with `--extra` for the five published pan-CAF signatures, same four runs. |
| `genesets/` | The five pan-CAF gene lists and the three related-program gene sets as used. |
| `null_diagnostic/` | `R2_null_diagnostic.csv`. |
| `figures/` | `Figure_1.pdf`, `Figure_3.pdf`, `Figure_4.pdf` and the source values for Figures 3 and 4. |
| `logs/` | Console output of the eight runs. |
| `MANIFEST_sha256.txt` | SHA-256 checksum of every file in this deposit. |

## Result files and where they are reported

File names are the same in every run folder.

| File | Content | Reported in |
|---|---|---|
| `R2_0_reproduction.csv` | Agreement of the recomputed primary profile with the first-revision profile (`revision_outputs/credibility/<run>/CR1_profile_WBP5.csv.gz`) | — |
| `R2_1_patient_weights.csv` | Cells and share of the pooled estimate per patient | Table S10, Section A |
| `R2_2_null_and_weighting.csv` | p-EMT signature enrichment under the expression-matched and the expression- and detection-matched null, for the pooled and the patient-balanced profile | Table S10, Section B; Fig. 4C |
| `R2_3_patient_level.csv` | Patient-level sign and Wilcoxon signed-rank tests | Table S10, Section C |
| `R2_4_downsample_all.csv`, `R2_4_downsample_summary.csv` | Capped down-sampling, 100 repeats | Table S10, Section D |
| `R2_5_program_sources.csv`, `R2_5_program_overlap.csv` | Gene sets, evaluable genes and pairwise overlap | Table S13 |
| `R2_6_program_comparison.csv` | Enrichment of intact and mutually exclusive gene sets | Table S14 |
| `R2_7_conditional_enrichment.csv`, `R2_8_module_partials.csv` | Conditional enrichment and module-level partial coefficients | Table S15 |
| `R2x_5_…` to `R2x_8_…` (in `caf_sensitivity/`) | The same four blocks for the pan-CAF signatures | Table S16 |
| `R2_CAF_derived_GSE103322.txt`, `…_table.csv` (in `specificity/gse103322/`) | The fibroblast-identity signature derived in GSE103322 and the statistics behind it | Table S13 |
| `null_diagnostic/R2_null_diagnostic.csv` | Detection rate of signature genes relative to expression-matched genes, and the expected mean of each null | Table S10, Section E |

`R2_run_info.txt` and `R2x_run_info.txt` record the arguments, package versions and the
SHA-256 checksums of the scripts as run. Those checksums match the scripts deposited in
`scripts/` (revisions b and c) and in `revision_outputs/scripts/`
(`CRED_profiles_20261021.py`, `wbp5_rev_common.py`). Paths in these files are the
authors' local working paths and are preserved as written.

## Script revisions

- **Revision a** (`R2_specificity_sensitivity_20261116.py`) was the first version run on
  the data. Its own reproduction check (block R2-0) showed that the routine used for the
  sensitivity profiles (down-sampling, conditional analysis) divided by sums of squares in
  both cohorts, whereas the primary profile of GSE198315 divides by the number of cells;
  the two differ for sparsely detected genes (maximum absolute difference 0.051 in the
  seed-42 run). The primary profile itself was not affected: it agrees with the first-revision
  profile to within 1e-15 in every run (`R2_0_reproduction.csv`).
- **Revision b** uses the normalization of the primary profile in each cohort for every
  profile. All pre-specified results reported in the manuscript come from revision b. The
  output of revision a is not deposited and is not used anywhere.
- **Revision c** adds the `--extra` option and nothing else to the computation; with
  `--extra` it runs blocks R2-5 to R2-8 for the added gene sets only and does not recompute
  the pre-specified output.

## Gene sets

- `genesets/galbo_pan_*.txt` — the five pan-cancer CAF subtype signatures of Galbo, Zang
  and Zheng, Clin Cancer Res 2021;27:2636–2647 (doi 10.1158/1078-0432.CCR-20-4226),
  Supplementary Table S3, as published. The head and neck data set used in that study is
  GSE103322, the discovery cohort of the present work, so these signatures are published
  but not independent of that cohort.
- `genesets/related_programs_used.gmt` — the three library gene sets as downloaded from
  Enrichr on 16 November 2026: MSigDB_Hallmark_2020 "Epithelial Mesenchymal Transition"
  and "TGF-beta Signaling", and GO_Biological_Process_2023 "Extracellular Matrix
  Organization (GO:0030198)".
- The p-EMT signature is that of Puram et al. (2017) and is deposited with the first
  revision.

## Null diagnostic

`R2_null_diagnostic_20261117.py` was written on 17 November 2026, after the pre-specified
sensitivity analyses had been run and read. It is descriptive: it tests nothing and
changes no result. It reads the four first-revision profiles
(`revision_outputs/credibility/<run>/CR1_profile_WBP5.csv.gz`).

## Figures

| File | Script | Input |
|---|---|---|
| `Figure_1.pdf` | `make_fig1_vector_20261117.py` | none (schematic) |
| `Figure_3.pdf`, `Figure_3_source.csv` | `make_fig3_leiden_20261116d.py` | the annotated GSE198315 object (seed 42) written by `GSE198315_reannotation_20261015.py` |
| `Figure_4.pdf`, `Figure_4B_source.csv`, `Figure_4C_source.csv` | `make_fig4_enrichment_20261117.py` | `revision_outputs/gse198315/D2_depth_adjusted_per_patient.csv` (first revision) and `specificity/<run>/R2_2_null_and_weighting.csv` |

Figures are vector PDF files written by Matplotlib. Figure 2 is unchanged from the first
revision. The figures numbered 3 and 4 in the first revision are Supplementary Figures S5
and S4 in this revision, and the former Figure 5 is the basis of the present Figure 4.

`make_full_from_annotated.py` rebuilds `<name>_full_annotated.h5ad` from
`<name>_annotated.h5ad` (the second output of step 2 of
`GSE198315_reannotation_20261015.py`), so that the full-matrix object can be recovered
without re-running the annotation.

## Logs

The console logs are saved as UTF-8; progress messages are partly in Korean.
