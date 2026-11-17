# Pre-analysis plan — second-round revision (CSBJ-D-26-00923R1)

Written 16 November 2026, before `R2_specificity_sensitivity_20261116.py` was run on either cohort.
The script was tested on simulated data only. The rules below are not to be changed after the results are seen;
any rule that is not met is reported as not met.

## Scope

Editor/reviewer point 1 (biological specificity of the WBP5–p-EMT association) and point 3
(patient weighting and the null model of the principal gene-set enrichment analysis).

## Fixed specifications (unchanged from the first round)

- Cells: malignant cells with an assignable patient identifier (GSE103322); marker-annotated epithelial cells from
  tumor-core and tumor-periphery samples (GSE198315; subsampling seeds 42, 2026 and 7).
- Profile: within-patient standardized ranks, first-order partial correlation adjusting for the technical covariate
  (log1p number of detected genes in GSE103322; log1p total counts in GSE198315).
- Primary signature: Puram p-EMT 100-gene list; secondary: its first 30 genes.
- One-sided permutation test of the set mean, 10,000 draws, seed 42.

## New specifications

| Element | Definition |
|---|---|
| Expression- and detection-matched null | 25 mean-expression bins, each split into quartiles of detection rate (100 strata); a stratum with fewer than 10 eligible genes falls back to its expression bin, and the number of fallbacks is reported |
| Patient-balanced profile | per-patient adjusted coefficients (patients with ≥ 20 cells), Fisher z-transformed and averaged with equal weight; genes evaluable in fewer than half of the patients are excluded |
| Patient-level test | per patient, set mean minus its matched expectation; sign test and Wilcoxon signed-rank test across patients |
| Capped down-sampling | every patient capped at the median cell count of the eligible patients; 100 repeats; 2,000 draws each |
| Related programs | Hallmark epithelial–mesenchymal transition and Hallmark TGF-β signaling (MSigDB_Hallmark_2020), GO:0030198 extracellular matrix organization (GO_Biological_Process_2023), and a fibroblast-identity signature derived in GSE103322 from cell-type labels only (top 100 genes by Cliff's δ, fibroblast versus malignant, detected in ≥ 25% of fibroblasts) |
| Disjoint sets | "p-EMT only" = p-EMT genes in no related program; "X only" = genes of program X not in the p-EMT list |
| Conditional analysis | enrichment of p-EMT-only genes after additionally adjusting for the module score of X-only genes, and the reverse; null draws exclude all program genes |

## Decision rules

### Point 3 — robustness of the principal finding (PEMT100)

The finding is described as robust to the null model and to patient weighting only if, in the discovery cohort
and in all three secondary-cohort subsamples:

1. permutation p ≤ 0.01 under the expression- and detection-matched null (pooled profile);
2. permutation p ≤ 0.01 for the patient-balanced profile under both nulls;
3. capped down-sampling: p < 0.05 in at least 95 of 100 repeats.

The patient-level test (positive excess in at least 8 of 10 patients, Wilcoxon p < 0.05) is reported for each
cohort as found and is not a condition of the claim, because with ten patients its power is limited.
If any of 1–3 fails, the manuscript states which specification the finding depends on.

### Point 1 — specificity

Evaluated in each cohort; for GSE198315 a criterion counts as met when it holds in at least two of three subsamples.

- (a) p-EMT-only genes are enriched (p ≤ 0.01);
- (b) p-EMT-only genes remain enriched after adjustment for the module of each related program (p ≤ 0.05 for every program);
- (c) for Hallmark EMT and for the fibroblast signature, the X-only genes are no longer enriched after adjustment for the p-EMT-only module (p > 0.05).

| Outcome | Wording adopted |
|---|---|
| (a), (b) and (c) met in both cohorts | the association is preferentially with the p-EMT program relative to the related programs tested |
| (a) and (b) met, (c) not met | the association includes a p-EMT component that is not explained by the related programs, within a broader EMT/mesenchymal co-expression pattern |
| (b) not met | the association cannot be separated from a broader EMT/mesenchymal program; title, Abstract and Conclusions are revised accordingly |
| cohorts disagree | reported per cohort; the weaker reading is adopted for the cross-cohort claim |

All programs, all four profile-by-null combinations and both directions of the conditional analysis are reported
whatever the outcome.
