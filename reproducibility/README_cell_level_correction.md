# Correction to `reproducibility/GSE103322_WBP5_cell_level_raw.csv`

This file replaces the version deposited previously. Three defects were present
in the earlier file. None of them affects any value reported in the manuscript
or the Supplementary Information; all three were introduced when the file was
written out, not in the analysis.

## 1. Eighteen rows carried `#NAME?` in `non_cancer_type` and `compartment`

The GSE103322 annotation string for this category is `-Fibroblast`. Because the
string begins with a hyphen, a spreadsheet application parses it as a formula
and writes back `#NAME?` on save. The category is the one tabulated in
Table S12 as **Fibroblast (excl.)** (n = 18 in the full file; n = 17 after the
pooled-plate exclusion described below).

In this version both columns carry the Table S12 label `Fibroblast (excl.)`.
The original GEO string is recorded here rather than in a column of its own,
because a column holding `-Fibroblast` would reintroduce the same failure mode.

A reviewer recomputing the cancer-versus-fibroblast contrast from the earlier
file obtained n = 1,422 rather than the n = 1,440 stated in Section 3.1.1.
The manuscript figure combines the two fibroblast rows of Table S12
(1,422 + 18); that combination is stated in the text and is reproducible from
this file by selecting `compartment in {"Fibroblast", "Fibroblast (excl.)"}`.

## 2. One hundred and seventy-three rows carried an empty `patient`

These cells use the identifier form `HNSCC_17_P10_G09_S273_comb`, with an
underscore after the cohort prefix, whereas the remaining cells of the same
patient use `HNSCC17_...`. The identifier parser that wrote the file matched
only the second form, so these cells — all T cells of patient 17 — lost their
patient assignment.

The analysis retained them. With these cells included, the
malignant-versus-pooled-non-malignant paired comparison of Section 3.1.1
returns the reported mean difference of 0.79 log₂ units, Cohen's d_z = 0.73 and
an exact Wilcoxon signed-rank p = 0.049; without them it returns 0.72. The
published values are therefore the correct ones, and `patient` is set to 17 for
these rows in this version.

## 3. A further 109 rows carry an empty `patient` by design

These are the pooled "combo" plate cells that Methods Section 2.2 excludes
(15 malignant, 94 non-malignant). They remain unassigned, and are now flagged
in a new `pooled_plate` column so that the exclusion is reproducible from the
file alone: filtering to `pooled_plate == FALSE` returns 2,200 malignant and
3,593 non-malignant cells, the counts stated in Methods Section 2.2.

## Other changes

`compartment` and `non_cancer_type` are normalized to the labels used in
Table S12, so that the table can be reproduced from this file without an
external mapping. Specifically, `0` and empty annotations become `Unannotated`,
`myocyte` becomes `Myocyte`, and `non_cancer_type` is left empty for malignant
cells rather than holding `0`. No numeric value in the file was changed.

## Values reproduced from this file

| Quantity | Reported | From this file |
|---|---|---|
| Table S12, all eleven compartment rows (n, mean, PCT+) | — | reproduced exactly |
| Cancer vs. all non-malignant (Cliff's δ, p) | +0.308, 1.3×10⁻¹⁰⁵ | +0.3075, 1.33×10⁻¹⁰⁵ |
| Cancer vs. fibroblast, both rows combined | n = 1,440, δ = −0.023, p = 0.23 | n = 1,440, δ = −0.0228, p = 0.23 |
| Cancer vs. endothelial | n = 260, δ = +0.010, p = 0.78 | n = 260, δ = +0.0104, p = 0.78 |
| Cancer vs. T / dendritic / B / mast / macrophage | +0.675 / +0.662 / +0.570 / +0.528 / +0.439 | +0.6750 / +0.6623 / +0.5702 / +0.5284 / +0.4394 |
| Patient-level paired difference | 0.79, d_z = 0.73, p = 0.049 | 0.7937, d_z = 0.7280, p = 0.048828 |
| Pooled-plate exclusion | 109 excluded; 2,200 / 3,593 retained | 109; 2,200 / 3,593 |

The script that produced this file from the earlier version is deposited
alongside it as `fix_cell_level.py`.

## Avoiding recurrence

The corruption arises only when a CSV holding a label that begins with `=`,
`+`, `-` or `@` is opened and saved in a spreadsheet application. No label in
this version begins with any of those characters, and the file should be edited
with a text editor or a script rather than a spreadsheet.
