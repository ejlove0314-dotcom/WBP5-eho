# Correction to the deposited GSE103322 cell-level files

Two cell-level files are deposited under `reproducibility/`:

    GSE103322_WBP5_cell_level_raw.csv          WBP5 only
    GSE103322_WBP5_EMT_markers_cell_level.csv  WBP5 plus ITGA5, LAMC2, VIM, KRT14, KRT17

Both were written out by the same export step and both carried the same three
defects. None of them affects any value reported in the manuscript or the
Supplementary Information; all three were introduced when the files were
written, not in the analysis.

**An earlier release repaired only the first of the two files.** The repair
script tested for a missing patient with `not field.strip()`, which is the
spelling used in `..._cell_level_raw.csv`; `..._EMT_markers_cell_level.csv`
writes the string `NA` instead, so the test never fired and the script reported
success without changing anything. Both files are repaired in this release, and
the script now treats an empty field, `NA`, `NaN`, `N/A` and `null` alike. Its
checks were also rewritten as post-conditions on the output, so the script is
idempotent: running it on an already-repaired file leaves the file byte for byte
unchanged and reports zero repairs.

## 1. Eighteen rows carried `#NAME?` in `non_cancer_type` and `compartment`

The GSE103322 annotation string for this category is `-Fibroblast`. Because the
string begins with a hyphen, a spreadsheet application parses it as a formula
and writes back `#NAME?` on save. The category is the one tabulated in
Table S12 as **Fibroblast (excluded)**; the short form **Fibroblast (excl.)** is
the label carried in these files and shown in the Fig. 2A legend, where the full
form does not fit. It holds n = 18 rows in the full file and n = 17 after the
pooled-plate exclusion described in point 3.

The original GEO string is recorded here rather than in a column of its own,
because a column holding `-Fibroblast` would reintroduce the same failure mode.

A reviewer recomputing the cancer-versus-fibroblast contrast from an
unrepaired file obtains n = 1,422 rather than the n = 1,440 stated in
Section 3.1.1, because the eighteen rows appear under `#NAME?` instead of under
a fibroblast label. The manuscript figure combines the two fibroblast rows of
Table S12 (1,422 + 18); that combination is stated in the text and is
reproducible from either repaired file by selecting
`compartment in {"Fibroblast", "Fibroblast (excl.)"}`.

## 2. One hundred and seventy-three rows had lost their patient identifier

These cells use the identifier form `HNSCC_17_P10_G09_S273_comb`, with an
underscore after the cohort prefix, whereas the remaining cells of the same
patient use `HNSCC17_...`. The identifier parser that wrote the files matched
only the second form, so these cells — all T cells of patient 17 — lost their
patient assignment.

The analysis retained them. With these cells included, the
malignant-versus-pooled-non-malignant paired comparison of Section 3.1.1
returns the reported mean difference of 0.79 log₂ units, Cohen's d_z = 0.73 and
an exact Wilcoxon signed-rank p = 0.049; without them it returns 0.72. The
published values are therefore the correct ones, and `patient` is set to 17 for
these rows in both files.

Treating the first underscore-delimited token of a cell name as the patient
identifier yields 21 apparent identifiers for the 17 patients actually present:
patients 25, 26 and 28 each appear under both an `HN…` and an `HNSCC…` prefix,
and patient 17 appears under the additional `HNSCC_17_…` form described above.

## 3. A further 109 rows carry no patient identifier by design

These are the pooled "combo" plate cells that Methods Section 2.2 excludes
(15 malignant, 94 non-malignant). They remain unassigned, and are flagged in a
`pooled_plate` column so that the exclusion is reproducible from the file alone:
filtering to `pooled_plate == FALSE` returns 2,200 malignant and 3,593
non-malignant cells, the counts stated in Methods Section 2.2.

## Other changes

`compartment` and `non_cancer_type` are normalized to the labels used in
Table S12 and in the Fig. 2A legend, so that the table can be reproduced from
either file without an external mapping. Specifically, `0` and empty annotations
become `Unannotated`, `myocyte` becomes `Myocyte`, and `non_cancer_type` is left
empty for malignant cells rather than holding `0`. A missing patient is left as
the file's own marker — empty in `..._cell_level_raw.csv`, `NA` in
`..._EMT_markers_cell_level.csv` — so that each file keeps the convention it was
deposited with. No numeric value in either file was changed.

After repair the two files agree on every column they share: compartment,
non-cancer type, patient and pooled-plate status are identical row for row.

## Values reproduced from these files

| Quantity | Reported | From these files |
|---|---|---|
| Table S12, all eleven compartment rows (n, mean, PCT+) | — | reproduced exactly |
| Cancer vs. all non-malignant (Cliff's δ, p) | +0.308, 1.3×10⁻¹⁰⁵ | +0.3075, 1.33×10⁻¹⁰⁵ |
| Cancer vs. fibroblast, both rows combined | n = 1,440, δ = −0.023, p = 0.23 | n = 1,440, δ = −0.0228, p = 0.23 |
| Cancer vs. endothelial | n = 260, δ = +0.010, p = 0.78 | n = 260, δ = +0.0104, p = 0.78 |
| Cancer vs. T / dendritic / B / mast / macrophage | +0.675 / +0.662 / +0.570 / +0.528 / +0.439 | +0.6750 / +0.6623 / +0.5702 / +0.5284 / +0.4394 |
| Patient-level paired difference | 0.79, d_z = 0.73, p = 0.049 | 0.7937, d_z = 0.7280, p = 0.048828 |
| Pooled-plate exclusion | 109 excluded; 2,200 / 3,593 retained | 109; 2,200 / 3,593 |
| Within-patient Spearman ρ, all five markers | unchanged by this repair | identical before and after |
| Leave-one-patient-out ranges, all five markers | unchanged by this repair | identical before and after |

The last two rows were checked by running the deposited
`verification/patient_level_correlation.py` and
`verification/verify_TableS6_and_pseudobulk.py` against the file before and
after repair: every printed coefficient, p-value, q-value and range is
identical, and the three saved CSV outputs match the copies already deposited
under `verification/`. The only difference in the printed output is that the
compartment table now names the eighteen rows `Fibroblast (excl.)` instead of
`#NAME?`, the 324 unannotated rows `Unannotated` instead of `0`, and the
nineteen myocytes `Myocyte` instead of `myocyte`.

The script that produced these files is deposited alongside them as
`fix_cell_level.py`.

## Avoiding recurrence

The corruption arises only when a CSV holding a label that begins with `=`,
`+`, `-` or `@` is opened and saved in a spreadsheet application. No label in
either file begins with any of those characters, and the files should be edited
with a text editor or a script rather than a spreadsheet.
