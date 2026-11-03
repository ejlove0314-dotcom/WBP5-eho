# WBP5-eho

Code and data supporting the manuscript "Single-cell analysis identifies WBP5
(TCEAL9) as a partial-EMT program-associated gene shared by malignant and
mesenchymal compartments in HNSCC."

## Repository layout

| Path | Contents |
|------|----------|
| `verification/` | The numerical audit carried out while the manuscript was under review (87 ledger rows, July-August 2026), together with `README_verification_ledger_scope.md`, which records which of its locations the revision moved and which two values it revised. The recorded values are current; its section and table numbers are those of the submitted version. |
| `reproducibility/` | Scripts and verified outputs for manuscript Tables 1 and 2, the corrected cell-level file and the note describing its correction. |
| `revision_outputs/` | Analyses added during revision: the GSE103322 patient-level and permutation outputs (`gse103322_puram/`), the GSE198315 re-annotation and compartment analyses (`gse198315*/`), the four credibility analyses (`credibility/`), subsample barcodes, structural outputs, reference files and the scripts that produced them. |
| `tables/` | Source tables underlying Supplementary Tables S3 and S4. |
| `timer3/` | Raw TIMER3 exports underlying Table S4 and Fig. S3. |
| `figures/` | Scripts that regenerate Figs. S1 and S3. |
| `gse198315_recalc/` | Depth-corrected recomputation of the GSE198315 correlations. |
| `docking_P0/` | AutoDock Vina docking at pocket P0 - the protocol reported in the manuscript. |
| `docking_P0/mtiopenscreen/` | MTiOpenScreen Drugs-lib screen at the same pocket (job V06821248773098). |
| `phase1_singlecell/` | GSE198315 single-cell analysis pipeline. |
| `pockets/`, `residues/` | DoGSiteScorer pocket detection: run parameters, descriptor table, grid maps and pocket-lining residues. |
| `deprecated/` | Superseded scripts, ligand files and reports. Nothing here underlies a reported result, with one exception: Section 4.2 of `SwissADME_Analysis_Report.docx` is the source of Supplementary Table S5, Section D. See `deprecated/README.md`. |

Root-level files are the AlphaFold model (`AF-Q9UHQ7-F1-model_v6.pdb`, and the
identical working copy `WBP5_AlphaFold.pdb` used as the docking receptor), the
REINVENT4 configuration files (`rl_wbp5.toml`, `rl_wbp5_fixed.toml`,
`reinvent_wbp5_sampling.toml`), the corrected ligand set
(`all_candidates_corrected.smi`), the individual structure files for
Candidate_2 through Candidate_5, and the conda environment specification
(`environment.yml`).

## Environment

```bash
conda env create -f environment.yml
conda activate wbp5
```

Package versions match those reported in Methods 2.11. AutoDock Vina v1.2.5,
PyRx v0.8 and REINVENT4 v4.7.15 are not installable via conda and must be
obtained separately; the web servers used (AlphaFold DB, DoGSiteScorer,
MTiOpenScreen, SwissADME, ProTox 3.0, TIMER3) require no local installation.

## Ligand structures

`all_candidates_corrected.smi` holds the structures used for every docking, ADMET and
toxicity result reported in the manuscript. Each string was read back from the
corresponding docked pose in `docking_P0/` and re-parsed with RDKit 2026.03.5.

The molecular weights below are RDKit values; its periodic table gives sulfur
32.067. Supplementary Tables S5 and S10 report the SwissADME values for the same
structures, held in `reproducibility/swissadme_final_20260820.csv`. For the six
compounds both report, they agree on Candidate_4 and are 0.01 Da lower on the other
five; the difference arises in the two tools' atomic-mass tables and rounding, not in
the structures.

| ID | Formula | MW (Da) | InChIKey |
|----|---------|---------|----------|
| Candidate_1 | C17H21N3O3S | 347.44 | JGUIGLSIQCDQRA-UHFFFAOYSA-N |
| Candidate_2 | C19H23N5O2 | 353.43 | CWFOXUOVUNELQQ-UHFFFAOYSA-N |
| Candidate_3 | C17H17N3O3S | 343.41 | DSVBRFOHTGCUGB-UHFFFAOYSA-N |
| Candidate_4 | C19H21N3O4 | 355.39 | MMBJDRLBTDFZJL-UHFFFAOYSA-N |
| Candidate_5 | C16H19N3O4S | 349.41 | UQFXIMOKNSXXKZ-UHFFFAOYSA-N |
| Pazopanib_Ref | C21H23N7O2S | 437.53 | CUIHSIWYWATEQL-UHFFFAOYSA-N |
| Dasatinib | C22H26ClN7O2S | 488.02 | ZBNZXTGUTAYRHI-UHFFFAOYSA-N |
| Verteporfin | C41H42N4O8 | 718.81 | CABKTHJNHVBKCC-ZSFNYQMMSA-N |
| Simvastatin | C25H38O5 | 418.57 | RYMZZMVNJRMUDD-HGQWONQESA-N |
| Dobutamine | C18H23NO3 | 301.39 | JRWZLRBJNMZMFE-CYBMUJFWSA-N |

Two entries in the earlier file `all_candidates.smi` were incorrect. That file, the
corresponding Candidate_1 and Pazopanib_Ref structure files, and the input and report
files derived from them have been moved to `deprecated/`; see `deprecated/README.md`
for the discrepancies and how they were resolved. The structure files for Candidate_2
through Candidate_5 were unaffected and remain at the repository root.

## Numerical verification

Every numerical claim in the manuscript has been traced to a deposited source
file or recomputed from one. See [`verification/VERIFICATION.md`](verification/VERIFICATION.md)
for the narrative report and [`verification/verification_ledger.csv`](verification/verification_ledger.csv)
for the machine-readable ledger (manuscript location, quantity, manuscript
value, source file, source value, status).

To reproduce the derived statistics that are not stored in any export - the two
chi-square tests in Section 3.1.2, the Benjamini-Hochberg q values in Table S4,
the docking margins in Supplementary Note S1 Sections 3.5 and 3.7, and the
depth-correction deltas in Table S3:

```bash
python verification/verify_recompute.py
```

Two scripts require the full cell-level matrix:

```bash
python verification/verify_TableS6_and_pseudobulk.py GSE103322_WBP5_EMT_markers_cell_level.csv
python verification/patient_level_correlation.py GSE103322_WBP5_EMT_markers_cell_level.csv
```

The first checks the compartment-wise positivity rates reported in Table S6. The
second recomputes the Table 2 coefficients at patient resolution, both
within-patient and by pseudobulk aggregation over malignant cells, addressing
the pseudoreplication limitation discussed in Section 4.4. Its outputs are
deposited alongside it as `WBP5_within_patient_rho.csv`,
`WBP5_pseudobulk_rho.csv` and `WBP5_patient_means.csv`.

## TIMER3 bulk-cohort analysis

`timer3/` contains the raw exports underlying Table S4 and Fig. S3: ten
`genecorr_table-{gene}_{adjusted|raw}.csv` files and nine
`genecorr_plot {gene} {cohort}.jpg` scatter plots.

Query parameters: TIMER3 Gene_Correlation module, gene of interest `TCEAL9`,
TCGA HNSC, accessed 2026-08-27, each marker run twice with Purity Adjustment
on and off.

> **Note on the exported plots.** In TIMER3 output panels the x-axis titles are
> assigned by facet position rather than by content, so they are transposed
> whenever the queried gene sorts alphabetically before "Purity" (ITGA5, LAMC2,
> KRT14, KRT17). The facet strip titles are correct. Fig. S3 was assembled from
> these files with corrected axis titles by
> [`figures/make_figS3_from_timer3.py`](figures/make_figS3_from_timer3.py);
> the unmodified exports are deposited here as received.

## Figure regeneration

`figures/make_figS1_revised.py` regenerates Fig. S1 (multi-resolution
correlation analysis). Panel C carries TIMER3 purity-adjusted coefficients for
the TCGA HNSC cohort; bars that did not reach q < 0.05 after purity adjustment
are hatched and marked ns.

`figures/make_figS3_from_timer3.py` assembles Fig. S3 from the TIMER3 exports
in `timer3/`.

## Docking protocol

All docking values reported in the manuscript derive from pocket P0
(grid center -23.835, -0.451, 23.910; box 25 x 25 x 25 A; exhaustiveness 8),
in five replicates per ligand using fixed seeds 42, 1042, 2042, 3042 and 4042.

`pockets/wbp5alphafoldpdb*_desc.txt` is the DoGSiteScorer descriptor table from which
the pocket geometry reported in Table S9 is taken: for P0, volume 229.62 A^3, surface
811.54 A^2, depth 24.90 A, enclosure 0.36 and Drug Score 0.837. The same table gives
the corresponding values for the seven other detected pockets, none of which exceeds a
Drug Score of 0.40. `pockets/PocXlsDescriptors.txt` explains each column, and
`pockets/dogsitescorer_run_parameters.txt` records the server settings used.

## Exploratory analyses not reported in the manuscript

`deprecated/boltz2/` contains a ligand-protein co-folding run performed with Boltz-2
during an early exploratory phase. It is retained for transparency but is not cited in
the manuscript. The run recorded pose-confidence metrics (confidence score, ipTM, pLDDT)
only; the affinity-prediction module was not used. Under those metrics the positive
control (EGFR-erlotinib) behaved as expected, whereas six chemically unrelated
negative-control ligands returned higher scores than the reference binder at WBP5.
Because pose-confidence metrics are not intended to rank binding strength, the run was
treated as inconclusive rather than as evidence about the target, and no result from it
informed any conclusion in the manuscript.

## Data availability

Source data are available from GEO under accessions GSE103322 and GSE198315.
The AlphaFold model is AF-Q9UHQ7-F1-model_v6 (UniProt Q9UHQ7).

## Data files and their use in the manuscript

| File | Used in |
|---|---|
| `reproducibility/swissadme_final_20260820.csv` | Table S5 Sections B, C and E; Table S10 |
| `deprecated/SwissADME_Analysis_Report.docx`, Section 4.2 | Table S5 Section D |
| `reproducibility/docking_replicates_P0_summary_v2.csv` | Table S8; Table S10 |
| `reproducibility/ProTox-3.0 - Prediction of TOXicity of chemicals *.xlsx` (six files: candidates 1-5 and pazopanib) | Table S5 Sections F1-F2 |
| `phase1_singlecell/gse198315_wbp5_marker_correlations.tsv` | Fig. 5B; Fig. S1 |
| `phase1_singlecell/make_gse198315_wbp5_marker_pipeline.py` | Script generating the above TSV |
| `tables/TableS3_gse198315_depth_correction.csv` | Table S3 |
| `tables/TableS4_timer3_hnsc_correlations.csv` | Table S4 |
| `gse198315_recalc/headline.tsv` | Section 3.2; Table S3 |
| `gse198315_recalc/correlations_full.tsv` | Fig. S1A |

### Note on superseded files

Files in `deprecated/` correspond to an earlier analysis run and are retained
for provenance only. With one exception they are not the basis of any value
reported in the manuscript: the Candidate_1 and pazopanib structures were both
corrected after the initial SwissADME submission, and the physicochemical,
drug-likeness and ADMET values in the manuscript derive from
`reproducibility/swissadme_final_20260820.csv`.

The exception is the CYP450 predictions in Supplementary Table S5, Section D.
That module returned no output when the corrected structures were re-submitted
in August 2026 — the five CYP columns of `swissadme_final_20260820.csv` read
`n/d` for all six compounds — so the four rows Table S5 reports come from
Section 4.2 of `deprecated/SwissADME_Analysis_Report.docx`. Methods Section
2.10 states this, and `deprecated/README.md` sets out how the molecular weights
establish which of the report's rows are usable.

## Running the deposited scripts

Two things are worth knowing before running anything here.

**Absolute paths.** Twenty-two of the fifty-five scripts carry absolute paths
from the machine each analysis was run on, as `C:/Temp/...`, `/home/...` or
`/mnt/...`; twelve of those are in `deprecated/` and are not the basis of any
reported result. The paths have deliberately not been rewritten: these files
are the record of what was executed, and editing them now would mean the
repository no longer holds the scripts that produced the deposited outputs.
Adjust the paths to your own layout before running a script, or read the
deposited output it produced, which is stored beside it.

**Checking the claims in this documentation.** `verification/check_deposit_claims.py`
resolves every path the README files here name, recomputes the two script counts
stated above from the files themselves, and lists the sentences in this
documentation that make a universal or exclusive claim about the deposit. It needs
only the Python standard library and takes no arguments; passing the manuscript,
Supplementary Information and response-letter `.docx` files extends the first and
third checks over the paper as well. It is a checking tool rather than one of the
fifty-five analysis scripts counted above, and it excludes itself from those counts,
so a reader who counts `.py` files in this repository will find fifty-six.

**Comment language.** Eighteen scripts carry comments in Korean, most of them
in `revision_outputs/scripts/`. The code, the variable names and the printed
output are in English throughout, and every value these scripts produce is
also recorded in the deposited CSV outputs and in `verification/`, so a reader
who does not read Korean can still follow what each script computes and check
its results.

## Correction records

Four deposited files were corrected after the first revision was prepared. Each
correction is documented next to the file it concerns, and in every case the
manuscript value was already correct:

| Note | Concerns |
|---|---|
| `reproducibility/README_cell_level_correction.md` | 18 rows whose annotation had been corrupted to `#NAME?` by spreadsheet round-tripping, 173 rows that had lost their patient identifier, and a new `pooled_plate` column that makes the Methods 2.2 exclusion reproducible from the file alone |
| `revision_outputs/gse198315/README_J1_correction.md` | One cell missing from the Fig. 4A verification export |
| `revision_outputs/structure/README_RSA_correction.md` | The sphere-sampling density used for the per-residue solvent accessibility |
| `verification/README_verification_ledger_scope.md` | The scope of the verification ledger and audit report, both of which use the section and table numbering of the submitted version |

## License

MIT. See [`LICENSE`](LICENSE).
