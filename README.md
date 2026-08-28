# WBP5-eho

Code and data supporting the manuscript "Single-cell-to-structure nomination of
WBP5 (TCEAL9) as a p-EMT-associated target in HNSCC."

## Repository layout

| Path | Contents |
|------|----------|
| `verification/` | Numerical verification of every value reported in the manuscript. **Start here.** |
| `reproducibility/` | Scripts and verified outputs for manuscript Tables 1-4. |
| `tables/` | Source tables underlying Supplementary Tables S3 and S4. |
| `timer3/` | Raw TIMER3 exports underlying Table S4 and Fig. S3. |
| `figures/` | Scripts that regenerate Figs. S1 and S3. |
| `gse198315_recalc/` | Depth-corrected recomputation of the GSE198315 correlations. |
| `docking_P0/` | AutoDock Vina docking at pocket P0 - the protocol reported in the manuscript. |
| `docking_P0/mtiopenscreen/` | MTiOpenScreen Drugs-lib screen at the same pocket (job V06821248773098). |
| `phase1_singlecell/` | GSE198315 single-cell analysis pipeline. |
| `pockets/`, `residues/` | DoGSiteScorer pocket detection: run parameters, descriptor table, grid maps and pocket-lining residues. |
| `deprecated/` | Superseded scripts, ligand files and reports. **Not used for any reported result.** |

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
the docking margins in Sections 3.5 and 3.7, and the depth-correction deltas in
Table S3:

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
the pocket geometry reported in Table 3 is taken: for P0, volume 229.62 A^3, surface
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
| `reproducibility/swissadme_final_20260820.csv` | Table 4; Table S5 Sections B-C |
| `reproducibility/docking_replicates_P0_summary_v2.csv` | Table 3; Table 4 |
| `reproducibility/ProTox-3.0 - ... candidate 1-5, pazopanib.xlsx` | Table S5 Sections F1-F2 |
| `phase1_singlecell/gse198315_wbp5_marker_correlations.tsv` | Fig. 5B; Fig. S1 |
| `phase1_singlecell/make_gse198315_wbp5_marker_pipeline.py` | Script generating the above TSV |
| `tables/TableS3_gse198315_depth_correction.csv` | Table S3 |
| `tables/TableS4_timer3_hnsc_correlations.csv` | Table S4 |
| `gse198315_recalc/headline.tsv` | Section 3.2; Table S3 |
| `gse198315_recalc/correlations_full.tsv` | Fig. S1A |

### Note on superseded files

Files in `deprecated/` correspond to an earlier analysis run and are retained
for provenance only. They are not the basis of any value reported in the
manuscript. In particular, the Candidate_1 structure was corrected after the
initial SwissADME submission (see Table S5, Section A footnote); all values in
the manuscript derive from `swissadme_final_20260820.csv`.
Molecular weights in this table were recomputed with RDKit 2026.03.5 and may
differ from the SwissADME values reported in the manuscript by 0.01 Da owing
to differing atomic mass tables.

## License

MIT. See [`LICENSE`](LICENSE).
