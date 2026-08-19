# WBP5-eho

Code and data supporting the manuscript "Single-cell-to-structure analysis identifies
WBP5 (TCEAL9) as a candidate p-EMT target in HNSCC."

## Repository layout

| Path | Contents |
|------|----------|
| `reproducibility/` | Scripts and verified outputs for manuscript Tables 1-4. **Start here.** |
| `docking_P0/` | AutoDock Vina docking at pocket P0 - the protocol reported in the manuscript. |
| `docking_P0/mtiopenscreen/` | MTiOpenScreen Drugs-lib screen at the same pocket (job V06821248773098). |
| `phase1_singlecell/` | GSE198315 single-cell analysis pipeline. |
| `pockets/`, `residues/` | DoGSiteScorer pocket detection: run parameters, descriptor table, grid maps and pocket-lining residues. |
| `deprecated/` | Superseded scripts, ligand files and reports. **Not used for any reported result.** |

Root-level files are the AlphaFold model (`AF-Q9UHQ7-F1-model_v6.pdb`, and the
identical working copy `WBP5_AlphaFold.pdb` used as the docking receptor), the
REINVENT4 configuration files (`rl_wbp5.toml`, `rl_wbp5_fixed.toml`,
`reinvent_wbp5_sampling.toml`), the corrected ligand set
(`all_candidates_corrected.smi`), and the individual structure files for
Candidate_2 through Candidate_5.

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
