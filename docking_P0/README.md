# Docking at pocket P0

AutoDock Vina docking of ten ligands (five de novo candidates, pazopanib, and four
approved comparators) at pocket P0 of the WBP5 AlphaFold model.

- Grid center: -23.835, -0.451, 23.910
- Box: 25 x 25 x 25 A
- Exhaustiveness: 8
- Five replicates per ligand (seeds 42, 1042, 2042, 3042, 4042)

`candidate1_pyridine_P0.csv` supersedes the Candidate_1 row of
`docking_replicates_P0_summary.csv`; the SMILES was corrected to the pyridin-4-yl
variant (MW 347.44) and the ligand re-docked. The merged, corrected summary is in
`reproducibility/docking_replicates_P0_summary_v2.csv`.
