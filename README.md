# WBP5-eho

Code and data supporting the manuscript "Single-cell-to-structure analysis identifies
WBP5 (TCEAL9) as a candidate p-EMT target in HNSCC."

## Repository layout

| Path | Contents |
|------|----------|
| `reproducibility/` | Scripts and verified outputs for manuscript Tables 1-4. **Start here.** |
| `docking_P0/` | AutoDock Vina docking at pocket P0 — the protocol reported in the manuscript. |
| `phase1_singlecell/` | GSE198315 single-cell analysis pipeline. |
| `pockets/`, `residues/` | DoGSiteScorer pocket detection output for the AlphaFold model. |

All docking values reported in the manuscript derive from pocket P0
(grid center -23.835, -0.451, 23.910; box 25 x 25 x 25 A).

Source data are available from GEO under accessions GSE103322 and GSE198315.

---
# WBP5/HNSCC Figures — Arial-unified

## Font policy
All figures use a sans-serif family with Arial-metric widths. On the
Linux generation environment (no real Arial installed) this resolves to
`Liberation Sans`, which is byte-equivalent in width to Arial. Running
the scripts on the Windows PC where Arial is installed will produce
identical layouts rendered with real Arial, because `_style.py` prefers
`Arial` when present.

## Files

| File | Source data | Script |
|------|-------------|--------|
| `fig_workflow.png` | Manuscript outline (text-based) | `_scripts/fig_workflow.py` |
| `fig_structure_pocket_view1.png` / `view2.png` | `WBP5_AlphaFold.pdb` | `_scripts/fig_structure.py` |
| `fig_heatmap_3panel.png` | Transcribed from original heatmap figure | `_scripts/fig_heatmap.py` |
| `fig_regional_prevalence.png` | `GSE198315 region_compartment breakdown heatmap ... .xlsx` | `_scripts/fig_prevalence.py` |
| `fig_stacked_composition.png` | `GSE198315_*distribution_summary.xlsx` + double-high breakdown | `_scripts/fig_stacked.py` |
| `fig_external_validation.png` | Transcribed from the original external-validation figure | `_scripts/fig_extvalid.py` |
| `fig_compartment_4panel.png` | Transcribed from the original compartment/triage figure | `_scripts/fig_compartment.py` |
| `fig_umap_gse198315_4panel.png` | `GSE198315/out/umap_coords.csv` + `obs.csv` | `_scripts/fig_umap_gse198315.py` |
| `fig_umap_overlay_4panel.png` | same + top-10% thresholds on expr_WBP5/expr_EGFR | `_scripts/fig_umap_overlay.py` |
| `fig_umap_gse103322_4panel.png` | `GSE103322/out/umap_coords.csv` + `obs.csv` + WBP5 row from matrix | `_scripts/fig_umap_gse103322.py` |

## Not included (require external tools)

| Missing | Reason | How to get it |
|---------|--------|---------------|
| Candidate molecule grids (`all_candidates_grid`, `top_vs_pazopanib`) | RDKit not available in the sandbox | Run `_scripts/fig_molecules_ARIAL.py` on the Windows PC |
| PyMOL ribbon screenshot (`wbp5alphafoldpdb*.png`) | PyMOL/ChimeraX not available | Open `WBP5_AlphaFold.pdb` in PyMOL, set `set cartoon_fancy_helices,1; set_cartoon_transparency,0`, export with Arial labels |

## Regenerating everything

```bash
cd _scripts
for f in fig_*.py; do python3 "$f"; done
```

## Font override

If Arial is installed system-wide, `_style.py` picks it up automatically.
Otherwise the rendered font is Liberation Sans, which matches Arial's
kerning and advance widths exactly.
## Note on the ligand SMILES files

Two entries in `all_candidates.smi` (deposited 13 April 2026) required correction
before the analyses reported in the manuscript. The original file is retained
unmodified for provenance; the corrected strings are provided separately as
`all_candidates_corrected.smi`. All docking, ADMET and toxicity results reported
in the manuscript were obtained with the corrected structures.

### Candidate_1

The deposited string

```
CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccnc2)cc1 Candidate_1
```

encodes a five-membered aromatic ring bearing a single nitrogen without an
explicit hydrogen, and therefore cannot be kekulized by RDKit or Open Babel.
The molecular weight (347 Da) and topological polar surface area (88.16 A^2)
recorded for this molecule in the REINVENT4 run correspond to the
pyridine-containing analogue (C17H21N3O3S; 347.44 Da; TPSA 88.16 A^2) rather
than to the pyrrole analogue (C16H21N3O3S; 335.43 Da; TPSA 91.06 A^2). All
analyses were therefore performed with

```
CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccncc2)cc1     pyridin-4-yl
InChIKey JGUIGLSIQCDQRA-UHFFFAOYSA-N
```

The position of the ring nitrogen (2-, 3- or 4-pyridyl) could not be
independently confirmed from the archived output, because the three positional
isomers share the same molecular weight, topological polar surface area and
aromatic ring count. The 4-pyridyl assignment follows the ligand definition used
in the original docking script.

### Pazopanib_Ref

The deposited string encoded a thiophene-carboxamide/piperidine analogue
(C21H24N6O3S2; 472.60 Da), not pazopanib. The reference compound was corrected to
the free base of pazopanib (PubChem CID 10113978):

```
CN(c1ccc2c(c1)nn(c2C)C)c1ccnc(n1)Nc1ccc(c(c1)S(=O)(=O)N)C
C21H23N7O2S; 437.53 Da; InChIKey CUIHSIWYWATEQL-UHFFFAOYSA-N
```

### Verification

Every entry in `all_candidates_corrected.smi` parses under RDKit 2026.03.5 and
matches the molecular formula and molecular weight reported in the manuscript:

| ID            | Formula      | MW (Da) | InChIKey                    |
|---------------|--------------|---------|-----------------------------|
| Candidate_1   | C17H21N3O3S  | 347.44  | JGUIGLSIQCDQRA-UHFFFAOYSA-N |
| Candidate_2   | C19H23N5O2   | 353.43  | CWFOXUOVUNELQQ-UHFFFAOYSA-N |
| Candidate_3   | C17H17N3O3S  | 343.41  | DSVBRFOHTGCUGB-UHFFFAOYSA-N |
| Candidate_4   | C19H21N3O4   | 355.39  | MMBJDRLBTDFZJL-UHFFFAOYSA-N |
| Candidate_5   | C16H19N3O4S  | 349.41  | UQFXIMOKNSXXKZ-UHFFFAOYSA-N |
| Pazopanib_Ref | C21H23N7O2S  | 437.53  | CUIHSIWYWATEQL-UHFFFAOYSA-N |

