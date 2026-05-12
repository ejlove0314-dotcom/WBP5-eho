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
