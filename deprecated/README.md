# Deprecated files

With one exception, nothing in this directory was used to produce any result
reported in the manuscript; these files are retained for provenance only.

The exception is `SwissADME_Analysis_Report.docx`. Section 4.2 of that report is the
source of the CYP450 predictions tabulated in Supplementary Table S5, Section D, for
Candidate_2 to Candidate_5. Every other section of the same report is superseded, and
the values it records for Candidate_1 and for pazopanib are not used anywhere. The
reasons are set out under the warning below.

## Contents

| Path | Why it is deprecated |
|------|----------------------|
| `run_docking.py` | Uses grid center `2.402, -3.010, 1.561`, which is pocket **P_3**, not the reported pocket P0. Superseded by `docking_P0/redock_P0_correct.py`. |
| `run_vina_docking.py`, `test_vina.py` | Earlier drafts of the docking driver. Superseded by `docking_P0/redock_P0_correct.py`. |
| `all_candidates.smi`, `Candidate_1.smi/.sdf/.pdb`, `Pazopanib_Ref.smi/.sdf/.pdb` | Contain the two incorrect structures described below. Superseded by `all_candidates_corrected.smi`. |
| `convert_simple.py`, `convert_simple2.py`, `convert_simple3.py`, `convert_smiles_to_sdf.py`, `convert_to_pdb.py` | Ad hoc format-conversion helpers used with the deprecated SMILES files. |
| `gen_img_simple.py`, `generate_images.py`, `generate_protein_image.py` | Draft figure-rendering scripts. Manuscript figures were produced separately. |
| `PyRx_docking_SMILES.txt`, `SwissADME_input_SMILES.txt` | Input lists for the earlier docking and ADMET submissions. Both carry the two incorrect structures described below. |
| `SwissADME_Analysis_Report.docx` | ADMET report dated 13 April 2026, generated from those inputs. **Section 4.2 is the source of Supplementary Table S5, Section D**; every other section is superseded. See the warning below. |
| `fda.sdf` | Partial local import of an approved-compound library (103 entries with ZINC identifiers). The local import was abandoned; the library-scale screen reported in the manuscript was run on the MTiOpenScreen server instead (Methods 2.7). |
| `boltz2/` | Exploratory Boltz-2 co-folding run; not reported (see root README). |

## Warning: `SwissADME_Analysis_Report.docx`

This report predates the structure corrections and **contradicts the manuscript on
three points**, listed below. Apart from Section 4.2, it is retained only to document
what was originally computed.

| Statement in the report | Status |
|-------------------------|--------|
| Candidate_1 has MW 335 | Superseded. That is the pyrrole analogue; the compound analyzed is the pyridin-4-yl variant, MW 347.44. |
| Pazopanib has MW 472.58 and TPSA 166.93, and fails the Veber, Egan and Muegge filters | Incorrect. Those values belong to the wrong structure. Pazopanib is C21H23N7O2S, MW 437.53. The corresponding claim was withdrawn from the manuscript. |
| Candidates #3 and #5 are the top picks | Superseded. After redocking at P0 and re-running ADMET on the corrected structures, Candidate_4 is the affinity lead and Candidate_2 the property lead. |

The ADMET and toxicity values reported in the manuscript were obtained after these
corrections and are tabulated in Supplementary Tables S5 and S10.

## Why Section 4.2 of this report is still the source of Table S5, Section D

The SwissADME CYP450 module returned no output when the corrected structures were
re-submitted in August 2026: `reproducibility/swissadme_final_20260820.csv` carries
`n/d` in all five CYP columns for all six compounds. The April 2026 session is
therefore the only run that produced CYP values, and Section 4.2 of this report is the
record of it.

Four of its six rows were computed on the structures the manuscript analyzes. Because
both runs are SwissADME, their molecular weights can be compared directly against the
August file: this report records 353.42, 343.40, 355.39 and 349.40 for Candidate_2,
Candidate_3, Candidate_4 and Candidate_5, which are the August values to the
hundredth, whereas its 335.42 for Candidate_1 and 472.58 for pazopanib belong to the
two superseded structures described above (347.43 and 437.52 in the August file).
Supplementary Table S5, Section D reports the four and gives `n.d.` for the other two,
which is why those two rows are blank there.

Molecular weights in the root `README.md` are RDKit values and differ from the
SwissADME values above by at most 0.01 Da; see the note there.

## Pocket assignment

Grid centers for the two pockets differ by 32.75 A:

| Pocket | Grid center | DoGSiteScorer Drug Score | Used in manuscript |
|--------|-------------|--------------------------|--------------------|
| P0 | -23.835, -0.451, 23.910 | 0.837 | yes |
| P_3 | 2.402, -3.010, 1.561 | 0.281 | no |

Every docking value in the manuscript was recomputed at P0. Any output produced by
`run_docking.py` refers to P_3 and should not be compared with the reported energies.

## Ligand structure corrections

### Candidate_1

The deposited string

```
CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccnc2)cc1 Candidate_1
```

encodes a five-membered aromatic ring bearing a single nitrogen without an explicit
hydrogen, and therefore cannot be kekulized by RDKit or Open Babel. The molecular weight
(347 Da) and topological polar surface area (88.16 A^2) recorded for this molecule in the
REINVENT4 run correspond to the pyridine-containing analogue (C17H21N3O3S; 347.44 Da;
TPSA 88.16 A^2) rather than to the pyrrole analogue (C16H21N3O3S; 335.43 Da;
TPSA 91.06 A^2). All reported analyses used

```
CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccncc2)cc1     pyridin-4-yl
InChIKey JGUIGLSIQCDQRA-UHFFFAOYSA-N
```

The position of the ring nitrogen (2-, 3- or 4-pyridyl) could not be independently
confirmed from the archived output, because the three positional isomers share the same
molecular weight, topological polar surface area and aromatic ring count. The 4-pyridyl
assignment follows the ligand definition used in the original docking script.

### Pazopanib_Ref

The deposited string

```
CC1=C(C=C(C=C1)NC2=NC=CC(=N2)N3CCC(CC3)NC(=O)C4=CC=CS4)S(=O)(=O)N
```

encodes a thiophene-carboxamide/piperidine analogue (C21H24N6O3S2; 472.60 Da), not
pazopanib. The reference compound was corrected to the free base of pazopanib
(PubChem CID 10113978):

```
Cc1ccc(Nc2nccc(N(C)c3ccc4c(C)n(C)nc4c3)n2)cc1S(N)(=O)=O
C21H23N7O2S; 437.53 Da; InChIKey CUIHSIWYWATEQL-UHFFFAOYSA-N
```

The corrected structure is the one present in
`docking_P0/Pazopanib_Ref_seed42_out.pdbqt` (see its `REMARK SMILES` record), so the
reported energy of -6.86 +/- 0.05 kcal/mol refers to pazopanib itself.
