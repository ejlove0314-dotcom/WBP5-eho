"""Relative solvent accessibility (RSA) of the AlphaFold WBP5/TCEAL9 model.

SASA is computed with the Shrake-Rupley algorithm as implemented in Biopython,
using the default 100-point sphere sampling and a 1.40 A probe.  RSA = SASA /
Tien et al. (2013) theoretical maximum accessible surface area for the residue
type.  The 100-point setting is the one reported in Supplementary Table S18.

Usage:  python compute_rsa_shrake_rupley.py MODEL.pdb OUT.csv [n_points]

An earlier version of this script used n_points=960.  The finer sampling shifts
individual RSA values by up to about 0.03 but does not alter the conclusion
drawn in Table S18: over the P0 region (residues 38-58) it gives a range of
0.290-0.624 with a median of 0.521, against 0.266-0.593 and a median of 0.516
at 100 points, so every P0 residue is solvent-exposed under either setting.
Pass 960 as the third argument to reproduce the earlier values.
"""
import csv, sys
from Bio.PDB import PDBParser, ShrakeRupley

# Tien et al. 2013 (PLoS ONE 8:e80635), theoretical Gly-X-Gly maxima
TIEN = {"ALA":129,"ARG":274,"ASN":195,"ASP":193,"CYS":167,"GLN":225,"GLU":223,
        "GLY":104,"HIS":224,"ILE":197,"LEU":201,"LYS":236,"MET":224,"PHE":240,
        "PRO":159,"SER":155,"THR":172,"TRP":285,"TYR":263,"VAL":174}
THREE2ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
             "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
             "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
POCKET = set(range(38, 59))          # P0 pocket span, residues 38-58

pdb, out = sys.argv[1], sys.argv[2]
n_points = int(sys.argv[3]) if len(sys.argv) > 3 else 100
s = PDBParser(QUIET=True).get_structure("wbp5", pdb)
ShrakeRupley(n_points=n_points).compute(s[0], level="R")

rows = []
for res in s[0].get_residues():
    rn = res.get_resname()
    if rn not in TIEN:
        continue
    num = res.get_id()[1]
    sasa = res.sasa
    plddt = sum(a.get_bfactor() for a in res) / len(res)
    rows.append({"residue_number": num,
                 "residue": rn,
                 "one_letter": THREE2ONE[rn],
                 "sasa_A2": round(sasa, 2),
                 "max_asa_A2_Tien2013": TIEN[rn],
                 "rsa": round(sasa / TIEN[rn], 4),
                 "mean_pLDDT": round(plddt, 2),
                 "in_P0_pocket_38_58": int(num in POCKET)})

with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

pk = sorted(r["rsa"] for r in rows if r["in_P0_pocket_38_58"])
mid = pk[len(pk)//2] if len(pk) % 2 else (pk[len(pk)//2-1] + pk[len(pk)//2]) / 2
print(f"residues: {len(rows)}  (span {rows[0]['residue_number']}-{rows[-1]['residue_number']})")
print(f"sphere sampling: {n_points} points")
print(f"P0 region (38-58): n={len(pk)}  RSA {min(pk):.4f}-{max(pk):.4f}  median {mid:.4f}")
