from rdkit import Chem
from rdkit.Chem import Draw
import os

output_dir = r"C:\Users\user\Desktop\WBP5 Small molecule test\images"
os.makedirs(output_dir, exist_ok=True)

mols_data = [
    ("Candidate_2 (-6.0)", "O=C(CCc1ccccc1)NNC(=O)C1CCN(c2ncccn2)CC1"),
    ("Candidate_3 (-5.6)", "CC1=CN(Cc2ccc(S(N)(=O)=O)cc2)C(=O)c2cccnc2C1"),
    ("Candidate_5 (-4.9)", "CN1CCN(S(=O)(=O)c2ccc(NC(=O)c3ccco3)cc2)CC1"),
    ("Candidate_1 (-4.9)", "CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccncc2)cc1"),
    ("Candidate_4 (-5.6)", "Cc1ccc(OCc2cc(C(=O)NCc3c(C)noc3C)no2)cc1C"),
    ("Pazopanib (-6.1)", "CC1=C(C=C(C=C1)NC2=NC=CC(=N2)N3CCC(CC3)NC(=O)C4=CC=CS4)S(=O)(=O)N"),
]

mols = []
legends = []
for name, smi in mols_data:
    print(f"Processing {name}...")
    mol = Chem.MolFromSmiles(smi)
    if mol:
        mols.append(mol)
        legends.append(name)

print("Creating grid image...")
img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(400,350), legends=legends)
img.save(os.path.join(output_dir, "all_candidates_grid.png"))

print("Creating top picks vs pazopanib...")
img2 = Draw.MolsToGridImage([mols[0],mols[1],mols[5]], molsPerRow=3, subImgSize=(450,400),
    legends=[legends[0],legends[1],legends[5]])
img2.save(os.path.join(output_dir, "top_vs_pazopanib.png"))

print("Done!")
