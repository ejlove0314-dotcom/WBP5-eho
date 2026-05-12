from rdkit import Chem
from rdkit.Chem import AllChem
import os

output_dir = r"C:\Users\user\Desktop\WBP5 Small molecule test"

candidates = [
    ("Candidate_3", "CC1=CN(Cc2ccc(S(N)(=O)=O)cc2)C(=O)c2cccnc2C1"),
    ("Candidate_5", "CN1CCN(S(=O)(=O)c2ccc(NC(=O)c3ccco3)cc2)CC1"),
    ("Candidate_2", "O=C(CCc1ccccc1)NNC(=O)C1CCN(c2ncccn2)CC1"),
    ("Candidate_1", "CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccnc2)cc1"),
    ("Candidate_4", "Cc1ccc(OCc2cc(C(=O)NCc3c(C)noc3C)no2)cc1C"),
    ("Pazopanib_Ref", "CC1=C(C=C(C=C1)NC2=NC=CC(=N2)N3CCC(CC3)NC(=O)C4=CC=CS4)S(=O)(=O)N"),
]

for name, smi in candidates:
    print(f"Processing {name}...")
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol.SetProp("_Name", name)
    path = os.path.join(output_dir, f"{name}.sdf")
    w = Chem.SDWriter(path)
    w.write(mol)
    w.close()
    print(f"  -> {name}.sdf done")

print("\nAll done!")
