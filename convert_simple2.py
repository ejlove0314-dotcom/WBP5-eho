from rdkit import Chem
from rdkit.Chem import AllChem
import os

output_dir = r"C:\Users\user\Desktop\WBP5 Small molecule test"

# Remaining candidates that failed or weren't processed
candidates = [
    ("Candidate_1", "CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccnc2)cc1"),
    ("Candidate_4", "Cc1ccc(OCc2cc(C(=O)NCc3c(C)noc3C)no2)cc1C"),
    ("Pazopanib_Ref", "CC1=C(C=C(C=C1)NC2=NC=CC(=N2)N3CCC(CC3)NC(=O)C4=CC=CS4)S(=O)(=O)N"),
]

for name, smi in candidates:
    print(f"Processing {name}...")
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        # Try sanitize=False then manual sanitize
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is not None:
            try:
                Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                print(f"  (used partial sanitize for {name})")
            except:
                print(f"  ERROR: Cannot sanitize {name}")
                continue
        else:
            print(f"  ERROR: Cannot parse {name}")
            continue
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    mol.SetProp("_Name", name)
    path = os.path.join(output_dir, f"{name}.sdf")
    w = Chem.SDWriter(path)
    w.write(mol)
    w.close()
    print(f"  -> {name}.sdf done")

print("\nAll done!")
