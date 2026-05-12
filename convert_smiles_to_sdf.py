"""Convert SMILES to SDF files for PyRx docking"""
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os

# Output directory
output_dir = r"C:\Users\user\Desktop\WBP5 Small molecule test"

# Candidates with SMILES (SwissADME ranking order)
candidates = {
    "Candidate_3": "CC1=CN(Cc2ccc(S(N)(=O)=O)cc2)C(=O)c2cccnc2C1",
    "Candidate_5": "CN1CCN(S(=O)(=O)c2ccc(NC(=O)c3ccco3)cc2)CC1",
    "Candidate_2": "O=C(CCc1ccccc1)NNC(=O)C1CCN(c2ncccn2)CC1",
    "Candidate_1": "CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccnc2)cc1",
    "Candidate_4": "Cc1ccc(OCc2cc(C(=O)NCc3c(C)noc3C)no2)cc1C",
    "Pazopanib_Ref": "CC1=C(C=C(C=C1)NC2=NC=CC(=N2)N3CCC(CC3)NC(=O)C4=CC=CS4)S(=O)(=O)N",
}

# Create individual SDF files
for name, smiles in candidates.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"ERROR: Failed to parse {name}: {smiles}")
        continue

    # Add hydrogens and generate 3D coordinates
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result == -1:
        print(f"WARNING: 3D embedding failed for {name}, trying random coords...")
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), useRandomCoords=True)

    # Optimize geometry with MMFF94
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except:
        print(f"WARNING: MMFF optimization failed for {name}, trying UFF...")
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except:
            print(f"WARNING: UFF also failed for {name}")

    # Set molecule name
    mol.SetProp("_Name", name)

    # Write individual SDF
    sdf_path = os.path.join(output_dir, f"{name}.sdf")
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()
    print(f"Created: {sdf_path}")

# Also create a combined SDF with all molecules
combined_path = os.path.join(output_dir, "all_candidates.sdf")
writer = Chem.SDWriter(combined_path)
for name, smiles in candidates.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        continue
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result == -1:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), useRandomCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except:
            pass
    mol.SetProp("_Name", name)
    writer.write(mol)
writer.close()
print(f"\nCreated combined: {combined_path}")
print("\nDone! All SDF files ready for PyRx.")
