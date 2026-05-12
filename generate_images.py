"""Generate molecular structure images using RDKit"""
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import os

output_dir = r"C:\Users\user\Desktop\WBP5 Small molecule test\images"
os.makedirs(output_dir, exist_ok=True)

# All candidates + Pazopanib
molecules = {
    "Candidate_2\n(TOP PICK, -6.0 kcal/mol)": "O=C(CCc1ccccc1)NNC(=O)C1CCN(c2ncccn2)CC1",
    "Candidate_3\n(TOP PICK, -5.6 kcal/mol)": "CC1=CN(Cc2ccc(S(N)(=O)=O)cc2)C(=O)c2cccnc2C1",
    "Candidate_5\n(Backup, -4.9 kcal/mol)": "CN1CCN(S(=O)(=O)c2ccc(NC(=O)c3ccco3)cc2)CC1",
    "Candidate_1\n(Backup, -4.9 kcal/mol)": "CC(C)NS(=O)(=O)c1ccc(CCC(=O)Nc2ccncc2)cc1",
    "Candidate_4\n(Deprioritized, -5.6 kcal/mol)": "Cc1ccc(OCc2cc(C(=O)NCc3c(C)noc3C)no2)cc1C",
    "Pazopanib\n(Reference, -6.1 kcal/mol)": "CC1=C(C=C(C=C1)NC2=NC=CC(=N2)N3CCC(CC3)NC(=O)C4=CC=CS4)S(=O)(=O)N",
}

# Generate individual images
mols = []
legends = []
for name, smi in molecules.items():
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"ERROR: {name}")
        continue
    AllChem.Compute2DCoords(mol)
    mols.append(mol)
    legends.append(name)

    # Individual image
    safe_name = name.split("\n")[0].replace(" ", "_")
    img_path = os.path.join(output_dir, f"{safe_name}.png")
    Draw.MolToFile(mol, img_path, size=(400, 350), legend=name.replace("\n", " "))
    print(f"Created: {img_path}")

# Grid image - all 6 molecules
print("\nCreating grid image...")
grid_path = os.path.join(output_dir, "all_candidates_grid.png")
img = Draw.MolsToGridImage(
    mols,
    molsPerRow=3,
    subImgSize=(450, 400),
    legends=legends,
    useSVG=False
)
img.save(grid_path)
print(f"Created: {grid_path}")

# Top 2 picks comparison
print("\nCreating top picks comparison...")
top_mols = mols[:2]  # Candidate_2 and Candidate_3
top_legends = legends[:2]
top_path = os.path.join(output_dir, "top_picks.png")
img2 = Draw.MolsToGridImage(
    top_mols,
    molsPerRow=2,
    subImgSize=(500, 450),
    legends=top_legends,
    useSVG=False
)
img2.save(top_path)
print(f"Created: {top_path}")

# Top 2 vs Pazopanib comparison
print("\nCreating comparison with Pazopanib...")
compare_mols = [mols[0], mols[1], mols[5]]  # Cand_2, Cand_3, Pazopanib
compare_legends = [legends[0], legends[1], legends[5]]
compare_path = os.path.join(output_dir, "top_vs_pazopanib.png")
img3 = Draw.MolsToGridImage(
    compare_mols,
    molsPerRow=3,
    subImgSize=(450, 400),
    legends=compare_legends,
    useSVG=False
)
img3.save(compare_path)
print(f"Created: {compare_path}")

print(f"\nAll images saved to: {output_dir}")
print("Done!")
