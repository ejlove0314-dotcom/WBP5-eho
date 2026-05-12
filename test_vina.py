import os
import subprocess

print("Step 1: Check files...")

pdb = r"C:\Users\user\Desktop\WBP5 Small molecule test\WBP5_AlphaFold.pdb"
print(f"PDB exists: {os.path.exists(pdb)}")

receptor = r"C:\Users\user\.mgltools\PyRx\Macromolecules\WBP5_AlphaFold_model1\WBP5_AlphaFold_model1.pdbqt"
print(f"Receptor exists: {os.path.exists(receptor)}")

lig_dir = r"C:\Users\user\.mgltools\PyRx\Ligands"
for f in os.listdir(lig_dir):
    if "Candidate" in f or "Pazopanib" in f:
        print(f"  Ligand: {f}")

print("\nStep 2: Test vina...")
result = subprocess.run(
    [r"C:\Program Files (x86)\PyRx\vina.exe", "--help"],
    capture_output=True, text=True
)
print(f"Vina stdout: {result.stdout[:200]}")
print(f"Vina stderr: {result.stderr[:200]}")

print("\nStep 3: Calculate pocket center...")
x, y, z = [], [], []
with open(pdb, 'r') as f:
    for line in f:
        if line.startswith("ATOM"):
            res_num = int(line[22:26].strip())
            if 64 <= res_num <= 74:
                x.append(float(line[30:38]))
                y.append(float(line[38:46]))
                z.append(float(line[46:54]))

print(f"Pocket atoms: {len(x)}")
if len(x) > 0:
    cx, cy, cz = sum(x)/len(x), sum(y)/len(y), sum(z)/len(z)
    print(f"Center: {cx:.3f}, {cy:.3f}, {cz:.3f}")

print("\nDone!")
