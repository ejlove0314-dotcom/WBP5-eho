import os
import subprocess

VINA = r"C:\Program Files (x86)\PyRx\vina.exe"
RECEPTOR = r"C:\Users\user\.mgltools\PyRx\Macromolecules\WBP5_AlphaFold_model1\WBP5_AlphaFold_model1.pdbqt"
LIG_DIR = r"C:\Users\user\.mgltools\PyRx\Ligands"
OUT_DIR = r"C:\Users\user\Desktop\WBP5 Small molecule test\docking_results"
os.makedirs(OUT_DIR, exist_ok=True)

# P_0 pocket center (residues 64-74)
cx, cy, cz = 2.402, -3.010, 1.561

ligands = [
    "Candidate_1.pdbqt",
    "Candidate_2.pdbqt",
    "Candidate_3.pdbqt",
    "Candidate_4.pdbqt",
    "Candidate_5.pdbqt",
    "Pazopanib_Ref.pdbqt",
]

results = []

for lig in ligands:
    lig_path = os.path.join(LIG_DIR, lig)
    name = lig.replace(".pdbqt", "")
    out_path = os.path.join(OUT_DIR, f"{name}_out.pdbqt")
    log_path = os.path.join(OUT_DIR, f"{name}_log.txt")

    print(f"\nDocking {name}...")

    cmd = [
        VINA,
        "--receptor", RECEPTOR,
        "--ligand", lig_path,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", "25",
        "--size_y", "25",
        "--size_z", "25",
        "--exhaustiveness", "8",
        "--num_modes", "9",
        "--out", out_path,
        "--log", log_path,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            log = f.read()
        for line in log.split('\n'):
            if line.strip().startswith('1'):
                parts = line.split()
                try:
                    affinity = float(parts[1])
                    results.append((name, affinity))
                    print(f"  Best affinity: {affinity} kcal/mol")
                except:
                    pass
                break
    else:
        print(f"  ERROR: {result.stderr[:300]}")

print("\n" + "=" * 50)
print("DOCKING RESULTS SUMMARY")
print("=" * 50)
print(f"{'Molecule':<20} {'Affinity (kcal/mol)'}")
print("-" * 40)
results.sort(key=lambda x: x[1])
for name, aff in results:
    tag = " <-- REF" if "Pazopanib" in name else ""
    print(f"{name:<20} {aff}{tag}")
print("\nMore negative = stronger binding")
