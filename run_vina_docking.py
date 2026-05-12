"""
AutoDock Vina docking script for WBP5 P_0 pocket
Runs docking for all 6 candidates against WBP5
"""
import os
import subprocess

# Paths
VINA = r'"C:\Program Files (x86)\PyRx\vina.exe"'
RECEPTOR = r"C:\Users\user\.mgltools\PyRx\Macromolecules\WBP5_AlphaFold_model1\WBP5_AlphaFold_model1.pdbqt"
LIGAND_DIR = r"C:\Users\user\.mgltools\PyRx\Ligands"
OUTPUT_DIR = r"C:\Users\user\Desktop\WBP5 Small molecule test\docking_results"

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Step 1: Calculate grid box center from P_0 pocket residues (64-74)
print("=" * 60)
print("Calculating grid box center from WBP5 P_0 pocket (res 64-74)")
print("=" * 60)

pdb_file = r"C:\Users\user\Desktop\WBP5 Small molecule test\WBP5_AlphaFold.pdb"
x_coords, y_coords, z_coords = [], [], []

with open(pdb_file, 'r') as f:
    for line in f:
        if line.startswith("ATOM"):
            res_num = int(line[22:26].strip())
            if 64 <= res_num <= 74:
                x_coords.append(float(line[30:38]))
                y_coords.append(float(line[38:46]))
                z_coords.append(float(line[46:54]))

center_x = sum(x_coords) / len(x_coords)
center_y = sum(y_coords) / len(y_coords)
center_z = sum(z_coords) / len(z_coords)

print(f"P_0 pocket center: X={center_x:.3f}, Y={center_y:.3f}, Z={center_z:.3f}")
print(f"Number of atoms in pocket: {len(x_coords)}")

# Grid box size (25 Angstrom cube)
size_x, size_y, size_z = 25, 25, 25

# Step 2: Define ligands to dock
ligands = [
    "Candidate_3.pdbqt",
    "Candidate_5.pdbqt",
    "Candidate_2.pdbqt",
    "Candidate_1.pdbqt",
    "Candidate_4.pdbqt",
    "Pazopanib_Ref.pdbqt",
]

# Step 3: Run Vina for each ligand
print("\n" + "=" * 60)
print("Running AutoDock Vina Docking")
print("=" * 60)

results = []

for ligand_file in ligands:
    ligand_path = os.path.join(LIGAND_DIR, ligand_file)
    if not os.path.exists(ligand_path):
        print(f"\nWARNING: {ligand_file} not found, skipping...")
        continue

    name = ligand_file.replace(".pdbqt", "")
    out_path = os.path.join(OUTPUT_DIR, f"{name}_out.pdbqt")
    log_path = os.path.join(OUTPUT_DIR, f"{name}_log.txt")

    print(f"\nDocking {name}...")

    cmd = (
        f'{VINA}'
        f' --receptor "{RECEPTOR}"'
        f' --ligand "{ligand_path}"'
        f' --center_x {center_x:.3f}'
        f' --center_y {center_y:.3f}'
        f' --center_z {center_z:.3f}'
        f' --size_x {size_x}'
        f' --size_y {size_y}'
        f' --size_z {size_z}'
        f' --exhaustiveness 8'
        f' --num_modes 9'
        f' --out "{out_path}"'
        f' --log "{log_path}"'
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Parse the log for binding affinity
    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            log_content = f.read()
        print(log_content[-500:])  # Print last part of log

        # Extract best binding affinity
        for line in log_content.split('\n'):
            line = line.strip()
            if line.startswith('   1'):
                parts = line.split()
                if len(parts) >= 2:
                    affinity = float(parts[1])
                    results.append((name, affinity))
                    print(f"  -> Best affinity: {affinity} kcal/mol")
                break
    else:
        print(f"  -> ERROR: No log file generated")
        if result.stderr:
            print(f"  -> STDERR: {result.stderr[:300]}")

# Step 4: Print summary
print("\n" + "=" * 60)
print("DOCKING RESULTS SUMMARY")
print("=" * 60)
print(f"{'Molecule':<20} {'Binding Affinity (kcal/mol)':<30}")
print("-" * 50)

# Sort by affinity (most negative = best)
results.sort(key=lambda x: x[1])
for name, affinity in results:
    marker = " <-- REFERENCE" if "Pazopanib" in name else ""
    print(f"{name:<20} {affinity:<30}{marker}")

print("\nMore negative = stronger binding")
print(f"\nResults saved to: {OUTPUT_DIR}")
print("Done!")
