"""Generate WBP5 protein structure image with pocket highlighted"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os

output_dir = r"C:\Users\user\Desktop\WBP5 Small molecule test\images"
os.makedirs(output_dir, exist_ok=True)

pdb_file = r"C:\Users\user\Desktop\WBP5 Small molecule test\WBP5_AlphaFold.pdb"

# Parse PDB - extract CA atoms
all_ca = []  # (x, y, z, res_num)
all_atoms = []  # all atoms for pocket
pocket_atoms = []

with open(pdb_file, 'r') as f:
    for line in f:
        if line.startswith("ATOM"):
            atom_name = line[12:16].strip()
            res_num = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if atom_name == "CA":
                all_ca.append((x, y, z, res_num))

            if 64 <= res_num <= 74:
                pocket_atoms.append((x, y, z))

print(f"Total CA atoms: {len(all_ca)}")
print(f"Pocket atoms (res 64-74): {len(pocket_atoms)}")

# Create figure
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot protein backbone (CA atoms)
ca_x = [a[0] for a in all_ca]
ca_y = [a[1] for a in all_ca]
ca_z = [a[2] for a in all_ca]

# Color by residue number for rainbow effect
res_nums = [a[3] for a in all_ca]
colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(all_ca)))

# Plot backbone as line
ax.plot(ca_x, ca_y, ca_z, color='#4A90D9', alpha=0.6, linewidth=1.5, label='WBP5 backbone')

# Plot CA atoms
ax.scatter(ca_x, ca_y, ca_z, c='#4A90D9', s=15, alpha=0.4)

# Highlight P_0 pocket in red
pocket_x = [a[0] for a in pocket_atoms]
pocket_y = [a[1] for a in pocket_atoms]
pocket_z = [a[2] for a in pocket_atoms]

ax.scatter(pocket_x, pocket_y, pocket_z, c='red', s=40, alpha=0.8,
           label='P_0 pocket (res 64-74)', edgecolors='darkred', linewidth=0.5)

# Mark pocket center
cx = sum(pocket_x) / len(pocket_x)
cy = sum(pocket_y) / len(pocket_y)
cz = sum(pocket_z) / len(pocket_z)
ax.scatter([cx], [cy], [cz], c='yellow', s=200, marker='*',
           edgecolors='red', linewidth=1.5, label=f'Pocket center ({cx:.1f}, {cy:.1f}, {cz:.1f})', zorder=5)

# Formatting
ax.set_xlabel('X (Angstrom)', fontsize=10, labelpad=10)
ax.set_ylabel('Y (Angstrom)', fontsize=10, labelpad=10)
ax.set_zlabel('Z (Angstrom)', fontsize=10, labelpad=10)
ax.set_title('WBP5 (Q9UHQ7) AlphaFold Structure\nwith P_0 Druggable Pocket (Drug Score 0.84)',
             fontsize=14, fontweight='bold', pad=20)
ax.legend(loc='upper left', fontsize=9, framealpha=0.9)

# Set viewing angle
ax.view_init(elev=20, azim=45)
ax.set_box_aspect([1,1,1])

# Add text annotation
ax.text2D(0.02, 0.02, 'P_0 pocket: HLSNEDMFREV (res 64-74)\nDrug Score: 0.84',
          transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
          bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
save_path = os.path.join(output_dir, "WBP5_structure_pocket.png")
plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
print(f"Saved: {save_path}")

# Create a second view (rotated)
ax.view_init(elev=10, azim=135)
save_path2 = os.path.join(output_dir, "WBP5_structure_pocket_view2.png")
plt.savefig(save_path2, dpi=200, bbox_inches='tight', facecolor='white')
print(f"Saved: {save_path2}")

print("Done!")
