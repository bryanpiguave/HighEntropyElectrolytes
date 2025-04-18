import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Atom index, type, coordinates
atom_lines = [
    (1, 'O2', [22.96, 22.875, 10.778]),
    (2, 'O3', [21.424, 22.692, 9.144]),
    (3, 'O1', [22.912, 24.437, 9.098]),
    (4, 'C2', [22.349, 21.602, 10.920]),
    (5, 'C3', [21.074, 21.757, 10.155]),
    (6, 'C1', [22.472, 23.427, 9.625]),
    (7, 'H', [22.181, 21.363, 11.974]),
    (8, 'H', [23.012, 20.848, 10.481]),
    (9, 'H', [20.73, 20.817, 9.716]),
    (10, 'H', [20.271, 22.186, 10.762])
]

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Plot atoms
for idx, atype, pos in atom_lines:
    x, y, z = pos
    ax.scatter(x, y, z, s=150, label=atype if idx == 1 else "", alpha=0.8)
    ax.text(x + 0.1, y + 0.1, z + 0.1, f"{idx}\n{atype}", fontsize=9, ha='center')

# Formatting
ax.set_title("EC Molecule: Atom Indices and Types", fontsize=14)
ax.set_xlabel("X (Å)")
ax.set_ylabel("Y (Å)")
ax.set_zlabel("Z (Å)")
ax.view_init(elev=20, azim=30)
ax.grid(False)
plt.tight_layout()

plt.show()
