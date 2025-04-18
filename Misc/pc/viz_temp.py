import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Atom index, LAMMPS type (int), chemical type (str), coordinates
atom_lines = [
    (1, 6, 'O1', [20.760000, 23.016001, 10.858000]),
    (2, 7, 'O2', [21.527000, 23.042999, 8.741000]),
    (3, 8, 'O3', [20.033001, 24.650000, 9.417000]),
    (4, 1, 'C1', [21.844999, 22.098000, 10.790000]),
    (5, 2, 'C2', [21.934999, 21.811001, 9.321000]),
    (6, 3, 'C3', [21.563000, 20.871000, 11.634000]),
    (7, 4, 'C4', [20.712000, 23.662001, 9.651000]),
    (8, 5, 'H',  [22.742001, 22.611000, 11.156000]),
    (9, 5, 'H',  [22.950001, 21.556000, 9.007000]),
    (10, 5, 'H', [21.247000, 21.028000, 8.984000]),
    (11, 5, 'H', [20.624001, 20.393000, 11.335000]),
    (12, 5, 'H', [22.372000, 20.139999, 11.554000]),
    (13, 5, 'H', [21.448000, 21.153999, 12.687000])
]

fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot atoms
for idx, lmp_type, chem_type, pos in atom_lines:
    x, y, z = pos
    ax.scatter(x, y, z, s=150, alpha=0.8)
    ax.text(x + 0.1, y + 0.1, z + 0.1, f"{idx}\n{chem_type}", fontsize=9, ha='center')

# Formatting
ax.set_title("Propylene Carbonate: Atom Indices and Types", fontsize=14)
ax.set_xlabel("X (Å)")
ax.set_ylabel("Y (Å)")
ax.set_zlabel("Z (Å)")
ax.view_init(elev=20, azim=30)
ax.grid(False)
plt.tight_layout()
plt.show()
