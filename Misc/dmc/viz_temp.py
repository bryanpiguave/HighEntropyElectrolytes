import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Atom index, LAMMPS type, chemical type, coordinates
atom_lines = [
    (1, 7, 'O2', [16.037001, 19.371000, 12.196000]),
    (2, 8, 'O3', [16.431000, 19.638000, 14.343000]),
    (3, 6, 'O1', [15.178000, 21.223000, 13.249000]),
    (4, 1, 'C1', [15.817000, 20.183001, 13.261000]),
    (5, 2, 'C2', [15.444000, 19.827999, 10.983000]),
    (6, 3, 'C3', [16.277000, 20.398001, 15.542000]),
    (7, 4, 'H',  [15.682000, 19.106001, 10.196000]),
    (8, 4, 'H',  [14.356000, 19.882000, 11.081000]),
    (9, 4, 'H',  [15.856000, 20.799000, 10.693000]),
    (10, 4, 'H', [16.721001, 21.392000, 15.428000]),
    (11, 4, 'H', [15.221000, 20.475000, 15.817000]),
    (12, 4, 'H', [16.804001, 19.875999, 16.344999]),
]

fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot atoms
for idx, lmp_type, chem_type, pos in atom_lines:
    x, y, z = pos
    ax.scatter(x, y, z, s=150, alpha=0.8)
    ax.text(x + 0.1, y + 0.1, z + 0.1, f"{idx}\n{chem_type}", fontsize=9, ha='center')

# Formatting
ax.set_title("Dimethyl Carbonate: Atom Indices and Types", fontsize=14)
ax.set_xlabel("X (Å)")
ax.set_ylabel("Y (Å)")
ax.set_zlabel("Z (Å)")
ax.view_init(elev=20, azim=30)
ax.grid(False)
plt.tight_layout()
plt.show()
