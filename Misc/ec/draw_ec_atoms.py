from rdkit import Chem
from rdkit.Chem import Draw

# Load EC molecule
mol = Chem.MolFromMolFile("ec.sdf", removeHs=False)
assert mol is not None, "Failed to load ec.sdf"

# 1. Draw molecule with indices
for atom in mol.GetAtoms():
    atom.SetProp("atomLabel", str(atom.GetIdx()))
img = Draw.MolToImage(mol, size=(500, 500), kekulize=True)
img.save("ec_indexed.png")
print("✅ Saved molecule as ec_indexed.png")

# 2. Print basic atom info
print("\nIndex\tElement\tNeighbors (by index)")
for atom in mol.GetAtoms():
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
    print(f"{atom.GetIdx():<6}\t{atom.GetSymbol():<6}\t{neighbors}")

# 3. Placeholder: assign atom types manually here
# -------------------------------------------------
# Update this dictionary once you’ve identified types from ec_indexed.png
atom_type_map = {
    0: "O2",
    1: "O3",
    2: "O1",
    3: "C2",
    4: "C3",
    5: "C1",
    6: "H1",
    7: "H2",
    8: "H3",
    9: "H4"
}
# -------------------------------------------------

# 4. Sanity check
print("\nAssigned Atom Types:")
for idx, atom_type in atom_type_map.items():
    print(f"Atom {idx} → {atom_type}")
