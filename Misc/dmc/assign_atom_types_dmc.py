from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

# === Manual Mapping for Dimethyl Carbonate (DMC) ===
atom_type_map = {
    0: "O2",
    1: "O3",
    2: "O1",
    3: "C1",
    4: "C2",
    5: "C3",
    6: "H",
    7: "H",
    8: "H",
    9: "H",
    # Extend for more H if needed
}

atom_params = {
    "O1": {"sigma": 2.96, "epsilon": 0.210, "charge": -0.6452},
    "O2": {"sigma": 3.00, "epsilon": 0.170, "charge": -0.4684},
    "O3": {"sigma": 3.00, "epsilon": 0.170, "charge": -0.4684},
    "C1": {"sigma": 3.75, "epsilon": 0.105, "charge": 1.0996},
    "C2": {"sigma": 3.50, "epsilon": 0.066, "charge": 0.0330},
    "C3": {"sigma": 3.50, "epsilon": 0.066, "charge": 0.0330},
    "H":  {"sigma": 2.50, "epsilon": 0.030, "charge": 0.1041},
}

# === Load Molecule ===
mol = Chem.MolFromMol2File("dmc.mol2", removeHs=False)
conf = mol.GetConformer()
atoms = mol.GetAtoms()

# === Write MOL2 with atom types + charges ===
with open("dmc_custom.mol2", "w") as f:
    f.write("@<TRIPOS>MOLECULE\nDMC\n{}\n0\n0\n0\n0\nSMALL\n".format(len(atoms)))
    f.write("\n@<TRIPOS>ATOM\n")
    for atom in atoms:
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        atype = atom_type_map.get(idx, "H")
        charge = atom_params[atype]["charge"]
        f.write("{:>6} {:<4} {:>9.4f} {:>9.4f} {:>9.4f} {:<4} 1 DMC {:>7.4f}\n".format(
            idx + 1, atype, pos.x, pos.y, pos.z, atype, charge))
    f.write("@<TRIPOS>BOND\n")

print("✅ Saved: dmc_custom.mol2")

# === Generate PNG with atom indices ===
AllChem.Compute2DCoords(mol)
drawer = rdMolDraw2D.MolDraw2DCairo(500, 400)
drawer.drawOptions().addAtomIndices = True
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
with open("dmc_atom_indices.png", "wb") as img_file:
    img_file.write(drawer.GetDrawingText())

print("🖼️  Saved: dmc_atom_indices.png")
