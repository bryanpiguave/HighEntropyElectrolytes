import os
import re

def short_tag(atom_name, mol_prefix):
    """
    Generate a short 4-character atom name tag.
    Uses:
    - First letter of mol_prefix (e.g., 'E' for EC)
    - Atom index or name-based number (e.g., '1')
    - Abbreviated core atom name (e.g., 'C' for C1, O for O2)
    Example: EC_C1 → E1C
    """
    idx = ''.join(re.findall(r'\d+', atom_name)) or "X"
    core = ''.join(re.findall(r'[A-Za-z]+', atom_name))[:2] or "X"
    tag = (mol_prefix[0].upper() + idx + core[0].upper())[:4]
    return tag

def tag_atom_names_short(input_file, output_file, mol_id):
    """
    Updates atom names in a .mol2 file with unique 4-character tags.
    This ensures compatibility with Packmol (PDB format limitation).
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()

    in_atoms = False
    new_lines = []

    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            new_lines.append(line)
            continue
        elif line.startswith("@<TRIPOS>"):
            in_atoms = False
            new_lines.append(line)
            continue

        if in_atoms:
            parts = line.split()
            if len(parts) >= 2:
                original_name = parts[1]
                parts[1] = short_tag(original_name, mol_id)
            new_lines.append(" ".join(parts) + "\n")
        else:
            new_lines.append(line)

    with open(output_file, 'w') as f:
        f.writelines(new_lines)

    print(f"Tagged: {input_file} → {output_file} using mol_id '{mol_id}'")

def main():
    # Define your molecule name-to-file mapping here
    mol2_files = {
        "EC": "ec.mol2",
        "DMC": "dmc.mol2",
        "PC": "pc.mol2"
    }

    for mol_id, infile in mol2_files.items():
        if not os.path.isfile(infile):
            print(f"❌ Skipping missing file: {infile}")
            continue
        outfile = f"{mol_id.lower()}_tagged.mol2"
        tag_atom_names_short(infile, outfile, mol_id)

if __name__ == "__main__":
    main()
