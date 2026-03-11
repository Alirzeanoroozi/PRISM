import numpy as np
import json
import os
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

CONTACT_DISTANCE_THRESHOLD = 4.0 # Å

CONTACT_DIR = "templates/contacts"
os.makedirs(CONTACT_DIR, exist_ok=True)

def get_contacts(template):
    protein = template[:4].lower()
    chain1 = template[4]
    chain2 = template[5]
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(template, f"templates/pdbs/{protein}.pdb")
    model = structure[0]

    chain1 = model[chain1]
    chain2 = model[chain2]

    def get_atoms(residue):
        for atom in residue:
            if atom.element == "H":
                continue
            yield atom

    contacts = []
    for res1 in chain1:
        if not is_aa(res1):
            continue
        coords1 = np.array([a.get_coord() for a in get_atoms(res1)])
        if len(coords1) == 0:
            continue

        for res2 in chain2:
            if not is_aa(res2):
                continue
            coords2 = np.array([a.get_coord() for a in get_atoms(res2)])
            if len(coords2) == 0:
                continue

            # Minimum distance between any atom pair
            dists = np.linalg.norm(coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :], axis=2)
            if np.min(dists) <= CONTACT_DISTANCE_THRESHOLD:
                contacts.append((res1.get_id()[1], res2.get_id()[1]))
    with open(f"{CONTACT_DIR}/{template}.json", "w") as f:
        json.dump(contacts, f)
    return contacts

def get_contacts_from_atom_lines(pdb_path, output_path, atom_lines_0, atom_lines_1):
    """
    Compute inter-chain contacts between two sets of ATOM lines and write
    residue pairs to output_path. Supports multi-chain receptor/ligand.
    """
    import numpy as np
    from Bio.PDB.Polypeptide import is_aa

    def parse_atom_line(line):
        if len(line) < 54:
            return None
        try:
            res_name = line[17:20].strip()
            chain = line[21]
            res_id = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atom_name = line[12:16].strip()
            return (chain, res_id, res_name, atom_name, np.array([x, y, z]))
        except (ValueError, IndexError):
            return None

    def get_atoms_by_residue(lines):
        """Group atoms by (chain, res_id)."""
        from collections import defaultdict
        res_atoms = defaultdict(list)
        for line in lines:
            parsed = parse_atom_line(line)
            if parsed is None:
                continue
            chain, res_id, res_name, atom_name, coord = parsed
            if atom_name == "H":
                continue
            key = (chain, res_id, res_name)
            res_atoms[key].append(coord)
        return res_atoms

    residues_0 = get_atoms_by_residue(atom_lines_0)
    residues_1 = get_atoms_by_residue(atom_lines_1)
    contacts = []

    for (c1, r1, name1), coords1 in residues_0.items():
        if not is_aa(name1, standard=True) and name1 not in ("HOH", "H2O", "WAT"):
            continue
        arr1 = np.array(coords1)
        if len(arr1) == 0:
            continue

        for (c2, r2, name2), coords2 in residues_1.items():
            if not is_aa(name2, standard=True) and name2 not in ("HOH", "H2O", "WAT"):
                continue
            arr2 = np.array(coords2)
            if len(arr2) == 0:
                continue

            dists = np.linalg.norm(arr1[:, np.newaxis, :] - arr2[np.newaxis, :, :], axis=2)
            if np.min(dists) <= CONTACT_DISTANCE_THRESHOLD:
                contacts.append((r1, r2))

    with open(output_path, "w") as f:
        for r1, r2 in contacts:
            f.write(f"{r1}\t{r2}\n")


if __name__ == "__main__":
    template = "1a28AB"
    print(get_contacts(template))