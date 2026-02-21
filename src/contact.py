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

if __name__ == "__main__":
    template = "1a28AB"
    print(get_contacts(template))