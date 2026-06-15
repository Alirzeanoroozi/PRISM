"""Residue contact maps for templates and refined structures."""

import os
import json
import numpy as np
from itertools import combinations
from collections import defaultdict
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

CONTACT_DISTANCE_THRESHOLD = 4.0  # Å

CONTACT_DIR = "templates/contacts"
os.makedirs(CONTACT_DIR, exist_ok=True)


def _residue_heavy_atoms(chain):
    for residue in chain:
        if not is_aa(residue, standard=True):
            continue
        coords = []
        for atom in residue:
            if atom.element == "H":
                continue
            coords.append(atom.get_coord())
        if coords:
            yield int(residue.id[1]), np.array(coords)


def _pairwise_contacts(coords1_map, coords2_map):
    contacts = []
    for r1, arr1 in coords1_map:
        for r2, arr2 in coords2_map:
            dists = np.linalg.norm(arr1[:, None, :] - arr2[None, :, :], axis=2)
            if dists.min() <= CONTACT_DISTANCE_THRESHOLD:
                contacts.append((r1, r2))
    return contacts


def get_contacts(template, templates_root="templates"):
    """Compute residue-residue contact pairs for all chain pairs of a template.

    Returns a dict keyed by `(chain_a, chain_b)` -> list of (res_a, res_b).
    Also writes a JSON next to the template. Backwards-compatible: a 2-chain
    template writes a flat list of pairs (matching legacy artifact format).
    """
    if len(template) < 6:
        raise ValueError(f"Template {template!r} must be PDBID + >=2 chain letters")
    protein = template[:4].lower()
    chain_ids = list(template[4:])

    pdb_path = os.path.join(templates_root, "pdbs", f"{protein}.pdb")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(template, pdb_path)
    model = structure[0]

    coords_per_chain = {}
    for cid in chain_ids:
        if cid in model.child_dict:
            coords_per_chain[cid] = list(_residue_heavy_atoms(model[cid]))
        else:
            coords_per_chain[cid] = []

    pair_contacts = {}
    for c1, c2 in combinations(chain_ids, 2):
        pair_contacts[(c1, c2)] = _pairwise_contacts(coords_per_chain[c1], coords_per_chain[c2])

    out_path = os.path.join(CONTACT_DIR, f"{template}.json")
    if len(chain_ids) == 2:
        c1, c2 = chain_ids
        payload = pair_contacts[(c1, c2)]
    else:
        payload = {f"{c1}{c2}": v for (c1, c2), v in pair_contacts.items()}
    with open(out_path, "w") as f:
        json.dump(payload, f)
    return pair_contacts


def get_contacts_from_atom_lines(pdb_path, output_path, atom_lines_0, atom_lines_1):
    """Compute inter-partner residue contacts from two sets of ATOM lines."""
    def parse_atom_line(line):
        if len(line) < 54:
            return None
        try:
            res_name = line[17:20].strip()
            chain = line[21]
            res_id = int(line[22:26].strip())
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            atom_name = line[12:16].strip()
            return (chain, res_id, res_name, atom_name, np.array([x, y, z]))
        except (ValueError, IndexError):
            return None

    def group_by_residue(lines):
        residues = defaultdict(list)
        for line in lines:
            parsed = parse_atom_line(line)
            if parsed is None:
                continue
            chain, res_id, res_name, atom_name, coord = parsed
            if atom_name.startswith("H"):
                continue
            residues[(chain, res_id, res_name)].append(coord)
        return residues

    res0 = group_by_residue(atom_lines_0)
    res1 = group_by_residue(atom_lines_1)
    contacts = []
    for (c1, r1, n1), coords1 in res0.items():
        if not is_aa(n1, standard=True):
            continue
        arr1 = np.asarray(coords1)
        for (c2, r2, n2), coords2 in res1.items():
            if not is_aa(n2, standard=True):
                continue
            arr2 = np.asarray(coords2)
            dists = np.linalg.norm(arr1[:, None, :] - arr2[None, :, :], axis=2)
            if dists.min() <= CONTACT_DISTANCE_THRESHOLD:
                contacts.append((c1, r1, c2, r2))
    with open(output_path, "w") as f:
        for c1, r1, c2, r2 in contacts:
            f.write(f"{c1}\t{r1}\t{c2}\t{r2}\n")
    return contacts


if __name__ == "__main__":
    print(get_contacts("1a28AB"))
