"""Hotspot prediction (HotPoint-style: burial + contact potential)."""

import os
import json
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

from .utils import distance_calculator, three2one, PAIR_POTENTIAL
from .sasa_utils import get_asa_flat

RELATIVE_ASA_THRESHOLD = 20.0
CONTACT_POTENTIAL_THRESHOLD = 18.0
RESIDUE_DISTANCE_THRESHOLD = 3
CONTACT_DISTANCE_THRESHOLD = 7.0

HOTSPOT_DIR = "templates/hotspots"
os.makedirs(HOTSPOT_DIR, exist_ok=True)


def hotspot_creator(template, templates_root="templates"):
    """Identify hotspot residues for a 2-chain template (e.g. `1a28AB`)."""
    if len(template) < 6:
        raise ValueError(f"Template id {template!r} must be PDBID + 2 chain letters")
    chain_id1, chain_id2 = template[4], template[5]

    asa_complex = get_asa_flat(template, templates_root)
    contact_potentials = get_contact_potentials(template, templates_root)

    hotspot_dict = {chain_id1: [], chain_id2: []}
    for key, rel_asa in asa_complex.items():
        chain = key.split("_")[2]
        if chain not in hotspot_dict:
            continue
        cp = contact_potentials.get(key, 0.0)
        if rel_asa <= RELATIVE_ASA_THRESHOLD and abs(cp) >= CONTACT_POTENTIAL_THRESHOLD:
            res_name, res_num, _ = key.split("_")
            hotspot_dict[chain].append((res_num, res_name))

    with open(f"{HOTSPOT_DIR}/{template}.json", "w") as f:
        json.dump(hotspot_dict, f, indent=4)
    return hotspot_dict


def get_contact_potentials(template, templates_root="templates"):
    contact_dict = contacting_residues(template, templates_root)
    contact_potentials = {}
    for res1, neighbours in contact_dict.items():
        total = 0.0
        for res2 in neighbours:
            aa1 = three2one(res1.split("_")[0])
            aa2 = three2one(res2.split("_")[0])
            if aa1 == "X" or aa2 == "X":
                continue
            key_pair = "-".join(sorted([aa1, aa2]))
            total += PAIR_POTENTIAL.get(key_pair, 0.0)
        contact_potentials[res1] = total
    return contact_potentials


def contacting_residues(template, templates_root="templates"):
    cms = _center_of_mass(template, templates_root)
    contact_dict = {res: [] for res in cms}
    keys = list(cms.keys())
    for i, res1 in enumerate(keys):
        coor1 = cms[res1]
        chain1 = res1.split("_")[2]
        num1 = int(res1.split("_")[1])
        for res2 in keys[i + 1:]:
            chain2 = res2.split("_")[2]
            num2 = int(res2.split("_")[1])
            if chain1 == chain2 and abs(num1 - num2) <= RESIDUE_DISTANCE_THRESHOLD:
                continue
            if distance_calculator(coor1, cms[res2]) <= CONTACT_DISTANCE_THRESHOLD:
                contact_dict[res1].append(res2)
                contact_dict[res2].append(res1)
    return contact_dict


def _center_of_mass(template, templates_root):
    all_atoms = _fetch_all_atoms_coordinates(template, templates_root)
    centers = {}
    for residue, atoms in all_atoms.items():
        if not atoms:
            continue
        coords = list(atoms.values())
        centers[residue] = [sum(c[i] for c in coords) / len(coords) for i in range(3)]
    return centers


def _fetch_all_atoms_coordinates(template, templates_root):
    pdb_id = template[:4].lower()
    chains_keep = set(template[4:])
    pdb_path = os.path.join(templates_root, "pdbs", f"{pdb_id}.pdb")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_path)
    model = structure[0]

    all_atoms = {}
    for chain in model:
        if chain.id not in chains_keep:
            continue
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            key = f"{residue.get_resname()}_{residue.id[1]}_{chain.id}"
            all_atoms[key] = {atom.get_name(): list(atom.get_coord()) for atom in residue}
    return all_atoms


if __name__ == "__main__":
    hotspot_creator("1a28AB")
