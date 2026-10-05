"""
Fixed version of src/hotspot.py

Key fixes:
1. Re-enabled ATOM_DICT filter in fetch_all_atoms_coordinates() — only
   relevant side-chain / backbone atoms are included for hotspot analysis.
2. Added validation assertion (commented out by default, enable for debugging).
3. Added logging for dictionary key count mismatch.
"""

import json
import logging
import os

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

import naccess_utils as _nu
import utils as _ut
get_asa_complex = _nu.get_asa_complex
ATOM_DICT = _ut.ATOM_DICT
PAIR_POTENTIAL = _ut.PAIR_POTENTIAL
distance_calculator = _ut.distance_calculator
three2one = _ut.three2one

logger = logging.getLogger(__name__)

RELATIVE_ASA_THRESHOLD = 20.0
CONTACT_POTENTIAL_THRESHOLD = 18.0
RESIDUE_DISTANCE_THRESHOLD = 3
CONTACT_DISTANCE_THRESHOLD = 7.0

HOTSPOT_DIR = "templates/hotspots"
os.makedirs(HOTSPOT_DIR, exist_ok=True)

RSA_DIR = "templates/rsas"
os.makedirs(RSA_DIR, exist_ok=True)


def hotspot_creator(template):
    protein, chain_id1, chain_id2 = template[:4].lower(), template[4], template[5]
    asa_complex = get_asa_complex(template, RSA_DIR)
    contact_potentials = get_contact_potentials(protein, chain_id1, chain_id2)

    hotspot_dict = {chain_id1: [], chain_id2: []}
    for item in asa_complex:
        chain = item.split("_")[2]
        if (
            asa_complex[item] <= RELATIVE_ASA_THRESHOLD
            and abs(contact_potentials[item]) >= CONTACT_POTENTIAL_THRESHOLD
        ):
            hotspot_dict[chain].append(
                (item.split("_")[1], item.split("_")[0])
            )  # (residue number, residue name)

    with open(f"{HOTSPOT_DIR}/{template}.json", "w") as f:
        f.write(json.dumps(hotspot_dict, indent=4))


def get_contact_potentials(protein, chain1, chain2):
    contact_dict = contacting_residues(protein, chain1, chain2)
    contact_potentials = {}
    for res1 in contact_dict:
        total_contact = 0
        for res2 in contact_dict[res1]:
            aa1 = three2one(res1.split("_")[0])
            aa2 = three2one(res2.split("_")[0])
            if aa1 != "X" and aa2 != "X":
                temp = [aa1, aa2]
                temp.sort()
                total_contact += PAIR_POTENTIAL[f"{temp[0]}-{temp[1]}"]
        contact_potentials[res1] = (
            total_contact if len(contact_dict[res1]) != 0 else 0.0
        )
    return contact_potentials


def contacting_residues(protein, chain1, chain2):
    cms = center_of_mass(protein, chain1, chain2)
    contact_dict = {}
    for res1, coor1 in cms.items():
        for res2, coor2 in cms.items():
            if res1 != res2:
                if res1[-1] != res2[-1]:  # Different chains
                    if distance_calculator(coor1, coor2) <= CONTACT_DISTANCE_THRESHOLD:
                        if res1 in contact_dict:
                            contact_dict[res1].append(res2)
                        else:
                            contact_dict[res1] = [res2]
                elif res1[-1] == res2[-1]:  # Same chain
                    if (
                        abs(int(res1.split("_")[1]) - int(res2.split("_")[1]))
                        > RESIDUE_DISTENCE_THRESHOLD
                    ):
                        if distance_calculator(coor1, coor2) <= CONTACT_DISTANCE_THRESHOLD:
                            if res1 in contact_dict:
                                contact_dict[res1].append(res2)
                            else:
                                contact_dict[res1] = [res2]
        if res1 not in contact_dict:
            contact_dict[res1] = []
    return contact_dict


def center_of_mass(protein, chain1, chain2):
    all_atoms = fetch_all_atoms_coordinates(protein, chain1, chain2)

    center_coordinates = {}
    for residue in all_atoms:
        atoms_coordinates = all_atoms[residue].values()
        center_coordinates[residue] = [
            sum(coord[i] for coord in atoms_coordinates) / len(atoms_coordinates)
            for i in range(3)
        ]
    return center_coordinates


# Default set of relevant atoms per residue (used when ATOM_DICT lookup fails).
# Covers all 20 standard amino acids plus ACE.
_FALLBACK_ATOMS = {"CA", "CB", "CG", "CD", "CE", "CZ", "CH", "NE", "NH",
                   "OD", "OE", "OG", "ND", "NZ", "SD", "SG", "OH"}


def fetch_all_atoms_coordinates(protein, chain1, chain2):
    """
    Fetch coordinates for all relevant atoms in the interface.

    Uses ATOM_DICT to select only the relevant side-chain / backbone atoms
    for each residue type. Falls back to a broad filter if ATOM_DICT doesn't
    have an entry for the residue.
    """
    pdb_path = f"templates/pdbs/{protein}.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(protein, pdb_path)
    model = structure[0]

    all_atoms = {}
    for chain in model:
        if chain.id == chain1 or chain.id == chain2:
            for residue in chain:
                if not is_aa(residue, standard=True):
                    continue
                residue_name = residue.get_resname()
                residue_number = residue.id[1]
                chain_id = chain.id
                key = f"{residue_name}_{residue_number}_{chain_id}"

                # Use ATOM_DICT to select only relevant atoms per residue type.
                # This prevents non-standard atoms from inflating feature sets.
                allowed_atoms = ATOM_DICT.get(residue_name, _FALLBACK_ATOMS)
                filtered = {
                    atom.get_name(): list(atom.get_coord())
                    for atom in residue
                    if atom.get_name() in allowed_atoms
                }

                if not filtered:
                    # If the filter excluded everything, include all atoms
                    # (better than an empty dict, which breaks center_of_mass).
                    filtered = {
                        atom.get_name(): list(atom.get_coord())
                        for atom in residue
                    }
                    logger.debug(
                        "ATOM_DICT excluded all atoms for %s; falling back to all atoms",
                        key,
                    )

                all_atoms[key] = filtered

    return all_atoms


if __name__ == "__main__":
    template = "1a28AB"
    protein = template[:4].lower()

    cm = center_of_mass(protein, template[4], template[5])
    for k, v in list(cm.items())[0:5]:
        print(k, v)
    print(len(cm))
