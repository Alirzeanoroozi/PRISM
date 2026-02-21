import os
from Bio.PDB import PDBParser
import json
from Bio.PDB.Polypeptide import is_aa

from .utils import vdw_radii_extended, distance_calculator, DEFAULT_VDW

# Nearby residue cutoff: Cα of non-interacting residue within this distance of an interacting residue (same chain)
NEARBY_CA_CUTOFF = 6.0  # Å

INTERFACE_DIR = "templates/interfaces"
os.makedirs(INTERFACE_DIR, exist_ok=True)

INTERFACE_LIST_DIR = "templates/interfaces_lists"
os.makedirs(INTERFACE_LIST_DIR, exist_ok=True)

BFACTOR_DIR = "templates/bfactor_pdbs"
os.makedirs(BFACTOR_DIR, exist_ok=True)

def _format_pdb_line(atom_serial, atom_name, res_name, chain_id, res_seq, x, y, z, occupancy=1.00, bfactor=0.00, element="C"):
    """Format a single ATOM line. atom_name should be 4 chars (e.g. ' CA ')."""
    return (
        f"ATOM  {atom_serial:5d} {atom_name:4s} {res_name:3s} {chain_id}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}           {element:2s}  \n"
    )

def write_ca_only_pdb(path, chain_id, chains_cas, interface_res):
    serial = 1
    with open(path, "w") as f:
        for res_seq in sorted(chains_cas[chain_id].keys()):
            if res_seq not in interface_res:
                continue
            res_name, (x, y, z) = chains_cas[chain_id][res_seq]
            f.write(_format_pdb_line(serial, " CA ", res_name, chain_id, res_seq, x, y, z, bfactor=1.00))
            serial += 1

def generate_interface(template: str):
    protein = template[:4].lower()
    chain1 = template[4]
    chain2 = template[5]
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(template, f"templates/pdbs/{protein}.pdb")
    model = structure[0]

    # Collect all atoms (and Cα) for the two chains
    # chain_data[c][res_seq] = list of (atom_name, coords, vdw)
    # chain_ca[c][res_seq] = (res_name, coords) for Cα
    chains_atoms = {chain1: {}, chain2: {}}
    chains_cas = {chain1: {}, chain2: {}}

    for chain in model:
        if chain.id not in (chain1, chain2):
            continue
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            res_seq = residue.id[1]
            res_name = residue.get_resname()
            vdw_map = vdw_radii_extended(res_name)
            atoms = []
            ca_coords = None
            for atom in residue:
                name = atom.get_name()
                coords = list(atom.get_coord())
                vdw = vdw_map.get(name, DEFAULT_VDW)
                atoms.append((name, coords, vdw))
                if name == "CA":
                    ca_coords = coords
            if atoms and ca_coords is not None:
                chains_atoms[chain.id][res_seq] = (res_name, atoms)
                chains_cas[chain.id][res_seq] = (res_name, ca_coords)

    # Find interacting residues: any atom from chain1 vs any atom from chain2 within vdW sum + 0.5
    interacting1 = set()
    interacting2 = set()
    for res_seq1, (res_name1, atoms1) in chains_atoms[chain1].items():
        for res_seq2, (res_name2, atoms2) in chains_atoms[chain2].items():
            for (name1, coords1, vdw1) in atoms1:
                if name1[0] == "H":
                    continue
                for (name2, coords2, vdw2) in atoms2:
                    if name2[0] == "H":
                        continue
                    if distance_calculator(coords1, coords2) <= vdw1 + vdw2 + 0.5:
                        interacting1.add(res_seq1)
                        interacting2.add(res_seq2)

    # Add nearby residues: same chain, Cα within 6 Å of an interacting residue
    def add_nearby(chain_id):
        temp_set = set()
        interact = interacting1 if chain_id == chain1 else interacting2
        ca_dict = chains_cas[chain_id]
        for res_seq, (res_name, ca_coords) in ca_dict.items():
            if res_seq in interact:
                continue
            for other_seq in interact:
                _, other_ca = chains_cas[chain_id][other_seq]
                if distance_calculator(ca_coords, other_ca) <= NEARBY_CA_CUTOFF:
                    temp_set.add(res_seq)
        interact.update(temp_set)

    add_nearby(chain1)
    add_nearby(chain2)

    # Write Cα-only PDB files (interface residues only, one file per chain)
    ca_path1 = os.path.join(INTERFACE_DIR, f"{template}_{chain1}_int.pdb")
    write_ca_only_pdb(ca_path1, chain1, chains_cas, interacting1)

    ca_path2 = os.path.join(INTERFACE_DIR, f"{template}_{chain2}_int.pdb")
    write_ca_only_pdb(ca_path2, chain2, chains_cas, interacting2)

    # Write full PDB with bFactor = 1 for interface, 0 for non-interface
    bfactor_path = os.path.join(BFACTOR_DIR, f"{template}_bfactor.pdb")
    lines_out = []
    serial = 1
    for chain in model:
        cid = chain.id
        interface_set = interacting1 if cid == chain1 else (interacting2 if cid == chain2 else set())
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            res_seq = residue.id[1]
            bfactor_val = 99.0 if res_seq in interface_set else 0.0
            for atom in residue:
                name = atom.get_name()
                coords = atom.get_coord()
                res_name = residue.get_resname()
                elem = "C" if name == "CA" else name[0]
                atom_name_4 = (" " + name).ljust(4)[:4]  # PDB atom name 4 chars
                occ = atom.get_occupancy()
                lines_out.append(
                    _format_pdb_line(
                        serial, atom_name_4, res_name, cid, res_seq,
                        coords[0], coords[1], coords[2],
                        occupancy=occ, bfactor=bfactor_val, element=elem
                    )
                )
                serial += 1
        break  # single model
    with open(bfactor_path, "w") as f:
        f.writelines(lines_out)
        f.write("END\n")

    interface_lists = {chain1:  list(interacting1), chain2: list(interacting2)}
    with open(os.path.join(INTERFACE_LIST_DIR, f"{template}.json"), "w") as f:
        f.write(json.dumps(interface_lists, indent=4))

    return {
        "ca_chain1": ca_path1,
        "ca_chain2": ca_path2,
        "bfactor_pdb": bfactor_path,
        "interface_residues_chain1": interacting1,
        "interface_residues_chain2": interacting2,
    }

if __name__ == "__main__":
    protein = "1a28"
    chain1 = "A"
    chain2 = "B"
    result = generate_interface(protein, chain1, chain2)
    print(result)
