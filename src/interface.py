"""Interface residue extraction for templates (heavy-atom vdW overlap + nearby Cα)."""

import os
import json
from itertools import combinations
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

from .utils import vdw_radii_extended, distance_calculator, DEFAULT_VDW

NEARBY_CA_CUTOFF = 6.0  # Å

INTERFACE_DIR = "templates/interfaces"
INTERFACE_LIST_DIR = "templates/interfaces_lists"
os.makedirs(INTERFACE_DIR, exist_ok=True)
os.makedirs(INTERFACE_LIST_DIR, exist_ok=True)


def _format_pdb_line(atom_serial, atom_name, res_name, chain_id, res_seq, x, y, z,
                    occupancy=1.00, bfactor=0.00, element="C"):
    return (
        f"ATOM  {atom_serial:5d} {atom_name:4s} {res_name:3s} {chain_id}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}           {element:2s}  \n"
    )


def _collect_chain_atoms(model, chain_ids):
    """Collect heavy atoms + Cα per residue for the requested chains."""
    chains_atoms = {cid: {} for cid in chain_ids}
    chains_cas = {cid: {} for cid in chain_ids}
    for chain in model:
        if chain.id not in chain_ids:
            continue
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            res_seq = int(residue.id[1])
            res_name = residue.get_resname()
            vdw_map = vdw_radii_extended(res_name)
            atoms = []
            ca_coords = None
            for atom in residue:
                name = atom.get_name()
                if name.startswith("H"):
                    continue
                coords = list(atom.get_coord())
                vdw = vdw_map.get(name, DEFAULT_VDW)
                atoms.append((name, coords, vdw))
                if name == "CA":
                    ca_coords = coords
            if atoms and ca_coords is not None:
                chains_atoms[chain.id][res_seq] = (res_name, atoms)
                chains_cas[chain.id][res_seq] = (res_name, ca_coords)
    return chains_atoms, chains_cas


def _find_interacting(atoms_a, atoms_b, tolerance=0.5):
    interacting_a, interacting_b = set(), set()
    for ra, (_, list_a) in atoms_a.items():
        for rb, (_, list_b) in atoms_b.items():
            for (_, ca, vdwa) in list_a:
                for (_, cb, vdwb) in list_b:
                    if distance_calculator(ca, cb) <= vdwa + vdwb + tolerance:
                        interacting_a.add(ra)
                        interacting_b.add(rb)
                        break
                else:
                    continue
                break
    return interacting_a, interacting_b


def _expand_by_neighbor(interacting, chain_cas, cutoff=NEARBY_CA_CUTOFF):
    extras = set()
    if not interacting:
        return interacting
    interact_cas = [chain_cas[r][1] for r in interacting]
    for res_seq, (_, ca) in chain_cas.items():
        if res_seq in interacting:
            continue
        for ic in interact_cas:
            if distance_calculator(ca, ic) <= cutoff:
                extras.add(res_seq)
                break
    interacting.update(extras)
    return interacting


def _write_ca_only(path, chain_id, chains_cas, interface_res):
    serial = 1
    with open(path, "w") as f:
        for res_seq in sorted(chains_cas[chain_id].keys()):
            if res_seq not in interface_res:
                continue
            res_name, (x, y, z) = chains_cas[chain_id][res_seq]
            f.write(_format_pdb_line(serial, " CA ", res_name, chain_id, res_seq, x, y, z, bfactor=1.00))
            serial += 1


def generate_interface(template, templates_root="templates"):
    """Generate interface artifacts for a template id like `1a28AB` (>=2 chains).

    For templates with more than 2 chains, generates per-pair artifacts where the
    interface PDB file for `chain` contains the union of contacts vs all other
    declared chains. Returns the interacting residue sets keyed by chain id.
    """
    if len(template) < 6:
        raise ValueError(f"Template {template!r} must be PDBID + >=2 chain letters")
    protein = template[:4].lower()
    chain_ids = list(template[4:])

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(template, os.path.join(templates_root, "pdbs", f"{protein}.pdb"))
    model = structure[0]
    chains_atoms, chains_cas = _collect_chain_atoms(model, set(chain_ids))

    interacting = {cid: set() for cid in chain_ids}
    for c1, c2 in combinations(chain_ids, 2):
        if not chains_atoms.get(c1) or not chains_atoms.get(c2):
            continue
        i1, i2 = _find_interacting(chains_atoms[c1], chains_atoms[c2])
        interacting[c1].update(i1)
        interacting[c2].update(i2)

    for cid in chain_ids:
        interacting[cid] = _expand_by_neighbor(interacting[cid], chains_cas[cid])

    for cid in chain_ids:
        path = os.path.join(INTERFACE_DIR, f"{template}_{cid}_int.pdb")
        _write_ca_only(path, cid, chains_cas, interacting[cid])

    with open(os.path.join(INTERFACE_LIST_DIR, f"{template}.json"), "w") as f:
        json.dump({cid: sorted(interacting[cid]) for cid in chain_ids}, f, indent=4)

    return interacting


if __name__ == "__main__":
    generate_interface("1a28AB")
