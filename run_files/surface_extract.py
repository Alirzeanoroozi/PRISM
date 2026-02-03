import os
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

from naccess_utils import get_asa_complex
from utils import distance_calculator

RSATHRESHOLD = 15.0
SCFFTHRESHOLD = 1.4

SURFACE_EXTRACTION_DIR = "surface_extraction"
os.makedirs(SURFACE_EXTRACTION_DIR, exist_ok=True)

def extract_surfaces(queries):
    for protein in queries:
        extract_surface(protein)

def extract_surface(protein):
    asa_complex = get_asa_complex(protein, SURFACE_EXTRACTION_DIR)
    rsa_residues = [key for key in asa_complex if asa_complex[key] > RSATHRESHOLD]
    if len(rsa_residues) == 0:
        return {}
    pdb_path = f"pdbs/{protein}.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(protein, pdb_path)

    # Build CA coords for rsa_residues: key "RESNAME_NUM_CHAIN" -> (chain, res_num)
    rsa_keys = set()
    for key in rsa_residues:
        parts = key.split("_")
        if len(parts) >= 3:
            rsa_keys.add((parts[2], parts[1]))  # (chain, res_num)

    rsa_ca_coords = []
    all_ca_data = []  # (chain_id, res_num_str, res_name, res_seq, coords)

    for model in structure:
        for chain in model:
            for residue in chain:
                if not is_aa(residue, standard=True):
                    continue
                hetflag, res_seq, icode = residue.id
                if hetflag != " ":
                    continue
                res_num_str = str(res_seq)
                res_name = residue.get_resname()
                if "CA" in residue:
                    ca = residue["CA"]
                    coords = list(ca.get_coord())
                    all_ca_data.append((chain.id, res_num_str, res_name, res_seq, coords))
                    if (chain.id, res_num_str) in rsa_keys:
                        rsa_ca_coords.append(coords)
        break  # use first model only

    # Find CAs within SCFFTHRESHOLD of any rsa_residue CA
    asa_ca_lines = []
    for chain_id, res_num_str, res_name, res_seq, coords in all_ca_data:
        for rsa_coords in rsa_ca_coords:
            if distance_calculator(coords, rsa_coords) <= SCFFTHRESHOLD:
                asa_ca_lines.append((res_name, chain_id, res_seq, coords))
                break

    # Write with sequential serial numbers
    asa_ca_lines = [
        _format_ca_pdb_line(i, res_name, chain_id, res_seq, coords[0], coords[1], coords[2])
        for i, (res_name, chain_id, res_seq, coords) in enumerate(asa_ca_lines, start=1)
    ]

    asa_path = f"{SURFACE_EXTRACTION_DIR}/{protein}.asa.pdb"
    with open(asa_path, "w") as f:
        for line in asa_ca_lines:
            f.write(line)
        f.write("END\n")

    return {i: line for i, line in enumerate(asa_ca_lines)}

def _format_ca_pdb_line(serial, res_name, chain_id, res_seq, x, y, z):
    """Format a CA ATOM line for PDB output."""
    return (
        f"ATOM  {serial:5d}  CA  {res_name:3s} {chain_id}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"
    )

if __name__ == "__main__":
    extract_surfaces(["1a28"])