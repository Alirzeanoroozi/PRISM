import os
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

from .sasa_utils import get_asa_complex
from .utils import distance_calculator

RSATHRESHOLD = 25.0
SCFFTHRESHOLD = 1.4

SURFACE_EXTRACTION_DIR = "processed/surface_extraction"
os.makedirs(SURFACE_EXTRACTION_DIR, exist_ok=True)

def extract_surfaces(queries):
    failed_count = 0
    for protein in queries:
        if not extract_surface(protein):
            failed_count += 1
    return failed_count

def extract_surface(protein):
    asa_complex = get_asa_complex(protein, "processed")
    rsa_residues = {chain: [res_num for res_num in asa_complex[chain] if asa_complex[chain][res_num] > RSATHRESHOLD] for chain in asa_complex}

    pdb_path = f"processed/pdbs/{protein[:4].lower()}.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(protein, pdb_path)
    model = structure[0]

    rsa_ca_coords = {chain.id: [] for chain in model}
    all_ca_data = {chain.id: [] for chain in model}  # {chain: [(res_num_str, res_name, res_seq, coords)]}

    for chain in model:
        if chain.id not in rsa_residues:
            continue
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            res_seq = residue.id[1]
            res_num_str = str(res_seq)
            res_name = residue.get_resname()
            if "CA" in residue:
                ca = residue["CA"]
                coords = list(ca.get_coord())
                all_ca_data[chain.id].append((res_num_str, res_name, res_seq, coords))
                if chain.id in rsa_residues and res_num_str in rsa_residues[chain.id]:
                    rsa_ca_coords[chain.id].append(coords)

    # Find CAs within SCFFTHRESHOLD of any rsa_residue CA
    asa_ca_lines = []
    for chain in model:
        if chain.id in rsa_residues:
            for res_num_str, res_name, res_seq, coords in all_ca_data[chain.id]:
                for rsa_coords in rsa_ca_coords[chain.id]:
                    if distance_calculator(coords, rsa_coords) <= SCFFTHRESHOLD:
                        asa_ca_lines.append((res_name, chain.id, res_seq, coords))
                        break
            
    # Write with sequential serial numbers
    asa_ca_lines = [
        _format_ca_pdb_line(i, res_name, chain_id, res_seq, coords[0], coords[1], coords[2])
        for i, (res_name, chain_id, res_seq, coords) in enumerate(asa_ca_lines, start=1)
    ]
    if len(asa_ca_lines) == 0:
        print(f"No ASA CA lines found for {protein}")
        return False
    
    with open(f"{SURFACE_EXTRACTION_DIR}/{protein}_asa.pdb", "w") as f:
        for line in asa_ca_lines:
            f.write(line)
        f.write("END\n")

def _format_ca_pdb_line(serial, res_name, chain_id, res_seq, x, y, z):
    """Format a CA ATOM line for PDB output."""
    return (
        f"ATOM  {serial:5d}  CA  {res_name:3s} {chain_id}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"
    )

if __name__ == "__main__":
    # extract_surfaces(["1fgnHL", "1fgnL", "1fgnH"])
    extract_surfaces(["1fgj", "1fgnL", "1fgnH"])
    