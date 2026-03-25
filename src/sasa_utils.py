import freesasa
from Bio.PDB import PDBParser

def get_asa_complex(pdb, save_directory: str):
    chains = [c.upper() for c in pdb[4:]]
    pdb_path = f"{save_directory}/pdbs/{pdb[:4].lower()}.pdb"
    try:
        parser = PDBParser()
        structure = parser.get_structure("target", pdb_path)
        result, _ = freesasa.calcBioPDB(structure)
        residue_areas = result.residueAreas()
    except Exception as e:
        raise RuntimeError(f"FreeSASA ASA calculation failed for {pdb_path}: {e}")

    relative_asa = {}
    for chain_id, chain_dict in residue_areas.items():
        if chain_id not in chains:
            continue
        relative_asa[chain_id] = {res_num: ra.relativeTotal * 100.0 for res_num, ra in chain_dict.items()}
    return relative_asa
