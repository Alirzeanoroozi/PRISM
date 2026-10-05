import os
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

from .naccess_utils import get_asa_complex_target
from .utils import distance_calculator

RSATHRESHOLD = 15.0
SCFFTHRESHOLD = 5.0

SURFACE_EXTRACTION_DIR = "processed/surface_extraction"
os.makedirs(SURFACE_EXTRACTION_DIR, exist_ok=True)
SURFACE_FAILURE_LOG = os.path.join(SURFACE_EXTRACTION_DIR, "failures.tsv")

def _resolve_scaffold_threshold(value=None):
    if value is not None:
        return float(value)
    return float(os.environ.get("PRISM_SCFF_THRESHOLD", SCFFTHRESHOLD))


def extract_surfaces(queries, scaffold_threshold=None):
    scaffold_threshold = _resolve_scaffold_threshold(scaffold_threshold)
    for protein in queries:
        # Skip chains that already have a valid (non-placeholder) surface PDB.
        # Surface extraction is shared across every pipeline run on a workspace;
        # without this guard, concurrent jobs race on NACCESS/FreeSASA fixed
        # filenames in the same directory. A 4-byte "END\n" file is only ever the
        # recorded placeholder for a failed/empty extraction, so it is re-tried.
        # Use the canonical (lowercase-pdb + uppercase-chains) name so the guard
        # matches the files both the precompute step and the aligner use, even
        # when the inputs.csv spelling differs in case.
        from .pdb_download import normalize_target_id
        try:
            canonical_protein = normalize_target_id(protein)
        except Exception:
            canonical_protein = protein
        asa_path = os.path.join(SURFACE_EXTRACTION_DIR, f"{canonical_protein}.asa.pdb")
        if os.path.exists(asa_path) and os.path.getsize(asa_path) > 4:
            print(f"Surface already present for {canonical_protein}, skipping extraction.")
            continue
        try:
            extract_surface(protein, scaffold_threshold=scaffold_threshold)
        except Exception as exc:
            # Preserve the pair in the batch while making the stage failure
            # explicit for later comparison/reporting.
            os.makedirs(SURFACE_EXTRACTION_DIR, exist_ok=True)
            with open(f"{SURFACE_EXTRACTION_DIR}/{canonical_protein}.asa.pdb", "w") as handle:
                handle.write("END\n")
            with open(SURFACE_FAILURE_LOG, "a") as handle:
                handle.write(f"{canonical_protein}\t{type(exc).__name__}: {exc}\n")

def extract_surface(protein, scaffold_threshold=None):
    scaffold_threshold = _resolve_scaffold_threshold(scaffold_threshold)
    print(f"Extracting surface for {protein}...")
    asa_complex = get_asa_complex_target(protein, SURFACE_EXTRACTION_DIR)
    rsa_residues = [key for key in asa_complex if asa_complex[key] > RSATHRESHOLD]
    if len(rsa_residues) == 0:
        # Keep the downstream alignment contract even when a structure has no
        # exposed residues above the RSA threshold.  An empty CA-only PDB is a
        # data outcome that can be recorded, whereas a missing file aborts the
        # entire multi-pair batch during parser startup.
        os.makedirs(SURFACE_EXTRACTION_DIR, exist_ok=True)
        with open(f"{SURFACE_EXTRACTION_DIR}/{protein}.asa.pdb", "w") as handle:
            handle.write("END\n")
        return {}
    # Materialized single-chain PDBs use the canonical (lowercase-pdb +
    # uppercase-chains) name, e.g. 2CV5A -> processed/pdbs/2cv5A.pdb.
    from .pdb_download import normalize_target_id
    canonical = normalize_target_id(protein)
    chain_pdb_path = f"processed/pdbs/{canonical}.pdb"
    pdb_path = chain_pdb_path if os.path.exists(chain_pdb_path) else f"processed/pdbs/{protein[:4].lower()}.pdb"
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
            if distance_calculator(coords, rsa_coords) <= scaffold_threshold:
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
