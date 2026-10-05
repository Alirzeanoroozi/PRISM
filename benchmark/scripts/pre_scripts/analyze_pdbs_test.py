import os
import pandas as pd
from .utils import STANDARD_AA
from Bio.PDB import PDBParser
from .pdb_download import download_pdb_file

TEMPLATES_PDBS_DIR = "benchmark/benchmark5.5/structures"

def analyse_pdb(filepath):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", filepath)

    num_models = len(structure)

    residues_seen = set()
    chain_set = set()
    first_model = next(structure.get_models())
    for chain in first_model:
        chain_id = chain.id
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname in STANDARD_AA:
                resseq = residue.id[1]
                residues_seen.add((chain_id, resseq))
                chain_set.add(chain_id)

    chain_names = ",".join(sorted(chain_set))
    num_chains = len(chain_set)
    total_amino_acids = len(residues_seen)

    # Get EXPDTA from file manually (since Biopython does not parse it)
    expdta = ""
    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("EXPDTA"):
                    expdta = line[10:80].strip()
                    break
    except Exception:
        pass

    # Structure type: NMR vs single
    exp_upper = expdta.upper()
    if "NMR" in exp_upper or "NUCLEAR MAGNETIC" in exp_upper or num_models > 1:
        structure_type = "NMR"
    else:
        structure_type = "single"

    return {
        "pdb_id": filepath.split("/")[-1].split(".")[0].lower(),
        "num_chains": num_chains,
        "chain_names": chain_names,
        "total_amino_acids": total_amino_acids,
        "structure_type": structure_type,
        "expdta": expdta,
        "num_models": num_models,
    }

def extract_pdb_id(filename):
    """
    Extract PDB ID from filenames like: 1A2K_I_b.pdb, 1A2K_l_U.pdb, etc.
    Returns the first 4 characters (e.g., '1A2K')
    """
    base_name = filename.split(".")[0]  # Remove .pdb extension
    pdb_id = base_name.split("_")[0].lower()  # Get first part before underscore
    return pdb_id

def run_analysis():
    # Get all .pdb files from the benchmark folder
    if not os.path.exists(TEMPLATES_PDBS_DIR):
        print(f"Directory {TEMPLATES_PDBS_DIR} does not exist!")
        return [], 0, 0
    
    pdb_files = sorted([f for f in os.listdir(TEMPLATES_PDBS_DIR) if f.endswith(".pdb")])
    print(f"Found {len(pdb_files)} PDB files in {TEMPLATES_PDBS_DIR}")

    results = []
    failed_templates = []
    for pdb_file in pdb_files:
        try:
            filepath = os.path.join(TEMPLATES_PDBS_DIR, pdb_file)
            info = analyse_pdb(filepath)
            # Override pdb_id with properly extracted ID
            info["pdb_id"] = extract_pdb_id(pdb_file)
            info["filename"] = pdb_file
            results.append(info)
        except Exception as e:
            print(f"Failed to analyse {pdb_file}: {e}")
            failed_templates.append(pdb_file)

    df = pd.DataFrame(results)
    df.to_csv("benchmark/benchmark5.5/benchmark_analysis.csv", index=False)
    filtered = df[(df["num_chains"] > 1) & (df["structure_type"] == "single")]
    filtered.to_csv("benchmark/benchmark5.5/templates_filtered.csv", index=False)

    templates_updatex_path = "processed/templates.txt"
    os.makedirs("processed", exist_ok=True)
    with open(templates_updatex_path, "w", encoding="utf-8") as f:
        for pdb_id in filtered["pdb_id"].unique():
            f.write(f"{pdb_id}\n")

    missing_templates_path = "processed/missing_templates.txt"
    with open(missing_templates_path, "w", encoding="utf-8") as f:
        for pdb_file in failed_templates:
            f.write(f"{pdb_file}\n")

    return results, len(filtered), len(failed_templates)

if __name__ == "__main__":
    results, filtered_count, failed_count = run_analysis()
    print(f"Analysed {len(results)} template files.")
    print(f"Filtered {filtered_count} templates")
    print(f"Failed {failed_count} templates")