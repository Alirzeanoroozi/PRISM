import os
import gzip
import shutil
import urllib.request
import pandas as pd

TARGET_DIR = "processed/pdbs"
os.makedirs(TARGET_DIR, exist_ok=True)

def download_pdb_file(pdb_name, pdb_dir):
    final_pdb = f"{pdb_dir}/{pdb_name}.pdb"
    gz_file = f"{pdb_dir}/{pdb_name}.ent.gz"

    try:
        url = f"https://files.pdbj.org/pub/pdb/data/structures/all/pdb/pdb{pdb_name[:4].lower()}.ent.gz"
        response = urllib.request.urlopen(url)

        with open(gz_file, "wb") as fh:
            fh.write(response.read())

        with gzip.open(gz_file, "rb") as f_in, open(final_pdb, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(gz_file)

        return True
    except Exception as e:
        print(f"PDB download failed for {pdb_name}: {e}")
        return False

def pdb_downloader():
    receptor_targets = []
    ligand_targets = []

    # Read pair list
    df = pd.read_csv("inputs.csv")
    for _, row in df.iterrows():
        if len(row["Receptor"]) == 5 and len(row["Ligand"]) == 5:
            receptor_targets.append(row["Receptor"])
            ligand_targets.append(row["Ligand"])
        else:
            print(f"Skipping target {row['Receptor']} or {row['Ligand']} because it is not a 5-letter PDB ID. We only accept single chain PDB IDs.")
            continue
    
    # Download and process PDBs
    for target in list(set(receptor_targets + ligand_targets)):
        if not os.path.exists(f"{TARGET_DIR}/{target[:4].lower()}.pdb"):
            if not download_pdb_file(target[:4].lower(), TARGET_DIR):
                print(f"Failed to download PDB {target}")
                continue
        else:
            print(f"PDB {target} already exists")

    return list(set(receptor_targets)), list(set(ligand_targets))
