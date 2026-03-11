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
        url = f"https://files.pdbj.org/pub/pdb/data/structures/all/pdb/pdb{pdb_name}.ent.gz"
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

def pdb_downloader(args):
    receptor_targets = []
    ligand_targets = []

    df = pd.read_csv(args.inputs_csv)

    for receptor, ligand in zip(df["Receptor"], df["Ligand"]):
        receptor_pdb_id = receptor[:4].lower()
        ligand_pdb_id = ligand[:4].lower()
        if not os.path.exists(f"{TARGET_DIR}/{receptor_pdb_id}.pdb") and not download_pdb_file(receptor_pdb_id, TARGET_DIR):
            print(f"Failed to download PDB {receptor_pdb_id}")
            continue
        if not os.path.exists(f"{TARGET_DIR}/{ligand_pdb_id}.pdb") and not download_pdb_file(ligand_pdb_id, TARGET_DIR):
            print(f"Failed to download PDB {ligand_pdb_id}")
            continue

        receptor_targets.append(receptor_pdb_id)
        ligand_targets.append(ligand_pdb_id)

    assert len(receptor_targets) == len(ligand_targets), "Number of receptor and ligand targets must be the same"
    return receptor_targets, ligand_targets
