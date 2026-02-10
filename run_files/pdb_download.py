import os
import subprocess
import gzip
import shutil
import urllib.request

TARGET_DIR = "processed/pdbs"
os.makedirs(TARGET_DIR, exist_ok=True)

def convert_mmcif_to_pdb(mmcif_filepath, pdb_dir):
    try:
        beem_exe = "external_tools/BeEM-master/BeEM"

        ret = subprocess.call([beem_exe, mmcif_filepath], cwd=pdb_dir)
        if ret != 0:
            raise RuntimeError("BeEM failed")

        merge_script = "run_files/merge_bundles.py"
        subprocess.call(["python", merge_script, "*-bundle*.pdb", pdb_dir])
        return True
    except Exception as e:
        print(f"Error converting mmCIF to PDB: {e}")
        return False

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

    try:
        mmcif_gz = f"{pdb_dir}/{pdb_name}.cif.gz"
        mmcif_file = f"{pdb_dir}/{pdb_name}.cif"

        mmcif_url = f"https://files.pdbj.org/pub/pdb/data/structures/all/mmCIF/{pdb_name[:4].lower()}.cif.gz"
        response = urllib.request.urlopen(mmcif_url)

        with open(mmcif_gz, "wb") as fh:
            fh.write(response.read())

        with gzip.open(mmcif_gz, "rb") as f_in, open(mmcif_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(mmcif_gz)

        return convert_mmcif_to_pdb(mmcif_file, final_pdb)

    except Exception as e:
        print(f"Error downloading {pdb_name}: {e}")
        return False

def pdb_downloader():
    targets = []

    # Read pair list
    with open("input_pair_list.txt", "r") as filehnd:
        for line in filehnd:
            line = line.strip().split()
            if len(line) == 2:
                if len(line[0]) == 5 and len(line[1]) == 5:
                    targets.append(line[0])
                    targets.append(line[1])
                else:
                    print(f"Skipping target {line[0]} or {line[1]} because it is not a 5-letter PDB ID."
                          f"We only accept single chain PDB IDs.")
                    continue
    
    # Download and process PDBs
    for target in list(set(targets)):
        if not os.path.exists(f"{TARGET_DIR}/{target[:4].lower()}.pdb"):
            if not download_pdb_file(target[:4].lower(), TARGET_DIR):
                print(f"Failed to download PDB {target}")
                continue
        else:
            print(f"PDB {target} already exists")

    return list(set(targets))
