import os
import subprocess
from threading import Lock

from .utils import standard_data

# Local change: serialize NACCESS because it writes fixed filenames in cwd.
_NACCESS_LOCK = Lock()

def get_asa_complex(template, save_directory):
    run_naccess(template, save_directory)
    chain_id1 = template[4]
    chain_id2 = template[5]
    
    relative_asa_complex = {}    
    with open(f"{save_directory}/{template}.rsa", 'r') as f:
        for rsaline in f.readlines():
            if rsaline.startswith("RES"):
                items = rsaline.split()
                residue_name = items[1]
                standard = standard_data(residue_name)
                if standard != -1:
                    chain = items[2]
                    if chain not in (chain_id1, chain_id2):
                        continue
                    residue_number = items[3]
                    absolute_asa = float(items[4])
                    relative_asa = absolute_asa * 100 / standard
                    relative_asa_complex[f"{residue_name}_{residue_number}_{chain}"] = relative_asa    
    
    return relative_asa_complex

def get_asa_complex_target(template, save_directory):
    run_naccess(template, save_directory, is_target=True)
    chain_id = template[4]
    
    relative_asa_complex = {}    
    with open(f"{save_directory}/{template}.rsa", 'r') as f:
        for rsaline in f.readlines():
            if rsaline.startswith("RES"):
                items = rsaline.split()
                residue_name = items[1]
                standard = standard_data(residue_name)
                if standard != -1:
                    chain = items[2]
                    if chain != chain_id:
                        continue
                    residue_number = items[3]
                    absolute_asa = float(items[4])
                    relative_asa = absolute_asa * 100 / standard
                    relative_asa_complex[f"{residue_name}_{residue_number}_{chain}"] = relative_asa    
    
    return relative_asa_complex

def run_naccess(template, save_directory, is_target=False):
    if is_target:
        pdb_path = f"processed/pdbs/{template[:4].lower()}.pdb"
    else:
        pdb_path = f"templates/pdbs/{template[:4].lower()}.pdb"
    # Local change: ensure target directory exists before moving NACCESS outputs.
    os.makedirs(save_directory, exist_ok=True)

    pdb_id = template[:4].lower()
    rsa_src = f"{pdb_id}.rsa"
    asa_src = f"{pdb_id}.asa"
    log_src = f"{pdb_id}.log"
    out_path = f".naccess_{template}.out"

    # Local change: use subprocess + lock for safer error handling and threading.
    # NACCESS uses fixed filenames in the cwd (e.g. accall.input), so concurrent
    # runs from template threads can clobber each other unless serialized.
    with _NACCESS_LOCK:
        with open(out_path, "w") as out_file:
            result = subprocess.run(
                ["external_tools/naccess/naccess", pdb_path],
                stdout=out_file,
                stderr=subprocess.STDOUT,
                check=False,
            )

        if result.returncode != 0 or not os.path.exists(rsa_src):
            output = ""
            if os.path.exists(out_path):
                with open(out_path, "r") as f:
                    output = f.read().strip()
            raise RuntimeError(
                f"NACCESS failed for {template} (pdb: {pdb_path}, exit code: {result.returncode}). "
                f"{output or 'No output captured.'}"
            )

        os.rename(rsa_src, f"{save_directory}/{template}.rsa")
        for tmp in (asa_src, log_src, out_path):
            if os.path.exists(tmp):
                os.remove(tmp)
