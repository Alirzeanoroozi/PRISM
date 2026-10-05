import os
import json
import subprocess
import sys
from threading import Lock

from .utils import standard_data
from .pdb_download import normalize_target_id, target_chain_ids

# Local change: serialize NACCESS because it writes fixed filenames in cwd.
_NACCESS_LOCK = Lock()


def _surface_backend():
    backend = os.environ.get("PRISM_SURFACE_BACKEND", "naccess").strip().lower()
    if backend not in {"naccess", "freesasa"}:
        raise ValueError(
            "PRISM_SURFACE_BACKEND must be 'naccess' or 'freesasa', "
            f"got {backend!r}"
        )
    return backend

def _split_chain_resnum(chain_token):
    """Split a NACCESS RES chain column into (chain, residue_number).

    For ordinary structures the column is just the chain (e.g. ``A``), but for
    structures with 4-digit residue numbers (>=1000) NACCESS fuses chain and
    number into a single token (e.g. ``D1230``). Return ``(chain, resnum)``
    where resnum is the trailing digits (or '' when none are fused).
    """
    import re as _re
    m = _re.match(r"^([A-Za-z]+)(\d*)$", (chain_token or ""))
    if not m:
        return chain_token, ""
    chain = m.group(1)
    resnum = m.group(2)
    return chain, resnum


def get_asa_complex(template, save_directory):
    if _surface_backend() == "freesasa":
        areas = run_freesasa(template, save_directory)
        return _relative_areas(areas, target_chain_ids(template))

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
                    chain, fused_resnum = _split_chain_resnum(items[2])
                    if chain not in (chain_id1, chain_id2):
                        continue
                    residue_number = fused_resnum or items[3]
                    absolute_asa = float(items[4])
                    relative_asa = absolute_asa * 100 / standard
                    relative_asa_complex[f"{residue_name}_{residue_number}_{chain}"] = relative_asa    
    
    return relative_asa_complex

def get_asa_complex_target(template, save_directory):
    if _surface_backend() == "freesasa":
        areas = run_freesasa(template, save_directory, is_target=True)
        chains = target_chain_ids(template)
        return _relative_areas(areas, chains)

    run_naccess(template, save_directory, is_target=True)
    chain_ids = set(target_chain_ids(template))
    
    relative_asa_complex = {}    
    with open(f"{save_directory}/{template}.rsa", 'r') as f:
        for rsaline in f.readlines():
            if rsaline.startswith("RES"):
                items = rsaline.split()
                residue_name = items[1]
                standard = standard_data(residue_name)
                if standard != -1:
                    chain, fused_resnum = _split_chain_resnum(items[2])
                    if chain_ids and chain not in chain_ids:
                        continue
                    residue_number = fused_resnum or items[3]
                    absolute_asa = float(items[4])
                    relative_asa = absolute_asa * 100 / standard
                    relative_asa_complex[f"{residue_name}_{residue_number}_{chain}"] = relative_asa    
    
    return relative_asa_complex


def _relative_areas(areas, allowed_chains):
    relative_asa = {}
    for chain, residues in areas.items():
        # Empty allowed_chains means "no chain filtering" (a bare PDB ID like
        # ``5zng`` from the CSV). This mirrors the NACCESS path which only
        # filters when chain_ids is non-empty.
        if allowed_chains and chain not in allowed_chains:
            continue
        for residue_number, residue in residues.items():
            residue_name = residue["residue_name"]
            standard = standard_data(residue_name)
            if standard == -1:
                continue
            absolute_asa = float(residue["total"])
            relative_asa[f"{residue_name}_{residue_number}_{chain}"] = (
                absolute_asa * 100 / standard
            )
    return relative_asa


def run_freesasa(template, save_directory, is_target=False):
    if is_target:
        target_id = normalize_target_id(template)
        chain_path = f"processed/pdbs/{target_id}.pdb"
        pdb_path = chain_path if os.path.exists(chain_path) else f"processed/pdbs/{target_id[:4]}.pdb"
    else:
        pdb_path = f"templates/pdbs/{template[:4].lower()}.pdb"

    os.makedirs(save_directory, exist_ok=True)
    output_path = os.path.join(save_directory, f"{template}.freesasa.json")
    python_bin = os.environ.get("PRISM_FREESASA_PYTHON", sys.executable)
    # freesasa_runner lives under the `src` package of the PRISM repo. The repo
    # root must be on PYTHONPATH so `python -m src.freesasa_runner` resolves even
    # when the pipeline runs from a different working directory (e.g. a batch
    # RUN_ROOT workspace).
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    _env = dict(os.environ)
    _existing = _env.get("PYTHONPATH", "")
    _env["PYTHONPATH"] = repo_root + (os.pathsep + _existing if _existing else "")
    result = subprocess.run(
        [python_bin, "-m", "src.freesasa_runner", pdb_path, output_path],
        capture_output=True,
        text=True,
        check=False,
        env=_env,
    )
    if result.returncode != 0 or not os.path.exists(output_path):
        raise RuntimeError(
            f"FreeSASA failed for {template} (pdb: {pdb_path}, "
            f"exit code: {result.returncode}). "
            f"{result.stderr.strip() or result.stdout.strip() or 'No output captured.'}"
        )
    with open(output_path, "r") as handle:
        return json.load(handle)

def run_naccess(template, save_directory, is_target=False):
    if is_target:
        target_id = normalize_target_id(template)
        chain_path = f"processed/pdbs/{target_id}.pdb"
        pdb_path = chain_path if os.path.exists(chain_path) else f"processed/pdbs/{target_id[:4]}.pdb"
    else:
        pdb_path = f"templates/pdbs/{template[:4].lower()}.pdb"
    # Local change: ensure target directory exists before moving NACCESS outputs.
    os.makedirs(save_directory, exist_ok=True)

    # NACCESS names outputs after the actual input filename. Chain-qualified
    # targets therefore produce e.g. ``4dn3HL.rsa`` rather than ``4dn3.rsa``.
    pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
    rsa_src = f"{pdb_id}.rsa"
    asa_src = f"{pdb_id}.asa"
    log_src = f"{pdb_id}.log"
    out_path = f".naccess_{template}.out"

    # Local change: use subprocess + lock for safer error handling and threading.
    # NACCESS uses fixed filenames in the cwd (e.g. accall.input), so concurrent
    # runs from template threads can clobber each other unless serialized.
    with _NACCESS_LOCK:
        with open(out_path, "w") as out_file:
            naccess_bin = os.environ.get(
                "PRISM_NACCESS_EXECUTABLE", "external_tools/naccess/naccess"
            )
            result = subprocess.run(
                [naccess_bin, pdb_path],
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
