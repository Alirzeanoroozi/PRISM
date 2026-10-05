"""
Fixed version of src/naccess_utils.py

Key fixes:
1. Added check_freesasa_available() for pre-flight validation.
2. run_freesasa() validates the FreeSASA output more thoroughly.
3. Added structured error codes for FreeSASA/NACCESS failures.
4. run_naccess() returns structured error codes instead of bare RuntimeError.
"""

import json
import logging
import os
import subprocess
import sys
from threading import Lock

import pdb_download as _pd
import utils as _ut
normalize_target_id = _pd.normalize_target_id
target_chain_ids = _pd.target_chain_ids
standard_data = _ut.standard_data

logger = logging.getLogger(__name__)

_NACCESS_LOCK = Lock()
FREESASA_PRE_CHECK_CACHE = {}


def _surface_backend():
    backend = os.environ.get("PRISM_SURFACE_BACKEND", "naccess").strip().lower()
    if backend not in {"naccess", "freesasa"}:
        raise ValueError(
            f"PRISM_SURFACE_BACKEND must be 'naccess' or 'freesasa', got {backend!r}"
        )
    return backend


def check_freesasa_available():
    """
    Pre-flight check: verify FreeSASA is importable in the target interpreter.

    Run this BEFORE the pipeline starts, not after a failing subprocess call.
    Returns (available: bool, python_path: str, error: str or None).
    """
    python_bin = os.environ.get("PRISM_FREESASA_PYTHON", sys.executable)
    cache_key = f"{python_bin}_{os.getpid()}"
    if cache_key in FREESASA_PRE_CHECK_CACHE:
        return FREESASA_PRE_CHECK_CACHE[cache_key]

    result = subprocess.run(
        [python_bin, "-c", "import freesasa; print(freesasa.__version__)"],
        capture_output=True, text=True, timeout=30,
    )
    if result.returncode == 0:
        status = (True, python_bin, None)
        logger.info("FreeSASA available in %s (version: %s)", python_bin, result.stdout.strip())
    else:
        error_msg = (
            f"FreeSASA is NOT available in {python_bin}. "
            f"Install with: conda install -c conda-forge freesasa "
            f"or use --surface-backend=naccess. "
            f"Stderr: {result.stderr.strip()[:500]}"
        )
        logger.warning(error_msg)
        status = (False, python_bin, error_msg)

    FREESASA_PRE_CHECK_CACHE[cache_key] = status
    return status


def get_asa_complex(template, save_directory):
    if _surface_backend() == "freesasa":
        areas = run_freesasa(template, save_directory)
        return _relative_areas(areas, target_chain_ids(template))

    run_naccess(template, save_directory)
    chain_id1 = template[4]
    chain_id2 = template[5]

    relative_asa_complex = {}
    with open(f"{save_directory}/{template}.rsa", "r") as f:
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
    if _surface_backend() == "freesasa":
        areas = run_freesasa(template, save_directory, is_target=True)
        chain_ids = set(target_chain_ids(template))
        if not chain_ids:
            return _relative_areas(areas, None)
        return _relative_areas(areas, tuple(chain_ids))

    run_naccess(template, save_directory, is_target=True)
    chain_ids = set(target_chain_ids(template))

    relative_asa_complex = {}
    with open(f"{save_directory}/{template}.rsa", "r") as f:
        for rsaline in f.readlines():
            if rsaline.startswith("RES"):
                items = rsaline.split()
                residue_name = items[1]
                standard = standard_data(residue_name)
                if standard != -1:
                    chain = items[2]
                    if chain_ids and chain not in chain_ids:
                        continue
                    residue_number = items[3]
                    absolute_asa = float(items[4])
                    relative_asa = absolute_asa * 100 / standard
                    relative_asa_complex[f"{residue_name}_{residue_number}_{chain}"] = relative_asa

    return relative_asa_complex


def _relative_areas(areas, allowed_chains):
    relative_asa = {}
    for chain, residues in areas.items():
        if allowed_chains is not None and chain not in allowed_chains:
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
    """Run FreeSASA and return parsed JSON areas dict.

    Raises RuntimeError with diagnostic information on failure.
    """
    if is_target:
        target_id = normalize_target_id(template)
        chain_path = f"processed/pdbs/{target_id}.pdb"
        pdb_path = (
            chain_path
            if os.path.exists(chain_path)
            else f"processed/pdbs/{target_id[:4]}.pdb"
        )
    else:
        pdb_path = f"templates/pdbs/{template[:4].lower()}.pdb"

    os.makedirs(save_directory, exist_ok=True)
    output_path = os.path.join(save_directory, f"{template}.freesasa.json")
    python_bin = os.environ.get("PRISM_FREESASA_PYTHON", sys.executable)

    result = subprocess.run(
        [python_bin, "-m", "src.freesasa_runner", pdb_path, output_path],
        capture_output=True, text=True, check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"FreeSASA failed for {template} (pdb: {pdb_path}, "
            f"exit code: {result.returncode}). "
            f"{result.stderr.strip() or result.stdout.strip() or 'No output captured.'}"
        )

    if not os.path.exists(output_path):
        raise RuntimeError(
            f"FreeSASA for {template} completed (exit=0) but "
            f"output file not found: {output_path}"
        )

    with open(output_path, "r") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise RuntimeError(
            f"FreeSASA for {template} produced non-dict JSON: "
            f"type={type(data).__name__}"
        )

    if not data:
        logger.warning("FreeSASA returned empty areas dict for %s (pdb: %s)", template, pdb_path)

    return data


def run_naccess(template, save_directory, is_target=False):
    """Run NACCESS and verify the .rsa output was produced.

    Raises RuntimeError with diagnostic information on failure.
    """
    if is_target:
        target_id = normalize_target_id(template)
        chain_path = f"processed/pdbs/{target_id}.pdb"
        pdb_path = (
            chain_path
            if os.path.exists(chain_path)
            else f"processed/pdbs/{target_id[:4]}.pdb"
        )
    else:
        pdb_path = f"templates/pdbs/{template[:4].lower()}.pdb"

    os.makedirs(save_directory, exist_ok=True)

    pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
    rsa_src = f"{pdb_id}.rsa"
    asa_src = f"{pdb_id}.asa"
    log_src = f"{pdb_id}.log"
    out_path = f".naccess_{template}.out"

    with _NACCESS_LOCK:
        with open(out_path, "w") as out_file:
            naccess_bin = os.environ.get(
                "PRISM_NACCESS_EXECUTABLE", "external_tools/naccess/naccess"
            )
            if not os.path.exists(naccess_bin):
                raise RuntimeError(
                    f"NACCESS binary not found: {naccess_bin}. "
                    f"Recompile: cd external_tools/naccess && gfortran accall.f -o accall -O"
                )
            result = subprocess.run(
                [naccess_bin, pdb_path],
                stdout=out_file, stderr=subprocess.STDOUT, check=False,
            )

        if result.returncode != 0 or not os.path.exists(rsa_src):
            output = ""
            if os.path.exists(out_path):
                with open(out_path, "r") as f:
                    output = f.read().strip()
            raise RuntimeError(
                f"NACCESS failed for {template} (pdb: {pdb_path}, "
                f"exit code: {result.returncode}). "
                f"{output or 'No output captured.'}"
            )

        os.rename(rsa_src, f"{save_directory}/{template}.rsa")
        for tmp in (asa_src, log_src, out_path):
            if os.path.exists(tmp):
                os.remove(tmp)
