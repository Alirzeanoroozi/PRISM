#!/usr/bin/env python3
"""
Helpers for frontier protein-DNA model benchmarking in PRISM-prescript.
"""

import csv
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1


DNA_RESIDUES = {
    "DA": "A",
    "DC": "C",
    "DG": "G",
    "DT": "T",
    "DI": "I",
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "I": "I",
}

PROTEIN_KIND = "protein"
DNA_KIND = "dna"


@dataclass
class AdapterAvailability:
    available: bool
    reason: str = ""
    binary_path: str = ""
    module_name: str = ""


def load_registry(registry_path):
    with open(registry_path, "r") as handle:
        return json.load(handle)


def parse_chain_list(value):
    if value in ("", None, "NA"):
        return []
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _is_dna_residue(residue):
    return residue.get_resname().strip().upper() in DNA_RESIDUES


def _protein_one_letter(resname):
    resname = resname.strip().upper()
    try:
        return seq1(resname, undef_code="X")
    except Exception:
        return {
            "MSE": "M",
            "SEC": "U",
            "PYL": "O",
        }.get(resname, "X")


def extract_chain_sequences(pdb_path, chain_ids=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("frontier_sequences", pdb_path)
    protein_entries = []
    dna_entries = []
    for model in structure:
        for chain in model:
            if chain_ids and chain.id not in chain_ids:
                continue
            protein_seq = []
            dna_seq = []
            for residue in chain:
                if is_aa(residue, standard=True):
                    protein_seq.append(_protein_one_letter(residue.get_resname()))
                elif _is_dna_residue(residue):
                    dna_seq.append(DNA_RESIDUES[residue.get_resname().strip().upper()])
            if protein_seq:
                protein_entries.append({"chain_id": chain.id, "kind": PROTEIN_KIND, "sequence": "".join(protein_seq)})
            if dna_seq:
                dna_entries.append({"chain_id": chain.id, "kind": DNA_KIND, "sequence": "".join(dna_seq)})
        break
    return protein_entries, dna_entries


def write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_csv(path, rows, fieldnames=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path, entries):
    lines = []
    for entry in entries:
        lines.append(f">{entry['kind']}|name={entry['chain_id']}")
        lines.append(entry["sequence"])
    write_text(path, "\n".join(lines) + ("\n" if lines else ""))


def write_chai_input(workspace_dir, manifest_row):
    input_dir = Path(workspace_dir) / "frontier_inputs" / "chai1"
    input_dir.mkdir(parents=True, exist_ok=True)
    protein_entries, dna_entries = extract_chain_sequences(manifest_row["protein_unbound_pdb"])
    dna_protein_entries, dna_dna_entries = extract_chain_sequences(manifest_row["dna_unbound_pdb"])
    combined = protein_entries + dna_entries + dna_protein_entries + dna_dna_entries
    fasta_path = input_dir / f"{manifest_row['pair_id']}.fasta"
    write_fasta(fasta_path, combined)
    meta_path = input_dir / f"{manifest_row['pair_id']}.json"
    write_text(meta_path, json.dumps({"pair_id": manifest_row["pair_id"], "entries": combined}, indent=2))
    return {"input_path": str(fasta_path), "metadata_path": str(meta_path), "input_dir": str(input_dir)}


def write_boltz_input(workspace_dir, manifest_row):
    input_dir = Path(workspace_dir) / "frontier_inputs" / "boltz2"
    input_dir.mkdir(parents=True, exist_ok=True)
    protein_entries, _ = extract_chain_sequences(manifest_row["protein_unbound_pdb"])
    _, dna_entries = extract_chain_sequences(manifest_row["dna_unbound_pdb"])
    lines = ["version: 1", "sequences:"]
    for entry in protein_entries:
        lines.extend(
            [
                "  - protein:",
                f"      id: {entry['chain_id']}",
                f"      sequence: {entry['sequence']}",
            ]
        )
    for entry in dna_entries:
        lines.extend(
            [
                "  - dna:",
                f"      id: {entry['chain_id']}",
                f"      sequence: {entry['sequence']}",
            ]
        )
    yaml_path = input_dir / f"{manifest_row['pair_id']}.yaml"
    write_text(yaml_path, "\n".join(lines) + "\n")
    meta_path = input_dir / f"{manifest_row['pair_id']}.json"
    write_text(meta_path, json.dumps({"pair_id": manifest_row["pair_id"], "protein_entries": protein_entries, "dna_entries": dna_entries}, indent=2))
    return {"input_path": str(yaml_path), "metadata_path": str(meta_path), "input_dir": str(input_dir)}


def write_alphafold3_input(workspace_dir, manifest_row):
    input_dir = Path(workspace_dir) / "frontier_inputs" / "alphafold3"
    input_dir.mkdir(parents=True, exist_ok=True)
    protein_entries, _ = extract_chain_sequences(manifest_row["protein_unbound_pdb"])
    _, dna_entries = extract_chain_sequences(manifest_row["dna_unbound_pdb"])
    sequences = []
    for entry in protein_entries:
        sequences.append({"protein": {"id": entry["chain_id"], "sequence": entry["sequence"], "description": f"PRISM {manifest_row['pair_id']} protein {entry['chain_id']}"}})
    for entry in dna_entries:
        sequences.append({"dna": {"id": entry["chain_id"], "sequence": entry["sequence"], "description": f"PRISM {manifest_row['pair_id']} dna {entry['chain_id']}"}})
    payload = {
        "name": manifest_row["pair_id"],
        "modelSeeds": [1],
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 1,
    }
    json_path = input_dir / f"{manifest_row['pair_id']}.json"
    write_text(json_path, json.dumps(payload, indent=2))
    return {"input_path": str(json_path), "input_dir": str(input_dir), "metadata_path": str(input_dir / f"{manifest_row['pair_id']}.json")}


def write_rosettafoldna_inputs(workspace_dir, manifest_row):
    input_dir = Path(workspace_dir) / "frontier_inputs" / "rosettafoldna"
    input_dir.mkdir(parents=True, exist_ok=True)
    protein_entries, _ = extract_chain_sequences(manifest_row["protein_unbound_pdb"])
    _, dna_entries = extract_chain_sequences(manifest_row["dna_unbound_pdb"])
    protein_fastas = []
    dna_fastas = []
    for entry in protein_entries:
        path = input_dir / f"{manifest_row['pair_id']}_{entry['chain_id']}.fa"
        write_fasta(path, [entry])
        protein_fastas.append(path)
    for entry in dna_entries:
        path = input_dir / f"{manifest_row['pair_id']}_{entry['chain_id']}.fa"
        write_fasta(path, [entry])
        dna_fastas.append(path)
    return {"input_dir": str(input_dir), "protein_fastas": [str(path) for path in protein_fastas], "dna_fastas": [str(path) for path in dna_fastas]}


def write_rosettafold_all_atom_input(workspace_dir, manifest_row):
    input_dir = Path(workspace_dir) / "frontier_inputs" / "rosettafold_all_atom"
    input_dir.mkdir(parents=True, exist_ok=True)
    protein_entries, _ = extract_chain_sequences(manifest_row["protein_unbound_pdb"])
    _, dna_entries = extract_chain_sequences(manifest_row["dna_unbound_pdb"])
    protein_inputs = {}
    na_inputs = {}
    for entry in protein_entries:
        fasta_path = input_dir / f"{manifest_row['pair_id']}_{entry['chain_id']}.fa"
        write_fasta(fasta_path, [entry])
        protein_inputs[entry["chain_id"]] = {"fasta_file": str(fasta_path)}
    for entry in dna_entries:
        fasta_path = input_dir / f"{manifest_row['pair_id']}_{entry['chain_id']}.fa"
        write_fasta(fasta_path, [entry])
        na_inputs[entry["chain_id"]] = {"fasta": str(fasta_path), "input_type": "dna"}
    config_lines = [
        "defaults:",
        "  - base",
        "",
        f"job_name: \"{manifest_row['pair_id']}\"",
        "",
        "protein_inputs:",
    ]
    for chain_id, payload in protein_inputs.items():
        config_lines.extend([f"  {chain_id}:", f"    fasta_file: {payload['fasta_file']}"])
    config_lines.append("")
    config_lines.append("na_inputs:")
    for chain_id, payload in na_inputs.items():
        config_lines.extend([f"  {chain_id}:", f"    fasta: {payload['fasta']}", f"    input_type: {payload['input_type']}"])
    config_lines.extend(["", "loader_params:", "  MAXCYCLE: 10", ""])
    config_path = input_dir / f"{manifest_row['pair_id']}.yaml"
    write_text(config_path, "\n".join(config_lines))
    return {"input_path": str(config_path), "input_dir": str(input_dir), "metadata_path": str(input_dir / f"{manifest_row['pair_id']}.json")}


def _python_module_available(python_executable, module_name):
    if not module_name:
        return False
    try:
        proc = subprocess.run(
            [python_executable, "-c", f"import {module_name}"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=30,
        )
    except Exception:
        return False
    return proc.returncode == 0


def _conda_executable():
    return shutil.which("conda") or os.environ.get("CONDA_EXE", "")


def _conda_which(conda_env, candidate):
    conda_exe = _conda_executable()
    if not conda_exe or not conda_env or not candidate:
        return ""
    try:
        proc = subprocess.run(
            [conda_exe, "run", "-n", conda_env, "which", candidate],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=120,
        )
    except Exception:
        return ""
    if proc.returncode != 0:
        return ""
    return proc.stdout.decode("utf-8", errors="ignore").strip()


def resolve_binary(tool, python_executable=sys.executable):
    env_var = tool.get("runner_env")
    if env_var:
        override = os.environ.get(env_var)
        if override:
            override_path = Path(override)
            if override_path.exists():
                return str(override_path)
            found = shutil.which(override)
            if found:
                return found

    conda_env = tool.get("conda_env", "")
    if conda_env:
        for candidate in tool.get("binary_candidates", []):
            found = _conda_which(conda_env, candidate)
            if found:
                return found

    for candidate in tool.get("binary_candidates", []):
        found = shutil.which(candidate)
        if found:
            return found

    return ""


def alphafold3_container_available():
    wrapper_path = Path("/opt/ohpc/pub/apps/alphafold/3.0.1/alphafold")
    sif_path = Path("/opt/ohpc/pub/apps/alphafold/3.0.1/alphafold3.sif")
    return wrapper_path.exists() or sif_path.exists()


def alphafold3_supported_gpu():
    try:
        proc = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=name,compute_cap",
                "--format=csv,noheader",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            timeout=30,
        )
    except Exception as exc:
        return False, f"nvidia-smi_unavailable:{exc}"

    if proc.returncode != 0:
        return False, "nvidia-smi_query_failed"

    lines = [line.strip() for line in proc.stdout.decode("utf-8", errors="ignore").splitlines() if line.strip()]
    if not lines:
        return False, "no_gpu_reported"

    allowed_markers = ("A100", "H100")
    for line in lines:
        if any(marker in line for marker in allowed_markers):
            return True, line
    return False, f"unsupported_gpu_family:{';'.join(lines)}"


def alphafold3_database_root():
    """Return the AF3 database root path, preferring PRISM_AF3_DB_DIR env var."""
    env_root = os.environ.get("PRISM_AF3_DB_DIR")
    if env_root:
        return Path(env_root)
    return Path("/datasets/alphafold3")


def alphafold3_database_available():
    """Check whether the AF3 database directory exists with key files."""
    db_root = alphafold3_database_root()
    if not db_root.is_dir():
        return False, f"missing_af3_database_root:{db_root}"
    marker = db_root / "bfd-first_non_consensus_sequences.fasta"
    if marker.is_file():
        return True, str(db_root)
    return False, f"missing_af3_database_root:{db_root}"


def alphafold3_model_root():
    """Return the AF3 model weights path, preferring PRISM_AF3_MODEL_DIR env var."""
    env_root = os.environ.get("PRISM_AF3_MODEL_DIR")
    if env_root:
        return Path(env_root)
    return Path("/datasets/alphafold3/models")


def boltz_runtime_env(binary_path):
    env = os.environ.copy()
    if not binary_path:
        return env
    binary_root = Path(binary_path).resolve().parent.parent
    site_packages = binary_root / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}" / "site-packages"
    lib_dirs = [
        site_packages / "nvidia" / "cu13" / "lib",
        site_packages / "nvidia" / "cuda_nvrtc" / "lib",
        site_packages / "nvidia" / "cudnn" / "lib",
        site_packages / "nvidia" / "nvshmem" / "lib",
    ]
    existing = env.get("LD_LIBRARY_PATH", "")
    merged = [str(path) for path in lib_dirs if path.exists()]
    if existing:
        merged.append(existing)
    if merged:
        env["LD_LIBRARY_PATH"] = ":".join(merged)
    return env


def probe_tool(tool, python_executable=sys.executable):
    if tool.get("tool_id") == "alphafold3":
        db_available, db_reason = alphafold3_database_available()
        if not db_available:
            return AdapterAvailability(available=False, reason=db_reason)
        if alphafold3_container_available():
            wrapper_path = Path("/opt/ohpc/pub/apps/alphafold/3.0.1/alphafold")
            sif_path = Path("/opt/ohpc/pub/apps/alphafold/3.0.1/alphafold3.sif")
            if wrapper_path.exists():
                return AdapterAvailability(available=True, binary_path=str(wrapper_path))
            if sif_path.exists():
                return AdapterAvailability(available=True, binary_path=str(sif_path))

    binary_path = resolve_binary(tool, python_executable=python_executable)
    module_name = tool.get("python_module", "")
    if binary_path:
        if module_name and _python_module_available(python_executable, module_name):
            return AdapterAvailability(available=True, binary_path=binary_path)
        if not module_name:
            return AdapterAvailability(available=True, binary_path=binary_path)
    if module_name and _python_module_available(python_executable, module_name):
        if tool.get("tool_id") == "rosettafold_all_atom":
            return AdapterAvailability(available=True, module_name=module_name)
        if tool.get("integration_kind") != "local_cli":
            return AdapterAvailability(available=True, module_name=module_name)
    if tool.get("integration_kind") == "reference_local":
        return AdapterAvailability(available=False, reason="tool_not_installed_or_missing_runner")
    return AdapterAvailability(available=False, reason="tool_not_installed")


def model_counts_by_kind(manifest_row):
    protein_entries, _ = extract_chain_sequences(manifest_row["protein_unbound_pdb"])
    _, dna_entries = extract_chain_sequences(manifest_row["dna_unbound_pdb"])
    return len(protein_entries), len(dna_entries)


def _classify_chain(chain):
    aa_count = 0
    dna_count = 0
    for residue in chain:
        if is_aa(residue, standard=True):
            aa_count += 1
        elif _is_dna_residue(residue):
            dna_count += 1
    if aa_count >= dna_count and aa_count > 0:
        return PROTEIN_KIND
    if dna_count > 0:
        return DNA_KIND
    return "unknown"


def normalize_structure_chains(model_path, output_path, protein_chain_ids, dna_chain_ids):
    parser = PDBParser(QUIET=True)
    model_path = Path(model_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = model_path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("frontier_model", str(model_path))
    model = next(structure.get_models())
    chains = list(model.get_chains())
    chain_kinds = [_classify_chain(chain) for chain in chains]

    protein_queue = list(protein_chain_ids)
    dna_queue = list(dna_chain_ids)
    assignments = {}
    for chain, kind in zip(chains, chain_kinds):
        if kind == PROTEIN_KIND and protein_queue:
            assignments[chain.id] = protein_queue.pop(0)
        elif kind == DNA_KIND and dna_queue:
            assignments[chain.id] = dna_queue.pop(0)
        elif protein_queue:
            assignments[chain.id] = protein_queue.pop(0)
        elif dna_queue:
            assignments[chain.id] = dna_queue.pop(0)
        else:
            assignments[chain.id] = chain.id

    original_ids = [chain.id for chain in chains]
    temp_ids = {chain_id: f"TMP{i:02d}" for i, chain_id in enumerate(original_ids)}
    for chain in chains:
        chain.id = temp_ids[chain.id]
    for chain in chains:
        for original_id, temp_id in temp_ids.items():
            if chain.id == temp_id:
                chain.id = assignments.get(original_id, original_id)
                break

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))
    return str(output_path)


def normalize_prediction_file(path, output_path, manifest_row):
    """Normalize prediction chain names to match manifest chain IDs.

    Wraps normalize_structure_chains with a manifest-row interface.
    """
    return normalize_structure_chains(
        path, output_path,
        manifest_row.get("protein_chain_ids", ""),
        manifest_row.get("dna_chain_ids", ""),
    )


def candidate_model_files(output_dir):
    output_dir = Path(output_dir)
    candidates = []
    for pattern in ("**/*.pdb", "**/*.cif", "**/*.mmcif"):
        candidates.extend(output_dir.glob(pattern))
    filtered = []
    for path in candidates:
        if path.is_file():
            filtered.append(path)
    priority = []
    for path in filtered:
        name = path.name.lower()
        score = 100
        if "ranked_0" in name or "model_00" in name or "best" in name:
            score = 0
        elif "model_0" in name or "seed_0" in name:
            score = 1
        elif path.suffix.lower() == ".pdb":
            score = 2
        priority.append((score, str(path)))
    return [Path(path) for _, path in sorted(priority, key=lambda item: (item[0], item[1]))]


def normalize_prediction_file(model_path, normalized_dir, manifest_row):
    normalized_dir = Path(normalized_dir)
    normalized_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{manifest_row['pair_id']}_{Path(model_path).stem}"
    normalized_path = normalized_dir / f"{stem}.pdb"
    return normalize_structure_chains(
        model_path,
        normalized_path,
        protein_chain_ids=parse_chain_list(manifest_row["protein_chain_ids"]),
        dna_chain_ids=parse_chain_list(manifest_row["dna_chain_ids"]),
    )
