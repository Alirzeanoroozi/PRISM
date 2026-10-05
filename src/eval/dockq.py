"""
Calculate DockQ between a model (predicted) PDB and a native (reference) PDB.

DockQ is a continuous quality measure (0–1) for protein–protein docking models,
combining fnat (fraction native contacts), LRMS (ligand RMSD), and iRMS (interface RMSD).
CAPRI interpretation: <0.23 Incorrect, 0.23–0.49 Acceptable, 0.49–0.80 Medium, ≥0.80 High.

Requires: pip install DockQ
Ref: https://github.com/wallnerlab/DockQ
"""
import json
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import uuid
from pathlib import Path


def resolve_dockq_executable():
    override = os.environ.get("DOCKQ_BIN")
    if override:
        return [override]
    override_python = os.environ.get("DOCKQ_PYTHON")
    if override_python:
        # Preserve the virtual-environment entry path exactly. Resolving its
        # symlink can discard the environment-local DockQ installation.
        return [override_python, "-m", "DockQ"]
    try:
        import DockQ  # noqa: F401
        return [sys.executable, "-m", "DockQ"]
    except Exception:
        pass
    env_candidate = Path(sys.executable).resolve().parent / "DockQ"
    if env_candidate.is_file():
        return [str(env_candidate)]
    found = shutil.which("DockQ")
    if found:
        return [found]
    return None


def _pdb_residues(pdb_path, chain_id):
    """Return ordered standard residue identities for one PDB chain."""
    residues = []
    with open(pdb_path, encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM") or len(line) < 27 or line[21].strip() != chain_id:
                continue
            key = (line[22:26].strip(), line[26].strip())
            if not residues or residues[-1][:2] != key:
                residues.append((*key, line[17:20].strip()))
    return residues


def validate_no_align_mapping(model_pdb, native_pdb, mapping):
    """Require a bijective chain map with matching residue IDs and names."""
    try:
        model_chains, native_chains = str(mapping).split(":", 1)
    except ValueError as exc:
        raise ValueError(f"invalid DockQ mapping: {mapping!r}") from exc
    if not model_chains or not native_chains or len(model_chains) != len(native_chains):
        raise ValueError(f"non-bijective DockQ mapping: {mapping!r}")
    if len(set(model_chains)) != len(model_chains) or len(set(native_chains)) != len(native_chains):
        raise ValueError(f"duplicate chain in DockQ mapping: {mapping!r}")
    for model_chain, native_chain in zip(model_chains, native_chains):
        model_residues = _pdb_residues(model_pdb, model_chain)
        native_residues = _pdb_residues(native_pdb, native_chain)
        if not model_residues or not native_residues:
            raise ValueError(
                f"mapped chain is absent from model/native PDB: {model_chain}:{native_chain}"
            )
        if len(model_residues) != len(native_residues):
            raise ValueError(f"residue correspondence differs for no-align mapping {model_chain}:{native_chain}")
        for model_residue, native_residue in zip(model_residues, native_residues):
            if model_residue[:2] != native_residue[:2]:
                raise ValueError(f"residue numbering differs for no-align mapping {model_chain}:{native_chain}")
            if model_residue[2] != native_residue[2]:
                raise ValueError(f"residue identity differs for no-align mapping {model_chain}:{native_chain}")


def calculate_dockq(
    model_pdb,
    native_pdb,
    mapping=None,
    work_dir=None,
    no_align=False,
    n_cpu=1,
):
    """
    Run DockQ to compare a model PDB against a native reference PDB.

    Parameters
    ----------
    model_pdb : str
        Path to model (predicted) PDB file.
    native_pdb : str
        Path to native (reference) PDB file.
    mapping : str, optional
        Chain mapping MODELCHAINS:NATIVECHAINS (e.g., "AB:HL"). Omit to let DockQ auto-detect.
    work_dir : str, optional
        Working directory for temp output. Default: system temp.
    n_cpu : int, optional
        Explicit CPU count passed to DockQ as ``--n_cpu``. Defaults to one
        to avoid unexpected parallelism when this low-level helper is called
        outside a resource-aware wrapper.

    Returns
    -------
    dict
        - dockq: DockQ score (0–1)
        - fnat: Fraction of native contacts
        - irmsd: Interface RMSD (Å)
        - lrmsd: Ligand RMSD (Å)
        - fnonnat: Fraction of non-native contacts
        - f1: F1 score (harmonic mean of precision and recall)
        - clashes: Number of clashing interfacial residues
        - interfaces: Per-interface results when multiple interfaces exist
    """
    model_pdb = Path(model_pdb).resolve()
    native_pdb = Path(native_pdb).resolve()
    if not model_pdb.exists():
        raise FileNotFoundError(f"Model PDB not found: {model_pdb}")
    if not native_pdb.exists():
        raise FileNotFoundError(f"Native PDB not found: {native_pdb}")

    if no_align:
        if mapping is None:
            raise ValueError("--no_align requires an explicit chain mapping")
        validate_no_align_mapping(model_pdb, native_pdb, mapping)

    work_dir = Path(work_dir) if work_dir else Path(tempfile.gettempdir())
    work_dir.mkdir(parents=True, exist_ok=True)
    json_file = work_dir / f"dockq_{os.getpid()}_{uuid.uuid4().hex}.json"

    dockq_bin = resolve_dockq_executable()
    if not dockq_bin:
        raise RuntimeError("DockQ not found. Install with: pip install DockQ")

    cmd = [
        *dockq_bin,
        str(model_pdb),
        str(native_pdb),
        "--json",
        str(json_file),
    ]
    if mapping is not None:
        cmd.extend(["--mapping", str(mapping)])
    if no_align:
        cmd.append("--no_align")
    if n_cpu is not None:
        try:
            n_cpu = int(n_cpu)
        except (TypeError, ValueError) as exc:
            raise ValueError("n_cpu must be a positive integer") from exc
        if n_cpu < 1:
            raise ValueError("n_cpu must be a positive integer")
        cmd.extend(["--n_cpu", str(n_cpu)])

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(work_dir),
            timeout=120,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"DockQ failed (exit {result.returncode}): {result.stderr or result.stdout}"
            )
    except subprocess.TimeoutExpired:
        raise RuntimeError("DockQ timed out")
    except Exception as e:
        raise RuntimeError(f"DockQ failed: {e}") from e

    if not json_file.exists():
        raise RuntimeError("DockQ did not produce JSON output")

    with open(json_file) as f:
        data = json.load(f)

    parsed = _parse_dockq_json(data)
    parsed.update(
        {
            "raw_dockq_json": str(json_file),
            "raw_dockq_json_sha256": hashlib.sha256(json_file.read_bytes()).hexdigest(),
            "dockq_argv": cmd,
            "dockq_no_align": bool(no_align),
            "dockq_n_cpu": n_cpu,
        }
    )
    return parsed


def _parse_dockq_json(data):
    """Extract a flat result dict from DockQ JSON (wallnerlab DockQ format)."""
    # DockQ JSON: best_result is keyed by interface and best_dockq is a sum.
    best_result = data.get("best_result", {})
    if isinstance(best_result, dict):
        interfaces = list(best_result.values())
    else:
        interfaces = best_result if isinstance(best_result, list) else []

    dockq_global = data.get("GlobalDockQ")
    if dockq_global is not None:
        try:
            dockq_global = float(dockq_global)
        except (TypeError, ValueError) as exc:
            raise ValueError("GlobalDockQ must be numeric") from exc
        if not 0.0 <= dockq_global <= 1.0:
            raise ValueError(f"GlobalDockQ is outside [0, 1]: {dockq_global}")

    fnat = irmsd = lrmsd = fnonnat = f1 = clashes = None
    if len(interfaces) == 1:
        first = interfaces[0]
        fnat = first.get("fnat")
        irmsd = first.get("iRMSD")
        lrmsd = first.get("LRMSD")
        fnonnat = first.get("fnonnat")
        f1 = first.get("F1")
        clashes = first.get("clashes", 0)

    return {
        "dockq": dockq_global,
        "dockq_global": dockq_global,
        "dockq_sum": data.get("best_dockq"),
        "dockq_json_status": "valid" if dockq_global is not None else "valid_unscored",
        "fnat": fnat,
        "irmsd": irmsd,
        "lrmsd": lrmsd,
        "fnonnat": fnonnat,
        "f1": f1,
        "clashes": clashes,
        "interfaces": interfaces,
    }


def dockq_to_capri_class(dockq):
    """
    Map DockQ score to CAPRI quality class.

    Returns
    -------
    str
        'Incorrect', 'Acceptable', 'Medium', or 'High'
    """
    if dockq is None:
        return "Unknown"
    if dockq < 0.23:
        return "Incorrect"
    if dockq < 0.49:
        return "Acceptable"
    if dockq < 0.80:
        return "Medium"
    return "High"


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Calculate DockQ between model and native PDB (benchmark metric)"
    )
    parser.add_argument("model_pdb", help="Model (predicted) PDB path")
    parser.add_argument("native_pdb", help="Native (reference) PDB path")
    parser.add_argument(
        "--mapping",
        default=1,
        help="Chain mapping MODELCHAINS:NATIVECHAINS (e.g., AB:HL)",
    )
    parser.add_argument(
        "--no-align",
        action="store_true",
        help="Pass --no_align to DockQ",
    )
    parser.add_argument(
        "--n-cpu",
        type=int,
        default=None,
        help="Explicit CPU count passed to DockQ as --n_cpu",
    )
    args = parser.parse_args()

    result = calculate_dockq(
        args.model_pdb,
        args.native_pdb,
        mapping=args.mapping,
        no_align=args.no_align,
        n_cpu=args.n_cpu,
    )
    print("DockQ:", result["dockq"])
    print("CAPRI class:", dockq_to_capri_class(result["dockq"]))
    print("fnat:", result["fnat"])
    print("iRMSD (Å):", result["irmsd"])
    print("LRMSD (Å):", result["lrmsd"])
    print("fnonnat:", result["fnonnat"])
    print("F1:", result["f1"])
    print("clashes:", result["clashes"])
