"""
Calculate DockQ between a model (predicted) PDB and a native (reference) PDB.

DockQ is a continuous quality measure (0–1) for protein–protein docking models,
combining fnat (fraction native contacts), LRMS (ligand RMSD), and iRMS (interface RMSD).
CAPRI interpretation: <0.23 Incorrect, 0.23–0.49 Acceptable, 0.49–0.80 Medium, ≥0.80 High.

Requires: pip install DockQ
Ref: https://github.com/wallnerlab/DockQ
"""
import json
import subprocess
import tempfile
from pathlib import Path


def calculate_dockq(model_pdb, native_pdb, mapping=None, work_dir=None):
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

    work_dir = Path(work_dir) if work_dir else Path(tempfile.gettempdir())
    work_dir.mkdir(parents=True, exist_ok=True)
    json_file = work_dir / "dockq_out.json"

    cmd = [
        "DockQ",
        str(model_pdb),
        str(native_pdb),
        "--json",
        str(json_file),
    ]
    if mapping is not None:
        cmd.extend(["--mapping", str(mapping)])

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
    except FileNotFoundError:
        raise RuntimeError(
            "DockQ not found. Install with: pip install DockQ"
        ) from None
    except subprocess.TimeoutExpired:
        raise RuntimeError("DockQ timed out")
    except Exception as e:
        raise RuntimeError(f"DockQ failed: {e}") from e

    if not json_file.exists():
        raise RuntimeError("DockQ did not produce JSON output")

    with open(json_file) as f:
        data = json.load(f)

    return _parse_dockq_json(data)


def _parse_dockq_json(data):
    """Extract a flat result dict from DockQ JSON (wallnerlab DockQ format)."""
    # DockQ JSON: best_result is dict keyed by (chain1, chain2), best_dockq is total
    best_result = data.get("best_result", {})
    if isinstance(best_result, dict):
        interfaces = list(best_result.values())
    else:
        interfaces = best_result if isinstance(best_result, list) else []

    dockq_total = data.get("best_dockq", data.get("GlobalDockQ"))
    if dockq_total is None and interfaces:
        scores = [i.get("DockQ") for i in interfaces if i.get("DockQ") is not None]
        dockq_total = max(scores) if scores else None

    fnat = irmsd = lrmsd = fnonnat = f1 = clashes = None
    if interfaces:
        first = interfaces[0]
        fnat = first.get("fnat")
        irmsd = first.get("iRMSD")
        lrmsd = first.get("LRMSD")
        fnonnat = first.get("fnonnat")
        f1 = first.get("F1")
        clashes = first.get("clashes", 0)

    return {
        "dockq": dockq_total,
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
        default=None,
        help="Chain mapping MODELCHAINS:NATIVECHAINS (e.g., AB:HL)",
    )
    args = parser.parse_args()

    result = calculate_dockq(
        args.model_pdb,
        args.native_pdb,
        mapping=args.mapping,
    )
    print("DockQ:", result["dockq"])
    print("CAPRI class:", dockq_to_capri_class(result["dockq"]))
    print("fnat:", result["fnat"])
    print("iRMSD (Å):", result["irmsd"])
    print("LRMSD (Å):", result["lrmsd"])
    print("fnonnat:", result["fnonnat"])
    print("F1:", result["f1"])
    print("clashes:", result["clashes"])
