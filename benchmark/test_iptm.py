"""
Calculate TM-score between an input PDB and a reference PDB using TMalign.

Note: ipTM (interface predicted TM-score) is an AlphaFold-internal confidence metric
produced during prediction—it cannot be computed from a PDB alone. This script
computes TM-score via TMalign, the standard post-hoc metric for benchmarking
predicted vs reference structures.
"""
import os
import re
import subprocess
import tempfile
from pathlib import Path


def calculate_tmscore(input_pdb, reference_pdb, tmalign_path=None, work_dir=None):
    """
    Run TMalign between input and reference PDBs and parse TM-score.

    Parameters
    ----------
    input_pdb : str
        Path to input (predicted) PDB file.
    reference_pdb : str
        Path to reference (native) PDB file.
    tmalign_path : str, optional
        Path to TMalign executable. Default: external_tools/TMalign (relative to PRISM root).
    work_dir : str, optional
        Working directory for temp output. Default: system temp.

    Returns
    -------
    dict
        - tmscore_ref: TM-score normalized by reference length (use this for benchmarking)
        - tmscore_query: TM-score normalized by query length
        - rmsd: RMSD of aligned residues
        - aligned_length: number of aligned residues
        - seq_id: sequence identity (n_identical/n_aligned)
    """
    input_pdb = Path(input_pdb).resolve()
    reference_pdb = Path(reference_pdb).resolve()
    if not input_pdb.exists():
        raise FileNotFoundError(f"Input PDB not found: {input_pdb}")
    if not reference_pdb.exists():
        raise FileNotFoundError(f"Reference PDB not found: {reference_pdb}")

    if tmalign_path is None:
        # Assume PRISM root is parent of benchmark/
        prism_root = Path(__file__).resolve().parent.parent
        tmalign_path = prism_root / "external_tools" / "TMalign"
    tmalign_path = Path(tmalign_path).resolve()
    if not tmalign_path.exists():
        raise FileNotFoundError(f"TMalign not found: {tmalign_path}")

    work_dir = Path(work_dir) if work_dir else Path(tempfile.gettempdir())
    work_dir.mkdir(parents=True, exist_ok=True)
    out_file = work_dir / "tmalign_out.tm"

    cmd = [
        str(tmalign_path),
        str(input_pdb),
        str(reference_pdb),
    ]
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(work_dir),
            timeout=60,
        )
        output = result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        raise RuntimeError("TMalign timed out")
    except Exception as e:
        raise RuntimeError(f"TMalign failed: {e}")

    return _parse_tmalign_output(output)


def _parse_tmalign_output(output):
    """Parse TMalign stdout/stderr for TM-score and RMSD."""
    tmscore_ref = None
    tmscore_query = None
    rmsd = None
    aligned_length = None
    seq_id = None

    lines = output.splitlines()
    for line in lines:
        if line.startswith("Aligned length"):
            # Format: Aligned length= 51, RMSD=   0.00, Seq_ID=n_identical/n_aligned= 1.000
            m = re.search(r"Aligned length=\s*(\d+)", line)
            if m:
                aligned_length = int(m.group(1))
            m = re.search(r"RMSD=\s*([\d.]+)", line)
            if m:
                rmsd = float(m.group(1))
            m = re.search(r"Seq_ID=n_identical/n_aligned=\s*([\d.]+)", line)
            if m:
                seq_id = float(m.group(1))
        elif "TM-score=" in line and "(if normalized" in line:
            # Format: TM-score= 0.20319 (if normalized by length of Chain_1, ...)
            try:
                tmscore = float(line.split("TM-score=")[1].split()[0])
                if "Chain_1" in line:
                    tmscore_query = tmscore
                elif "Chain_2" in line:
                    tmscore_ref = tmscore
            except (ValueError, IndexError):
                pass

    # TMalign prints Chain_1 first, Chain_2 second; assign by order if labels missing
    if tmscore_ref is None and tmscore_query is None:
        tm_scores = []
        for line in lines:
            if "TM-score=" in line and "(if normalized" in line:
                try:
                    val = float(line.split("TM-score=")[1].split()[0])
                    tm_scores.append(val)
                except (ValueError, IndexError):
                    pass
        if len(tm_scores) >= 2:
            tmscore_query, tmscore_ref = tm_scores[0], tm_scores[1]
        elif len(tm_scores) == 1:
            tmscore_ref = tm_scores[0]

    return {
        "tmscore_ref": tmscore_ref,
        "tmscore_query": tmscore_query,
        "rmsd": rmsd,
        "aligned_length": aligned_length,
        "seq_id": seq_id,
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Calculate TM-score between input and reference PDB (benchmark metric)"
    )
    parser.add_argument("input_pdb", help="Input (predicted) PDB path")
    parser.add_argument("reference_pdb", help="Reference (native) PDB path")
    parser.add_argument("--tmalign", default=None, help="Path to TMalign executable")
    args = parser.parse_args()
    result = calculate_tmscore(args.input_pdb, args.reference_pdb, tmalign_path=args.tmalign)
    print("TM-score (normalized by reference):", result["tmscore_ref"])
    print("TM-score (normalized by query):", result["tmscore_query"])
    print("RMSD:", result["rmsd"])
    print("Aligned length:", result["aligned_length"])
    print("Seq ID:", result["seq_id"])
