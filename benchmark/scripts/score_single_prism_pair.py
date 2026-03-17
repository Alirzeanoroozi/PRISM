#!/usr/bin/env python3
"""
Score one model PDB or a folder of model PDBs against a native/reference PDB
with the same benchmark metrics used by the PRISM evaluation pipeline.

Single-file mode:
    python benchmark/scripts/score_single_prism_pair.py model.pdb native.pdb

Folder mode:
    python benchmark/scripts/score_single_prism_pair.py models_dir native.pdb --out-csv results.csv
"""

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd, timeout_sec, env=None):
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_sec,
            env=env,
        )
    except subprocess.TimeoutExpired:
        return {
            "ok": False,
            "stdout": "",
            "stderr": f"timed out after {timeout_sec}s",
            "returncode": None,
        }

    return {
        "ok": proc.returncode == 0,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
        "returncode": proc.returncode,
    }


def maybe_float(value):
    try:
        return float(value)
    except Exception:
        return None


def infer_chain_order(pdb_path):
    chains = []
    seen = set()
    with open(pdb_path) as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain_id = line[21].strip() or "_"
            if chain_id not in seen:
                seen.add(chain_id)
                chains.append(chain_id)
    if len(chains) != 2:
        raise ValueError(
            f"Automatic chain inference requires exactly 2 chains in {pdb_path}, found {len(chains)}: {chains}"
        )
    return chains[0], chains[1]


def resolve_model_inputs(path_arg):
    input_path = Path(path_arg)
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        return sorted(p for p in input_path.iterdir() if p.is_file() and p.suffix.lower() == ".pdb")
    raise FileNotFoundError(f"Input path not found: {input_path}")


def parse_dockq_short(stdout):
    values = {}
    tokens = stdout.strip().split()
    for i, token in enumerate(tokens[:-1]):
        if token in {"DockQ", "iRMSD", "LRMSD", "fnat", "fnonnat", "F1", "clashes"}:
            values[token] = maybe_float(tokens[i + 1])
    return values


def score_one(
    model_pdb,
    native_pdb,
    score_python,
    timeout_sec,
    dockq_no_align,
    model_receptor=None,
    model_ligand=None,
    native_receptor=None,
    native_ligand=None,
):
    model_receptor = model_receptor or infer_chain_order(model_pdb)[0]
    model_ligand = model_ligand or infer_chain_order(model_pdb)[1]
    native_receptor = native_receptor or infer_chain_order(native_pdb)[0]
    native_ligand = native_ligand or infer_chain_order(native_pdb)[1]

    repo_root = Path.cwd()
    irmsd_script = repo_root / "benchmark/scripts/irmsd.py"

    result = {
        "model_pdb": str(model_pdb.resolve()),
        "native_pdb": str(native_pdb.resolve()),
        "model_receptor": model_receptor,
        "model_ligand": model_ligand,
        "native_receptor": native_receptor,
        "native_ligand": native_ligand,
        "irmsd": None,
        "dockq": None,
        "dockq_irmsd": None,
        "dockq_lrmsd": None,
        "dockq_fnat": None,
        "dockq_fnonnat": None,
        "dockq_f1": None,
        "dockq_clashes": None,
        "error_irmsd": "",
        "error_dockq": "",
    }

    irmsd_cmd = [
        str(score_python),
        str(irmsd_script),
        str(model_pdb),
        model_receptor,
        model_ligand,
        str(native_pdb),
        native_receptor,
        native_ligand,
    ]
    irmsd_run = run_cmd(irmsd_cmd, timeout_sec)
    if irmsd_run["ok"]:
        result["irmsd"] = maybe_float(irmsd_run["stdout"])
    else:
        result["error_irmsd"] = irmsd_run["stderr"] or irmsd_run["stdout"] or "iRMSD failed"

    mapping = f"{model_receptor}{model_ligand}:{native_receptor}{native_ligand}"
    env = dict(os.environ)
    env["PATH"] = f"{score_python.parent}:{env.get('PATH', '')}"
    dockq_cmd = [
        str(score_python),
        "-m",
        "DockQ",
        str(model_pdb),
        str(native_pdb),
        "--mapping",
        mapping,
        "--short",
    ]
    if dockq_no_align:
        dockq_cmd.append("--no_align")
    dockq_run = run_cmd(dockq_cmd, timeout_sec, env=env)
    if dockq_run["ok"]:
        detail = parse_dockq_short(dockq_run["stdout"])
        result["dockq"] = detail.get("DockQ")
        result["dockq_irmsd"] = detail.get("iRMSD")
        result["dockq_lrmsd"] = detail.get("LRMSD")
        result["dockq_fnat"] = detail.get("fnat")
        result["dockq_fnonnat"] = detail.get("fnonnat")
        result["dockq_f1"] = detail.get("F1")
        result["dockq_clashes"] = detail.get("clashes")
    else:
        result["error_dockq"] = dockq_run["stderr"] or dockq_run["stdout"] or "DockQ failed"

    return result


def write_csv(rows, out_csv):
    fieldnames = [
        "model_pdb",
        "native_pdb",
        "model_receptor",
        "model_ligand",
        "native_receptor",
        "native_ligand",
        "irmsd",
        "dockq",
        "dockq_irmsd",
        "dockq_lrmsd",
        "dockq_fnat",
        "dockq_fnonnat",
        "dockq_f1",
        "dockq_clashes",
        "error_irmsd",
        "error_dockq",
    ]
    with open(out_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Score one model PDB or a folder of model PDBs with iRMSD and DockQ")
    parser.add_argument("model_input", help="Model PDB file or directory containing model PDB files")
    parser.add_argument("native_pdb", help="Native/reference PDB")
    parser.add_argument("--model-receptor", default=None, help="Model receptor chain(s), e.g. L")
    parser.add_argument("--model-ligand", default=None, help="Model ligand chain(s), e.g. H")
    parser.add_argument("--native-receptor", default=None, help="Native receptor chain(s), e.g. A")
    parser.add_argument("--native-ligand", default=None, help="Native ligand chain(s), e.g. B")
    parser.add_argument(
        "--score-python",
        default=sys.executable,
        help="Python interpreter used to run benchmark scoring helpers",
    )
    parser.add_argument(
        "--timeout-sec",
        type=int,
        default=120,
        help="Timeout for each scoring subprocess",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print machine-readable JSON for single-file mode",
    )
    parser.add_argument(
        "--out-csv",
        default=None,
        help="Write results to CSV. Required for folder mode; optional for single-file mode.",
    )
    parser.add_argument(
        "--dockq-no-align",
        action="store_true",
        help="Pass --no_align to DockQ for faster scoring when residue numbering/mapping is already trusted",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    score_python = Path(args.score_python)
    if not score_python.is_absolute():
        score_python = repo_root / score_python

    native_pdb = Path(args.native_pdb)
    if not native_pdb.is_absolute():
        native_pdb = repo_root / native_pdb

    model_inputs = resolve_model_inputs(args.model_input)
    if not model_inputs:
        raise SystemExit("No .pdb files found in input")

    rows = []
    for model_pdb in model_inputs:
        if not model_pdb.is_absolute():
            model_pdb = repo_root / model_pdb
        rows.append(
            score_one(
                model_pdb=model_pdb,
                native_pdb=native_pdb,
                score_python=score_python,
                timeout_sec=args.timeout_sec,
                dockq_no_align=args.dockq_no_align,
                model_receptor=args.model_receptor,
                model_ligand=args.model_ligand,
                native_receptor=args.native_receptor,
                native_ligand=args.native_ligand,
            )
        )

    if args.out_csv:
        out_csv = Path(args.out_csv)
        if not out_csv.is_absolute():
            out_csv = repo_root / out_csv
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        write_csv(rows, out_csv)

    if len(rows) == 1 and args.json:
        print(json.dumps(rows[0], indent=2, sort_keys=True))
        return 0

    if len(rows) == 1 and not args.out_csv:
        row = rows[0]
        print(f"model: {row['model_pdb']}")
        print(f"native: {row['native_pdb']}")
        print(
            "chains: "
            f"model {row['model_receptor']}/{row['model_ligand']} -> "
            f"native {row['native_receptor']}/{row['native_ligand']}"
        )
        print(f"iRMSD: {row['irmsd'] if row['irmsd'] is not None else 'NA'}")
        print(f"DockQ: {row['dockq'] if row['dockq'] is not None else 'NA'}")
        if row["error_irmsd"] or row["error_dockq"]:
            print("errors:")
            if row["error_irmsd"]:
                print(f"  irmsd: {row['error_irmsd']}")
            if row["error_dockq"]:
                print(f"  dockq: {row['error_dockq']}")
        return 0

    if args.out_csv:
        print(f"wrote {len(rows)} result rows to {out_csv}")
    else:
        print(json.dumps(rows, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
