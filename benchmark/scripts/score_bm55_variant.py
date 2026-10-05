#!/usr/bin/env python3
"""
Score one BM5.5 pipeline variant against T_Rigid, T_medium, and T_difficult.

The BM5.5 parallel pipeline produces model filenames like:
    {template}_{receptor}_{ligand}_o{ori}_L_{template}_{receptor}_{ligand}_o{ori}_R_rosetta_0001_0001.pdb

This script:
1. Scans the variant's rosetta_refinement/ (or fiberdock_refinement/) for model PDBs
2. Extracts (template, receptor, ligand) from the filename
3. Maps the native complex (template[:4]) to benchmark CSV rows
4. Scores each model with DockQ and iRMSD using the benchmark's native PDBs
5. Writes per_prediction.csv and pair_summary.csv for each benchmark set
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Regex for BM5.5 combined Rosetta refinement filenames
# Example: 1ahwAF_1fgnHL_1tfhA_o1_L_1ahwAF_1fgnHL_1tfhA_o1_R_rosetta_0001_0001.pdb
ROSETTA_MODEL_RE = re.compile(
    r"^(?P<template>[A-Za-z0-9]+)_"
    r"(?P<receptor>[A-Za-z0-9]+)_"
    r"(?P<ligand>[A-Za-z0-9]+)_"
    r"o(?P<orientation>\d+)_[LR]_"
    r"[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_o\d+_[LR]_"
    r"rosetta_\d+_\d+\.pdb$"
)

# Regex for FiberDock refinement filenames (if any)
# Example: 1ahwAF_1fgnHL_1tfhA_o1_L_1ahwAF_1fgnHL_1tfhA_o1_R.ref.pdb
FIBERDOCK_MODEL_RE = re.compile(
    r"^(?P<template>[A-Za-z0-9]+)_"
    r"(?P<receptor>[A-Za-z0-9]+)_"
    r"(?P<ligand>[A-Za-z0-9]+)_"
    r"o(?P<orientation>\d+)_[LR]_"
    r"[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+_o\d+_[LR]\.ref\.pdb$"
)


def load_benchmark_set(csv_path: Path) -> Dict[str, List[dict]]:
    """Load a benchmark CSV and index by native complex PDB ID (first 4 chars)."""
    rows = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            complex_raw = row.get("Complex", row.get("PDB ID 1", "")).strip()
            pdb4 = complex_raw[:4].lower()
            rows.append({
                "complex": complex_raw,
                "pdb_id_1": row.get("PDB ID 1", "").strip(),
                "pdb_id_2": row.get("PDB ID 2", "").strip(),
                "pdb4": pdb4,
                "row_index": len(rows),
            })
    # Index by pdb4
    index: Dict[str, List[dict]] = defaultdict(list)
    for r in rows:
        index[r["pdb4"]].append(r)
    return dict(index)


def find_model_pdbs(variant_root: Path) -> List[Tuple[Path, str, str, str, str]]:
    """Find all model PDBs and extract (path, template, receptor, ligand, native_pdb4)."""
    models = []

    # Rosetta refinement dir
    for pdb_path in variant_root.rglob("processed/rosetta_refinement/*_0001_0001.pdb"):
        m = ROSETTA_MODEL_RE.match(pdb_path.name)
        if m:
            template = m.group("template")
            receptor = m.group("receptor")
            ligand = m.group("ligand")
            native_pdb4 = template[:4].lower()
            models.append((pdb_path, template, receptor, ligand, native_pdb4))

    # FiberDock refinement dir
    for pdb_path in variant_root.rglob("processed/fiberdock_refinement/*.ref.pdb"):
        m = FIBERDOCK_MODEL_RE.match(pdb_path.name)
        if m:
            template = m.group("template")
            receptor = m.group("receptor")
            ligand = m.group("ligand")
            native_pdb4 = template[:4].lower()
            models.append((pdb_path, template, receptor, ligand, native_pdb4))

    return models


def run_dockq(
    score_python: str,
    model_pdb: str,
    native_pdb: str,
    model_receptor: str,
    model_ligand: str,
    native_receptor: Optional[str],
    native_ligand: Optional[str],
    timeout_sec: int = 120,
) -> dict:
    """Run DockQ on a model-native pair."""
    if native_receptor and native_ligand:
        mapping = f"{model_receptor}{model_ligand}:{native_receptor}{native_ligand}"
    else:
        mapping = f"{model_receptor}{model_ligand}:{model_receptor}{model_ligand}"

    cmd = [
        score_python, "-m", "DockQ",
        model_pdb, native_pdb,
        "--mapping", mapping,
        "--short",
    ]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)
        stdout = proc.stdout.strip()
        stderr = proc.stderr.strip()

        result = {"dockq": None, "fnat": None, "irmsd": None, "lrmsd": None, "error": ""}
        if proc.returncode == 0 and stdout:
            # Parse DockQ short output: DockQ 0.85 iRMSD 1.2 LRMSD 3.4 fnat 0.7 ...
            parts = stdout.split()
            for i, p in enumerate(parts):
                if p == "DockQ" and i + 1 < len(parts):
                    result["dockq"] = float(parts[i + 1])
                elif p == "iRMSD" and i + 1 < len(parts):
                    result["irmsd"] = float(parts[i + 1])
                elif p == "LRMSD" and i + 1 < len(parts):
                    result["lrmsd"] = float(parts[i + 1])
                elif p == "fnat" and i + 1 < len(parts):
                    result["fnat"] = float(parts[i + 1])
        else:
            result["error"] = stderr[:200] if stderr else f"exit={proc.returncode}"
        return result
    except subprocess.TimeoutExpired:
        return {"dockq": None, "fnat": None, "irmsd": None, "lrmsd": None, "error": "timeout"}
    except Exception as e:
        return {"dockq": None, "fnat": None, "irmsd": None, "lrmsd": None, "error": str(e)[:200]}


PER_PREDICTION_FIELDS = [
    "model_pdb", "template", "receptor", "ligand",
    "model_receptor_chain", "model_ligand_chain",
    "native_complex", "native_pdb",
    "dockq", "fnat", "irmsd", "lrmsd",
    "mapping", "error",
]


def _detect_chains(pdb_path: str) -> List[str]:
    """Detect protein chains (ATOM records) in a PDB file."""
    chains = []
    seen = set()
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM") and len(line) > 21:
                    ch = line[21]
                    if ch not in seen:
                        seen.add(ch)
                        chains.append(ch)
    except Exception:
        pass
    return chains


def score_variant(
    variant_root: Path,
    benchmark_sets: Dict[str, Path],
    native_dir: Path,
    score_python: str,
    out_dir: Path,
    variant_name: str,
    verbose_every: int = 200,
    timeout_sec: int = 120,
) -> None:
    """Score all models in a variant against benchmark sets."""
    models = find_model_pdbs(variant_root)
    print(f"[{variant_name}] Found {len(models)} model PDBs")

    if not models:
        print(f"[{variant_name}] No models found — skipping")
        return

    # Index models by native complex
    models_by_native: Dict[str, List] = defaultdict(list)
    for model_info in models:
        models_by_native[model_info[4]].append(model_info)

    for set_name, csv_path in benchmark_sets.items():
        out_subdir = out_dir / f"{variant_name}_{set_name}"
        out_subdir.mkdir(parents=True, exist_ok=True)

        bench_index = load_benchmark_set(csv_path)
        if not bench_index:
            print(f"[{variant_name}] [{set_name}] No benchmark rows loaded — skipping")
            continue

        native_pdb_dir = native_dir / f"native_bound_complexes_{set_name}"
        os.makedirs(native_pdb_dir, exist_ok=True)

        # Download missing natives
        all_native_pdb4s = set()
        for rows in bench_index.values():
            for r in rows:
                all_native_pdb4s.add(r["pdb4"])

        for pdb4 in sorted(all_native_pdb4s):
            npdb = native_pdb_dir / f"{pdb4}.pdb"
            if not npdb.exists():
                url = f"https://files.rcsb.org/download/{pdb4}.pdb"
                try:
                    subprocess.run(["curl", "-s", "-o", str(npdb), url], timeout=60, check=True)
                except Exception:
                    print(f"  Could not download {pdb4}.pdb")

        # Score models
        results = []
        total = len(models)
        for idx, (model_path, template, receptor, ligand, native_pdb4) in enumerate(models, 1):
            if native_pdb4 not in bench_index:
                continue

            npdb = native_pdb_dir / f"{native_pdb4}.pdb"
            if not npdb.exists():
                continue

            # Detect actual chains in the model PDB
            model_chains = _detect_chains(str(model_path))
            if not model_chains:
                continue

            # The model has chains A, B (or sometimes renamed).
            # Rosetta refinement produces partners with chains assigned by template order.
            # For DockQ: if model has 2 chains, use first as receptor, second as ligand.
            if len(model_chains) >= 2:
                rec_chain, lig_chain = model_chains[0], model_chains[1]
            else:
                continue

            # For this model, try each benchmark row that matches this native
            for bench_row in bench_index[native_pdb4]:
                dockq_result = run_dockq(
                    score_python,
                    str(model_path),
                    str(npdb),
                    model_receptor=rec_chain,
                    model_ligand=lig_chain,
                    native_receptor=None,
                    native_ligand=None,
                    timeout_sec=timeout_sec,
                )

                mapping = f"{rec_chain}:{bench_row['pdb_id_1']}_{lig_chain}:{bench_row['pdb_id_2']}"
                results.append({
                    "model_pdb": str(model_path),
                    "template": template,
                    "receptor": receptor,
                    "ligand": ligand,
                    "model_receptor_chain": rec_chain,
                    "model_ligand_chain": lig_chain,
                    "native_complex": bench_row["complex"],
                    "native_pdb": str(npdb),
                    "dockq": dockq_result.get("dockq", ""),
                    "fnat": dockq_result.get("fnat", ""),
                    "irmsd": dockq_result.get("irmsd", ""),
                    "lrmsd": dockq_result.get("lrmsd", ""),
                    "mapping": mapping,
                    "error": dockq_result.get("error", ""),
                })

            if idx % verbose_every == 0 or idx == total:
                matched = sum(1 for r in results if r["dockq"] is not None and r["dockq"] != "")
                print(f"[{variant_name}] [{set_name}] Scored {idx}/{total} models ({matched} with DockQ)")

        # Write per_prediction.csv
        pp_csv = out_subdir / "per_prediction.csv"
        with open(pp_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=PER_PREDICTION_FIELDS)
            writer.writeheader()
            writer.writerows(results)
        print(f"[{variant_name}] [{set_name}] Wrote {pp_csv} ({len(results)} rows)")


def main():
    parser = argparse.ArgumentParser(description="Score BM5.5 pipeline variants against benchmarks")
    parser.add_argument("--variant-root", required=True, help="Root directory of one pipeline variant run")
    parser.add_argument("--variant-name", required=True, help="Name for this variant (e.g. naccess_gt_external_rosetta)")
    parser.add_argument("--benchmark-dir", default="benchmark/data", help="Directory with T_Rigid.csv, T_medium.csv, T_difficult.csv")
    parser.add_argument("--native-dir", default="benchmark/prism_processed_results", help="Directory for native bound complex PDBs")
    parser.add_argument("--out-dir", default="benchmark/prism_processed_results", help="Output directory for scoring CSVs")
    parser.add_argument("--score-python", default=sys.executable, help="Python with DockQ installed")
    parser.add_argument("--sets", default="rigid,medium,difficult", help="Comma-separated benchmark sets")
    parser.add_argument("--verbose-every", type=int, default=200)
    parser.add_argument("--timeout-sec", type=int, default=120)
    args = parser.parse_args()

    variant_root = Path(args.variant_root)
    if not variant_root.exists():
        parser.error(f"Variant root not found: {variant_root}")

    benchmark_dir = Path(args.benchmark_dir)
    native_dir = Path(args.native_dir) if Path(args.native_dir).is_absolute() else Path.cwd() / args.native_dir
    out_dir = Path.cwd() / args.out_dir

    set_map = {}
    for s in args.sets.split(","):
        s = s.strip()
        if s == "rigid":
            set_map["t_rigid"] = benchmark_dir / "T_Rigid.csv"
        elif s == "medium":
            set_map["t_medium"] = benchmark_dir / "T_medium.csv"
        elif s == "difficult":
            set_map["t_difficult"] = benchmark_dir / "T_difficult.csv"

    score_variant(
        variant_root=variant_root,
        benchmark_sets=set_map,
        native_dir=native_dir,
        score_python=args.score_python,
        out_dir=out_dir,
        variant_name=args.variant_name,
        verbose_every=args.verbose_every,
        timeout_sec=args.timeout_sec,
    )


if __name__ == "__main__":
    main()
