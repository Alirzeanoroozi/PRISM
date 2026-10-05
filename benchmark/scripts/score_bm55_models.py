#!/usr/bin/env python3
"""
Score BM5.5 pipeline variant models against benchmark CSVs using DockQ.

Handles the combined Rosetta output filenames from the BM5.5 parallel pipeline
by detecting the TER-split chain groups in each model PDB.
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.PDB import PDBParser


def infer_chain_order_from_ter(pdb_path: str) -> Tuple[str, str]:
    """Split model PDB on TER markers to infer receptor/ligand chain groups."""
    chains = []
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('TER') and chains:
                break
            if line.startswith('ATOM') and len(line) > 21:
                ch = line[21]
                if ch not in seen:
                    seen.add(ch)
                    chains.append(ch)
    if not chains:
        raise ValueError(f"No ATOM chains in {pdb_path}")
    mid = len(chains) // 2 + (len(chains) % 2)
    return ''.join(chains[:mid]), ''.join(chains[mid:])


def load_benchmark_set(csv_path: Path) -> Dict[str, List[dict]]:
    """Load a benchmark CSV; index by native complex PDB ID."""
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            complex_raw = row.get("Complex", "").strip()
            if ":" not in complex_raw or len(complex_raw) < 6:
                continue
            pdb4 = complex_raw[:4].lower()
            rows.append({
                "complex": complex_raw,
                "pdb_id_1": row.get("PDB ID 1", "").strip(),
                "pdb_id_2": row.get("PDB ID 2", "").strip(),
                "benchmark_set": csv_path.stem,
            })
    index = defaultdict(list)
    for r in rows:
        index[r["complex"][:4].lower()].append(r)
    return dict(index)


def run_dockq(score_python: str, model_pdb: str, native_pdb: str,
              mapping: str, timeout: int = 120) -> dict:
    cmd = [score_python, "-m", "DockQ", model_pdb, native_pdb,
           "--mapping", mapping, "--short"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        stdout = proc.stdout.strip()
        stderr = proc.stderr.strip()
        result = {"dockq": "", "fnat": "", "irmsd": "", "lrmsd": "", "error": ""}
        if proc.returncode == 0 and stdout:
            parts = stdout.split()
            for i, p in enumerate(parts):
                if p == "DockQ" and i + 1 < len(parts):
                    result["dockq"] = parts[i+1]
                elif p == "fnat" and i + 1 < len(parts):
                    result["fnat"] = parts[i+1]
                elif p == "iRMSD" and i + 1 < len(parts):
                    result["irmsd"] = parts[i+1]
                elif p == "LRMSD" and i + 1 < len(parts):
                    result["lrmsd"] = parts[i+1]
        else:
            result["error"] = (stderr[:200] if stderr else f"exit={proc.returncode}")
        return result
    except subprocess.TimeoutExpired:
        return {"dockq": "", "fnat": "", "irmsd": "", "lrmsd": "", "error": "timeout"}
    except Exception as e:
        return {"dockq": "", "fnat": "", "irmsd": "", "lrmsd": "", "error": str(e)[:200]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant-root", required=True)
    parser.add_argument("--variant-name", required=True)
    parser.add_argument("--benchmark-dir", default="benchmark/data")
    parser.add_argument("--native-dir", default="benchmark/prism_processed_results")
    parser.add_argument("--out-dir", default="benchmark/prism_processed_results")
    parser.add_argument("--score-python", default="/scratch/tmp/prism-dockq-env/bin/python")
    parser.add_argument("--sets", default="rigid,medium,difficult")
    parser.add_argument("--verbose-every", type=int, default=200)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()

    variant_root = Path(args.variant_root)
    benchmark_dir = Path(args.benchmark_dir)
    native_base = Path(args.native_dir)
    out_dir = Path(args.out_dir)

    set_map = {}
    for s in args.sets.split(","):
        s = s.strip()
        csv_path = benchmark_dir / f"T_{s.capitalize()}.csv"
        if csv_path.exists():
            set_map[s] = csv_path

    # Discover model PDBs — handle multiple naming conventions
    models = []
    # External Rosetta: "rosetta_refinement/*_0001_0001.pdb"
    models.extend(variant_root.rglob("processed/rosetta_refinement/*_0001_0001.pdb"))
    # PyRosetta: "pyrosetta_refinement/structures/*_rosetta.pdb"
    models.extend(variant_root.rglob("processed/pyrosetta_refinement/structures/*_rosetta.pdb"))
    # FiberDock: "fiberdock_refinement/*/*.ref.pdb" (one subdir per pair)
    for pdb_path in variant_root.rglob("processed/fiberdock_refinement/*/*.ref.pdb"):
        models.append(pdb_path)

    if not models:
        print(f"[{args.variant_name}] No model PDBs found")
        return
    print(f"[{args.variant_name}] Found {len(models)} model PDBs")

    set_dir_map = {
        "rigid": "t_rigid",
        "medium": "t_medium",
        "difficult": "t_difficult",
    }

    for set_name, csv_path in set_map.items():
        set_dir = set_dir_map.get(set_name, set_name)
        bench_index = load_benchmark_set(csv_path)
        native_pdb_dir = native_base / f"native_bound_complexes_{set_dir}"
        os.makedirs(native_pdb_dir, exist_ok=True)

        # Ensure natives exist
        for pdb4 in sorted(bench_index.keys()):
            npdb = native_pdb_dir / f"{pdb4}.pdb"
            if not npdb.exists():
                try:
                    subprocess.run(["curl", "-s", "-o", str(npdb),
                                    f"https://files.rcsb.org/download/{pdb4}.pdb"],
                                   timeout=60, check=True)
                except Exception:
                    pass

        out_subdir = out_dir / f"{args.variant_name}_{set_name}"
        out_subdir.mkdir(parents=True, exist_ok=True)
        results = []
        total = 0

        for pdb_path in models:
            # Determine native_pdb4 differently for FiberDock vs Rosetta filenames
            if "fiberdock" in str(pdb_path):
                # FiberDock files: parent dir contains template info (e.g. 1a4lAC_1rrpAB_1yrgB_o2)
                native_pdb4 = pdb_path.parent.name.split("_")[0][:4].lower()
            else:
                native_pdb4 = Path(pdb_path.name).stem[:4].lower()
            if native_pdb4 not in bench_index:
                continue

            # Detect model chains by TER splitting
            try:
                model_rec, model_lig = infer_chain_order_from_ter(str(pdb_path))
            except (ValueError, OSError):
                continue

            # For each benchmark row matching this native complex
            for bench_row in bench_index[native_pdb4]:
                native_complex = bench_row["complex"]
                npdb = native_pdb_dir / f"{native_pdb4}.pdb"
                if not npdb.exists():
                    continue

                # Native chains from complex string: e.g., 1AHW_A:B -> A, B
                native_str = native_complex.split("_", 1)[1] if "_" in native_complex else native_complex
                native_r = native_str.split(":")[0]
                native_l = native_str.split(":")[1]

                mapping = f"{model_rec}{model_lig}:{native_r}{native_l}"
                result = run_dockq(args.score_python, str(pdb_path), str(npdb), mapping, args.timeout)

                row = {
                    "model_pdb": str(pdb_path),
                    "native_complex": native_complex,
                    "native_pdb": str(npdb),
                    "mapping": mapping,
                    "dockq": result["dockq"],
                    "fnat": result["fnat"],
                    "irmsd": result["irmsd"],
                    "lrmsd": result["lrmsd"],
                    "error": result["error"],
                }
                results.append(row)
                total += 1

            if total % args.verbose_every == 0 and total > 0:
                scored = sum(1 for r in results if r["dockq"])
                print(f"  [{set_name}] Processed {total} models ({scored} scored)")

        # Write CSV
        pp_csv = out_subdir / "per_prediction.csv"
        fields = ["model_pdb", "native_complex", "native_pdb", "mapping",
                   "dockq", "fnat", "irmsd", "lrmsd", "error"]
        with open(pp_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(results)
        scored = sum(1 for r in results if r["dockq"])
        print(f"[{args.variant_name}] [{set_name}] {pp_csv}: {len(results)} rows ({scored} with DockQ)")


if __name__ == "__main__":
    main()
