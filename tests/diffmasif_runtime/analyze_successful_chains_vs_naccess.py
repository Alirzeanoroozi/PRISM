import argparse
import csv
import json
import os
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional


CURRENT_DIR = Path(__file__).resolve().parent
REPLACEMENT_DIR = CURRENT_DIR.parent / "diffmasif_replacement"
if str(REPLACEMENT_DIR) not in sys.path:
    sys.path.insert(0, str(REPLACEMENT_DIR))
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from compare_with_naccess import parse_rsa_file, summarize_comparison
from export_ply_vertices_to_csv import parse_ply_vertices
from export_residue_accessibility import build_proxy_accessibility_from_points


def successful_chain_stems(status_dir: Path) -> List[str]:
    stems: List[str] = []
    for path in sorted(status_dir.glob("*.status")):
        lines = path.read_text().splitlines()
        if lines and lines[0] == "ok":
            stems.append(path.stem)
    return stems


def run_naccess_chain_only(pdb_path: Path, chain_id: str) -> Dict[str, float]:
    naccess_exec = Path("external_tools/naccess/naccess").resolve()
    if not naccess_exec.exists():
        raise FileNotFoundError(f"Naccess executable not found: {naccess_exec}")

    with tempfile.TemporaryDirectory(prefix="masif_chain_naccess_", dir=".") as tmpdir:
        tmpdir_path = Path(tmpdir)
        local_pdb = tmpdir_path / pdb_path.name
        local_pdb.write_bytes(pdb_path.read_bytes())
        import subprocess

        result = subprocess.run(
            [str(naccess_exec), local_pdb.name],
            cwd=str(tmpdir_path),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            check=False,
        )
        rsa_candidates = [
            tmpdir_path / f"{local_pdb.stem}.rsa",
            tmpdir_path / f"{pdb_path.stem}.rsa",
            tmpdir_path / f"{pdb_path.stem.lower()}.rsa",
            tmpdir_path / f"{pdb_path.stem[:4].lower()}.rsa",
        ]
        rsa_path = next((candidate for candidate in rsa_candidates if candidate.exists()), None)
        if result.returncode != 0 or rsa_path is None:
            raise RuntimeError(
                f"Naccess failed for {pdb_path.name} (exit code {result.returncode}).\n{result.stdout}"
            )
        return parse_rsa_file(str(rsa_path), allowed_chains=[chain_id])


def compare_chain(stem: str, masif_root: Path, assignment_cutoff: float) -> Dict[str, object]:
    pdb_id, chain_id = stem.split("_")
    benchmark_pdb = masif_root / "data/masif_site/data_preparation/01-benchmark_pdbs" / f"{stem}.pdb"
    surface_ply = masif_root / "data/masif_site/data_preparation/01-benchmark_surfaces" / f"{stem}.ply"
    if not benchmark_pdb.exists():
        raise FileNotFoundError(f"Benchmark PDB not found: {benchmark_pdb}")
    if not surface_ply.exists():
        raise FileNotFoundError(f"Surface PLY not found: {surface_ply}")

    raw_points = parse_ply_vertices(str(surface_ply))
    points = [
        {
            "x": float(point["x"]),
            "y": float(point["y"]),
            "z": float(point["z"]),
        }
        for point in raw_points
    ]
    proxy_stats = build_proxy_accessibility_from_points(
        pdb_path=str(benchmark_pdb),
        points=points,
        allowed_chains=[chain_id],
        assignment_cutoff=assignment_cutoff,
    )
    naccess = run_naccess_chain_only(benchmark_pdb, chain_id=chain_id)
    summary = summarize_comparison(naccess, proxy_stats)
    return {
        "stem": stem,
        "pdb_id": pdb_id,
        "chain_id": chain_id,
        **summary,
        **proxy_stats["_meta"],
    }


def aggregate(rows: List[Dict[str, object]]) -> Dict[str, Optional[float]]:
    ok_rows = [row for row in rows if row["status"] == "ok"]
    if not ok_rows:
        return {
            "successful_chains": 0,
            "failed_chains": len(rows),
            "avg_buried_threshold_agreement_at_20": None,
            "avg_surface_threshold_agreement_at_15": None,
            "avg_mean_absolute_error": None,
        }

    def avg(key: str) -> float:
        values = [float(row[key]) for row in ok_rows if row[key] is not None]
        return sum(values) / len(values) if values else None

    return {
        "successful_chains": len(ok_rows),
        "failed_chains": len(rows) - len(ok_rows),
        "avg_buried_threshold_agreement_at_20": avg("buried_threshold_agreement_at_20"),
        "avg_surface_threshold_agreement_at_15": avg("surface_threshold_agreement_at_15"),
        "avg_mean_absolute_error": avg("mean_absolute_error"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--masif-root", default="/scratch/rshadi25/GitHub/masif")
    parser.add_argument("--status-dir", default="tests/diffmasif_runtime/slurm_logs")
    parser.add_argument("--assignment-cutoff", type=float, default=4.0)
    parser.add_argument("--max-chains", type=int, default=None)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--summary-json",
        default="tests/diffmasif_replacement/output/masif_successful_chains_vs_naccess.json",
    )
    parser.add_argument(
        "--summary-csv",
        default="tests/diffmasif_replacement/output/masif_successful_chains_vs_naccess.csv",
    )
    args = parser.parse_args()

    masif_root = Path(args.masif_root).resolve()
    stems = successful_chain_stems(Path(args.status_dir))
    if args.max_chains is not None:
        stems = stems[: args.max_chains]

    rows: List[Dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_map = {
            executor.submit(compare_chain, stem, masif_root, args.assignment_cutoff): stem
            for stem in stems
        }
        for future in as_completed(future_map):
            stem = future_map[future]
            try:
                row = future.result()
                row["status"] = "ok"
                row["error"] = ""
            except Exception as exc:
                row = {
                    "stem": stem,
                    "pdb_id": stem.split("_")[0],
                    "chain_id": stem.split("_")[1],
                    "status": "failed",
                    "error": str(exc),
                    "buried_threshold_agreement_at_20": None,
                    "surface_threshold_agreement_at_15": None,
                    "mean_absolute_error": None,
                    "common_residue_count": None,
                    "naccess_only_residue_count": None,
                    "proxy_only_residue_count": None,
                    "assignment_cutoff": args.assignment_cutoff,
                    "surface_point_count": None,
                    "assigned_surface_point_count": None,
                    "unassigned_surface_point_count": None,
                }
            rows.append(row)

    rows.sort(key=lambda item: item["stem"])
    summary_csv = Path(args.summary_csv)
    summary_json = Path(args.summary_json)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "stem",
        "pdb_id",
        "chain_id",
        "status",
        "buried_threshold_agreement_at_20",
        "surface_threshold_agreement_at_15",
        "mean_absolute_error",
        "common_residue_count",
        "naccess_only_residue_count",
        "naccess_only_residues",
        "proxy_only_residue_count",
        "proxy_only_residues",
        "assignment_cutoff",
        "surface_point_count",
        "assigned_surface_point_count",
        "unassigned_surface_point_count",
        "error",
    ]
    with summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "masif_root": str(masif_root),
        "status_dir": str(Path(args.status_dir).resolve()),
        "chain_count": len(stems),
        "aggregate": aggregate(rows),
        "chains": rows,
    }
    with summary_json.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)

    print(summary_json)
    print(json.dumps(payload["aggregate"], sort_keys=True))


if __name__ == "__main__":
    main()
