import argparse
import csv
import json
import os
import subprocess
import sys
from typing import Dict, List, Tuple


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
REPLACEMENT_DIR = os.path.join(os.path.dirname(CURRENT_DIR), "diffmasif_replacement")
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from export_ply_vertices_to_csv import parse_ply_vertices


def resolve_target_structure(structure_id: str) -> Tuple[str, List[str]]:
    if len(structure_id) < 5:
        raise ValueError(f"Target id must include a chain: {structure_id}")
    pdb_path = os.path.abspath(f"processed/pdbs/{structure_id[:4].lower()}.pdb")
    chain_id = structure_id[4]
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"Target PDB not found: {pdb_path}")
    return pdb_path, [chain_id]


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def run_command(args: List[str], cwd: str = None) -> None:
    result = subprocess.run(
        args,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(args)}\n{result.stdout}")


def write_surface_csv(ply_path: str, output_csv: str) -> None:
    rows = parse_ply_vertices(ply_path)
    ensure_parent_dir(output_csv)
    with open(output_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["x", "y", "z"])
        writer.writeheader()
        writer.writerows(rows)


def summarize_successes(rows: List[Dict[str, object]]) -> Dict[str, object]:
    successes = [row for row in rows if row["status"] == "ok"]
    if not successes:
        return {
            "successful_runs": 0,
            "failed_runs": len(rows),
            "average_buried_threshold_agreement_at_20": None,
            "average_surface_threshold_agreement_at_15": None,
            "average_mean_absolute_error": None,
        }

    def avg(key: str) -> float:
        values = [float(row[key]) for row in successes]
        return sum(values) / len(values)

    return {
        "successful_runs": len(successes),
        "failed_runs": len(rows) - len(successes),
        "average_buried_threshold_agreement_at_20": avg("buried_threshold_agreement_at_20"),
        "average_surface_threshold_agreement_at_15": avg("surface_threshold_agreement_at_15"),
        "average_mean_absolute_error": avg("mean_absolute_error"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--masif-root", required=True)
    parser.add_argument("--structure-id", action="append", required=True)
    parser.add_argument("--output-dir", default="tests/diffmasif_runtime/output/batch")
    parser.add_argument(
        "--summary-json",
        default="tests/diffmasif_replacement/output/masif_batch_summary.json",
    )
    parser.add_argument(
        "--summary-csv",
        default="tests/diffmasif_replacement/output/masif_batch_summary.csv",
    )
    parser.add_argument("--assignment-cutoff", type=float, default=4.0)
    parser.add_argument("--reuse-existing", action="store_true")
    args = parser.parse_args()

    masif_root = os.path.abspath(args.masif_root)
    if not os.path.isdir(masif_root):
        raise FileNotFoundError(f"MaSIF root not found: {masif_root}")

    rows: List[Dict[str, object]] = []
    for structure_id in args.structure_id:
        structure_id = structure_id.upper()
        pdb_path, chains = resolve_target_structure(structure_id)
        chain_id = chains[0]
        masif_pair_id = f"{structure_id[:4]}_{chain_id}"
        surface_ply = os.path.join(
            masif_root,
            "data/masif_site/data_preparation/01-benchmark_surfaces",
            f"{masif_pair_id}.ply",
        )
        surface_csv = os.path.abspath(
            os.path.join(args.output_dir, f"{masif_pair_id}.surface_points.csv")
        )
        comparison_json = os.path.abspath(
            os.path.join(
                "tests/diffmasif_replacement/output",
                f"{masif_pair_id}.masif_surface.comparison.json",
            )
        )

        row: Dict[str, object] = {
            "structure_id": structure_id,
            "masif_pair_id": masif_pair_id,
            "status": "ok",
            "surface_csv": surface_csv,
            "comparison_json": comparison_json,
            "error": "",
        }
        try:
            if not (
                args.reuse_existing
                and os.path.exists(surface_csv)
                and os.path.exists(comparison_json)
            ):
                run_command(
                    ["./data_prepare_one.sh", "--file", pdb_path, masif_pair_id],
                    cwd=os.path.join(masif_root, "data/masif_site"),
                )
                if not os.path.exists(surface_ply):
                    raise FileNotFoundError(f"Expected MaSIF surface not found: {surface_ply}")
                write_surface_csv(surface_ply, surface_csv)
                run_command(
                    [
                        "python3",
                        os.path.join(REPLACEMENT_DIR, "compare_with_naccess.py"),
                        "--structure-kind",
                        "target",
                        "--structure-id",
                        structure_id,
                        "--surface-points",
                        surface_csv,
                        "--assignment-cutoff",
                        str(args.assignment_cutoff),
                        "--output",
                        comparison_json,
                    ]
                )
            with open(comparison_json, "r") as handle:
                data = json.load(handle)
            row.update(data["summary"])
            row.update(data["proxy_meta"])
        except Exception as exc:
            row["status"] = "failed"
            row["error"] = str(exc)
            row.setdefault("buried_threshold_agreement_at_20", "")
            row.setdefault("surface_threshold_agreement_at_15", "")
            row.setdefault("mean_absolute_error", "")
            row.setdefault("common_residue_count", "")
            row.setdefault("assigned_surface_point_count", "")
            row.setdefault("surface_point_count", "")
            row.setdefault("unassigned_surface_point_count", "")
        rows.append(row)

    ensure_parent_dir(args.summary_csv)
    ensure_parent_dir(args.summary_json)

    preferred_fieldnames = [
        "structure_id",
        "masif_pair_id",
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
        "assigned_surface_point_count",
        "surface_point_count",
        "unassigned_surface_point_count",
        "surface_csv",
        "comparison_json",
        "error",
    ]
    discovered = set()
    for row in rows:
        discovered.update(row.keys())
    fieldnames = preferred_fieldnames + sorted(discovered - set(preferred_fieldnames))
    with open(args.summary_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "masif_root": masif_root,
        "assignment_cutoff": args.assignment_cutoff,
        "targets": rows,
        "aggregate": summarize_successes(rows),
    }
    with open(args.summary_json, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)

    print(args.summary_json)
    print(json.dumps(payload["aggregate"], sort_keys=True))


if __name__ == "__main__":
    main()
