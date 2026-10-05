import argparse
import csv
import json
import os
import subprocess
from typing import List


def run_command(args: List[str]) -> None:
    result = subprocess.run(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(args)}\n{result.stdout}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure-kind", choices=["target", "template"], default="target")
    parser.add_argument("--structure-id", action="append", required=True)
    parser.add_argument("--output-dir", default="tests/diffmasif_replacement/output")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    summaries = []
    for structure_id in args.structure_id:
        if args.structure_kind == "target":
            pdb_path = f"processed/pdbs/{structure_id[:4].lower()}.pdb"
            chains = [structure_id[4]]
        else:
            pdb_path = f"templates/pdbs/{structure_id[:4].lower()}.pdb"
            chains = [structure_id[4], structure_id[5]]

        surface_csv = os.path.join(args.output_dir, f"{structure_id}.geometric_surface.csv")
        comparison_json = os.path.join(args.output_dir, f"{structure_id}.comparison.json")

        run_command(
            [
                "python3",
                "tests/diffmasif_replacement/generate_geometric_surface_points.py",
                "--pdb",
                pdb_path,
                "--output",
                surface_csv,
                "--chains",
                *chains,
            ]
        )
        run_command(
            [
                "python3",
                "tests/diffmasif_replacement/compare_with_naccess.py",
                "--structure-kind",
                args.structure_kind,
                "--structure-id",
                structure_id,
                "--surface-points",
                surface_csv,
                "--output",
                comparison_json,
            ]
        )
        with open(comparison_json, "r") as handle:
            data = json.load(handle)
        summaries.append(
            {
                "structure_id": structure_id,
                **data["summary"],
                **data["proxy_meta"],
            }
        )

    summary_csv = os.path.join(args.output_dir, f"{args.structure_kind}_summary.csv")
    with open(summary_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(summaries[0].keys()))
        writer.writeheader()
        writer.writerows(summaries)
    print(summary_csv)


if __name__ == "__main__":
    main()
