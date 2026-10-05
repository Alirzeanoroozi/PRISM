#!/usr/bin/env python3
"""Create identical ten-pair input batches for current and legacy PRISM runs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.pdb_download import normalize_target_id


BENCHMARK_FILES = (
    ("rigid", "T_Rigid.csv"),
    ("medium", "T_medium.csv"),
    ("difficult", "T_difficult.csv"),
)


def load_benchmark_rows(data_dir: Path) -> list[dict[str, str]]:
    rows = []
    for benchmark_set, filename in BENCHMARK_FILES:
        with (data_dir / filename).open(newline="") as handle:
            for row_number, row in enumerate(csv.DictReader(handle), start=2):
                receptor = normalize_target_id(row.get("PDB ID 1", ""))
                ligand = normalize_target_id(row.get("PDB ID 2", ""))
                rows.append(
                    {
                        "pair_id": f"{benchmark_set}_{row_number:04d}",
                        "benchmark_set": benchmark_set,
                        "source_row": str(row_number),
                        "complex": (row.get("Complex") or "").strip(),
                        "pdb_id_1_raw": (row.get("PDB ID 1") or "").strip(),
                        "pdb_id_2_raw": (row.get("PDB ID 2") or "").strip(),
                        "Receptor": receptor,
                        "Ligand": ligand,
                    }
                )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def prepare_batches(data_dir: Path, output_root: Path, batch_size: int = 10) -> list[Path]:
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    rows = load_benchmark_rows(data_dir)
    manifest_fields = [
        "pair_id", "benchmark_set", "source_row", "complex",
        "pdb_id_1_raw", "pdb_id_2_raw", "Receptor", "Ligand",
    ]
    write_csv(output_root / "shared_manifest.csv", rows, manifest_fields)

    batch_dirs = []
    for start in range(0, len(rows), batch_size):
        batch_rows = rows[start : start + batch_size]
        batch_dir = output_root / f"batch_{start // batch_size + 1:04d}"
        write_csv(batch_dir / "inputs.csv", batch_rows, manifest_fields)
        with (batch_dir / "pair_list").open("w") as handle:
            for row in batch_rows:
                handle.write(f"{row['Receptor']} {row['Ligand']}\n")
        batch_dirs.append(batch_dir)
    return batch_dirs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "benchmark/data")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=10)
    args = parser.parse_args()
    batches = prepare_batches(args.data_dir, args.output_root, args.batch_size)
    print(f"wrote {sum(1 for _ in (args.output_root / 'shared_manifest.csv').open()) - 1} pairs")
    print(f"wrote {len(batches)} batches of at most {args.batch_size} pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
