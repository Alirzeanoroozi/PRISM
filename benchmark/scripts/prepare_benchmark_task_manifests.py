#!/usr/bin/env python3
"""Prepare isolated ten-task (or final partial) manifests for the 257 rows.

The existing benchmark CSVs are the sole row source.  Each generated array
task owns one ``dataset_row_id`` and writes source/validation artifacts below
its task-local directory.  No pair is represented by a normalized PDB-pair
key.
"""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark.scripts.build_investigation_source_manifest import _read_rows
from benchmark.scripts.prepare_investigation_tasks import FIELDS, build_task_rows


def prepare_manifests(
    repo_root: str | Path,
    output_dir: str | Path,
    *,
    python_bin: str = "python3",
    chunk_size: int = 10,
) -> list[dict[str, str]]:
    if not 1 <= chunk_size <= 10:
        raise ValueError("chunk_size must be between 1 and 10")
    root = Path(repo_root).resolve()
    output = Path(output_dir).resolve()
    rows = _read_rows(root, limit=None)
    row_ids = [row["dataset_row_id"] for row in rows]
    resolved_python = Path(python_bin).expanduser()
    if not resolved_python.is_absolute():
        resolved_python = Path(shutil.which(python_bin) or python_bin)
    resolved_python = resolved_python.resolve()
    runner_path = root / "benchmark/scripts/isolated_kutem_runner.py"
    template_path = root / "benchmark/jobs/isolated_kutem_array.sbatch"
    input_paths = [
        str(root / "benchmark/data/T_Rigid.csv"),
        str(root / "benchmark/data/T_medium.csv"),
        str(root / "benchmark/data/T_difficult.csv"),
        str(root / "benchmark/originals/benchmark5.5.tgz"),
        str(root / "benchmark/scripts/build_investigation_source_manifest.py"),
        str(root / "benchmark/scripts/stage_curated_sources.py"),
        str(root / "benchmark/scripts/run_source_gate_task.py"),
        str(runner_path),
        str(template_path),
        str(resolved_python),
    ]
    command_template = (
        f"{shlex.quote(python_bin)} {shlex.quote(str(root / 'benchmark/scripts/run_source_gate_task.py'))} "
        f"--repo-root {shlex.quote(str(root))} --output-dir \"$PRISM_TASK_DIR/source\" --strict "
        "--dataset-row-id {dataset_row_id}"
    )
    records: list[dict[str, str]] = []
    for batch_number, start in enumerate(range(0, len(row_ids), chunk_size), start=1):
        batch_ids = row_ids[start : start + chunk_size]
        batch_dir = output / f"batch_{batch_number:03d}"
        manifest_path = batch_dir / "task_manifest.csv"
        batch_dir.mkdir(parents=True, exist_ok=True)
        task_rows = build_task_rows(
            batch_ids,
            command_template,
            input_paths=input_paths,
            output_paths=[
                "source/source_manifest.tsv",
                "source/structure_validation.tsv",
                "source/staged/staged_sources.tsv",
                "source/source_gate_summary.json",
            ],
            scientific_retry_id=f"source-gate-batch-{batch_number:03d}",
        )
        with manifest_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(task_rows)
        records.append(
            {
                "batch_id": f"source-gate-{batch_number:03d}",
                "task_manifest": str(manifest_path),
                "array_size": str(len(batch_ids)),
                "first_dataset_row_id": batch_ids[0],
                "last_dataset_row_id": batch_ids[-1],
                "task_count": str(len(batch_ids)),
            }
        )
    index_path = output / "batch_index.tsv"
    with index_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    return records


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--python-bin", default="python3")
    parser.add_argument("--chunk-size", type=int, default=10)
    args = parser.parse_args(argv)
    try:
        records = prepare_manifests(
            args.repo_root,
            args.output_dir,
            python_bin=args.python_bin,
            chunk_size=args.chunk_size,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"wrote {len(records)} isolated benchmark batch manifests to {Path(args.output_dir).resolve()}")
    print(f"dataset_rows={sum(int(record['task_count']) for record in records)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
