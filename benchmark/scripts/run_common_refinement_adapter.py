#!/usr/bin/env python3
"""Adapt compact candidate manifests to the validated common refinement worker.

GTalign retains assembled two-chain models rather than separate transformed
halves.  This adapter splits those models into run-scoped, resumable halves,
then invokes the same FiberDock/external-Rosetta/DockQ worker used by the
validated common refinement run.  Split inputs and one-row CSVs are removed
after the worker returns; the checkpoint retains source hashes and paths.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def split_combined_model(row: dict[str, str], split_root: Path) -> tuple[Path, Path]:
    model = Path(row["model"]).resolve()
    left_chain = "".join(str(row.get("chain_left", "")).split())
    right_chain = "".join(str(row.get("chain_right", "")).split())
    if not model.is_file() or not left_chain or not right_chain:
        raise ValueError("GTalign row needs an existing model and chain_left/chain_right")
    split_root.mkdir(parents=True, exist_ok=True)
    token = str(row.get("source_row", "0"))
    left = split_root / f"gtalign_source_{token}_L.pdb"
    right = split_root / f"gtalign_source_{token}_R.pdb"
    outputs = ((left, set(left_chain)), (right, set(right_chain)))
    handles = [(path.open("w", encoding="ascii"), chains) for path, chains in outputs]
    try:
        with model.open(encoding="ascii", errors="replace") as source:
            for line in source:
                if not line.startswith(("ATOM", "HETATM")) or len(line) <= 21:
                    continue
                chain = line[21].strip() or "_"
                for handle, chains in handles:
                    if chain in chains:
                        handle.write(line)
        for handle, _ in handles:
            handle.write("TER\nEND\n")
    finally:
        for handle, _ in handles:
            handle.close()
    if not all(path.is_file() and path.stat().st_size > 0 for path, _ in outputs):
        raise ValueError(f"model did not contain both requested chains: {model}")
    return left, right


def read_row(path: Path, index: int) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if index < 0 or index >= len(rows):
        raise IndexError(f"candidate index {index} outside 0..{len(rows) - 1}")
    return dict(rows[index])


def run_one(args: argparse.Namespace, index: int) -> int:
    row = read_row(args.selected_csv, index)
    temporary_root = args.comparison_root / "adapter_inputs" / f"task_{index:08d}"
    temporary_root.mkdir(parents=True, exist_ok=True)
    one_row = temporary_root / "candidate.csv"
    try:
        if row.get("model"):
            left, right = split_combined_model(row, temporary_root / "split")
            row["left"] = str(left)
            row["right"] = str(right)
        if not row.get("left") or not row.get("right"):
            raise ValueError("candidate row needs left/right or model")
        with one_row.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=sorted(row))
            writer.writeheader()
            writer.writerow(row)
        command = [
            sys.executable,
            str(args.worker_script.resolve()),
            "--comparison-root", str(args.comparison_root.resolve()),
            "--selected-csv", str(one_row.resolve()),
            "--index", "0",
            "--pipeline-repo", str(args.pipeline_repo.resolve()),
            "--fiberdock-dir", str(args.fiberdock_dir.resolve()),
            "--dockq-python", str(args.dockq_python),
            "--rosetta-prepack", str(args.rosetta_prepack),
            "--rosetta-dock", str(args.rosetta_dock),
            "--rosetta-db", str(args.rosetta_db),
        ]
        completed = subprocess.run(command, check=False)
        return completed.returncode
    finally:
        # The worker has checkpointed source hashes and its stage paths.  These
        # adapter-only copies are reproducibly regenerated on resume.
        for path in sorted(temporary_root.glob("*"), reverse=True):
            if path.is_file():
                path.unlink()
        for path in sorted((temporary_root / "split").glob("*") if (temporary_root / "split").is_dir() else [], reverse=True):
            if path.is_file():
                path.unlink()
        for directory in (temporary_root / "split", temporary_root):
            try:
                directory.rmdir()
            except OSError:
                pass


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-root", type=Path, required=True)
    parser.add_argument("--selected-csv", type=Path, required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True, help="exclusive candidate index")
    parser.add_argument("--worker-script", type=Path, required=True)
    parser.add_argument("--pipeline-repo", type=Path, required=True)
    parser.add_argument("--fiberdock-dir", type=Path, required=True)
    parser.add_argument("--dockq-python", required=True)
    parser.add_argument("--rosetta-prepack", required=True)
    parser.add_argument("--rosetta-dock", required=True)
    parser.add_argument("--rosetta-db", required=True)
    args = parser.parse_args()
    if args.start < 0 or args.end < args.start:
        parser.error("invalid [start,end) range")
    codes = {str(index): run_one(args, index) for index in range(args.start, args.end)}
    shard_id = os.environ.get("SLURM_ARRAY_TASK_ID", f"{args.start}_{args.end}")
    status_path = args.comparison_root / "shard_status" / f"shard_{shard_id}.json"
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status_path.write_text(
        json.dumps(
            {
                "schema_version": "prism-common-refinement-shard/v1",
                "start": args.start,
                "end": args.end,
                "candidate_count": max(0, args.end - args.start),
                "worker_return_codes": codes,
                "finished_at": datetime.now(timezone.utc).isoformat(),
                "failure_tolerant": True,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
