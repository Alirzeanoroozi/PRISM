#!/usr/bin/env python3
"""Collect pair- and batch-level execution status without dropping failures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def available_keys(path: Path) -> tuple[set[str], set[tuple[str, str]]]:
    """Read current pair IDs or legacy receptor/ligand-only availability."""
    pair_ids: set[str] = set()
    receptor_ligand: set[tuple[str, str]] = set()
    if not path.exists():
        return pair_ids, receptor_ligand
    for row in read_rows(path):
        if row.get("pair_id"):
            pair_ids.add(row["pair_id"])
        if row.get("Receptor") and row.get("Ligand"):
            receptor_ligand.add((row["Receptor"].strip().lower(), row["Ligand"].strip().lower()))
    return pair_ids, receptor_ligand


def stage_counts(run_root: Path, pipeline: str) -> dict[str, int]:
    if pipeline == "current":
        roots = {
            "pdb": run_root / "processed/pdbs",
            "surface": run_root / "processed/surface_extraction",
            "alignment": run_root / "processed/alignment",
            "transformation": run_root / "processed/transformation",
            "refinement": run_root / "processed/rosetta_refinement/structures",
        }
    else:
        job = run_root / "jobs" / run_root.name
        roots = {
            "pdb": job / "pdb",
            "surface": job / "surfaceExtract",
            "alignment": job / "alignment",
            "transformation": job / "transformation",
            "refinement": job / "fiberdock",
        }
    return {name: len(list(path.glob("*"))) if path.exists() else 0 for name, path in roots.items()}


def reason_from_log(log: Path) -> str:
    if not log.exists():
        return "no_log"
    text = log.read_text(errors="replace")
    for marker in ("Traceback", "RuntimeError:", "OperationalError:", "command not found", "No such file"):
        if marker in text:
            line = next((line.strip() for line in reversed(text.splitlines()) if marker in line), marker)
            return line[:240]
    return "completed_marker_missing"


def collect(batch_root: Path, run_root: Path, output: Path, current_root: list[Path] | None = None, legacy_root: list[Path] | None = None) -> None:
    rows = []
    for batch_dir in sorted(batch_root.glob("batch_*")):
        batch_rows = read_rows(batch_dir / "inputs.csv")
        for pipeline in ("current", "legacy"):
            roots = current_root if pipeline == "current" and current_root else legacy_root if pipeline == "legacy" and legacy_root else [run_root / pipeline]
            candidates = [root for root in roots if (root / batch_dir.name).exists()]
            pipeline_base = next(
                (root for root in candidates if (root / batch_dir.name / "status/completed").exists()),
                candidates[0] if candidates else roots[0],
            )
            pipeline_dir = pipeline_base / batch_dir.name
            available_ids, available_pairs = available_keys(pipeline_dir / "status/available_inputs.csv")
            unavailable_rows = read_rows(pipeline_dir / "status/unavailable_inputs.csv") if (pipeline_dir / "status/unavailable_inputs.csv").exists() else []
            unavailable = {row["pair_id"]: row.get("missing_pdb_ids", "") for row in unavailable_rows}
            counts = stage_counts(pipeline_dir, pipeline)
            completed = (pipeline_dir / "status/completed").exists()
            log = pipeline_dir / "logs/pipeline.out"
            if completed:
                batch_status = "completed"
                batch_reason = "completed_marker"
            elif log.exists():
                batch_reason = reason_from_log(log)
                if batch_reason == "completed_marker_missing":
                    batch_status = "running_or_incomplete"
                else:
                    batch_status = "failed_execution"
            else:
                batch_status = "not_started"
                batch_reason = "no_log"
            for row in batch_rows:
                if row["pair_id"] in unavailable:
                    status, reason = "input_missing", unavailable[row["pair_id"]]
                elif row["pair_id"] not in available_ids and (row.get("Receptor", "").lower(), row.get("Ligand", "").lower()) not in available_pairs:
                    status, reason = "batch_input_unverified", "missing_available_input_record"
                else:
                    status, reason = batch_status, batch_reason
                rows.append({
                    **row,
                    "pipeline": pipeline,
                    "status": status,
                    "reason": reason,
                    **{f"count_{key}": value for key, value in counts.items()},
                })
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["pair_id", "pipeline", "status", "reason"]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} pipeline-pair rows to {output}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-root", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, default=None)
    parser.add_argument("--current-root", type=Path, action="append", default=[])
    parser.add_argument("--legacy-root", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.run_root is None and (not args.current_root or not args.legacy_root):
        parser.error("provide --run-root or both --current-root and --legacy-root")
    collect(args.batch_root, args.run_root, args.output, args.current_root, args.legacy_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
