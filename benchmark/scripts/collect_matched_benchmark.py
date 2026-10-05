#!/usr/bin/env python3
"""Collect matched PRISM task evidence into fail-closed summaries.

The collector never turns a missing task or invalid score into a successful
prediction.  It accepts task-local ``exit.json`` records plus optional
``scores_global.tsv``/``scores_interfaces.tsv`` files emitted by the scoring
stage and writes stable pair-level summaries.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


SUMMARY_FIELDS = (
    "dataset_row_id", "arm", "attempted_tasks", "completed_tasks", "failed_tasks",
    "attempted_predictions", "scoreable_predictions", "completed_predictions", "failed_predictions", "prediction_success_rate",
    "best_GlobalDockQ_at_20", "mean_GlobalDockQ", "variance_GlobalDockQ",
    "best_iRMSD_at_20", "mean_iRMSD", "variance_iRMSD", "interface_rows",
    "best_TM_score_at_20", "mean_TM_score", "runtime_seconds", "throughput_predictions_per_second",
)
FAILURE_FIELDS = ("task_id", "dataset_row_id", "arm", "experiment", "status", "reason", "task_root")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _number(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _variance(values: list[float]) -> float | None:
    if not values:
        return None
    mean = sum(values) / len(values)
    return sum((value - mean) ** 2 for value in values) / len(values)


def _task_root(raw: str, manifest_path: Path) -> Path:
    path = Path(raw)
    return path.resolve() if path.is_absolute() else (manifest_path.parent / path).resolve()


def _pdb_is_valid(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        ca_by_chain: dict[str, int] = {}
        with path.open(errors="replace") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")):
                    chain = line[21].strip() or "_"
                    if line[12:16].strip() == "CA":
                        ca_by_chain[chain] = ca_by_chain.get(chain, 0) + 1
        # A refined complex must retain at least two structurally meaningful
        # partner chains; tiny one-residue remnants indicate chain/output
        # corruption even when the file is non-empty.
        return len(ca_by_chain) >= 2 and min(ca_by_chain.values()) >= 3
    except OSError:
        return False


def _artifact_failure(task: dict[str, str], root: Path, exit_data: dict[str, Any]) -> str | None:
    """Return a deterministic failure reason for declared task artifacts."""
    output_pdb = exit_data.get("output_pdb") or task.get("output_pdb")
    if output_pdb:
        output = Path(str(output_pdb))
        if not output.is_absolute():
            output = root / output
        if not _pdb_is_valid(output):
            return "invalid_or_missing_output_pdb"
    if task.get("experiment") == "strict_scoring":
        if not (root / "scores_global.tsv").is_file():
            return "missing_scores_global"
    return None


def collect(task_manifest: str | Path, output_dir: str | Path) -> dict[str, Path]:
    manifest_path = Path(task_manifest).resolve()
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    tasks = _read_tsv(manifest_path)
    if not tasks:
        raise ValueError(f"task manifest is empty or missing: {manifest_path}")
    seen: set[str] = set()
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    failures: list[dict[str, str]] = []
    for task in tasks:
        task_id = task.get("task_id", "")
        if not task_id or task_id in seen:
            raise ValueError(f"duplicate or empty task_id: {task_id!r}")
        seen.add(task_id)
        root = _task_root(task.get("output_root", ""), manifest_path)
        exit_path = root / "exit.json"
        exit_data: dict[str, Any] = {}
        if exit_path.is_file():
            try:
                exit_data = json.loads(exit_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as exc:
                failures.append({"task_id": task_id, "dataset_row_id": task.get("dataset_row_id", ""), "arm": task.get("arm", ""), "experiment": task.get("experiment", ""), "status": "failed", "reason": f"invalid_exit_json:{exc}", "task_root": str(root)})
        else:
            failures.append({"task_id": task_id, "dataset_row_id": task.get("dataset_row_id", ""), "arm": task.get("arm", ""), "experiment": task.get("experiment", ""), "status": "failed", "reason": "missing_exit_record", "task_root": str(root)})
        rc = exit_data.get("return_code")
        scientific_status = str(exit_data.get("scientific_status", ""))
        completed = rc == 0 and scientific_status not in {"failed", "blocked", "unavailable", "timeout", "completed_no_predictions"}
        if not completed and exit_path.is_file():
            failures.append({"task_id": task_id, "dataset_row_id": task.get("dataset_row_id", ""), "arm": task.get("arm", ""), "experiment": task.get("experiment", ""), "status": "failed", "reason": scientific_status or f"return_code:{rc}", "task_root": str(root)})
        artifact_failure = _artifact_failure(task, root, exit_data) if completed else None
        if artifact_failure:
            completed = False
            failures.append({"task_id": task_id, "dataset_row_id": task.get("dataset_row_id", ""), "arm": task.get("arm", ""), "experiment": task.get("experiment", ""), "status": "failed", "reason": artifact_failure, "task_root": str(root)})
        global_rows = _read_tsv(root / "scores_global.tsv")
        interface_rows = _read_tsv(root / "scores_interfaces.tsv")
        for row in global_rows:
            row["task_id"] = task_id
            row["dataset_row_id"] = task.get("dataset_row_id", "")
            row["arm"] = task.get("arm", "")
        key = (task.get("dataset_row_id", ""), task.get("arm", ""))
        groups.setdefault(key, []).append({"task": task, "root": root, "exit": exit_data, "completed": completed, "global": global_rows, "interfaces": interface_rows})

    summaries: list[dict[str, Any]] = []
    all_global: list[dict[str, Any]] = []
    all_interfaces: list[dict[str, Any]] = []
    for (dataset_row_id, arm), records in sorted(groups.items()):
        completed_tasks = sum(bool(record["completed"]) for record in records)
        failed_tasks = len(records) - completed_tasks
        global_rows = [row for record in records for row in record["global"]]
        interface_count = sum(len(record["interfaces"]) for record in records)
        dockq = [value for row in global_rows if (value := _number(row.get("GlobalDockQ"))) is not None]
        irmsd = [value for row in global_rows if (value := _number(row.get("iRMSD"))) is not None]
        tm_scores = [value for row in global_rows for raw in (row.get("TM-score"), row.get("TM_score"), row.get("tm_score")) if (value := _number(raw)) is not None]
        runtimes = [value for record in records if (value := _number(record["exit"].get("elapsed_seconds"))) is not None]
        runtime = sum(runtimes) if runtimes else None
        prediction_count = len(global_rows)
        scoreable = len(dockq)
        summaries.append({
            "dataset_row_id": dataset_row_id, "arm": arm,
            "attempted_tasks": len(records), "completed_tasks": completed_tasks, "failed_tasks": failed_tasks,
            "attempted_predictions": prediction_count, "scoreable_predictions": scoreable,
            "completed_predictions": scoreable, "failed_predictions": prediction_count - scoreable,
            "prediction_success_rate": scoreable / prediction_count if prediction_count else None,
            "best_GlobalDockQ_at_20": max(dockq) if dockq else 0.0,
            "mean_GlobalDockQ": sum(dockq) / len(dockq) if dockq else None,
            "variance_GlobalDockQ": _variance(dockq),
            "best_iRMSD_at_20": min(irmsd) if irmsd else None,
            "mean_iRMSD": sum(irmsd) / len(irmsd) if irmsd else None,
            "variance_iRMSD": _variance(irmsd), "interface_rows": interface_count,
            "best_TM_score_at_20": max(tm_scores) if tm_scores else None,
            "mean_TM_score": sum(tm_scores) / len(tm_scores) if tm_scores else None,
            "runtime_seconds": runtime,
            "throughput_predictions_per_second": prediction_count / runtime if runtime and prediction_count else None,
        })
        all_global.extend(global_rows)
        all_interfaces.extend({**row, "dataset_row_id": dataset_row_id, "arm": arm} for record in records for row in record["interfaces"])

    def write_rows(path: Path, fields: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    write_rows(output / "pair_summary.tsv", SUMMARY_FIELDS, summaries)
    write_rows(output / "failures.tsv", FAILURE_FIELDS, failures)
    write_rows(output / "scores_global.tsv", tuple(sorted({key for row in all_global for key in row})), all_global)
    write_rows(output / "scores_interfaces.tsv", tuple(sorted({key for row in all_interfaces for key in row})), all_interfaces)
    report = output / "comparison_report.md"
    total_tasks = len(tasks)
    total_predictions = sum(int(row["attempted_predictions"]) for row in summaries)
    total_scoreable = sum(int(row["completed_predictions"]) for row in summaries)
    report.write_text(
        "# Matched PRISM benchmark collection\n\n"
        "## Confirmed findings\n\n"
        f"- Tasks represented: {total_tasks}; task failures recorded: {len(failures)}.\n"
        f"- Prediction rows represented: {total_predictions}; completed/scoreable rows: {total_scoreable}; failed rows: {total_predictions - total_scoreable}.\n"
        "- iRMSD best values use the minimum; missing structural metrics remain null.\n\n"
        "## Likely explanations\n\n"
        "- None inferred by the collector; stage-specific evidence must be reviewed from failures.tsv and task logs.\n\n"
        "## Unresolved issues\n\n"
        "- A task is not scientific success unless its exit record, outputs, and score hashes reconcile.\n"
        "- This report is observational until the source/evaluator gates and shared-scoreability conditions are satisfied.\n",
        encoding="utf-8",
    )
    return {"pair_summary": output / "pair_summary.tsv", "failures": output / "failures.tsv", "scores_global": output / "scores_global.tsv", "scores_interfaces": output / "scores_interfaces.tsv", "report": report}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        paths = collect(args.task_manifest, args.output_dir)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    print(json.dumps({key: str(value) for key, value in paths.items()}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
