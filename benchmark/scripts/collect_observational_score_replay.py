#!/usr/bin/env python3
"""Collect score shards and compare them with the prior observational report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def numeric(rows: list[dict[str, str]], field: str) -> list[float]:
    values = []
    for row in rows:
        try:
            value = float(row.get(field, ""))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
    return values


def metrics(rows: list[dict[str, str]]) -> dict[str, object]:
    result = {"model_rows": len(rows), "ready_rows": sum(row.get("status") == "ready" for row in rows)}
    for field in ("dockq", "irmsd"):
        values = numeric(rows, field)
        result.update({
            f"{field}_n": len(values),
            f"{field}_mean": statistics.mean(values) if values else None,
            f"{field}_median": statistics.median(values) if values else None,
            f"{field}_best": (min(values) if field == "irmsd" else max(values)) if values else None,
        })
    result["score_errors"] = sum(bool(row.get("score_error")) for row in rows)
    result["pairs_with_score"] = len({row.get("pair_id") for row in rows if row.get("pair_id") and (row.get("dockq") or row.get("irmsd"))})
    return result


def validate_tasks(replay_root: Path) -> tuple[list[dict[str, str]], list[Path]]:
    """Fail closed unless the prepared shards and successful task outputs agree."""
    manifest_path = replay_root / "replay_manifest.json"
    if not manifest_path.is_file():
        raise ValueError(f"missing replay manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    shard_records = manifest.get("shards", [])
    expected_ids = {int(record["shard"]) for record in shard_records}
    if expected_ids != set(range(1, int(manifest.get("shard_count", 0)) + 1)):
        raise ValueError("replay manifest has incomplete or duplicate shard IDs")
    tasks_root = replay_root / "tasks"
    task_dirs = sorted(tasks_root.glob("task-*"))
    actual_ids = set()
    for task_dir in task_dirs:
        try:
            actual_ids.add(int(task_dir.name.removeprefix("task-")))
        except ValueError as exc:
            raise ValueError(f"invalid task directory: {task_dir}") from exc
    if actual_ids != expected_ids:
        raise ValueError(f"task IDs {sorted(actual_ids)} do not match manifest {sorted(expected_ids)}")

    rows: list[dict[str, str]] = []
    scored_paths: list[Path] = []
    record_by_id = {int(record["shard"]): record for record in shard_records}
    for task_id in sorted(expected_ids):
        task_dir = tasks_root / f"task-{task_id}"
        exit_path = task_dir / "exit.json"
        scored_path = task_dir / "scored_models.csv"
        if not exit_path.is_file() or not scored_path.is_file():
            raise ValueError(f"missing terminal task artifacts for task {task_id}")
        exit_record = json.loads(exit_path.read_text())
        if int(exit_record.get("array_task_id", -1)) != task_id:
            raise ValueError(f"task {task_id} has mismatched array_task_id")
        if int(exit_record.get("return_code", 1)) != 0 or exit_record.get("scientific_status") != "completed":
            raise ValueError(f"task {task_id} is not a successful completed task")
        expected_shard = (replay_root / "shards" / f"shard_{task_id:02d}.csv").resolve()
        if Path(exit_record.get("input", "")).resolve() != expected_shard:
            raise ValueError(f"task {task_id} input does not match declared shard")
        if Path(exit_record.get("output", "")).resolve() != scored_path.resolve():
            raise ValueError(f"task {task_id} output does not match task directory")
        declared_shard = record_by_id[task_id]
        if sha256(expected_shard) != declared_shard.get("sha256"):
            raise ValueError(f"shard hash changed after preparation: {expected_shard}")
        output_hash = exit_record.get("output_sha256")
        if not output_hash or sha256(scored_path) != output_hash:
            raise ValueError(f"task {task_id} output hash is missing or changed")
        shard_rows = read_rows(scored_path)
        if len(shard_rows) != int(declared_shard.get("rows", -1)):
            raise ValueError(f"task {task_id} row count does not match prepared shard")
        for row in shard_rows:
            raw_path_value = row.get("dockq_raw_json_path", "")
            raw_hash = row.get("dockq_raw_json_sha256", "")
            if raw_path_value or raw_hash:
                raw_path = Path(raw_path_value).resolve()
                if not raw_path.is_relative_to(task_dir.resolve()):
                    raise ValueError(f"task {task_id} raw DockQ JSON escapes task directory")
                if not raw_path.is_file() or not raw_hash or sha256(raw_path) != raw_hash:
                    raise ValueError(f"task {task_id} raw DockQ JSON is missing or changed")
        rows.extend(shard_rows)
        scored_paths.append(scored_path)

    model_manifest = Path(manifest["model_manifest"])
    expected_models = read_rows(model_manifest)
    expected_paths = {row.get("model_path") for row in expected_models if row.get("model_path")}
    actual_paths = {row.get("model_path") for row in rows if row.get("model_path")}
    if len(actual_paths) != len([row for row in rows if row.get("model_path")]):
        raise ValueError("duplicate model_path rows in replay outputs")
    if actual_paths != expected_paths:
        raise ValueError("replay output model coverage differs from prepared model manifest")
    return rows, scored_paths


def collect(replay_root: Path, previous_root: Path, output_root: Path) -> dict:
    rows, shard_paths = validate_tasks(replay_root)
    rows.sort(key=lambda row: (row.get("pipeline", ""), row.get("model_path", "")))
    output_root.mkdir(parents=True, exist_ok=False)
    combined = output_root / "scored_models.csv"
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with combined.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    by_pipeline = {}
    for pipeline in sorted({row.get("pipeline", "") for row in rows}):
        by_pipeline[pipeline] = metrics([row for row in rows if row.get("pipeline") == pipeline])
    previous_rows = read_rows(previous_root / "scored_models.csv")
    previous_by_pipeline = {
        pipeline: metrics([row for row in previous_rows if row.get("pipeline") == pipeline])
        for pipeline in sorted({row.get("pipeline", "") for row in previous_rows})
    }
    current_pairs = {
        row.get("pair_id") for row in rows if row.get("pair_id") and (row.get("dockq") or row.get("irmsd"))
    }
    previous_pairs = {
        row.get("pair_id") for row in previous_rows if row.get("pair_id") and (row.get("dockq") or row.get("irmsd"))
    }
    summary = {
        "schema_version": "observational-score-replay-summary/v1",
        "status": "observational_replay",
        "rows": len(rows),
        "shards": len(shard_paths),
        "pipeline_metrics": by_pipeline,
        "previous_report": str((previous_root / "FINAL_REPORT.md").resolve()),
        "previous_scored_models_sha256": sha256(previous_root / "scored_models.csv"),
        "previous_pipeline_metrics": previous_by_pipeline,
        "scoreable_pair_intersection": sorted(current_pairs & previous_pairs),
        "scoreable_pair_intersection_count": len(current_pairs & previous_pairs),
        "scoreable_pair_current_only": sorted(current_pairs - previous_pairs),
        "scoreable_pair_previous_only": sorted(previous_pairs - current_pairs),
        "combined_scores_sha256": sha256(combined),
    }
    (output_root / "replay_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    lines = [
        "# Observational score replay versus previous report",
        "",
        "This is a fresh fail-closed re-score of existing model outputs; it is not a new model-generation run.",
        "",
        f"- Replay shards collected: {len(shard_paths)}",
        f"- Replay model rows: {len(rows)}",
        f"- Scoreable pair intersection with previous replay: {summary['scoreable_pair_intersection_count']}",
        "",
        "| pipeline | replay rows | replay DockQ n | replay DockQ mean | prior rows | prior DockQ n | prior DockQ mean |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for pipeline in sorted(set(by_pipeline) | set(previous_by_pipeline)):
        now, old = by_pipeline.get(pipeline, {}), previous_by_pipeline.get(pipeline, {})
        lines.append(
            f"| {pipeline} | {now.get('model_rows', 0)} | {now.get('dockq_n', 0)} | {now.get('dockq_mean')} | "
            f"{old.get('model_rows', 0)} | {old.get('dockq_n', 0)} | {old.get('dockq_mean')} |"
        )
    lines += ["", "The previous aggregate means are retained as observational context only.", ""]
    (output_root / "COMPARISON.md").write_text("\n".join(lines))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay-root", type=Path, required=True)
    parser.add_argument("--previous-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(collect(args.replay_root, args.previous_root, args.output_root), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
