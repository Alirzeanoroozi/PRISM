#!/usr/bin/env python3
"""Collect compact PRISM case/stage results without copying raw artifacts.

The collector reads a submission manifest, corrected per-aligner manifests,
per-case task ledgers, and the selected attempt run summaries.  It writes only
small CSV/JSON evidence files to a new output root.  It intentionally does not
walk or copy alignment JSON/PDB trees.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "prism-compact-results-v1"
DEFAULT_ALIGNERS = ("tmalign", "multiprot", "usalign")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def write_json(path: Path, payload: Any) -> None:
    atomic_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    temporary.replace(path)


def json_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def number(value: Any, default: Any = "") -> Any:
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else default


def read_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else None


def stage_rows(summary: dict[str, Any], pipeline: str, index: int, case_id: str) -> list[dict[str, Any]]:
    events = summary.get("stage_events") or []
    opened: dict[str, str] = {}
    rows: list[dict[str, Any]] = []
    for event in events:
        stage = str(event.get("stage", ""))
        event_name = str(event.get("event", ""))
        timestamp = str(event.get("timestamp", ""))
        if event_name == "started":
            opened[stage] = timestamp
            rows.append({
                "pipeline": pipeline,
                "case_index": index,
                "case_id": case_id,
                "stage": stage,
                "status": "started",
                "return_code": "",
                "started_at": timestamp,
                "completed_at": "",
                "wall_seconds": "",
                "detail": event.get("detail", ""),
            })
        elif event_name in {"completed", "failed", "skipped"}:
            started = opened.get(stage, "")
            elapsed = ""
            if started and timestamp:
                try:
                    elapsed = (datetime.fromisoformat(timestamp) - datetime.fromisoformat(started)).total_seconds()
                except ValueError:
                    elapsed = ""
            rows.append({
                "pipeline": pipeline,
                "case_index": index,
                "case_id": case_id,
                "stage": stage,
                "status": event_name,
                "return_code": event.get("return_code", ""),
                "started_at": started,
                "completed_at": timestamp,
                "wall_seconds": elapsed,
                "detail": event.get("detail", ""),
            })
            opened.pop(stage, None)
    return rows


def stage_statuses(summary: dict[str, Any]) -> dict[str, str]:
    statuses: dict[str, str] = {}
    for event in summary.get("stage_events") or []:
        stage = str(event.get("stage", ""))
        if stage:
            statuses[stage] = str(event.get("event", ""))
    return statuses


def find_manifest(submission: dict[str, Any], run_root: Path, aligner: str) -> Path:
    manifests = submission.get("manifests") or {}
    value = manifests.get(aligner)
    if value:
        candidate = Path(value)
        if candidate.is_file():
            return candidate
    candidate = run_root / aligner / "batch_manifest_corrected.json"
    if candidate.is_file():
        return candidate
    raise FileNotFoundError(f"no corrected batch manifest for {aligner}")


def status_row(
    pipeline: str,
    entry: dict[str, Any],
    task: dict[str, Any],
    summary: dict[str, Any] | None,
    task_path: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    index = int(entry["index"])
    case_id = str(entry["case_id"])
    attempt_root = Path(str(task.get("attempt_root", ""))) if task.get("attempt_root") else None
    summary_path = attempt_root / "run_summary.json" if attempt_root else None
    summary_hash_actual = sha256_file(summary_path) if summary_path and summary_path.is_file() else ""
    errors = task.get("validation_errors") or []
    effective = summary or {}
    observed_records = number(effective.get("alignment_records"), number(task.get("alignment_records")))
    successful = number(effective.get("successful_alignment_records"), number(task.get("successful_alignment_records")))
    failed = number(effective.get("failed_alignment_records"), number(task.get("failed_alignment_records")))
    transformed = number(effective.get("transformed_model_pairs"), number(task.get("transformed_model_pairs")))
    conservation = ""
    if all(isinstance(value, (int, float)) for value in (observed_records, successful, failed)):
        conservation = str(successful + failed == observed_records).lower()
    summary_hash_claimed = str(task.get("summary_sha256", ""))
    summary_present = bool(summary_path and summary_path.is_file())
    summary_hash_match = (
        str(summary_hash_actual == summary_hash_claimed).lower()
        if summary_hash_claimed
        else ""
    )
    stage = stage_statuses(effective)
    row = {
        "pipeline": pipeline,
        "dataset": entry.get("dataset", "bm55_full"),
        "split": entry.get("split", ""),
        "case_index": index,
        "case_id": case_id,
        "receptor": entry.get("receptor", ""),
        "ligand": entry.get("ligand", ""),
        "status": task.get("status", "missing"),
        "evidence_status": "complete_summary" if summary_present else "summary_missing",
        "attempt": task.get("attempt", ""),
        "job_id": task.get("job_id", effective.get("job_id", "")),
        "return_code": effective.get("exit_code", task.get("return_code", "")),
        "elapsed_seconds": task.get("elapsed_seconds", ""),
        "wall_seconds": effective.get("wall_seconds", ""),
        "expected_alignment_records": entry.get("expected_alignment_records", ""),
        "alignment_records": observed_records,
        "successful_alignment_records": successful,
        "failed_alignment_records": failed,
        "transformed_model_pairs": transformed,
        "conservation_success_plus_failed_equals_records": conservation,
        "template_count": effective.get("template_count", ""),
        "template_sha256": effective.get("template_sha256", ""),
        "rank": effective.get("rank", ""),
        "rank_method": effective.get("rank_method", ""),
        "top_k": effective.get("top_k", ""),
        "host": effective.get("host", ""),
        "run_id": effective.get("run_id", ""),
        "summary_sha256_claimed": summary_hash_claimed,
        "summary_sha256_actual": summary_hash_actual,
        "summary_sha256_match": summary_hash_match,
        "validation_errors_json": json_text(errors),
        "stage_status_json": json_text(stage),
        "stage_wall_seconds_json": json_text({
            stage_name: next((r["wall_seconds"] for r in stage_rows(effective, pipeline, index, case_id)
                              if r["stage"] == stage_name and r["status"] in {"completed", "failed", "skipped"}), "")
            for stage_name in stage
        }),
        "command_json": json_text(effective.get("command", [])),
        "task_status_path": str(task_path),
        "attempt_root": str(attempt_root or ""),
        "summary_path": str(summary_path or ""),
        "pipeline_log": effective.get("pipeline_log", ""),
        "updated_at": task.get("updated_at", ""),
    }
    return row, stage_rows(effective, pipeline, index, case_id), {
        "summary_present": summary_present,
        "summary_hash_match": summary_hash_match == "true" or not summary_hash_claimed,
        "conservation": conservation in {"true", ""},
    }


def compact_fieldnames() -> list[str]:
    return [
        "pipeline", "dataset", "split", "case_index", "case_id", "receptor", "ligand",
        "status", "evidence_status", "attempt", "job_id", "return_code", "elapsed_seconds",
        "wall_seconds", "expected_alignment_records", "alignment_records",
        "successful_alignment_records", "failed_alignment_records", "transformed_model_pairs",
        "conservation_success_plus_failed_equals_records", "template_count", "template_sha256",
        "rank", "rank_method", "top_k", "host", "run_id", "summary_sha256_claimed",
        "summary_sha256_actual", "summary_sha256_match", "validation_errors_json",
        "stage_status_json", "stage_wall_seconds_json", "command_json", "task_status_path",
        "attempt_root", "summary_path", "pipeline_log", "updated_at",
    ]


def collect_pipeline(
    pipeline: str,
    batch_manifest_path: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    manifest = read_json(batch_manifest_path)
    if not manifest:
        raise RuntimeError(f"invalid batch manifest: {batch_manifest_path}")
    case_rows: list[dict[str, Any]] = []
    timing_rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    checks: list[dict[str, Any]] = []
    entries = sorted(manifest.get("entries", []), key=lambda item: int(item["index"]))
    for entry in entries:
        task_path = Path(entry["case_root"]) / "task_status.json"
        task = read_json(task_path) or {"status": "missing", "index": entry["index"], "case_id": entry["case_id"]}
        attempt_root = Path(str(task.get("attempt_root", ""))) if task.get("attempt_root") else None
        summary_path = attempt_root / "run_summary.json" if attempt_root else None
        summary = read_json(summary_path) if summary_path else None
        row, stages, row_checks = status_row(pipeline, entry, task, summary, task_path)
        case_rows.append(row)
        timing_rows.extend(stages)
        checks.append(row_checks)
        if row["status"] != "completed" or row["evidence_status"] != "complete_summary" or row["validation_errors_json"] != "[]" or row["conservation_success_plus_failed_equals_records"] != "true":
            failures.append({
                "pipeline": pipeline,
                "case_index": row["case_index"],
                "case_id": row["case_id"],
                "status": row["status"],
                "evidence_status": row["evidence_status"],
                "validation_errors_json": row["validation_errors_json"],
                "conservation": row["conservation_success_plus_failed_equals_records"],
                "task_status_path": row["task_status_path"],
                "summary_path": row["summary_path"],
            })
    return case_rows, timing_rows, {
        "batch_manifest": str(batch_manifest_path),
        "batch_manifest_sha256": sha256_file(batch_manifest_path),
        "aligner": manifest.get("aligner", pipeline),
        "dataset": manifest.get("dataset", ""),
        "template_count": manifest.get("template_count", ""),
        "template_sha256": manifest.get("template_sha256", ""),
        "expected_alignment_records": manifest.get("expected_alignment_records", ""),
        "entry_count": len(entries),
        "entries": entries,
        "row_checks": checks,
    }, failures


def comparison_row(pipeline: str, rows: list[dict[str, Any]], manifest: dict[str, Any]) -> dict[str, Any]:
    def total(field: str) -> float:
        return sum(float(row[field]) for row in rows if isinstance(row.get(field), (int, float)))

    completed = [row for row in rows if row["status"] == "completed"]
    records = total("alignment_records")
    success = total("successful_alignment_records")
    failed = total("failed_alignment_records")
    return {
        "pipeline": pipeline,
        "case_count": len(rows),
        "completed_cases": len(completed),
        "failed_cases": sum(row["status"] == "failed" for row in rows),
        "not_started_cases": sum(row["status"] == "not_started" for row in rows),
        "summary_missing_cases": sum(row["evidence_status"] != "complete_summary" for row in rows),
        "expected_alignment_records": manifest.get("expected_alignment_records", ""),
        "observed_alignment_records": int(records) if records.is_integer() else records,
        "successful_alignment_records": int(success) if success.is_integer() else success,
        "failed_alignment_records": int(failed) if failed.is_integer() else failed,
        "alignment_success_rate_completed": (success / records) if records else "",
        "transformed_model_pairs": int(total("transformed_model_pairs")),
        "wall_seconds_sum_completed": sum(float(row["wall_seconds"]) for row in completed if isinstance(row.get("wall_seconds"), (int, float))),
        "wall_seconds_max_completed": max((float(row["wall_seconds"]) for row in completed if isinstance(row.get("wall_seconds"), (int, float))), default=""),
        "template_count": manifest.get("template_count", ""),
        "template_sha256": manifest.get("template_sha256", ""),
        "rank": rows[0].get("rank", "") if rows else "",
        "rank_method": rows[0].get("rank_method", "") if rows else "",
    }


def inventory_rows(run_root: Path, output_root: Path, manifests: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    known = [
        ("submission_manifest", run_root / "full_submission_manifest.json", "KEEP", "submission provenance"),
        ("compact_output", output_root, "KEEP", "generated compact evidence"),
    ]
    for pipeline, manifest in manifests.items():
        manifest_path = Path(manifest["batch_manifest"])
        pipeline_root = run_root / pipeline
        known.extend([
            (f"{pipeline}_batch_manifest", manifest_path, "KEEP", "input manifest and template contract"),
            (f"{pipeline}_case_ledgers", pipeline_root / "cases", "KEEP", "case completion and failure accounting"),
            (f"{pipeline}_attempt_summaries", pipeline_root / "cases/*/attempts/*/run_summary.json", "KEEP", "timing and command provenance"),
            (f"{pipeline}_aggregate_copy", pipeline_root / "aggregate_corrected", "REVIEW", "partial raw-artifact copy from timed-out aggregator"),
            (f"{pipeline}_raw_alignment", pipeline_root / "cases/*/attempts/*/processed/alignment_*", "REVIEW", "needed for future score/mapping analysis"),
            (f"{pipeline}_transformed_models", pipeline_root / "cases/*/attempts/*/processed/transformation", "REVIEW", "needed for DockQ/PyRosetta evaluation"),
        ])
    for role, path, classification, reason in known:
        exact = path.is_file() or path.is_dir()
        size = path.stat().st_size if path.is_file() else ""
        rows.append({
            "role": role,
            "path": str(path),
            "exists": str(exact).lower(),
            "kind": "file" if path.is_file() else "directory_or_pattern",
            "size_bytes": size,
            "size_status": "exact_file" if path.is_file() else "not_measured_recursive",
            "classification": classification,
            "reason": reason,
        })
    return rows


def cleanup_rows(run_root: Path, output_root: Path, manifests: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows = [
        {"path": str(output_root), "classification": "KEEP", "reason": "compact evidence package", "risk": "low", "action": "preserve"},
        {"path": str(run_root / "full_submission_manifest.json"), "classification": "PROTECTED", "reason": "submission provenance and scheduler IDs", "risk": "high", "action": "never delete"},
    ]
    for pipeline, manifest in manifests.items():
        root = run_root / pipeline
        rows.extend([
            {"path": manifest["batch_manifest"], "classification": "PROTECTED", "reason": "input/template contract and case foreign keys", "risk": "high", "action": "never delete"},
            {"path": str(root / "cases/*/task_status.json"), "classification": "KEEP", "reason": "durable case status ledger", "risk": "high", "action": "preserve"},
            {"path": str(root / "cases/*/attempts/*/run_summary.json"), "classification": "KEEP", "reason": "timings, counts, command, environment provenance", "risk": "high", "action": "preserve"},
            {"path": str(root / "cases/*/controller_logs"), "classification": "KEEP", "reason": "failure/retry evidence", "risk": "medium", "action": "preserve until review"},
            {"path": str(root / "aggregate_corrected"), "classification": "REVIEW", "reason": "partial duplicate tree created by timed-out aggregator", "risk": "medium", "action": "do not delete before review"},
            {"path": str(root / "cases/*/attempts/*/processed/transformation"), "classification": "REVIEW", "reason": "required for future DockQ/PyRosetta scoring", "risk": "high", "action": "retain if scoring remains planned"},
            {"path": str(root / "cases/*/attempts/*/processed/alignment_*"), "classification": "REVIEW", "reason": "raw mappings/scores may be required for candidate-level analysis", "risk": "high", "action": "retain until candidate metrics are extracted"},
            {"path": str(root / "cases/*/attempts/*/pipeline.console.log"), "classification": "KEEP", "reason": "diagnostic evidence", "risk": "medium", "action": "preserve until final validation"},
        ])
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submission-manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--aligner", action="append", dest="aligners", choices=DEFAULT_ALIGNERS)
    args = parser.parse_args()
    started = time.perf_counter()
    submission_path = args.submission_manifest.resolve()
    run_root = submission_path.parent.resolve()
    output_root = args.output_root.resolve()
    submission = read_json(submission_path)
    if not submission:
        raise SystemExit(f"invalid submission manifest: {submission_path}")
    aligners = args.aligners or list(DEFAULT_ALIGNERS)
    manifests: dict[str, dict[str, Any]] = {}
    all_case_rows: list[dict[str, Any]] = []
    all_timing_rows: list[dict[str, Any]] = []
    all_failures: list[dict[str, Any]] = []
    errors: list[str] = []
    for pipeline in aligners:
        try:
            path = find_manifest(submission, run_root, pipeline)
            case_rows, timing_rows, manifest, failures = collect_pipeline(pipeline, path)
            manifests[pipeline] = manifest
            all_case_rows.extend(case_rows)
            all_timing_rows.extend(timing_rows)
            all_failures.extend(failures)
        except (FileNotFoundError, OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
            errors.append(f"{pipeline}: {type(exc).__name__}: {exc}")

    output_root.mkdir(parents=True, exist_ok=True)
    write_csv(output_root / "case_results.csv", all_case_rows, compact_fieldnames())
    write_csv(output_root / "stage_timings.csv", all_timing_rows, [
        "pipeline", "case_index", "case_id", "stage", "status", "return_code", "started_at", "completed_at", "wall_seconds", "detail"
    ])
    write_csv(output_root / "failure_catalog.csv", all_failures, [
        "pipeline", "case_index", "case_id", "status", "evidence_status", "validation_errors_json", "conservation", "task_status_path", "summary_path"
    ])
    comparison = [comparison_row(pipeline, [row for row in all_case_rows if row["pipeline"] == pipeline], manifest) for pipeline, manifest in manifests.items()]
    write_csv(output_root / "comparison_summary.csv", comparison, list(comparison[0]) if comparison else ["pipeline"])
    inventory = inventory_rows(run_root, output_root, manifests)
    write_csv(output_root / "artifact_inventory.csv", inventory, ["role", "path", "exists", "kind", "size_bytes", "size_status", "classification", "reason"])
    cleanup = cleanup_rows(run_root, output_root, manifests)
    write_csv(output_root / "cleanup_manifest.csv", cleanup, ["path", "classification", "reason", "risk", "action"])

    validation: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": utc_now(),
        "source_submission_manifest": str(submission_path),
        "source_submission_manifest_sha256": sha256_file(submission_path),
        "errors": errors,
        "checks": {},
    }
    case_keys = [(row["pipeline"], row["case_index"]) for row in all_case_rows]
    validation["checks"]["unique_case_keys"] = len(case_keys) == len(set(case_keys))
    validation["checks"]["all_expected_case_counts"] = all(
        manifest.get("entry_count") == 257 for manifest in manifests.values()
    )
    validation["checks"]["all_case_ledgers_present"] = all(
        row["status"] != "missing" for row in all_case_rows
    )
    validation["checks"]["completed_cases_have_summaries"] = all(
        row["status"] != "completed" or row["evidence_status"] == "complete_summary" for row in all_case_rows
    )
    validation["checks"]["summary_hashes_match_when_claimed"] = all(
        row["summary_sha256_match"] in {"true", ""} for row in all_case_rows
    )
    validation["checks"]["per_case_count_conservation"] = all(
        row["conservation_success_plus_failed_equals_records"] in {"true", ""} for row in all_case_rows
    )
    validation["checks"]["template_contract_consistent"] = len({
        (manifest.get("template_count"), manifest.get("template_sha256")) for manifest in manifests.values()
    }) <= 1
    validation["pipeline_checks"] = {}
    for pipeline, manifest in manifests.items():
        rows = [row for row in all_case_rows if row["pipeline"] == pipeline]
        completed = [row for row in rows if row["status"] == "completed"]
        observed = sum(float(row["alignment_records"]) for row in rows if isinstance(row.get("alignment_records"), (int, float)))
        expected = manifest.get("expected_alignment_records")
        validation["pipeline_checks"][pipeline] = {
            "case_count": len(rows),
            "completed_count": len(completed),
            "observed_alignment_records": int(observed) if observed.is_integer() else observed,
            "expected_alignment_records": expected,
            "completed_total_matches_expected": (observed == expected) if len(completed) == len(rows) else "not_applicable_incomplete",
            "status_counts": dict(Counter(row["status"] for row in rows)),
        }
    validation["checks"]["no_raw_artifacts_copied"] = not any(
        path.name.endswith((".pdb", ".jsonl", ".parquet")) for path in output_root.rglob("*") if path.is_file()
    )
    validation["overall"] = bool(not errors and all(validation["checks"].values()))
    write_json(output_root / "validation.json", validation)
    write_json(output_root / "schema.json", {
        "schema_version": SCHEMA_VERSION,
        "grain": "one row per pipeline and expected case index",
        "case_results_fields": compact_fieldnames(),
        "stage_timings_grain": "one row per observed stage event",
        "known_limits": [
            "Candidate-level alignment scores and mapping rows are not extracted.",
            "DockQ and PyRosetta are not present in this source run because refinement was disabled.",
            "Recursive raw-artifact sizes are intentionally not measured by this collector.",
        ],
    })
    write_json(output_root / "collection_manifest.json", {
        "schema_version": SCHEMA_VERSION,
        "status": "completed" if not errors else "completed_with_errors",
        "generated_at": utc_now(),
        "collector": str(Path(__file__).resolve()),
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "source_submission_manifest": str(submission_path),
        "source_submission_manifest_sha256": sha256_file(submission_path),
        "source_run_root": str(run_root),
        "aligners": aligners,
        "output_root": str(output_root),
        "case_rows": len(all_case_rows),
        "stage_rows": len(all_timing_rows),
        "failure_rows": len(all_failures),
        "elapsed_seconds": time.perf_counter() - started,
        "copy_policy": "summary-only; no raw alignment JSON or transformed PDB copied",
    })
    print(json.dumps({
        "status": "completed" if not errors else "completed_with_errors",
        "output_root": str(output_root),
        "case_rows": len(all_case_rows),
        "stage_rows": len(all_timing_rows),
        "failure_rows": len(all_failures),
        "errors": errors,
        "validation_overall": validation["overall"],
    }, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
