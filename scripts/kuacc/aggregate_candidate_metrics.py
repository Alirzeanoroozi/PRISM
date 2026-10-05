#!/usr/bin/env python3
"""Aggregate PRISM candidate-audit JSONL into compact analysis tables.

All alignment-audit records are streamed once. Rejected records are reduced to
case/orientation/status/reason summary statistics; only non-threshold-rejected
records are retained individually for future transformed-model analysis.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from aggregate_compact_results import DEFAULT_ALIGNERS, find_manifest, read_json, sha256_file


SCHEMA_VERSION = "prism-candidate-metrics-v1"
NUMERIC_FIELDS = (
    "match_count_left", "match_count_right", "match_coverage_left", "match_coverage_right",
    "tm_score_left", "tm_score_right", "clash_count", "contact_count", "rosetta_interaction_score",
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]] | None, fieldnames: list[str], producer=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        if rows is not None:
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fieldnames})
        elif producer is not None:
            producer(writer)
    temporary.replace(path)


def numeric(value: Any) -> float | None:
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def candidate_fields() -> list[str]:
    return [
        "pipeline", "dataset", "split", "case_index", "case_id", "template", "orientation",
        "chain_left", "chain_right", "query_left", "query_right", "status", "error_reason",
        "match_count_left", "match_count_right", "match_coverage_left", "match_coverage_right",
        "tm_score_left", "tm_score_right", "clash_count", "contact_count",
        "rosetta_interaction_score", "source_pipeline",
    ]


def summary_fields() -> list[str]:
    fields = ["pipeline", "dataset", "split", "case_index", "case_id", "orientation", "status", "error_reason", "record_count"]
    for name in NUMERIC_FIELDS:
        fields.extend([f"{name}_min", f"{name}_mean", f"{name}_max"])
    return fields


class Aggregate:
    def __init__(self, pipeline: str, entry: dict[str, Any], item: dict[str, Any]) -> None:
        self.pipeline = pipeline
        self.entry = entry
        self.orientation = item.get("orientation", "")
        self.status = item.get("status", "")
        self.error_reason = item.get("error_reason") or ""
        self.count = 0
        self.stats: dict[str, list[float]] = defaultdict(list)

    def add(self, item: dict[str, Any]) -> None:
        self.count += 1
        for name in NUMERIC_FIELDS:
            value = numeric(item.get(name))
            if value is not None:
                self.stats[name].append(float(value))

    def row(self) -> dict[str, Any]:
        row: dict[str, Any] = {
            "pipeline": self.pipeline,
            "dataset": self.entry.get("dataset", "bm55_full"),
            "split": self.entry.get("split", ""),
            "case_index": self.entry.get("index", ""),
            "case_id": self.entry.get("case_id", ""),
            "orientation": self.orientation,
            "status": self.status,
            "error_reason": self.error_reason,
            "record_count": self.count,
        }
        for name in NUMERIC_FIELDS:
            values = self.stats[name]
            row[f"{name}_min"] = min(values) if values else ""
            row[f"{name}_mean"] = sum(values) / len(values) if values else ""
            row[f"{name}_max"] = max(values) if values else ""
        return row


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
    output_root.mkdir(parents=True, exist_ok=True)
    atomic_json(output_root / "candidate_aggregation_status.json", {"status": "running", "started_at": now(), "aligners": aligners})

    generated_rows: list[dict[str, Any]] = []
    aggregates: dict[tuple[Any, ...], Aggregate] = {}
    status_counts: Counter[str] = Counter()
    pipeline_status_counts: dict[str, Counter[str]] = defaultdict(Counter)
    threshold_payloads: dict[str, dict[str, Any]] = {}
    audit_file_count: Counter[str] = Counter()
    audit_record_count: Counter[str] = Counter()
    errors: list[str] = []
    case_count: Counter[str] = Counter()

    for pipeline in aligners:
        try:
            manifest_path = find_manifest(submission, run_root, pipeline)
            manifest = read_json(manifest_path)
            if not manifest:
                raise RuntimeError(f"invalid batch manifest: {manifest_path}")
            entries = sorted(manifest.get("entries", []), key=lambda item: int(item["index"]))
            for entry in entries:
                case_count[pipeline] += 1
                task_path = Path(entry["case_root"]) / "task_status.json"
                task = read_json(task_path) or {}
                attempt_root = Path(str(task.get("attempt_root", ""))) if task.get("attempt_root") else None
                if not attempt_root:
                    continue
                audit_files = sorted((attempt_root / "processed" / "candidate_audit").glob("*.jsonl"))
                if not audit_files:
                    continue
                audit_file_count[pipeline] += len(audit_files)
                for audit_path in audit_files:
                    with audit_path.open(encoding="utf-8") as handle:
                        for line_number, line in enumerate(handle, 1):
                            try:
                                item = json.loads(line)
                            except json.JSONDecodeError as exc:
                                errors.append(f"{audit_path}:{line_number}: {exc}")
                                continue
                            if not isinstance(item, dict):
                                errors.append(f"{audit_path}:{line_number}: non-object record")
                                continue
                            status = str(item.get("status", ""))
                            reason = item.get("error_reason") or ""
                            status_counts[f"{pipeline}:{status}"] += 1
                            pipeline_status_counts[pipeline][status] += 1
                            audit_record_count[pipeline] += 1
                            thresholds = (item.get("metadata") or {}).get("transformation_thresholds")
                            if isinstance(thresholds, dict):
                                threshold_payloads.setdefault(pipeline, thresholds)
                            key = (pipeline, entry["index"], item.get("orientation", ""), status, reason)
                            bucket = aggregates.get(key)
                            if bucket is None:
                                bucket = Aggregate(pipeline, entry, item)
                                aggregates[key] = bucket
                            bucket.add(item)
                            if status != "alignment_threshold_rejected":
                                generated_rows.append({
                                    "pipeline": pipeline,
                                    "dataset": entry.get("dataset", manifest.get("dataset", "bm55_full")),
                                    "split": entry.get("split", ""),
                                    "case_index": entry.get("index", ""),
                                    "case_id": entry.get("case_id", ""),
                                    "template": item.get("template", ""),
                                    "orientation": item.get("orientation", ""),
                                    "chain_left": item.get("chain_left", ""),
                                    "chain_right": item.get("chain_right", ""),
                                    "query_left": item.get("query_left", ""),
                                    "query_right": item.get("query_right", ""),
                                    "status": status,
                                    "error_reason": reason,
                                    **{name: item.get(name, "") for name in NUMERIC_FIELDS},
                                    "source_pipeline": item.get("source_pipeline", ""),
                                })
        except (FileNotFoundError, OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
            errors.append(f"{pipeline}: {type(exc).__name__}: {exc}")

    summary_rows = [aggregates[key].row() for key in sorted(aggregates, key=lambda value: tuple(str(part) for part in value))]
    write_csv(output_root / "candidate_generated.csv", generated_rows, candidate_fields())
    write_csv(output_root / "candidate_rejection_summary.csv", summary_rows, summary_fields())
    atomic_json(output_root / "candidate_thresholds.json", threshold_payloads)
    validation = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": now(),
        "source_submission_manifest": str(submission_path),
        "source_submission_manifest_sha256": sha256_file(submission_path),
        "errors": errors,
        "pipeline_checks": {
            pipeline: {
                "case_count": case_count[pipeline],
                "audit_file_count": audit_file_count[pipeline],
                "audit_record_count": audit_record_count[pipeline],
                "status_counts": dict(pipeline_status_counts[pipeline]),
                "candidate_rows_retained": sum(1 for row in generated_rows if row["pipeline"] == pipeline),
            }
            for pipeline in aligners
        },
        "retention_policy": {
            "retained_individually": "records whose status is not alignment_threshold_rejected",
            "aggregated_only": "alignment_threshold_rejected records, grouped by case/orientation/status/reason",
            "raw_candidate_audit_files_copied": False,
        },
    }
    validation["overall"] = not errors and all(
        validation["pipeline_checks"][pipeline]["audit_file_count"] > 0
        for pipeline in aligners if validation["pipeline_checks"][pipeline]["case_count"] > 0
    )
    atomic_json(output_root / "candidate_validation.json", validation)
    atomic_json(output_root / "candidate_collection_manifest.json", {
        "schema_version": SCHEMA_VERSION,
        "status": "completed" if not errors else "completed_with_errors",
        "generated_at": now(),
        "collector": str(Path(__file__).resolve()),
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "source_run_root": str(run_root),
        "aligners": aligners,
        "output_root": str(output_root),
        "candidate_rows_retained": len(generated_rows),
        "summary_rows": len(summary_rows),
        "audit_records": dict(audit_record_count),
        "elapsed_seconds": time.perf_counter() - started,
    })
    atomic_json(output_root / "candidate_aggregation_status.json", {
        "status": "completed" if not errors else "completed_with_errors",
        "updated_at": now(),
        "candidate_rows_retained": len(generated_rows),
        "summary_rows": len(summary_rows),
        "audit_records": dict(audit_record_count),
        "errors": errors,
    })
    print(json.dumps({
        "status": "completed" if not errors else "completed_with_errors",
        "candidate_rows_retained": len(generated_rows),
        "summary_rows": len(summary_rows),
        "audit_records": dict(audit_record_count),
        "errors": errors,
        "validation_overall": validation["overall"],
    }, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
