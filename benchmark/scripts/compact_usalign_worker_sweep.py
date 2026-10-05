#!/usr/bin/env python3
"""Compact a completed USalign worker sweep before scratch cleanup."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


SUMMARY_FIELDS = [
    "configuration", "worker_count", "expected_records", "success_records",
    "failure_records", "records_per_second", "wall_seconds",
    "child_cpu_seconds", "child_system_seconds", "child_max_rss_kb", "status",
]
RECORD_FIELDS = [
    "configuration", "worker_count", "template_id", "chain",
    "execution_status", "status", "return_code", "match_count",
    "tm_score_query", "tm_score_ref", "tm_score", "tm_score_contract",
    "elapsed_seconds", "raw_output_sha256",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compact(input_root: Path, output_root: Path) -> dict[str, object]:
    output_root.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []
    record_rows: list[dict[str, object]] = []
    for summary_path in sorted(input_root.glob("*_w*/summary.json")):
        summary = json.loads(summary_path.read_text())
        summary_rows.append({field: summary.get(field) for field in SUMMARY_FIELDS})
        configuration = summary["configuration"]
        workers = summary["worker_count"]
        records_dir = summary_path.parent / "records"
        for record_path in sorted(records_dir.glob("*.json")):
            record = json.loads(record_path.read_text())
            record_rows.append({
                "configuration": configuration,
                "worker_count": workers,
                **{field: record.get(field) for field in RECORD_FIELDS[2:]},
            })

    with (output_root / "worker_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(summary_rows)
    with (output_root / "alignment_record_ledger.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECORD_FIELDS)
        writer.writeheader()
        writer.writerows(record_rows)

    expected = [int(row["expected_records"]) for row in summary_rows]
    observed = [
        sum(row["configuration"] == summary["configuration"] and row["worker_count"] == summary["worker_count"] for row in record_rows)
        for summary in summary_rows
    ]
    validation = {
        "schema_version": "prism-usalign-worker-sweep-compact/v1",
        "configuration_count": len(summary_rows),
        "record_count": len(record_rows),
        "all_execution_failures_zero": all(int(row["failure_records"]) == 0 for row in summary_rows),
        "summary_expected_records": expected,
        "summary_observed_record_counts": observed,
        "counts_match": expected == observed,
    }
    left_root = input_root / "default_w16" / "records"
    right_root = input_root / "fast_w16" / "records"
    if left_root.is_dir() and right_root.is_dir():
        left = {
            path.name: json.loads(path.read_text())
            for path in left_root.glob("*.json")
        }
        right = {
            path.name: json.loads(path.read_text())
            for path in right_root.glob("*.json")
        }
        shared = sorted(set(left) & set(right))
        score_fields = ("tm_score_query", "tm_score_ref", "tm_score")
        comparison = {
            "left": "default_w16",
            "right": "fast_w16",
            "left_records": len(left),
            "right_records": len(right),
            "shared_records": len(shared),
            "missing_from_left": sorted(set(right) - set(left)),
            "missing_from_right": sorted(set(left) - set(right)),
            "status_differences": 0,
            "mapping_differences": 0,
            "transform_differences": 0,
            "score_differences": 0,
            "max_absolute_score_delta": 0.0,
        }
        for name in shared:
            a, b = left[name], right[name]
            comparison["status_differences"] += int(a.get("status") != b.get("status"))
            comparison["mapping_differences"] += int(a.get("match_dict") != b.get("match_dict"))
            comparison["transform_differences"] += int(
                a.get("translation") != b.get("translation")
                or a.get("rotation_mat") != b.get("rotation_mat")
            )
            deltas = [abs(float(a.get(field, 0.0)) - float(b.get(field, 0.0))) for field in score_fields]
            comparison["score_differences"] += int(any(delta != 0.0 for delta in deltas))
            comparison["max_absolute_score_delta"] = max(
                comparison["max_absolute_score_delta"], *deltas
            )
        validation["default_vs_fast_w16"] = comparison
        (output_root / "default_vs_fast_w16.json").write_text(
            json.dumps(comparison, indent=2, sort_keys=True) + "\n"
        )
    (output_root / "validation.json").write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    return validation


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    validation = compact(args.input_root, args.output_root)
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0 if validation["all_execution_failures_zero"] and validation["counts_match"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
