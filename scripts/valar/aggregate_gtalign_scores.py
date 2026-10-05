#!/usr/bin/env python3
"""Aggregate completed corrected GTalign score chunks with no-drop checks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()
    manifest_path = Path(args.input_manifest).resolve()
    output_root = Path(args.output_root).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    task_dir = output_root / "tasks"
    tasks = []
    for path in sorted(task_dir.glob("task_*.json")):
        tasks.append(json.loads(path.read_text(encoding="utf-8")))
    chunk_size = int(manifest.get("chunk_size", 400))
    expected = math.ceil(len(manifest["entries"]) / chunk_size) if manifest["entries"] else 0
    starts = {int(task["start"]) for task in tasks}
    expected_starts = {index * chunk_size for index in range(expected)}
    records = [record for task in tasks for record in task.get("records", [])]
    indices = [int(record["index"]) for record in records]
    duplicate_indices = sorted(index for index, count in Counter(indices).items() if count > 1)
    missing_indices = sorted(set(range(len(manifest["entries"]))) - set(indices))
    records.sort(key=lambda record: int(record["index"]))
    columns = sorted({key for record in records for key in record})
    csv_path = output_root / "dockq_irmsd.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    statuses = Counter(str(record.get("status", "UNKNOWN")) for record in records)
    complete = not duplicate_indices and not missing_indices and expected_starts == starts
    summary = {
        "schema_version": "prism-gtalign-score-summary-20260920",
        "status": "completed" if complete else "incomplete",
        "stage": manifest["stage"],
        "candidate_count": len(manifest["entries"]),
        "missing_refinement_count": len(manifest.get("missing", [])),
        "task_count": len(tasks),
        "expected_task_count": expected,
        "chunk_size": chunk_size,
        "rows": len(records),
        "status_counts": dict(sorted(statuses.items())),
        "duplicate_indices": duplicate_indices,
        "missing_indices": missing_indices,
        "manifest_sha256": sha256(manifest_path),
        "csv": str(csv_path),
        "ranking": False,
        "prodigy": False,
    }
    (output_root / "score_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
