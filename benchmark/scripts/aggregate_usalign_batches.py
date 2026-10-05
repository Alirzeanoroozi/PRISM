#!/usr/bin/env python3
"""Validate and merge resumable USalign batch summaries without duplication."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_union(path: Path, rows: list[dict[str, Any]], delimiter: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row}) or ["status"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return sha256_file(path)


def aggregate(run_root: Path, output_dir: Path, expected_batches: int) -> dict[str, Any]:
    candidate_rows: list[dict[str, str]] = []
    score_rows: list[dict[str, str]] = []
    interface_rows: list[dict[str, str]] = []
    batch_rows: list[dict[str, Any]] = []
    invalid: list[dict[str, Any]] = []
    for number in range(1, expected_batches + 1):
        batch_root = run_root / "current" / f"batch_{number:04d}"
        status_path = batch_root / "status" / "transformed_dockq_status.json"
        candidate_path = batch_root / "status" / "candidate_generated.csv"
        score_path = batch_root / "status" / "transformed_dockq.tsv"
        interface_path = batch_root / "status" / "transformed_dockq_interfaces.tsv"
        record: dict[str, Any] = {"batch": number, "batch_root": str(batch_root)}
        try:
            summary = json.loads(status_path.read_text(encoding="utf-8"))
            if summary.get("status") != "validated_compacted":
                raise ValueError(f"unexpected status {summary.get('status')}")
            candidates = read_table(candidate_path)
            scores = read_table(score_path)
            interfaces = read_table(interface_path)
            if int(summary.get("candidate_rows", -1)) != len(scores):
                raise ValueError(f"score row mismatch: summary={summary.get('candidate_rows')} actual={len(scores)}")
            candidate_rows.extend(candidates)
            score_rows.extend(scores)
            interface_rows.extend(interfaces)
            record.update(status="validated", candidate_rows=len(candidates), score_rows=len(scores), interface_rows=len(interfaces))
        except (OSError, json.JSONDecodeError, csv.Error, ValueError) as exc:
            record.update(status="invalid", error=f"{type(exc).__name__}: {exc}")
            invalid.append(record)
        batch_rows.append(record)

    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_hash = write_union(output_dir / "candidate_generated.tsv", candidate_rows, "\t")
    score_hash = write_union(output_dir / "transformed_dockq.tsv", score_rows, "\t")
    interface_hash = write_union(output_dir / "transformed_dockq_interfaces.tsv", interface_rows, "\t")
    batch_hash = write_union(output_dir / "batch_validation.tsv", batch_rows, "\t")
    complete = len(batch_rows) == expected_batches and not invalid and all(row.get("status") == "validated" for row in batch_rows)
    result = {
        "schema_version": "prism-usalign-batch-aggregation/v1",
        "status": "validated_compacted" if complete else "incomplete",
        "expected_batches": expected_batches,
        "validated_batches": sum(row.get("status") == "validated" for row in batch_rows),
        "invalid_batches": invalid,
        "candidate_rows": len(candidate_rows),
        "score_rows": len(score_rows),
        "interface_rows": len(interface_rows),
        "candidate_generated_path": str((output_dir / "candidate_generated.tsv").resolve()),
        "transformed_dockq_path": str((output_dir / "transformed_dockq.tsv").resolve()),
        "transformed_interfaces_path": str((output_dir / "transformed_dockq_interfaces.tsv").resolve()),
        "candidate_generated_sha256": candidate_hash,
        "transformed_dockq_sha256": score_hash,
        "transformed_interfaces_sha256": interface_hash,
        "batch_validation_sha256": batch_hash,
        "cleanup_eligible_for_raw_batch_outputs": complete,
        "cleanup_note": "Raw per-batch alignment/transformation outputs require downstream refinement validation before deletion.",
    }
    (output_dir / "aggregation_status.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-batches", type=int, required=True)
    args = parser.parse_args()
    result = aggregate(args.run_root.resolve(), args.output_dir.resolve(), args.expected_batches)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "validated_compacted" else 2


if __name__ == "__main__":
    raise SystemExit(main())
