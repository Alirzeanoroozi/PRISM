#!/usr/bin/env python3
"""Classify score failures without retrying or replacing scientific outcomes."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def failure_class(error: str) -> str:
    text = (error or "").lower()
    if "buffer has wrong number of dimensions" in text:
        return "dockq_runtime_buffer_dimensions"
    if "missing native" in text:
        return "missing_native_asset"
    if "no requested cross interfaces" in text:
        return "missing_cross_interface"
    if "chain contract" in text or "transformation half" in text:
        return "model_integrity_rejected"
    return "unclassified_score_failure"


def adjudicate(input_paths: list[Path], output: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted(input_paths):
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("score_status") != "score_failed":
                    continue
                error = row.get("score_error", "")
                rows.append(
                    {
                        "source_scores_tsv": str(path),
                        "pair_id": row.get("pair_id", ""),
                        "benchmark_set": row.get("benchmark_set", ""),
                        "source_model_path": row.get("source_model_path", ""),
                        "source_model_sha256": row.get("source_model_sha256", ""),
                        "native_pdb_sha256": row.get("native_pdb_sha256", ""),
                        "failure_class": failure_class(error),
                        "failure_detail_sha256": hashlib.sha256(error.encode("utf-8")).hexdigest(),
                        "terminal_status": "score_failed",
                        "retry_authorized": "false",
                    }
                )
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["source_scores_tsv", "pair_id", "benchmark_set", "source_model_path", "source_model_sha256", "native_pdb_sha256", "failure_class", "failure_detail_sha256", "terminal_status", "retry_authorized"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(f"adjudicated={len(adjudicate(args.input, args.output))}")
