#!/usr/bin/env python3
"""Join current and legacy batch outcomes by the shared pair manifest."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def pipeline_rows(status_csv: Path, pipeline: str) -> dict[str, dict[str, str]]:
    return {row["pair_id"]: row for row in read_csv(status_csv) if row.get("pipeline") == pipeline}


def classify(current: dict[str, str], legacy: dict[str, str]) -> tuple[str, str]:
    cs, ls = current.get("status", "missing"), legacy.get("status", "missing")
    if cs == "input_missing" or ls == "input_missing":
        return "input_missing", "benchmark structure unavailable during login-side staging"
    if cs == "completed" and ls == "completed":
        return "both_completed", "compare candidate counts and shared DockQ/iRMSD metrics"
    if cs == "completed" and ls != "completed":
        if "execution" in ls:
            return "current_only_execution_failure", "legacy execution or dependency failure"
        return "current_only", "legacy batch incomplete or produced no completed marker"
    if ls == "completed" and cs != "completed":
        if "execution" in cs:
            return "legacy_only_execution_failure", "current execution, alignment, or refinement failure"
        return "legacy_only", "current batch incomplete or produced no completed marker"
    if "execution" in cs or "execution" in ls:
        return "both_or_execution_failure", "execution failure before comparable final outputs"
    return "pending_or_unverified", "one or both batches have not completed"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--status", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = read_csv(args.manifest)
    current = pipeline_rows(args.status, "current")
    legacy = pipeline_rows(args.status, "legacy")
    rows = []
    for pair in manifest:
        c, l = current.get(pair["pair_id"], {}), legacy.get(pair["pair_id"], {})
        comparison, explanation = classify(c, l)
        rows.append({
            **pair,
            "current_status": c.get("status", "missing"),
            "current_reason": c.get("reason", ""),
            "legacy_status": l.get("status", "missing"),
            "legacy_reason": l.get("reason", ""),
            "comparison_class": comparison,
            "likely_explanation": explanation,
            "current_model_count": c.get("count_refinement", ""),
            "legacy_model_count": l.get("count_refinement", ""),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["pair_id", "comparison_class"]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    print(f"wrote {len(rows)} pairwise rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
