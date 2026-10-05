#!/usr/bin/env python3
"""Evaluate ranked candidate tables per independent benchmark row."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median


def _group_id(row: dict[str, str]) -> str:
    value = row.get("dataset_row_id") or row.get("native_complex_id") or ""
    if not value:
        raise ValueError("ranked rows require dataset_row_id or native_complex_id")
    return value


def summarize(rows: list[dict[str, str]]) -> dict[str, object]:
    """Return per-row top-1 and oracle summaries for labeled candidates."""
    eligible = [
        row for row in rows
        if row.get("baseline_rank") not in (None, "")
        and row.get("label_status", "labeled") == "labeled"
        and row.get("dockq") not in (None, "")
        and row.get("native_like") in {"0", "1"}
    ]
    groups: dict[str, list[dict[str, str]]] = {}
    for row in eligible:
        groups.setdefault(_group_id(row), []).append(row)

    score_versions = sorted({
        row.get("baseline_score_version", "").strip()
        for row in eligible
        if row.get("baseline_score_version", "").strip()
    })
    if len(score_versions) > 1:
        raise ValueError(
            "ranked table mixes baseline score versions: " + ", ".join(score_versions)
        )

    per_group = []
    for group_id, group_rows in sorted(groups.items()):
        top1 = [row for row in group_rows if row["baseline_rank"] == "1"]
        if len(top1) != 1:
            raise ValueError(f"{group_id}: expected one baseline_rank=1 row, found {len(top1)}")
        top1_dockq = float(top1[0]["dockq"])
        best_dockq = max(float(row["dockq"]) for row in group_rows)
        native_like_count = sum(row["native_like"] == "1" for row in group_rows)
        per_group.append({
            "group_id": group_id,
            "candidate_count": len(group_rows),
            "native_like_count": native_like_count,
            "top1_native_like": int(top1[0]["native_like"]),
            "top1_dockq": top1_dockq,
            "best_dockq": best_dockq,
            "dockq_regret": best_dockq - top1_dockq,
        })

    count = len(per_group)
    return {
        "input_row_count": len(rows),
        "eligible_labeled_row_count": len(eligible),
        "ranking_group_count": count,
        "baseline_score_version": score_versions[0] if score_versions else "unversioned",
        "oracle_native_like_group_count": sum(row["native_like_count"] > 0 for row in per_group),
        "top1_native_like_group_count": sum(row["top1_native_like"] for row in per_group),
        "oracle_native_like_rate": (
            sum(row["native_like_count"] > 0 for row in per_group) / count if count else 0.0
        ),
        "top1_native_like_rate": (
            sum(row["top1_native_like"] for row in per_group) / count if count else 0.0
        ),
        "median_top1_dockq": median(row["top1_dockq"] for row in per_group) if count else 0.0,
        "median_best_dockq": median(row["best_dockq"] for row in per_group) if count else 0.0,
        "median_dockq_regret": median(row["dockq_regret"] for row in per_group) if count else 0.0,
        "groups": per_group,
    }


def evaluate(input_path: Path, output_path: Path) -> dict[str, object]:
    with input_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    summary = summarize(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(evaluate(args.input, args.output), indent=2))


if __name__ == "__main__":
    main()
