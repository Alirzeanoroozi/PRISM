#!/usr/bin/env python3
"""Rank a PRISM candidate CSV with the deterministic biological baseline."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.candidate_ranker import BASELINE_SCORE_VERSION, rank_candidates


def _ranking_group(row: dict[str, str]) -> tuple[str, str]:
    """Return the most durable available benchmark-row grouping key."""
    for field in ("dataset_row_id", "native_complex_id"):
        value = row.get(field, "").strip()
        if value:
            return field, value
    return "query_pair", f"{row.get('query_left', '')}|{row.get('query_right', '')}"


def rank_csv(input_path: Path, output_path: Path) -> int:
    with input_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    grouped: dict[tuple[str, str], list[dict[str, str | int]]] = {}
    for source_index, row in enumerate(rows):
        grouped.setdefault(_ranking_group(row), []).append({**row, "_ranking_source_index": source_index})

    rank_by_source_index: dict[int, tuple[float, int]] = {}
    ranked_count = 0
    for group_rows in grouped.values():
        ranked_group = rank_candidates(group_rows)
        ranked_count += len(ranked_group)
        for rank, ranked_row in enumerate(ranked_group, start=1):
            source_index = int(ranked_row["_ranking_source_index"])
            rank_by_source_index[source_index] = (float(ranked_row["baseline_score"]), rank)

    output_rows = []
    for source_index, row in enumerate(rows):
        group_field, group_value = _ranking_group(row)
        score_and_rank = rank_by_source_index.get(source_index)
        output = dict(row)
        output["ranking_group"] = f"{group_field}={group_value}"
        output["baseline_score"] = "" if score_and_rank is None else f"{score_and_rank[0]:.8f}"
        output["baseline_rank"] = "" if score_and_rank is None else str(score_and_rank[1])
        output["baseline_score_version"] = (
            "" if score_and_rank is None else BASELINE_SCORE_VERSION
        )
        output_rows.append(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        fields = list(rows[0]) if rows else ["status"]
        for field in (
            "ranking_group", "baseline_score", "baseline_rank", "baseline_score_version"
        ):
            if field not in fields:
                fields.append(field)
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    return ranked_count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(f"ranked {rank_csv(args.input, args.output)} eligible candidates")


if __name__ == "__main__":
    main()
