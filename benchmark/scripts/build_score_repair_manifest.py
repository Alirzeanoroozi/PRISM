#!/usr/bin/env python3
"""Extract staged rows matching explicit failed canonical score identities."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_table(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix == ".tsv" else ","
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def identity(row: dict[str, str]) -> tuple[str, str]:
    return row.get("dataset_row_id", ""), row.get("source_model_sha256", "")


def build(stage_path: Path, scores_path: Path, output_path: Path) -> int:
    stages = read_table(stage_path)
    failed = [row for row in read_table(scores_path) if row.get("score_status") == "score_failed"]
    failed_keys = [identity(row) for row in failed]
    if any(not all(key) for key in failed_keys):
        raise ValueError("failed score row lacks dataset_row_id or source_model_sha256")
    if len(failed_keys) != len(set(failed_keys)):
        raise ValueError("duplicate failed score identity")

    stage_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for row in stages:
        key = identity(row)
        if not all(key):
            continue
        if key in stage_by_key:
            raise ValueError(f"duplicate stage identity: {key}")
        stage_by_key[key] = row
    missing = sorted(set(failed_keys) - set(stage_by_key))
    if missing:
        raise ValueError(f"failed identities absent from stage manifest: {missing[:5]}")

    selected = [stage_by_key[key] for key in failed_keys]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(stages[0]) if stages else []
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(selected)
    return len(selected)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-manifest", type=Path, required=True)
    parser.add_argument("--failed-scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rows", type=int)
    args = parser.parse_args()
    count = build(args.stage_manifest, args.failed_scores, args.output)
    if args.expected_rows is not None and count != args.expected_rows:
        raise SystemExit(f"expected {args.expected_rows} repair rows; found {count}")
    print(f"repair_rows={count} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
