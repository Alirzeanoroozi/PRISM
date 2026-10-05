#!/usr/bin/env python3
"""Build an exact CSV retry manifest for failed bijective score rows."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def key(row: dict[str, str]) -> tuple[str, str]:
    return row.get("dataset_row_id", ""), row.get("source_model_sha256", "")


def read_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def prepare(
    base_scores: Path,
    stage_manifest: Path,
    output: Path,
    *,
    score_status: str = "score_failed",
    required_field: str = "",
    required_value: str = "",
) -> list[dict[str, str]]:
    failed = [
        row for row in read_rows(base_scores, "\t")
        if row.get("score_status") == score_status
        and (not required_field or row.get(required_field) == required_value)
    ]
    failed_keys = {key(row) for row in failed}
    if len(failed_keys) != len(failed) or any(not all(item) for item in failed_keys):
        raise ValueError("failed score rows lack unique dataset_row_id + source_model_sha256 keys")

    stages = read_rows(stage_manifest, ",")
    by_key: dict[tuple[str, str], dict[str, str]] = {}
    for row in stages:
        row_key = key(row)
        if row_key in by_key:
            raise ValueError(f"duplicate stage key: {row_key}")
        by_key[row_key] = row
    missing = sorted(failed_keys - set(by_key))
    if missing:
        raise ValueError(f"failed rows absent from stage manifest: {missing[:3]}")

    retry_rows = [by_key[key(row)] for row in failed]
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(stages[0]) if stages else []
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(retry_rows)
    return retry_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-scores", type=Path, required=True)
    parser.add_argument("--stage-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--score-status", default="score_failed")
    parser.add_argument("--required-field", default="")
    parser.add_argument("--required-value", default="")
    args = parser.parse_args()
    if bool(args.required_field) != bool(args.required_value):
        parser.error("--required-field and --required-value must be supplied together")
    rows = prepare(
        args.base_scores,
        args.stage_manifest,
        args.output,
        score_status=args.score_status,
        required_field=args.required_field,
        required_value=args.required_value,
    )
    print(f"retry_rows={len(rows)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
