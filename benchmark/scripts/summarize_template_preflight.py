#!/usr/bin/env python3
"""Summarize a template asset TSV without changing the source asset tree."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def summarize(path: Path) -> dict[str, object]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped[row["template_id"]].append(row)
    fully = sum(all(row["fully_resolvable"] == "1" for row in rows) for rows in grouped.values())
    valid = sum(all(row["valid"] == "1" for row in rows) for rows in grouped.values())
    missing = sorted(template_id for template_id, rows in grouped.items() if not all(row["fully_resolvable"] == "1" for row in rows))
    return {
        "asset_tsv": str(path),
        "asset_row_count": sum(len(rows) for rows in grouped.values()),
        "unique_template_count": len(grouped),
        "valid_template_count": valid,
        "fully_resolvable_template_count": fully,
        "missing_template_count": len(missing),
        "missing_template_ids": missing,
        "confirmatory_asset_coverage_status": "incomplete" if missing else "complete",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summarize(args.input), indent=2, sort_keys=True) + "\n", encoding="utf-8")
