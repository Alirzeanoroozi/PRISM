#!/usr/bin/env python3
"""Reduce template-asset preflight rows to one auditable row per template."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDS = (
    "template_id",
    "asset_count",
    "asset_types",
    "existing_asset_count",
    "format_valid_asset_count",
    "valid",
    "fully_resolvable",
    "missing",
    "asset_sha256s",
)


def build_coverage(input_tsv: Path, output_tsv: Path) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with input_tsv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            template_id = (row.get("template_id") or "").strip()
            if template_id:
                grouped[template_id].append(row)

    rows: list[dict[str, str]] = []
    for template_id in sorted(grouped):
        assets = sorted(grouped[template_id], key=lambda row: (row.get("asset_type", ""), row.get("asset_path", "")))
        missing = sorted({item for row in assets for item in (row.get("missing", "") or "").split(";") if item})
        hashes = sorted({row.get("sha256", "") for row in assets if row.get("sha256")})
        rows.append(
            {
                "template_id": template_id,
                "asset_count": str(len(assets)),
                "asset_types": ";".join(row.get("asset_type", "") for row in assets),
                "existing_asset_count": str(sum(row.get("exists") == "True" for row in assets)),
                "format_valid_asset_count": str(sum(row.get("format_valid") == "True" for row in assets)),
                "valid": str(int(all(row.get("valid") == "1" for row in assets))),
                "fully_resolvable": str(int(all(row.get("fully_resolvable") == "1" for row in assets))),
                "missing": ";".join(missing),
                "asset_sha256s": ";".join(hashes),
            }
        )

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(f"wrote {len(build_coverage(args.input, args.output))} template coverage rows")
