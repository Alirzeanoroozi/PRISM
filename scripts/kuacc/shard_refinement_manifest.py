#!/usr/bin/env python3
"""Shard a validated refinement manifest into Slurm-array-sized CSVs.

The worker reads one shard per array window instead of reparsing the complete
88K-row manifest for every task.  ``manifest_index`` is stable across shards
and is retained in each checkpoint for no-drop reconciliation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--shard-size", type=int, default=1000)
    args = parser.parse_args()
    if args.shard_size < 1 or args.shard_size > 1000:
        raise SystemExit("--shard-size must be between 1 and 1000 (KUACC MaxArraySize=1001)")

    with args.input_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit("input manifest has no rows")
    required = {"pipeline", "case_id", "left", "right", "native_pdb"}
    missing = sorted(required - set(rows[0]))
    if missing:
        raise SystemExit(f"input manifest missing columns: {missing}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = sorted(set().union(*(set(row) for row in rows)) | {"manifest_index"})
    shards: list[dict[str, object]] = []
    for shard_number, start in enumerate(range(0, len(rows), args.shard_size)):
        shard_rows = []
        for manifest_index, row in enumerate(rows[start : start + args.shard_size], start=start):
            shard_rows.append({**row, "manifest_index": str(manifest_index)})
        shard_path = args.output_dir / f"shard_{shard_number:04d}.csv"
        with shard_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(shard_rows)
        shards.append(
            {
                "shard_number": shard_number,
                "path": str(shard_path.resolve()),
                "first_manifest_index": start,
                "last_manifest_index": start + len(shard_rows) - 1,
                "rows": len(shard_rows),
                "sha256": sha256(shard_path),
                "array": f"0-{len(shard_rows) - 1}",
            }
        )

    manifest = {
        "status": "ready",
        "input_csv": str(args.input_csv.resolve()),
        "input_csv_sha256": sha256(args.input_csv),
        "selected_count": len(rows),
        "shard_size": args.shard_size,
        "shard_count": len(shards),
        "array_max_index": max(item["rows"] for item in shards) - 1,
        "shards": shards,
        "status_contract": "one checkpoint is expected for every manifest_index; completed_with_stage_failures remains visible and resumable",
    }
    output_manifest = args.output_dir / "shard_manifest.json"
    output_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: manifest[key] for key in ("status", "selected_count", "shard_count", "shard_size", "input_csv_sha256")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
