#!/usr/bin/env python3
"""Merge independently hashed PRISM verification baseline shards."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({field for row in rows for field in row})
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def merge_baseline_shards(shard_dirs: list[Path], output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to mix baseline output into non-empty directory: {output_dir}")
    claims: list[dict[str, str]] = []
    artifacts: list[dict[str, str]] = []
    for shard in shard_dirs:
        claims_path = shard / "claims.tsv"
        artifacts_path = shard / "artifact_manifest.tsv"
        if not claims_path.is_file() or not artifacts_path.is_file():
            raise FileNotFoundError(f"incomplete baseline shard: {shard}")
        claims.extend(_read_tsv(claims_path))
        artifacts.extend(_read_tsv(artifacts_path))
    claim_ids = [row.get("claim_id", "") for row in claims]
    if len(set(claim_ids)) != len(claim_ids):
        raise ValueError("duplicate claim_id across baseline shards")
    artifact_keys = [(row.get("run_root", ""), row.get("relative_path", "")) for row in artifacts]
    if len(set(artifact_keys)) != len(artifact_keys):
        raise ValueError("duplicate artifact path across baseline shards")
    output_dir.mkdir(parents=True, exist_ok=False)
    _write_tsv(output_dir / "claims.tsv", sorted(claims, key=lambda row: row["claim_id"]))
    _write_tsv(
        output_dir / "artifact_manifest.tsv",
        sorted(artifacts, key=lambda row: (row["run_root"], row["relative_path"])),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard", action="append", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    merge_baseline_shards(args.shard, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
