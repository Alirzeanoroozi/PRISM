#!/usr/bin/env python3
"""Prepare isolated shards for a fail-closed replay of existing model outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark.scripts.score_comparison_models import build_manifest, read_csv, write_rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare(
    manifest: Path,
    current_roots: list[Path],
    legacy_root: Path,
    output_root: Path,
    shards: int,
    existing_model_manifest: Path | None = None,
) -> dict:
    if shards < 1:
        raise ValueError("shards must be positive")
    output_root.mkdir(parents=True, exist_ok=False)
    model_manifest = output_root / "model_manifest.csv"
    if existing_model_manifest is None:
        build_manifest(manifest, current_roots, legacy_root, model_manifest)
    else:
        shutil.copy2(existing_model_manifest, model_manifest)
    rows = sorted(read_csv(model_manifest), key=lambda row: (row.get("pipeline", ""), row.get("model_path", "")))
    shard_rows = [[] for _ in range(shards)]
    for index, row in enumerate(rows):
        shard_rows[index % shards].append(row)
    shard_dir = output_root / "shards"
    shard_dir.mkdir()
    records = []
    for index, shard in enumerate(shard_rows, start=1):
        path = shard_dir / f"shard_{index:02d}.csv"
        write_rows(path, shard or [{"status": ""}])
        records.append({"shard": index, "path": str(path), "rows": len(shard), "sha256": sha256(path)})
    value = {
        "schema_version": "observational-score-replay/v1",
        "status": "prepared_observational_replay",
        "source_manifest": str(manifest.resolve()),
        "source_manifest_sha256": sha256(manifest),
        "current_roots": [str(path.resolve()) for path in current_roots],
        "legacy_root": str(legacy_root.resolve()),
        "model_manifest": str(model_manifest.resolve()),
        "model_manifest_sha256": sha256(model_manifest),
        "existing_model_manifest": str(existing_model_manifest.resolve()) if existing_model_manifest else None,
        "shards": records,
        "shard_count": shards,
    }
    (output_root / "replay_manifest.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--current-root", type=Path, action="append", required=True)
    parser.add_argument("--legacy-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--shards", type=int, default=10)
    parser.add_argument("--existing-model-manifest", type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.manifest, args.current_root, args.legacy_root, args.output_root, args.shards, args.existing_model_manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
