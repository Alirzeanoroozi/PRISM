#!/usr/bin/env python3
"""Safely remove only validated common-refinement scratch directories."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_aggregate(root: Path) -> dict[str, Any]:
    status_path = root / "aggregate" / "aggregation_status.json"
    if not status_path.is_file():
        raise RuntimeError(f"missing validated aggregation status: {status_path}")
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "complete" or not status.get("cleanup_eligible"):
        raise RuntimeError("aggregation is not cleanup-eligible")
    for key in ("results_path", "interfaces_path"):
        path = Path(str(status[key])).resolve()
        if not path.is_file() or path.parent != (root / "aggregate").resolve():
            raise RuntimeError(f"retained aggregate path is unsafe: {path}")
        hash_key = "results_sha256" if key == "results_path" else "interfaces_sha256"
        if sha256_file(path) != status.get(hash_key):
            raise RuntimeError(f"retained aggregate hash mismatch: {path}")
    return status


def _planned(root: Path) -> list[Path]:
    planned: list[Path] = []
    for relative in ("results", "adapter_inputs"):
        directory = root / relative
        if not directory.is_dir():
            continue
        for child in sorted(directory.iterdir()):
            if child.is_dir() or child.is_file():
                planned.append(child)
    return planned


def cleanup(root: Path, *, apply: bool) -> dict[str, Any]:
    root = root.resolve()
    if not root.name.endswith("_common_refinement_19855"):
        raise RuntimeError(f"refusing unexpected cleanup root: {root}")
    aggregate_status = _validate_aggregate(root)
    planned = _planned(root)
    record: dict[str, Any] = {
        "schema_version": "prism-common-refinement-cleanup/v1",
        "status": "applied" if apply else "dry_run",
        "root": str(root),
        "aggregate_status_path": str((root / "aggregate" / "aggregation_status.json").resolve()),
        "aggregate_results_sha256": aggregate_status.get("results_sha256"),
        "planned_paths": [str(path.relative_to(root)) for path in planned],
        "planned_count": len(planned),
        "deleted_paths": [],
        "finished_at": datetime.now(timezone.utc).isoformat(),
    }
    if apply:
        for path in planned:
            if path.is_dir():
                shutil.rmtree(path)
            elif path.is_file():
                path.unlink()
            record["deleted_paths"].append(str(path.relative_to(root)))
    manifest = root / "cleanup_manifest.json"
    manifest.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--apply", action="store_true", help="delete the validated run-scoped scratch")
    args = parser.parse_args()
    print(json.dumps(cleanup(args.root, apply=args.apply), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
