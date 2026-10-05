#!/usr/bin/env python3
"""Create a deterministic one-candidate-per-array-task refinement manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def pairs(root: Path) -> list[tuple[Path, Path]]:
    transformation = root / "processed" / "transformation"
    result = []
    for left in sorted(transformation.glob("*_L.pdb")):
        right = left.with_name(left.name[:-6] + "_R.pdb")
        if right.is_file():
            result.append((left, right))
    return result


def key(left: Path, right: Path) -> str:
    return hashlib.sha256(f"{left.name}\n{right.name}".encode("utf-8")).hexdigest()[:24]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--output", default="")
    args = parser.parse_args()
    root = Path(args.run_root).resolve()
    summary = json.loads((root / "run_summary.json").read_text(encoding="utf-8"))
    if summary.get("status") != "completed":
        raise SystemExit("parent run must be completed before batching refinement")
    if summary.get("rank") is True or summary.get("rank_method") == "prodigy":
        raise SystemExit("batch refinement requires ranking and PRODIGY to be disabled")
    candidates = pairs(root)
    entries = [
        {
            "index": index,
            "key": key(left, right),
            "left": str(left),
            "right": str(right),
        }
        for index, (left, right) in enumerate(candidates)
    ]
    output = Path(args.output).resolve() if args.output else root / "downstream" / "refinement_batch_manifest.json"
    atomic_json(
        output,
        {
            "status": "ready",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "run_root": str(root),
            "parent_summary_sha256": hashlib.sha256((root / "run_summary.json").read_bytes()).hexdigest(),
            "ranking": False,
            "prodigy": False,
            "task_unit": "one_transformed_pair_per_array_task",
            "task_count": len(entries),
            "entries": entries,
        },
    )
    print(json.dumps({"status": "ready", "task_count": len(entries), "manifest": str(output)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
