#!/usr/bin/env python3
"""Delete only frozen upstream alignment directories, with resumable ledgers.

This utility is intentionally conservative.  It consumes the already-created
alignment cleanup plan, validates every target against the expected source
root, removes one exact directory at a time, and records both the attempt and
the final state.  The directory-level ledger is the authoritative skip list
for later rounds; it is deliberately independent of per-file logging so a
large directory can be handled without an expensive recursive Python walk.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Iterable


SOURCE_ROOT = Path(
    "/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919"
)
PLAN_PATH = Path(
    "/scratch/users/rshadi25/valar-remote-runs/cleanup-manifests/"
    "20260928-prism-intermediate-dry-run/alignment_cleanup_plan.json"
)
LEDGER_DIR = PLAN_PATH.parent
EVENTS_PATH = LEDGER_DIR / "alignment_cleanup_events.jsonl"
LEDGER_PATH = LEDGER_DIR / "deleted_alignment_directories.tsv"
SKIPLIST_PATH = LEDGER_DIR / "alignment_cleanup_skiplist.tsv"


def _event(event: dict[str, Any], handle: Any) -> None:
    event = {"timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **event}
    handle.write(json.dumps(event, sort_keys=True) + "\n")
    handle.flush()


def _path_from_value(value: Any) -> str | None:
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        for key in ("path", "directory", "alignment_directory", "target"):
            candidate = value.get(key)
            if isinstance(candidate, str):
                return candidate
    return None


def _planned_paths(payload: Any) -> list[Path]:
    """Read common planner output shapes without guessing arbitrary strings."""
    values: Iterable[Any]
    if isinstance(payload, list):
        values = payload
    elif isinstance(payload, dict):
        values = []
        for key in ("targets", "target_directories", "directories", "entries"):
            if isinstance(payload.get(key), list):
                values = payload[key]
                break
        else:
            raise ValueError("plan has no recognized target list")
    else:
        raise ValueError("plan must be a JSON object or list")

    result: list[Path] = []
    for value in values:
        path = _path_from_value(value)
        if path is None:
            raise ValueError(f"plan target has no path: {value!r}")
        result.append(Path(path))
    if not result:
        raise ValueError("plan contains no target directories")
    return result


def _validate(path: Path) -> None:
    resolved_root = SOURCE_ROOT.resolve(strict=False)
    resolved = path.resolve(strict=False)
    try:
        resolved.relative_to(resolved_root)
    except ValueError as exc:
        raise ValueError(f"target is outside source root: {path}") from exc
    parts = resolved.relative_to(resolved_root).parts
    if "processed" not in parts:
        raise ValueError(f"target is not under processed/: {path}")
    if not resolved.name.startswith("alignment_"):
        raise ValueError(f"target is not alignment_*: {path}")
    if resolved.is_symlink():
        raise ValueError(f"refusing symlink target: {path}")


def _identity(path: Path) -> tuple[str, str, str]:
    relative = path.resolve(strict=False).relative_to(SOURCE_ROOT.resolve(strict=False))
    parts = relative.parts
    pipeline = parts[0] if parts else "unknown"
    case = next((part for part in parts if re.fullmatch(r"case[_-].+", part)), "unknown")
    attempt = next((part for part in parts if re.fullmatch(r"attempt[_-].+", part)), "unknown")
    return pipeline, case, attempt


def _completed_paths() -> set[str]:
    completed: set[str] = set()
    if not LEDGER_PATH.exists():
        return completed
    for line in LEDGER_PATH.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) >= 7 and fields[4] in {"deleted", "already_absent"}:
            completed.add(fields[3])
    return completed


def _refresh_skiplist() -> None:
    """Materialize only completed paths for the next processing round."""
    paths: set[str] = set()
    if LEDGER_PATH.exists():
        for line in LEDGER_PATH.read_text(encoding="utf-8", errors="replace").splitlines():
            fields = line.split("\t")
            if len(fields) >= 7 and fields[4] in {"deleted", "already_absent"}:
                paths.add(fields[3])
    temporary = SKIPLIST_PATH.with_suffix(".tmp")
    temporary.write_text("\n".join(sorted(paths)) + ("\n" if paths else ""), encoding="utf-8")
    temporary.replace(SKIPLIST_PATH)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=PLAN_PATH)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    plan_path = args.plan
    payload = json.loads(plan_path.read_text(encoding="utf-8"))
    paths = _planned_paths(payload)
    for path in paths:
        _validate(path)
    if len({str(path) for path in paths}) != len(paths):
        raise ValueError("plan contains duplicate target paths")

    if args.dry_run:
        print(f"validated {len(paths)} target directories", flush=True)
        return 0

    _refresh_skiplist()
    done = _completed_paths()
    remaining = [path for path in paths if str(path) not in done or path.exists()]
    with EVENTS_PATH.open("a", encoding="utf-8") as events, LEDGER_PATH.open("a", encoding="utf-8") as ledger:
        _event({"event": "fast_cleanup_started", "target_count": len(paths), "remaining_count": len(remaining)}, events)
        processed = 0
        for path in remaining:
            pipeline, case, attempt = _identity(path)
            if not path.exists() and not path.is_symlink():
                status = "already_absent"
                return_code = 0
            else:
                _event({"event": "deleting_directory", "alignment_directory": str(path), "pipeline": pipeline}, events)
                completed_process = subprocess.run(["rm", "-rf", "--", str(path)], check=False)
                return_code = completed_process.returncode
                status = "deleted" if return_code == 0 and not path.exists() else "failed"

            record = {
                "pipeline": pipeline,
                "case": case,
                "attempt": attempt,
                "alignment_directory": str(path),
                "status": status,
                "return_code": return_code,
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
            ledger.write("\t".join(str(record[key]) for key in (
                "pipeline", "case", "attempt", "alignment_directory", "status", "return_code", "timestamp"
            )) + "\n")
            ledger.flush()
            if status in {"deleted", "already_absent"}:
                with SKIPLIST_PATH.open("a", encoding="utf-8") as skiplist:
                    skiplist.write(str(path) + "\n")
                    skiplist.flush()
            _event({"event": "directory_result", **record}, events)
            processed += 1
            if status == "failed":
                print(f"FAILED {path}", file=sys.stderr, flush=True)
                return 2
            if processed % 25 == 0:
                print(f"processed {processed}/{len(remaining)}", flush=True)

        _event({"event": "fast_cleanup_finished", "target_count": len(paths), "processed": processed}, events)
    print(f"completed {processed}/{len(remaining)} remaining target directories; {len(paths) - len(remaining)} already logged", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
