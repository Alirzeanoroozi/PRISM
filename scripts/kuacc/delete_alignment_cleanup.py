#!/usr/bin/env python3
"""Delete only exact upstream processed/alignment_* directories with a ledger."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path


def target_directories(source_root: Path) -> list[tuple[str, Path]]:
    targets: list[tuple[str, Path]] = []
    for pipeline in ("tmalign", "multiprot", "usalign"):
        cases = source_root / pipeline / "cases"
        if not cases.is_dir():
            continue
        for case_dir in cases.iterdir():
            attempts = case_dir / "attempts"
            if not attempts.is_dir():
                continue
            for attempt in attempts.iterdir():
                processed = attempt / "processed"
                if not processed.is_dir():
                    continue
                for child in processed.iterdir():
                    if child.is_dir() and not child.is_symlink() and child.name.startswith("alignment_"):
                        targets.append((pipeline, child))
    return sorted(targets, key=lambda item: str(item[1]))


def file_inventory(directory: Path) -> tuple[int, int, list[str]]:
    total = 0
    count = 0
    files: list[str] = []
    for root, _, names in os.walk(directory, followlinks=False):
        for name in names:
            path = Path(root) / name
            try:
                total += path.stat().st_size
            except OSError:
                continue
            count += 1
            files.append(str(path))
    return total, count, files


def delete_tree_with_log(directory: Path, filelog) -> tuple[int, int]:
    """Delete one exact tree while logging each file before unlinking it."""
    total = 0
    count = 0
    for entry in list(os.scandir(directory)):
        path = Path(entry.path)
        if entry.is_symlink():
            raise RuntimeError(f"refusing to delete symlink inside target: {path}")
        if entry.is_dir(follow_symlinks=False):
            child_bytes, child_count = delete_tree_with_log(path, filelog)
            total += child_bytes
            count += child_count
            path.rmdir()
            continue
        try:
            size = entry.stat(follow_symlinks=False).st_size
        except OSError:
            size = 0
        filelog.write(f"{directory}\t{path}\tDELETED\n")
        filelog.flush()
        path.unlink()
        total += size
        count += 1
    return total, count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    args.log_dir.mkdir(parents=True, exist_ok=True)
    targets = target_directories(args.source_root.resolve())
    plan_path = args.log_dir / "alignment_cleanup_plan.json"
    skiplist_path = args.log_dir / "alignment_cleanup_skiplist.tsv"
    filelog_path = args.log_dir / "deleted_aligned_files.tsv"
    eventlog_path = args.log_dir / "alignment_cleanup_events.jsonl"

    plan = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_root": str(args.source_root.resolve()),
        "exact_pattern": "*/cases/*/attempts/*/processed/alignment_*",
        "target_directory_count": len(targets),
        "execute": args.execute,
        "protected": [
            "processed/transformation",
            "processed/pdbs",
            "processed/rosetta_refinement",
            "processed/fiberdock_refinement",
            "processed/candidate_audit",
            "all other source and large-refinement roots",
        ],
        "targets": [str(path) for _, path in targets],
    }
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    mode = "a" if args.execute else "w"
    with skiplist_path.open(mode, encoding="utf-8") as skiplist, filelog_path.open(mode, encoding="utf-8") as filelog, eventlog_path.open(mode, encoding="utf-8") as events:
        if not args.execute:
            skiplist.write("pipeline\talignment_directory\tstatus\tfile_count\tbytes\n")
            filelog.write("alignment_directory\tfile_path\tstatus\n")
        total_bytes = 0
        total_files = 0
        for pipeline, directory in targets:
            if not args.execute:
                size, count, files = file_inventory(directory)
                total_bytes += size
                total_files += count
                skiplist.write(f"{pipeline}\t{directory}\tPLANNED\t{count}\t{size}\n")
                for path in files:
                    filelog.write(f"{directory}\t{path}\tPLANNED\n")
                continue

            timestamp = datetime.now(timezone.utc).isoformat()
            event = {
                "timestamp": timestamp,
                "pipeline": pipeline,
                "alignment_directory": str(directory),
                "file_count": None,
                "bytes": None,
                "status": "deleting",
            }
            events.write(json.dumps(event, sort_keys=True) + "\n")
            events.flush()
            try:
                size, count = delete_tree_with_log(directory, filelog)
                directory.rmdir()
            except OSError as exc:
                event["status"] = "delete_failed"
                event["error"] = repr(exc)
                events.write(json.dumps(event, sort_keys=True) + "\n")
                events.flush()
                continue
            except RuntimeError as exc:
                event["status"] = "delete_refused"
                event["error"] = repr(exc)
                events.write(json.dumps(event, sort_keys=True) + "\n")
                events.flush()
                continue
            event["file_count"] = count
            event["bytes"] = size
            total_bytes += size
            total_files += count
            skiplist.write(f"{pipeline}\t{directory}\tDELETED\t{count}\t{size}\n")
            event["status"] = "deleted"
            events.write(json.dumps(event, sort_keys=True) + "\n")
            events.flush()

    print(json.dumps({"execute": args.execute, "target_directories": len(targets), "files": total_files, "bytes": total_bytes}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
