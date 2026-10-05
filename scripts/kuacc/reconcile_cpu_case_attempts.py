#!/usr/bin/env python3
"""Correct validation-only pilot states after the full-panel count contract changes."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
try:
    from valar_agent.prism_cpu_batch import (
        build_batch_plan,
        reclassify_task_status,
        sha256_file,
    )
except ModuleNotFoundError:
    sys.path.insert(0, str(SCRIPT_DIR))
    from prism_cpu_batch import build_batch_plan, reclassify_task_status, sha256_file


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def ensure_source_link(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        if destination.resolve() != source.resolve():
            raise RuntimeError(f"source link points elsewhere: {destination}")
        return
    if destination.exists():
        raise RuntimeError(f"refusing to replace existing source path: {destination}")
    destination.symlink_to(source.resolve())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-manifest", required=True)
    parser.add_argument("--output-manifest", required=True)
    args = parser.parse_args()

    old_path = Path(args.batch_manifest).resolve()
    output_path = Path(args.output_manifest).resolve()
    old = json.loads(old_path.read_text(encoding="utf-8"))
    if output_path.exists():
        raise SystemExit(f"refusing to overwrite existing corrected manifest: {output_path}")
    plan = build_batch_plan(
        dataset_dir=Path(old["dataset_dir"]),
        template_list=Path(old["template_list"]),
        output_root=old_path.parent,
        aligner=old["aligner"],
        expected_template_count=old["template_count"],
    )
    corrected = {
        **old,
        **plan,
        "status": "ready_reconciled",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "supersedes_manifest": str(old_path),
        "superseded_manifest_sha256": sha256_file(old_path),
    }

    for entry in corrected["entries"]:
        source_dir = Path(entry["case_dataset_dir"]) / "sources"
        for query_name, query_source in zip(entry["query_source_names"], entry["query_sources"]):
            ensure_source_link(Path(query_source), source_dir / query_name)

    reused = 0
    for entry in corrected["entries"]:
        status_path = Path(entry["case_root"]) / "task_status.json"
        if not status_path.is_file():
            continue
        status = json.loads(status_path.read_text(encoding="utf-8"))
        attempt_root = Path(status.get("attempt_root", "")) if status.get("attempt_root") else None
        summary_path = attempt_root / "run_summary.json" if attempt_root else None
        if not summary_path or not summary_path.is_file():
            continue
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        expected = {
            "status": "completed",
            "template_count": corrected["template_count"],
            "template_sha256": corrected["template_sha256"],
            "rank": False,
            "rank_method": "baseline",
            "alignment_records": entry["expected_alignment_records"],
        }
        updated = reclassify_task_status(status, summary, expected)
        if updated is not None:
            updated["reconciled_at"] = datetime.now(timezone.utc).isoformat()
            atomic_json(status_path, updated)
            reused += 1

    corrected["reused_completed_attempts"] = reused
    atomic_json(output_path, corrected)
    print(json.dumps({"status": corrected["status"], "manifest": str(output_path), "reused": reused}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
