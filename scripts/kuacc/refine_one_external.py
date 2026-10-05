#!/usr/bin/env python3
"""Run one external-Rosetta candidate in an isolated batch work directory."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def append_event(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            handle.write(json.dumps(payload, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--pipeline-repo", required=True)
    parser.add_argument("--pipeline-python", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()

    root = Path(args.run_root).resolve()
    repo = Path(args.pipeline_repo).resolve()
    manifest_path = Path(args.manifest).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if Path(manifest.get("run_root", "")).resolve() != root:
        raise SystemExit("batch manifest run_root does not match RUN_ROOT")
    summary_path = root / "run_summary.json"
    if not summary_path.is_file() or sha256(summary_path) != manifest.get("parent_summary_sha256"):
        raise SystemExit("batch manifest parent summary hash does not match run root")
    active_manifest = root / "downstream" / "manifest.json"
    if active_manifest.is_file():
        active = json.loads(active_manifest.read_text(encoding="utf-8"))
        current_job = os.environ.get("SLURM_ARRAY_JOB_ID") or os.environ.get("SLURM_JOB_ID", "")
        if active.get("status") == "running" and active.get("job_id") != current_job:
            raise SystemExit(f"refusing overlap with active downstream job {active.get('job_id', '')}")
    entry = manifest["entries"][args.index]
    downstream = root / "downstream"
    checkpoint = downstream / "checkpoints" / "refinement" / f"{entry['key']}.json"
    if checkpoint.is_file():
        previous = json.loads(checkpoint.read_text(encoding="utf-8"))
        if previous.get("status") in {"completed", "terminal_no_output"}:
            print(json.dumps({"status": "resumed", "checkpoint": str(checkpoint)}, sort_keys=True))
            return 0

    workdir = downstream / "batch_work" / entry["key"]
    workdir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    status = "terminal_no_output"
    error = ""
    result = {"totalscore": "-", "interaction_score": "-", "structure": "-"}
    previous_cwd = Path.cwd()
    try:
        os.chdir(workdir)
        sys.path.insert(0, str(root))
        sys.path.insert(0, str(repo))
        from src.rosetta_refinement import calculate_energy

        totalscore, interaction_score, structure = calculate_energy(entry["left"], entry["right"])
        result = {"totalscore": totalscore, "interaction_score": interaction_score, "structure": "-"}
        if structure != "-" and Path(structure).is_file():
            status = "completed"
            result["structure"] = str(Path(structure).resolve())
    except Exception as exc:  # checkpoint failures for retry
        status = "failed"
        error = f"{type(exc).__name__}: {exc}"
    finally:
        os.chdir(previous_cwd)

    record = {
        "status": status,
        "batch_index": args.index,
        "batch_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", os.environ.get("SLURM_JOB_ID", "")),
        "batch_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
        "left": entry["left"],
        "right": entry["right"],
        "result": result,
        "error": error,
        "workdir": str(workdir),
        "started_at": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - started,
        "input_key": hashlib.sha256(f"{entry['left']}\n{entry['right']}".encode("utf-8")).hexdigest(),
    }
    atomic_json(checkpoint, record)
    append_event(downstream / "events.jsonl", {"stage": "refinement_batch", **record})
    print(json.dumps(record, indent=2, sort_keys=True))
    return 0 if status in {"completed", "terminal_no_output"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
