#!/usr/bin/env python3
"""Run or resume one full-template BM55 CPU case."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
try:
    from valar_agent.prism_cpu_batch import task_is_resumable, validate_case_summary, sha256_file
except ModuleNotFoundError:
    sys.path.insert(0, str(SCRIPT_DIR))
    from prism_cpu_batch import task_is_resumable, validate_case_summary, sha256_file


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def append_event(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def next_attempt(attempts_dir: Path) -> tuple[int, Path]:
    attempts_dir.mkdir(parents=True, exist_ok=True)
    numbers = []
    for path in attempts_dir.glob("attempt_*"):
        try:
            numbers.append(int(path.name.split("_", 1)[1]))
        except (IndexError, ValueError):
            continue
    number = max(numbers, default=0) + 1
    return number, attempts_dir / f"attempt_{number:04d}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-manifest", required=True)
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()

    manifest_path = Path(args.batch_manifest).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    entries = manifest["entries"]
    if args.index < 0 or args.index >= len(entries):
        raise SystemExit(f"index {args.index} is outside manifest task range")
    entry = entries[args.index]
    case_root = Path(entry["case_root"]).resolve()
    status_path = case_root / "task_status.json"
    events_path = case_root / "events.jsonl"
    status = json.loads(status_path.read_text(encoding="utf-8")) if status_path.is_file() else {}
    previous_summary = None
    if status.get("attempt_root"):
        summary_path = Path(status["attempt_root"]) / "run_summary.json"
        if summary_path.is_file():
            previous_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if task_is_resumable(status, previous_summary):
        append_event(events_path, {"event": "resumed", "at": now(), "attempt_root": status["attempt_root"]})
        print(json.dumps({"status": "resumed", "attempt_root": status["attempt_root"]}, sort_keys=True))
        return 0

    attempt_number, attempt_root = next_attempt(Path(entry["attempts_dir"]).resolve())
    started = time.perf_counter()
    atomic_json(
        status_path,
        {
            "status": "running",
            "index": args.index,
            "case_id": entry["case_id"],
            "attempt": attempt_number,
            "attempt_root": str(attempt_root),
            "job_id": os.environ.get("SLURM_ARRAY_JOB_ID", os.environ.get("SLURM_JOB_ID", "")),
            "started_at": now(),
        },
    )
    append_event(events_path, {"event": "started", "at": now(), "attempt": attempt_number})

    pipeline_python = Path(manifest["pipeline_python"]).resolve()
    helper = Path(manifest["helper"]).resolve()
    command = [
        str(pipeline_python),
        str(helper),
        "--aligner", manifest["aligner"],
        "--dataset-dir", entry["case_dataset_dir"],
        "--template-list", manifest["template_list"],
        "--template-limit", str(manifest["template_count"]),
        "--run-root", str(attempt_root),
        "--scenario", f"bm55_full_{entry['case_id']}",
        "--pipeline-repo", manifest["pipeline_repo"],
        "--pipeline-python", manifest["pipeline_python"],
        "--timed-runner", manifest["timed_runner"],
        "--usalign-wrapper", manifest["usalign_wrapper"],
        "--rank", "false",
        "--rank-method", "baseline",
        "--top-k", "5",
    ]
    environment = dict(os.environ)
    environment.update({
        "PYTHONNOUSERSITE": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    attempt_root.mkdir(parents=True, exist_ok=True)
    controller_log = case_root / "controller_logs" / f"attempt_{attempt_number:04d}.console.log"
    controller_log.parent.mkdir(parents=True, exist_ok=True)
    with controller_log.open("w", encoding="utf-8") as handle:
        process = subprocess.run(command, cwd=attempt_root, env=environment, stdout=handle, stderr=subprocess.STDOUT, check=False)
    elapsed = time.perf_counter() - started
    summary_path = attempt_root / "run_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.is_file() else {}
    expected = {
        "status": "completed",
        "template_count": manifest["template_count"],
        "template_sha256": manifest["template_sha256"],
        "rank": False,
        "rank_method": "baseline",
        "alignment_records": entry["expected_alignment_records"],
    }
    errors = validate_case_summary(summary, expected)
    final_status = "completed" if process.returncode == 0 and not errors else "failed"
    result = {
        "status": final_status,
        "index": args.index,
        "case_id": entry["case_id"],
        "attempt": attempt_number,
        "attempt_root": str(attempt_root),
        "return_code": process.returncode,
        "validation_errors": errors,
        "job_id": os.environ.get("SLURM_ARRAY_JOB_ID", os.environ.get("SLURM_JOB_ID", "")),
        "elapsed_seconds": elapsed,
        "summary_sha256": sha256_file(summary_path) if summary_path.is_file() else "",
        "alignment_records": summary.get("alignment_records", 0),
        "successful_alignment_records": summary.get("successful_alignment_records", 0),
        "failed_alignment_records": summary.get("failed_alignment_records", 0),
        "transformed_model_pairs": summary.get("transformed_model_pairs", 0),
        "updated_at": now(),
    }
    atomic_json(status_path, result)
    append_event(events_path, {"event": final_status, "at": now(), **result})
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if final_status == "completed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
