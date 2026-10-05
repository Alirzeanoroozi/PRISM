#!/usr/bin/env python3
"""Run one manifest task in an isolated KUTEM Slurm-array directory.

The manifest is intentionally separate from benchmark/evaluator manifests. It
is a small execution contract with one row for each array index 1--10:

    array_index,task_id,command,config_paths,input_paths,output_paths,scientific_retry_id

The fixed KUTEM resource profile uses an array upper bound of ten.  Individual
submission manifests may contain one through ten rows; the final partial
batch is submitted with ``--array=1-N`` while retaining the same resource
contract.

The *_paths fields are semicolon-separated. Config and input paths are
relative to the manifest directory. Relative output paths are required to stay
inside the task directory. ``command`` is recorded verbatim and executed by
``bash -lc`` with the task directory as its working directory.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
from typing import Any, Iterable, Mapping


ARRAY_START = 1
ARRAY_END = 10
TEMPLATE_NAME = "isolated_kutem_array.sbatch"
PROFILE: dict[str, Any] = {
    "partition": "kutem",
    "account": "kutem",
    "qos": "kutem",
    "nodes": 1,
    "ntasks": 1,
    "cpus_per_task": 2,
    "memory": "2G",
    "time": "00:05:00",
    "array": "1-10",
}
REQUIRED_COLUMNS = {
    "array_index",
    "task_id",
    "command",
    "config_paths",
    "input_paths",
    "output_paths",
}
SLURM_ID_KEYS = (
    "SLURM_JOB_ID",
    "SLURM_ARRAY_JOB_ID",
    "SLURM_ARRAY_TASK_ID",
    "SLURM_ARRAY_TASK_COUNT",
    "SLURM_STEP_ID",
    "SLURM_JOB_NAME",
    "SLURM_CLUSTER_NAME",
    "SLURM_NODELIST",
    "SLURM_JOB_NODELIST",
    "SLURM_SUBMIT_DIR",
)
SLURM_RESOURCE_KEYS = (
    "SLURM_PARTITION",
    "SLURM_ACCOUNT",
    "SLURM_QOS",
    "SLURM_NNODES",
    "SLURM_JOB_NUM_NODES",
    "SLURM_NTASKS",
    "SLURM_NTASKS_PER_NODE",
    "SLURM_CPUS_PER_TASK",
    "SLURM_JOB_CPUS_PER_NODE",
    "SLURM_MEM_PER_NODE",
    "SLURM_MEM_PER_CPU",
    "SLURM_TIMELIMIT",
    "SLURM_JOB_TIME_LIMIT",
    "SLURM_RESTART_COUNT",
)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def parse_path_list(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(";") if item.strip()]


def validate_array_size(array_size: int) -> int:
    if array_size < ARRAY_START or array_size > ARRAY_END:
        raise ValueError(f"array_size must be between {ARRAY_START} and {ARRAY_END}, got {array_size}")
    return array_size


def infer_array_size(rows: list[dict[str, str]] | None = None, explicit: int | None = None) -> int:
    if explicit is not None:
        return validate_array_size(explicit)
    for variable in ("PRISM_ARRAY_SIZE", "SLURM_ARRAY_TASK_COUNT"):
        raw = os.environ.get(variable)
        if raw:
            return validate_array_size(int(raw))
    if rows:
        return validate_array_size(max(int(row["array_index"]) for row in rows))
    return ARRAY_END


def load_manifest(path: Path, *, array_size: int | None = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS - fieldnames
        if missing:
            raise ValueError(f"manifest is missing required columns: {', '.join(sorted(missing))}")
        rows = [{key: (value or "") for key, value in row.items() if key is not None} for row in reader]

    expected_size = infer_array_size(rows, array_size)
    expected = set(range(ARRAY_START, expected_size + 1))
    by_index: dict[int, dict[str, str]] = {}
    task_ids: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        try:
            index = int(row["array_index"])
        except ValueError as exc:
            raise ValueError(f"manifest row {row_number} has a non-integer array_index") from exc
        if index in by_index:
            raise ValueError(f"manifest contains duplicate array_index={index}")
        if index not in expected:
            raise ValueError(f"array_index={index} is outside the required 1-{expected_size} array")
        task_id = row["task_id"].strip()
        if not task_id:
            raise ValueError(f"manifest row {row_number} has an empty task_id")
        if task_id in task_ids:
            raise ValueError(f"manifest contains duplicate task_id={task_id!r}")
        if not row["command"]:
            raise ValueError(f"manifest row {row_number} has an empty command")
        by_index[index] = row
        task_ids.add(task_id)

    if set(by_index) != expected:
        missing = sorted(expected - set(by_index))
        raise ValueError(f"manifest must contain exactly one row for every array index 1-{expected_size}; missing {missing}")
    return [by_index[index] for index in range(ARRAY_START, expected_size + 1)]


def task_directory(run_root: Path, array_index: int) -> Path:
    if array_index not in range(ARRAY_START, ARRAY_END + 1):
        raise ValueError(f"array_index must be in 1-10, got {array_index}")
    return (run_root / "tasks" / f"task_{array_index:04d}").resolve()


def resolve_manifest_path(raw_path: str, manifest_dir: Path) -> Path:
    path = Path(raw_path).expanduser()
    return path.resolve() if path.is_absolute() else (manifest_dir / path).resolve()


def resolve_output_path(raw_path: str, task_dir: Path) -> Path:
    path = Path(raw_path).expanduser()
    if path.is_absolute():
        raise ValueError(f"output path must be relative to the task directory: {raw_path!r}")
    resolved = (task_dir / path).resolve()
    try:
        resolved.relative_to(task_dir.resolve())
    except ValueError as exc:
        raise ValueError(f"output path escapes the task directory: {raw_path!r}") from exc
    return resolved


def file_records(paths: Iterable[tuple[str, Path]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for raw_path, path in paths:
        exists = path.is_file()
        records.append(
            {
                "path": raw_path,
                "resolved_path": str(path),
                "exists": exists,
                "size_bytes": path.stat().st_size if exists else None,
                "sha256": sha256_file(path) if exists else None,
            }
        )
    return records


def aggregate_hash(records: list[dict[str, Any]]) -> str:
    return sha256_bytes(canonical_json(records))


def observed_slurm_metadata(env: Mapping[str, str] | None = None) -> dict[str, Any]:
    source = os.environ if env is None else env
    ids = {key: source.get(key) for key in SLURM_ID_KEYS}
    resources = {key: source.get(key) for key in SLURM_RESOURCE_KEYS}
    return {"ids": ids, "resources": resources}


def scheduler_retry_id(env: Mapping[str, str] | None = None) -> str:
    source = os.environ if env is None else env
    explicit = source.get("PRISM_SCHEDULER_RETRY_ID")
    if explicit:
        return explicit
    job_id = source.get("SLURM_ARRAY_JOB_ID") or source.get("SLURM_JOB_ID")
    restart = source.get("SLURM_RESTART_COUNT") or "0"
    return f"{job_id}:restart:{restart}" if job_id else "unscheduled:0"


def requested_submission(manifest: Path, run_root: Path, template: Path, *, array_size: int | None = None) -> list[str]:
    repo_root = Path(__file__).resolve().parents[2]
    rows = load_manifest(manifest, array_size=array_size)
    expected_size = infer_array_size(rows, array_size)
    values = {
        "MANIFEST": str(manifest),
        "RUN_ROOT": str(run_root),
        "PRISM_REPO_ROOT": str(repo_root),
        "PRISM_RUNNER_PYTHON": os.environ.get("PRISM_RUNNER_PYTHON", "python3"),
        "PRISM_SCHEDULER_RETRY_ID": f"submission-{run_root.name}",
    }
    if any("," in value or "\n" in value for value in values.values()):
        raise ValueError("manifest, run-root, and repository paths must not contain commas or newlines")
    export = "ALL," + ",".join(f"{key}={value}" for key, value in values.items())
    return ["sbatch", f"--array=1-{expected_size}", f"--export={export}", str(template)]


def dry_run_plan(manifest: Path, run_root: Path, template: Path, *, array_size: int | None = None) -> dict[str, Any]:
    rows = load_manifest(manifest, array_size=array_size)
    manifest_dir = manifest.parent
    task_plans = []
    for row in rows:
        index = int(row["array_index"])
        task_dir = task_directory(run_root, index)
        config_paths = parse_path_list(row.get("config_paths"))
        input_paths = parse_path_list(row.get("input_paths"))
        output_paths = parse_path_list(row.get("output_paths"))
        task_plans.append(
            {
                "array_index": index,
                "task_id": row["task_id"],
                "task_directory": str(task_dir),
                "command_sha256": sha256_bytes(row["command"].encode("utf-8")),
                "config_paths": [str(resolve_manifest_path(item, manifest_dir)) for item in config_paths],
                "input_paths": [str(resolve_manifest_path(item, manifest_dir)) for item in input_paths],
                "output_paths": [str(resolve_output_path(item, task_dir)) for item in output_paths],
                "scientific_retry_id": row.get("scientific_retry_id", "0") or "0",
            }
        )
    return {
        "mode": "dry-run",
        "submits_job": False,
        "profile": PROFILE,
        "manifest": str(manifest),
        "manifest_sha256": sha256_file(manifest),
        "run_root": str(run_root),
        "array_size": len(rows),
        "submission_command": requested_submission(manifest, run_root, template, array_size=len(rows)),
        "tasks": task_plans,
    }


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def rss_kb() -> dict[str, Any]:
    try:
        value = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    except (AttributeError, OSError):
        return {"value": None, "unit": "KB", "source": None}
    return {"value": value, "unit": "KB", "source": "resource.getrusage(RUSAGE_CHILDREN)"}


def exit_record_base(
    *,
    row: Mapping[str, str],
    index: int,
    manifest: Path,
    run_root: Path,
    task_dir: Path,
    started_at: str,
    env: Mapping[str, str],
) -> dict[str, Any]:
    manifest_dir = manifest.parent
    config_paths = parse_path_list(row.get("config_paths"))
    input_paths = parse_path_list(row.get("input_paths"))
    output_paths = parse_path_list(row.get("output_paths"))
    config_files = [(item, resolve_manifest_path(item, manifest_dir)) for item in config_paths]
    input_files = [(item, resolve_manifest_path(item, manifest_dir)) for item in input_paths]
    output_files = [(item, resolve_output_path(item, task_dir)) for item in output_paths]
    config_records = file_records(config_files)
    input_records = file_records(input_files)
    output_records = file_records(output_files)
    config_hash = aggregate_hash(config_records)
    input_hash = aggregate_hash(input_records)
    output_hash = aggregate_hash(output_records)
    return {
        "schema_version": "isolated-kutem-exit/v1",
        "timestamps": {"started_at": started_at},
        "identity": {
            "array_index": index,
            "task_id": row["task_id"],
            "manifest": str(manifest),
            "manifest_sha256": sha256_file(manifest),
            "run_root": str(run_root),
            "task_directory": str(task_dir),
        },
        "command": {
            "shell": "bash -lc",
            "text": row["command"],
            "sha256": sha256_bytes(row["command"].encode("utf-8")),
        },
        "hashes": {
            "command_sha256": sha256_bytes(row["command"].encode("utf-8")),
            "config_sha256": config_hash,
            "input_sha256": input_hash,
            "output_sha256": output_hash,
        },
        "runtime": {
            "runner_script": str(Path(__file__).resolve()),
            "runner_script_sha256": sha256_file(Path(__file__).resolve()),
            "python_executable": sys.executable,
            "python_version": sys.version,
        },
        "artifacts": {
            "config": {"files": config_records, "aggregate_sha256": config_hash},
            "inputs": {"files": input_records, "aggregate_sha256": input_hash},
            "outputs": {"files": output_records, "aggregate_sha256": output_hash},
        },
        "termination": {"return_code": None, "signal": None},
        "resource_usage": {"rss": {"value": None, "unit": "KB", "source": None}},
        "slurm": {"requested": PROFILE, "observed": observed_slurm_metadata(env)},
        "retry_ids": {
            "scientific_retry_id": row.get("scientific_retry_id", "0") or "0",
            "scheduler_retry_id": scheduler_retry_id(env),
        },
        "scientific_result": {
            "pair_success": None,
            "status": "unknown",
            "reason": "runner_does_not_infer_pair_success_from_batch_completion",
        },
        "execution": {"status": "started", "stdout": str(task_dir / "stdout.log"), "stderr": str(task_dir / "stderr.log")},
    }


def execute_task(
    *,
    manifest: Path,
    run_root: Path,
    array_index: int,
    array_size: int | None = None,
    allow_existing_task_dir: bool = False,
    env: Mapping[str, str] | None = None,
) -> tuple[int, Path]:
    manifest = manifest.resolve()
    run_root = run_root.resolve()
    rows = load_manifest(manifest, array_size=array_size)
    row = next(item for item in rows if int(item["array_index"]) == array_index)
    task_dir = task_directory(run_root, array_index)
    if task_dir.exists() and not allow_existing_task_dir:
        raise FileExistsError(f"task directory already exists; use a new run root or --allow-existing-task-dir: {task_dir}")
    if task_dir.exists() and allow_existing_task_dir:
        stale = [
            task_dir / "exit.json",
            task_dir / "stdout.log",
            task_dir / "stderr.log",
        ]
        stale.extend(resolve_output_path(item, task_dir) for item in parse_path_list(row.get("output_paths")))
        existing = [str(path) for path in stale if path.exists()]
        if existing:
            raise FileExistsError(
                "existing task directory contains prior outputs/logs; use a fresh run root: "
                + ", ".join(existing)
            )
    task_dir.mkdir(parents=True, exist_ok=allow_existing_task_dir)
    started_at = utc_now()
    child_env = dict(os.environ if env is None else env)
    child_env.update(
        {
            "PRISM_TASK_DIR": str(task_dir),
            "PRISM_ARRAY_INDEX": str(array_index),
            "PRISM_TASK_ID": row["task_id"],
            "PRISM_SCIENTIFIC_RETRY_ID": row.get("scientific_retry_id", "0") or "0",
            "PRISM_SCHEDULER_RETRY_ID": scheduler_retry_id(child_env),
        }
    )
    record = exit_record_base(
        row=row,
        index=array_index,
        manifest=manifest,
        run_root=run_root,
        task_dir=task_dir,
        started_at=started_at,
        env=child_env,
    )
    exit_path = task_dir / "exit.json"
    stdout_path = task_dir / "stdout.log"
    stderr_path = task_dir / "stderr.log"
    stdout_path.touch()
    stderr_path.touch()
    return_code: int | None = None
    error: str | None = None
    try:
        missing = [
            item["path"]
            for category in ("config", "inputs")
            for item in record["artifacts"][category]["files"]
            if not item["exists"]
        ]
        if missing:
            raise FileNotFoundError(f"missing config/input files: {', '.join(missing)}")
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
            process = subprocess.Popen(
                ["bash", "-lc", row["command"]],
                cwd=task_dir,
                env=child_env,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
            return_code = process.wait()
    except (OSError, ValueError, FileNotFoundError) as exc:
        error = str(exc)
        return_code = None

    signal_name = None
    if return_code is not None and return_code < 0:
        try:
            signal_name = signal.Signals(-return_code).name
        except ValueError:
            signal_name = f"SIG{-return_code}"
    finished_at = utc_now()
    output_files = [(item, resolve_output_path(item, task_dir)) for item in parse_path_list(row.get("output_paths"))]
    output_records = file_records(output_files)
    output_hash = aggregate_hash(output_records)
    record["artifacts"]["outputs"] = {"files": output_records, "aggregate_sha256": output_hash}
    record["hashes"]["output_sha256"] = output_hash
    log_files = [("stdout.log", stdout_path), ("stderr.log", stderr_path)]
    log_records = file_records(log_files)
    record["artifacts"]["logs"] = {"files": log_records, "aggregate_sha256": aggregate_hash(log_records)}
    record["timestamps"].update({"finished_at": finished_at, "record_written_at": finished_at})
    record["termination"] = {"return_code": return_code, "signal": signal_name}
    record["resource_usage"] = {"rss": rss_kb()}
    record["execution"] = {
        "status": "completed" if return_code == 0 else "failed",
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }
    if error is not None:
        record["execution"]["error"] = error
    write_json(exit_path, record)
    return (0 if return_code == 0 else 1), exit_path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--template", type=Path, default=Path(__file__).resolve().parents[1] / "jobs" / TEMPLATE_NAME)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--dry-run", action="store_true", help="validate and print a plan; never execute or submit")
    modes.add_argument("--submit", action="store_true", help="submit the exact array template via sbatch")
    modes.add_argument("--run-task", action="store_true", help="execute one task, normally used by the Slurm template")
    parser.add_argument("--array-index", type=int, help="task index for --run-task; defaults to SLURM_ARRAY_TASK_ID")
    parser.add_argument("--array-size", type=int, help="manifest array size, one through ten; defaults to Slurm or manifest size")
    parser.add_argument("--allow-existing-task-dir", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    manifest = args.manifest.resolve()
    run_root = args.run_root.resolve()
    template = args.template.resolve()
    if args.dry_run:
        print(json.dumps(dry_run_plan(manifest, run_root, template, array_size=args.array_size), indent=2, sort_keys=True))
        return 0
    if args.submit:
        completed = subprocess.run(
            requested_submission(manifest, run_root, template, array_size=args.array_size),
            check=False,
        )
        return completed.returncode
    index = args.array_index
    if index is None:
        raw_index = os.environ.get("SLURM_ARRAY_TASK_ID")
        if raw_index is None:
            raise SystemExit("--array-index or SLURM_ARRAY_TASK_ID is required with --run-task")
        index = int(raw_index)
    status, exit_path = execute_task(
        manifest=manifest,
        run_root=run_root,
        array_index=index,
        array_size=args.array_size,
        allow_existing_task_dir=args.allow_existing_task_dir,
    )
    print(exit_path)
    return status


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileExistsError, ValueError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
