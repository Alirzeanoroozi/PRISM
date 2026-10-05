#!/usr/bin/env python3
"""Measure USalign worker/configuration scaling on a frozen template pilot.

The sweep is deliberately alignment-only.  It uses the shared PRISM parser,
retains both explicit TM-score normalizations, and writes compact per-run
metrics plus a bounded record ledger.  It does not perform transformation,
filtering, ranking, refinement, or DockQ evaluation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import tempfile
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in os.sys.path:
    os.sys.path.insert(0, str(REPO_ROOT))

from src.alignment import iter_bounded_results, parse_tmalign  # noqa: E402


def build_usalign_command(
    executable: str,
    query: str,
    reference: str,
    matrix_path: str,
    *,
    fast: bool = False,
) -> list[str]:
    command = [executable, query, reference]
    if fast:
        command.append("-fast")
    command.extend(["-outfmt", "-1", "-m", matrix_path])
    return command


def parse_worker_list(value: str) -> list[int]:
    workers = sorted({int(item.strip()) for item in value.split(",") if item.strip()})
    if not workers or any(worker < 1 for worker in workers):
        raise ValueError("workers must contain positive integers")
    return workers


def load_tasks(
    template_list: Path,
    interface_root: Path,
    query: Path,
    *,
    template_limit: int | None = None,
) -> list[dict[str, str]]:
    tasks: list[dict[str, str]] = []
    selected_templates = 0
    for raw_line in template_list.read_text().splitlines():
        template = raw_line.strip()
        if not template:
            continue
        if template_limit is not None and selected_templates >= template_limit:
            break
        selected_templates += 1
        if len(template) <= 4:
            raise ValueError(f"template entry has no chain suffix: {template!r}")
        for chain in template[4:]:
            reference = interface_root / f"{template}_{chain}_int.pdb"
            tasks.append(
                {
                    "query": str(query),
                    "template_id": template,
                    "chain": chain,
                    "reference": str(reference),
                }
            )
    return tasks


def _raw_sha256(matrix_path: Path, stdout: str) -> str:
    digest = hashlib.sha256()
    digest.update(matrix_path.read_bytes())
    digest.update(b"\0")
    digest.update(stdout.encode())
    return digest.hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run_one(
    task: dict[str, str],
    executable: str,
    output_dir: Path,
    scratch_root: Path,
    *,
    fast: bool,
) -> dict[str, Any]:
    started = time.perf_counter()
    template = task["template_id"]
    chain = task["chain"]
    result_path = output_dir / f"pilot_{template}_{chain}.json"
    if not Path(task["reference"]).is_file() or not Path(task["query"]).is_file():
        return {
            "template_id": template,
            "chain": chain,
            "status": "failure",
            "failure_reason": "query_or_interface_missing",
            "elapsed_seconds": time.perf_counter() - started,
        }

    try:
        with tempfile.TemporaryDirectory(prefix="usalign-", dir=scratch_root) as scratch:
            scratch_path = Path(scratch)
            matrix_path = scratch_path / "matrix.out"
            tm_path = scratch_path / "stdout.tm"
            command = build_usalign_command(
                executable,
                task["query"],
                task["reference"],
                str(matrix_path),
                fast=fast,
            )
            completed = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            tm_path.write_text(completed.stdout)
            if completed.returncode != 0 or not matrix_path.is_file():
                raise RuntimeError(
                    f"return_code={completed.returncode}: {(completed.stderr or '').strip()[:300]}"
                )
            parsed = parse_tmalign(
                task["query"],
                task["reference"],
                task["query"].split("/")[-1].removesuffix(".asa.pdb"),
                template,
                chain,
                str(matrix_path),
                str(tm_path),
                str(output_dir),
                raw_output_sha256=_raw_sha256(matrix_path, completed.stdout),
                return_code=completed.returncode,
                aligner_name="USalign",
            )
            # parse_tmalign writes the durable per-record JSON.  Reload it so
            # the compact ledger reflects the corrected parser contract.
            parsed_path = output_dir / (
                f"{task['query'].split('/')[-1].removesuffix('.asa.pdb')}_{template}_{chain}.json"
            )
            parsed = json.loads(parsed_path.read_text())
            parsed["template_id"] = template
            parsed["chain"] = chain
            parsed["elapsed_seconds"] = time.perf_counter() - started
            parsed["fast"] = fast
            parsed["execution_status"] = "success"
            parsed["status"] = parsed.get("status", "failure")
            result_path.write_text(json.dumps(parsed, sort_keys=True) + "\n")
            if parsed_path != result_path and parsed_path.exists():
                parsed_path.unlink()
            return parsed
    except (OSError, RuntimeError, ValueError, IndexError, json.JSONDecodeError) as exc:
        return {
            "template_id": template,
            "chain": chain,
            "execution_status": "failure",
            "status": "failure",
            "failure_reason": str(exc),
            "elapsed_seconds": time.perf_counter() - started,
        }


def run_configuration(
    tasks: list[dict[str, str]],
    executable: Path,
    output_root: Path,
    workers: int,
    *,
    fast: bool,
) -> dict[str, Any]:
    label = f"{'fast' if fast else 'default'}_w{workers}"
    config_root = output_root / label
    records_root = config_root / "records"
    scratch_root = config_root / "scratch"
    records_root.mkdir(parents=True, exist_ok=True)
    scratch_root.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    child_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    records: list[dict[str, Any]] = []

    def worker(task: dict[str, str]) -> dict[str, Any]:
        return _run_one(task, str(executable), records_root, scratch_root, fast=fast)

    # Use the shared bounded dispatcher so the sweep exercises the requested
    # worker count without materializing all subprocesses at once.
    for _task, record in iter_bounded_results(tasks, worker, workers):
        records.append(record)
    child_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    wall_seconds = time.perf_counter() - started
    success = sum(record.get("execution_status") == "success" for record in records)
    failures = len(records) - success
    summary = {
        "schema_version": "prism-usalign-worker-sweep/v1",
        "configuration": "fast" if fast else "default",
        "worker_count": workers,
        "expected_records": len(tasks),
        "success_records": success,
        "failure_records": failures,
        "records_per_second": success / wall_seconds if wall_seconds else 0.0,
        "wall_seconds": wall_seconds,
        "child_cpu_seconds": max(0.0, child_after.ru_utime - child_before.ru_utime),
        "child_system_seconds": max(0.0, child_after.ru_stime - child_before.ru_stime),
        "child_max_rss_kb": max(child_before.ru_maxrss, child_after.ru_maxrss),
        "status": "completed" if failures == 0 else "completed_with_failures",
    }
    (config_root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template-list", type=Path, required=True)
    parser.add_argument("--interface-root", type=Path, required=True)
    parser.add_argument("--query", type=Path, required=True)
    parser.add_argument("--usalign", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--template-limit", type=int)
    parser.add_argument("--workers", default="1,2,4,8,16")
    parser.add_argument("--configs", default="default,fast")
    args = parser.parse_args()

    precreated = {
        item.name for item in args.output_root.iterdir()
    } if args.output_root.exists() else set()
    unexpected = precreated - {"logs"}
    if unexpected:
        raise SystemExit(f"output root is not empty: {args.output_root}")
    args.output_root.mkdir(parents=True, exist_ok=True)
    tasks = load_tasks(
        args.template_list,
        args.interface_root,
        args.query,
        template_limit=args.template_limit,
    )
    workers = parse_worker_list(args.workers)
    configs = {item.strip() for item in args.configs.split(",") if item.strip()}
    unknown = configs - {"default", "fast"}
    if unknown:
        raise SystemExit(f"unknown configurations: {sorted(unknown)}")
    selected_template_ids = list(dict.fromkeys(task["template_id"] for task in tasks))
    selected_template_sha256 = hashlib.sha256(
        ("\n".join(selected_template_ids) + "\n").encode()
    ).hexdigest()
    version_probe = subprocess.run(
        [str(args.usalign)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, check=False, timeout=5,
    )
    manifest = {
        "schema_version": "prism-usalign-worker-sweep/v1",
        "template_list": str(args.template_list.resolve()),
        "interface_root": str(args.interface_root.resolve()),
        "query": str(args.query.resolve()),
        "usalign": str(args.usalign.resolve()),
        "usalign_sha256": sha256_file(args.usalign.resolve()),
        "usalign_version_probe": version_probe.stdout.splitlines()[:3],
        "template_list_sha256": sha256_file(args.template_list.resolve()),
        "selected_template_sha256": selected_template_sha256,
        "selected_template_count": len(selected_template_ids),
        "query_sha256": sha256_file(args.query.resolve()),
        "workers": workers,
        "configurations": sorted(configs),
        "task_count": len(tasks),
        "template_limit": args.template_limit,
    }
    (args.output_root / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    summaries = []
    for fast in (False, True):
        if ("fast" if fast else "default") in configs:
            for worker_count in workers:
                summaries.append(
                    run_configuration(
                        tasks,
                        args.usalign.resolve(),
                        args.output_root,
                        worker_count,
                        fast=fast,
                    )
                )
    (args.output_root / "summary.json").write_text(
        json.dumps({"schema_version": "prism-usalign-worker-sweep/v1", "runs": summaries}, indent=2, sort_keys=True)
        + "\n"
    )
    return 0 if all(item["failure_records"] == 0 for item in summaries) else 2


if __name__ == "__main__":
    raise SystemExit(main())
