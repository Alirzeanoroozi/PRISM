#!/usr/bin/env python3
"""Submit one isolated KUTEM array per source-gate manifest batch."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path


FIELDS = (
    "batch_id",
    "manifest",
    "run_root",
    "array_size",
    "command",
    "return_code",
    "job_id",
    "stdout",
    "stderr",
)


def submit(index_path: str | Path, output_dir: str | Path, *, python_bin: str, runner_script: str | Path, template: str | Path, actually_submit: bool = True) -> Path:
    index = Path(index_path).resolve()
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(index.open(newline="", encoding="utf-8"), delimiter="\t"))
    records = []
    for row in rows:
        manifest = Path(row["task_manifest"]).resolve()
        run_root = output / "runs" / row["batch_id"]
        command = [
            python_bin,
            str(Path(runner_script).resolve()),
            "--manifest",
            str(manifest),
            "--run-root",
            str(run_root),
            "--template",
            str(Path(template).resolve()),
            "--array-size",
            row["array_size"],
            "--submit",
        ]
        if actually_submit:
            completed = subprocess.run(
                command,
                check=False,
                text=True,
                capture_output=True,
                env={**os.environ, "PRISM_RUNNER_PYTHON": python_bin},
            )
            stdout = completed.stdout.strip()
            stderr = completed.stderr.strip()
            match = re.search(r"Submitted batch job (\d+)", stdout)
            job_id = match.group(1) if match else ""
            return_code = str(completed.returncode)
        else:
            stdout = "dry-run"
            stderr = ""
            job_id = ""
            return_code = "0"
        records.append({
            "batch_id": row["batch_id"],
            "manifest": str(manifest),
            "run_root": str(run_root),
            "array_size": row["array_size"],
            "command": " ".join(command),
            "return_code": return_code,
            "job_id": job_id,
            "stdout": stdout,
            "stderr": stderr,
        })
    submission_path = output / "array_submissions.tsv"
    with submission_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    if actually_submit and any(record["return_code"] != "0" for record in records):
        raise RuntimeError("one or more source-gate arrays failed submission; see array_submissions.tsv")
    return submission_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--python-bin", default="python3")
    parser.add_argument("--runner-script", type=Path, default=Path("benchmark/scripts/isolated_kutem_runner.py"))
    parser.add_argument("--template", type=Path, default=Path("benchmark/jobs/isolated_kutem_array.sbatch"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    try:
        path = submit(
            args.index,
            args.output_dir,
            python_bin=args.python_bin,
            runner_script=args.runner_script,
            template=args.template,
            actually_submit=not args.dry_run,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    print(f"wrote source array submission manifest: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
