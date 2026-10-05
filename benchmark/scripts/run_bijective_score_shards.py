#!/usr/bin/env python3
"""Run strict BM5.5 scoring in internally parallel deterministic shards."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def merge_tsv(paths: list[Path], output: Path) -> int:
    rows: list[dict[str, str]] = []
    fields: list[str] = []
    for path in paths:
        for row in read_tsv(path):
            rows.append(row)
            for field in row:
                if field not in fields:
                    fields.append(field)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-manifest", type=Path, required=True)
    parser.add_argument("--native-root", type=Path, required=True)
    parser.add_argument("--score-python", type=Path, required=True)
    parser.add_argument("--irmsd-script", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("--workers must be positive")
    if args.output_root.exists() and any(args.output_root.iterdir()):
        parser.error(f"refusing to mix scoring shards into non-empty directory: {args.output_root}")

    scorer = Path(__file__).with_name("score_bijective_benchmark_models.py")
    commands: list[list[str]] = []
    model_paths: list[Path] = []
    interface_paths: list[Path] = []
    for shard in range(args.workers):
        shard_root = args.output_root / "shards" / f"shard_{shard:03d}"
        model_output = shard_root / "scores_models.tsv"
        interface_output = shard_root / "scores_interfaces.tsv"
        model_paths.append(model_output)
        interface_paths.append(interface_output)
        commands.append(
            [
                sys.executable,
                str(scorer),
                "--stage-manifest", str(args.stage_manifest.resolve()),
                "--native-root", str(args.native_root.resolve()),
                # Preserve the venv entry path. Path.resolve() follows its
                # Python symlink into the base interpreter and loses DockQ's
                # site-packages in clean Slurm shells.
                "--score-python", str(args.score_python.absolute()),
                "--irmsd-script", str(args.irmsd_script.resolve()),
                "--output", str(model_output.resolve()),
                "--interfaces-output", str(interface_output.resolve()),
                "--raw-json-dir", str((shard_root / "raw-json").resolve()),
                "--timeout", str(args.timeout),
                "--shard-count", str(args.workers),
                "--shard-index", str(shard),
            ]
        )

    def run(command: list[str]) -> dict[str, object]:
        result = subprocess.run(command, capture_output=True, text=True)
        return {
            "argv": command,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }

    started = datetime.now(timezone.utc).isoformat()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(run, commands))
    args.output_root.mkdir(parents=True, exist_ok=True)
    run_manifest = {
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "workers": args.workers,
        "stage_manifest": str(args.stage_manifest.resolve()),
        "results": results,
    }
    (args.output_root / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2) + "\n", encoding="utf-8")
    failures = [result for result in results if result["returncode"] != 0]
    if failures:
        print(f"failed_shards={len(failures)} manifest={args.output_root / 'run_manifest.json'}")
        return 2
    model_count = merge_tsv(model_paths, args.output_root / "scores_models.tsv")
    interface_count = merge_tsv(interface_paths, args.output_root / "scores_interfaces.tsv")
    print(f"models={model_count} interfaces={interface_count} workers={args.workers} output={args.output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
