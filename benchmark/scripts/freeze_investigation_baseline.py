#!/usr/bin/env python3
"""Freeze hashes and pair-level evidence for the observational July baseline.

This command never copies, rewrites, or interprets the existing comparison
results.  It records deterministic hashes for explicitly supplied files and
directories, together with the status label that prevents the baseline from
being used as a confirmatory estimate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable


FIELDS = (
    "kind",
    "path",
    "relative_path",
    "exists",
    "size_bytes",
    "sha256",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_files(path: Path) -> Iterable[Path]:
    if path.is_file():
        yield path
        return
    if path.is_dir():
        yield from sorted((item for item in path.rglob("*") if item.is_file()), key=lambda item: item.as_posix())


def freeze(paths: Iterable[str | Path], output_dir: str | Path, repo_root: str | Path) -> tuple[Path, Path]:
    output = Path(output_dir).resolve()
    root = Path(repo_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    seen: set[Path] = set()
    for raw in paths:
        path = Path(raw).expanduser().resolve()
        if not path.exists():
            records.append({
                "kind": "missing",
                "path": str(path),
                "relative_path": "",
                "exists": False,
                "size_bytes": "",
                "sha256": "",
            })
            continue
        for file_path in iter_files(path):
            if file_path in seen:
                continue
            seen.add(file_path)
            try:
                relative = file_path.relative_to(root).as_posix()
            except ValueError:
                relative = "external:" + str(file_path)
            records.append({
                "kind": "file",
                "path": str(file_path),
                "relative_path": relative,
                "exists": True,
                "size_bytes": file_path.stat().st_size,
                "sha256": sha256_file(file_path),
            })
    records.sort(key=lambda row: (str(row["path"]), str(row["kind"])))
    manifest = output / "baseline_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    summary = {
        "schema_version": "observational-baseline/v1",
        "status": "nonstandardized_observational_baseline",
        "confirmatory_use": False,
        "repo_root": str(root),
        "record_count": len(records),
        "existing_file_count": sum(bool(row["exists"]) for row in records),
        "missing_path_count": sum(not bool(row["exists"]) for row in records),
        "manifest": str(manifest),
    }
    summary_path = output / "baseline_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest, summary_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--path", action="append", required=True, help="file or directory to hash; repeat")
    args = parser.parse_args(argv)
    manifest, summary = freeze(args.path, args.output_dir, args.repo_root)
    print(f"wrote baseline manifest: {manifest}")
    print(f"wrote baseline summary: {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
