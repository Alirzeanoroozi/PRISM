#!/usr/bin/env python3
"""Build a conservative, non-destructive cleanup manifest.

The script never deletes or moves files.  Candidates must be supplied
explicitly so benchmark inputs, references, and validated evidence cannot be
selected by a broad wildcard.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


PROTECTED_PARTS = {
    "benchmark/data",
    "references",
    "working_version",
    "runtime_manifest.json",
    "environment.yaml",
    "environment.yml",
}


def is_protected(relative: Path) -> bool:
    value = relative.as_posix()
    return any(value == part or value.startswith(part + "/") for part in PROTECTED_PARTS)


def size_bytes(path: Path) -> int:
    if path.is_file() or path.is_symlink():
        return path.lstat().st_size
    return sum(item.lstat().st_size for item in path.rglob("*") if item.is_file() or item.is_symlink())


def tracked(root: Path, relative: Path) -> bool:
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--error-unmatch", str(relative)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def build(root: Path, candidates: list[str]) -> list[dict[str, str]]:
    rows = []
    for candidate in candidates:
        path = (root / candidate).resolve()
        try:
            relative = path.relative_to(root.resolve())
        except ValueError:
            rows.append({"path": str(path), "kind": "outside-repository", "proposed_action": "preserve"})
            continue
        exists = path.exists() or path.is_symlink()
        protected = is_protected(relative)
        is_tracked = tracked(root, relative) if exists else False
        if protected or is_tracked:
            action = "preserve"
            reason = "protected benchmark/source/reference/environment or tracked path"
        elif relative.as_posix() in {".pytest_cache", "__pycache__"}:
            action = "delete-after-postcheck"
            reason = "reproducible generated cache"
        else:
            action = "review-or-quarantine"
            reason = "untracked prior-run artifact; retain until evidence coverage is proven"
        rows.append(
            {
                "path": str(relative),
                "kind": "directory" if path.is_dir() else "file",
                "exists": str(exists).lower(),
                "size_bytes": str(size_bytes(path) if exists else 0),
                "tracked": str(is_tracked).lower(),
                "proposed_action": action,
                "reason": reason,
                "evidence_replacement": "runtime_manifest.json, final exit.json, or retained report when applicable",
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("candidate", nargs="+")
    args = parser.parse_args()
    root = args.repo_root.resolve()
    rows = build(root, args.candidate)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ("path", "kind", "exists", "size_bytes", "tracked", "proposed_action", "reason", "evidence_replacement")
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote cleanup manifest: {args.output} ({len(rows)} candidates)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
