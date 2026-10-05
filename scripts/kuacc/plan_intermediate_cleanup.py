#!/usr/bin/env python3
"""Create a conservative, non-destructive KUACC cleanup manifest."""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter
from pathlib import Path


TERMINAL = {"completed", "completed_with_stage_failures", "failed"}
KEEP_SUFFIXES = {".sbatch", ".slurm", ".jsonl"}
KEEP_NAMES = {"AGENTS.md", "README", "README.md", "LICENSE", "LICENSE.md"}


def under(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def active_protected_dirs(refinement_root: Path) -> set[Path]:
    protected: set[Path] = set()
    checkpoint_dir = refinement_root / "checkpoints"
    for checkpoint_path in checkpoint_dir.glob("*.json"):
        try:
            checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            continue
        if str(checkpoint.get("status", "")) in TERMINAL:
            continue
        candidate = checkpoint.get("candidate") or {}
        attempt_root = candidate.get("attempt_root")
        if attempt_root:
            protected.add(Path(str(attempt_root)).resolve())
        protected.add((refinement_root / "results" / str(checkpoint.get("key", ""))).resolve())
    return protected


def classify(path: Path, source_root: Path, refinement_root: Path, protected: set[Path]) -> tuple[str, str, str]:
    resolved = path.resolve()
    if any(under(resolved, directory) for directory in protected):
        return "KEEP", "active or nonterminal case", "do_not_delete"
    if path.name in KEEP_NAMES or path.suffix in KEEP_SUFFIXES:
        return "KEEP", "reproducibility/control file", "preserve"

    parts = {part.lower() for part in resolved.parts}
    if "checkpoints" in parts or "manifests" in parts or "manifest" in parts or "scripts" in parts:
        return "KEEP", "checkpoint/manifest/script", "preserve"
    if under(resolved, source_root):
        if "transformation" in parts:
            return "KEEP", "transformed structure", "preserve"
        if "alignment" in parts or "aligned" in parts or path.suffix.lower() in {".aln", ".ali"}:
            return "DELETE_CANDIDATE", "alignment/intermediate artifact outside protected case", "review_then_delete"
        return "REVIEW", "source-run file with unresolved provenance", "manual_review"
    if under(resolved, refinement_root):
        if path.suffix.lower() in {".pdb", ".json", ".csv", ".tsv", ".txt"}:
            return "KEEP", "transformed/refined/result/provenance artifact", "preserve"
        return "REVIEW", "refinement-run file with unresolved provenance", "manual_review"
    return "REVIEW", "outside declared roots", "manual_review"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--refinement-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    source_root = args.source_root.resolve()
    refinement_root = args.refinement_root.resolve()
    protected = active_protected_dirs(refinement_root)
    rows: list[dict[str, object]] = []
    for root in (source_root, refinement_root):
        if not root.is_dir():
            continue
        for directory, _, filenames in os.walk(root, followlinks=False):
            for name in filenames:
                path = Path(directory) / name
                try:
                    size = path.stat().st_size
                except OSError:
                    continue
                classification, reason, action = classify(path, source_root, refinement_root, protected)
                rows.append(
                    {
                        "path": str(path),
                        "size": size,
                        "classification": classification,
                        "reason": reason,
                        "suggested_action": action,
                    }
                )
    rows.sort(key=lambda row: (-int(row["size"]), str(row["path"])))
    summary = {
        "source_root": str(source_root),
        "refinement_root": str(refinement_root),
        "protected_active_directories": sorted(str(path) for path in protected),
        "file_count": len(rows),
        "class_counts": dict(Counter(str(row["classification"]) for row in rows)),
        "bytes_by_class": {
            classification: sum(int(row["size"]) for row in rows if row["classification"] == classification)
            for classification in ("KEEP", "REVIEW", "DELETE_CANDIDATE")
        },
        "files": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in ("file_count", "class_counts", "bytes_by_class")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
