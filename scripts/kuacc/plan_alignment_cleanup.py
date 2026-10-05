#!/usr/bin/env python3
"""Plan cleanup of exact upstream alignment/intermediate directories only."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def tree_size(root: Path) -> tuple[int, int]:
    total = 0
    files = 0
    for directory, _, names in os.walk(root, followlinks=False):
        for name in names:
            try:
                total += (Path(directory) / name).stat().st_size
                files += 1
            except OSError:
                pass
    return total, files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    for pipeline in ("tmalign", "multiprot", "usalign"):
        cases = args.source_root / pipeline / "cases"
        if not cases.is_dir():
            continue
        for processed in cases.glob("*/attempts/*/processed"):
            if not processed.is_dir():
                continue
            for child in processed.iterdir():
                if not child.is_dir():
                    continue
                if child.name.startswith("alignment_"):
                    classification = "DELETE_CANDIDATE"
                    reason = "upstream alignment output; transformed files are in a separate transformation directory"
                    action = "review_then_delete_directory"
                elif child.name == "__pycache__":
                    classification = "DELETE_CANDIDATE"
                    reason = "recomputable Python cache"
                    action = "review_then_delete_directory"
                elif child.name == "surface_extraction":
                    classification = "REVIEW"
                    reason = "intermediate refinement input; retained until explicitly approved"
                    action = "manual_review"
                else:
                    continue
                size, files = tree_size(child)
                rows.append(
                    {
                        "path": str(child),
                        "pipeline": pipeline,
                        "size": size,
                        "file_count": files,
                        "classification": classification,
                        "reason": reason,
                        "suggested_action": action,
                    }
                )
    rows.sort(key=lambda row: (-int(row["size"]), str(row["path"])))
    summary = {
        "source_root": str(args.source_root.resolve()),
        "target_directory_count": len(rows),
        "target_file_count": sum(int(row["file_count"]) for row in rows),
        "bytes_by_class": {
            classification: sum(int(row["size"]) for row in rows if row["classification"] == classification)
            for classification in ("DELETE_CANDIDATE", "REVIEW")
        },
        "directories": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in ("target_directory_count", "target_file_count", "bytes_by_class")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
