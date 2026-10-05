#!/usr/bin/env python3
"""Prepare a VALAR fallback manifest and exact source-file transfer list."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--stage-root", required=True, type=Path)
    parser.add_argument("--files-list", required=True, type=Path)
    parser.add_argument("--path-map", required=True, type=Path)
    parser.add_argument("--start-index", required=True, type=int)
    args = parser.parse_args()
    if args.start_index < 0:
        raise SystemExit("--start-index must be non-negative")

    with args.input_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected: list[dict[str, str]] = []
    source_paths: set[str] = set()
    for index, row in enumerate(rows):
        if index < args.start_index:
            continue
        row = dict(row)
        row["origin_manifest_index"] = str(index)
        for field in ("left", "right", "native_pdb"):
            source = str(row[field])
            if not source.startswith("/"):
                raise SystemExit(f"non-absolute {field} at index {index}: {source}")
            source_paths.add(source)
            row[field] = str((args.stage_root / source.lstrip("/")).resolve())
        selected.append(row)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted(set().union(*(set(row) for row in selected)))
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(selected)

    args.files_list.parent.mkdir(parents=True, exist_ok=True)
    args.files_list.write_text(
        "".join(f"{source.lstrip('/')}\n" for source in sorted(source_paths)),
        encoding="utf-8",
    )
    with args.path_map.open("w", encoding="utf-8") as handle:
        handle.write("source\tstaged\n")
        for source in sorted(source_paths):
            staged = args.stage_root / source.lstrip("/")
            handle.write(f"{source}\t{staged.resolve()}\n")

    print(
        {
            "input_rows": len(rows),
            "selected_rows": len(selected),
            "first_origin_manifest_index": args.start_index,
            "last_origin_manifest_index": args.start_index + len(selected) - 1,
            "unique_source_files": len(source_paths),
            "output_csv": str(args.output_csv.resolve()),
            "files_list_sha256": hashlib.sha256(args.files_list.read_bytes()).hexdigest(),
        }
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
