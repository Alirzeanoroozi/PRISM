#!/usr/bin/env python3
"""Prepare the non-overlapping GTalign cross-partition repair manifest."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SELECTED = ROOT / "selected_candidates.csv"
OUTPUT = ROOT / "repair" / "cross_partition_candidates.csv"
MANIFEST = ROOT / "repair" / "manifest.json"
WRAPPER = ROOT / "worker" / "run_corrected_chain_partition_refinement.py"


def clean(value: str) -> str:
    return "".join(str(value).split())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    with SELECTED.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = []
    patterns = Counter()
    for row in rows:
        source_left = clean(row.get("chain_left", ""))
        source_right = clean(row.get("chain_right", ""))
        native_left = clean(row.get("native_receptor_chains", ""))
        native_right = clean(row.get("native_ligand_chains", ""))
        pattern = (len(source_left), len(native_left), len(source_right), len(native_right))
        if len(source_left) + len(source_right) != len(native_left) + len(native_right):
            continue
        if len(source_left) == len(native_left) and len(source_right) == len(native_right):
            continue
        selected.append(row)
        patterns["%d,%d,%d,%d" % pattern] += 1
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else []
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(selected)
    manifest = {
        "schema_version": "prism-gtalign-cross-partition-repair/v1",
        "source_csv": str(SELECTED),
        "source_csv_sha256": sha256(SELECTED),
        "repair_csv": str(OUTPUT),
        "repair_csv_sha256": sha256(OUTPUT),
        "repair_worker": str(WRAPPER),
        "repair_worker_sha256": sha256(WRAPPER),
        "selected_count": len(selected),
        "patterns": dict(patterns),
        "selection_rule": "equal total source/native chains with unequal source-side/native-side partitions",
        "scope": "repair only; no source/canonical input duplication",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
