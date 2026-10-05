#!/usr/bin/env python3
"""Prepare a compact USalign common-refinement manifest from scored candidates.

The manifest references the validated common-transformation halves in place;
it does not copy structures.  Generated candidates with a score status that
can support a downstream refinement consumer are selected.  All other rows
are retained in a rejection ledger with an explicit reason.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any


REFINABLE_SCORE_STATUSES = {"scored", "scored_cross_only", "valid_unscored"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(source: Path) -> list[dict[str, str]]:
    with source.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def prepare(source: Path, output: Path, rejected_output: Path) -> dict[str, Any]:
    selected: list[dict[str, str]] = []
    rejected: list[dict[str, str]] = []
    rows = read_rows(source)
    for index, row in enumerate(rows):
        reason = ""
        if row.get("status") != "generated":
            reason = f"candidate_status:{row.get('status', '') or 'missing'}"
        elif row.get("score_status") not in REFINABLE_SCORE_STATUSES:
            reason = f"score_status:{row.get('score_status', '') or 'missing'}"

        left = Path(row.get("transformed_left", ""))
        right = Path(row.get("transformed_right", ""))
        native = Path(row.get("native_pdb", ""))
        if not reason and (not left.is_file() or left.stat().st_size == 0):
            reason = "transformed_left_missing"
        if not reason and (not right.is_file() or right.stat().st_size == 0):
            reason = "transformed_right_missing"
        if not reason and (not native.is_file() or native.stat().st_size == 0):
            reason = "native_pdb_missing"
        required = ("native_receptor_chains", "native_ligand_chains")
        if not reason and any(not row.get(key, "").strip() for key in required):
            reason = "native_chain_mapping_missing"

        if reason:
            rejected.append({"source_row": str(index), "reason": reason, **row})
            continue

        selected.append(
            {
                "pipeline": "usalign",
                "dataset": row.get("dataset", "bm55_full"),
                "split": row.get("split", ""),
                "case_id": row.get("case_id", row.get("pair_id", "")),
                "template": row.get("template", ""),
                "orientation": row.get("orientation", ""),
                "query_left": row.get("query_left", ""),
                "query_right": row.get("query_right", ""),
                "left": str(left.resolve()),
                "right": str(right.resolve()),
                "native_pdb": str(native.resolve()),
                "native_receptor_chains": row.get("native_receptor_chains", ""),
                "native_ligand_chains": row.get("native_ligand_chains", ""),
                "source_score_status": row.get("score_status", ""),
                "source_score_sha256": row.get("raw_dockq_json_sha256", ""),
                "source_candidate_index": row.get("candidate_index", ""),
                "source_row": str(index),
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    rejected_output.parent.mkdir(parents=True, exist_ok=True)
    selected_fields = sorted({key for row in selected for key in row}) or ["pipeline"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=selected_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)
    rejected_fields = sorted({key for row in rejected for key in row}) or ["reason"]
    with rejected_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rejected_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rejected)

    return {
        "schema_version": "prism-usalign-refinement-manifest/v1",
        "source": str(source.resolve()),
        "source_sha256": sha256_file(source),
        "source_rows": len(rows),
        "selected_count": len(selected),
        "rejected_count": len(rejected),
        "selected_csv": str(output.resolve()),
        "selected_csv_sha256": sha256_file(output),
        "rejected_csv": str(rejected_output.resolve()),
        "rejected_csv_sha256": sha256_file(rejected_output),
        "score_status_counts": dict(sorted(Counter(row.get("score_status", "") for row in rows).items())),
        "rejection_reason_counts": dict(sorted(Counter(row["reason"] for row in rejected).items())),
        "retention": "manifest only; transformed halves remain at source paths until common refinement and paired DockQ validate",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rejected", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.source.resolve(), args.output.resolve(), args.rejected.resolve())
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
