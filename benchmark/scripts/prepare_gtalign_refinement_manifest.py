#!/usr/bin/env python3
"""Prepare a compact GTalign refinement manifest without copying structures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare(source: Path, output: Path, rejected_output: Path) -> dict[str, Any]:
    selected: list[dict[str, str]] = []
    rejected: list[dict[str, str]] = []
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for index, row in enumerate(rows):
        reason = ""
        raw_path = Path(row.get("raw_dockq_json", ""))
        payload: dict[str, Any] = {}
        if not raw_path.is_file():
            reason = "raw_dockq_json_missing"
        else:
            try:
                payload = json.loads(raw_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                reason = f"raw_dockq_json_invalid:{type(exc).__name__}"
        model = Path(str(payload.get("model", ""))) if payload else Path()
        native = Path(str(payload.get("native", ""))) if payload else Path()
        if not reason and (not model.is_file() or model.stat().st_size == 0):
            reason = "transformed_model_missing"
        if not reason and (not native.is_file() or native.stat().st_size == 0):
            reason = "native_pdb_missing"
        if reason:
            rejected.append({"source_row": str(index), "reason": reason, **row})
            continue
        selected.append(
            {
                "pipeline": "gtalign",
                "dataset": "bm55_full",
                "split": row.get("split", ""),
                "case_id": row.get("case_id", ""),
                "template": row.get("template", ""),
                "orientation": row.get("orientation", ""),
                "query_left": row.get("query_left", ""),
                "query_right": row.get("query_right", ""),
                "chain_left": row.get("chain_left", ""),
                "chain_right": row.get("chain_right", ""),
                "model": str(model.resolve()),
                "native_pdb": str(native.resolve()),
                "native_receptor_chains": row.get("native_receptor_chains", ""),
                "native_ligand_chains": row.get("native_ligand_chains", ""),
                "source_score_status": row.get("source_status", row.get("score_status", "")),
                "source_score_sha256": row.get("raw_dockq_json_sha256", ""),
                "source_row": str(index),
            }
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    rejected_output.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in selected for key in row})
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)
    rejected_fields = sorted({key for row in rejected for key in row}) or ["reason"]
    with rejected_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rejected_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rejected)
    manifest = {
        "schema_version": "prism-gtalign-refinement-manifest/v1",
        "source": str(source.resolve()),
        "source_sha256": sha256_file(source),
        "selected_count": len(selected),
        "rejected_count": len(rejected),
        "selected_csv": str(output.resolve()),
        "selected_csv_sha256": sha256_file(output),
        "rejected_csv": str(rejected_output.resolve()),
        "rejected_csv_sha256": sha256_file(rejected_output),
        "source_status_counts": {
            value: sum(row.get("source_score_status") == value for row in selected)
            for value in sorted({row.get("source_score_status", "") for row in selected})
        },
        "retention": "manifest only; model/native inputs remain at source paths until refinement and scoring validate",
    }
    return manifest


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
