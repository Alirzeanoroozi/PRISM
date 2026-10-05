#!/usr/bin/env python3
"""Execute one isolated source/staging gate unit for one dataset row."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark.scripts.build_investigation_source_manifest import (
    SOURCE_FIELDS,
    VALIDATION_FIELDS,
    _write_tsv,
    build_manifests,
)
from benchmark.scripts.stage_curated_sources import stage_sources


def run_task(repo_root: str | Path, output_dir: str | Path, dataset_row_id: str, *, strict: bool = False) -> dict[str, object]:
    root = Path(repo_root).resolve()
    output = Path(output_dir).resolve()
    source_rows, validation_rows = build_manifests(root, dataset_row_ids=[dataset_row_id])
    output.mkdir(parents=True, exist_ok=True)
    source_path = output / "source_manifest.tsv"
    validation_path = output / "structure_validation.tsv"
    _write_tsv(source_path, SOURCE_FIELDS, source_rows)
    _write_tsv(validation_path, VALIDATION_FIELDS, validation_rows)
    staged_path, staging_failures = stage_sources(source_path, output / "staged", repo_root=root, strict=False)
    pipeline_source_rows = [row for row in source_rows if row["source_scope"] == "pipeline"]
    pipeline_validation_rows = [row for row in validation_rows if row["source_scope"] == "pipeline"]
    source_failures = [
        row for row in pipeline_source_rows
        if row["resolution_status"] != "resolved" or row["candidate_status"] != "unique" or not row["archive_prefix"]
    ]
    validation_failures = [
        row for row in pipeline_validation_rows
        if row["parse_status"] != "ok" or row["chain_set_status"] not in {"ok", "not_declared"}
    ]
    role_set = {row["source_role"] for row in pipeline_source_rows}
    role_failure = role_set != {"pipeline_receptor", "pipeline_ligand", "native_receptor", "native_ligand"}
    summary = {
        "schema_version": "source-gate-task/v1",
        "dataset_row_id": dataset_row_id,
        "source_manifest": str(source_path),
        "structure_validation": str(validation_path),
        "staged_sources": str(staged_path),
        "pipeline_role_count": len(pipeline_source_rows),
        "source_failures": len(source_failures),
        "validation_failures": len(validation_failures),
        "staging_failures": staging_failures,
        "role_contract_failure": role_failure,
        "status": "success" if not (source_failures or validation_failures or staging_failures or role_failure) else "source_gate_failed",
        "failure_ids": [row["source_role"] for row in source_failures] + [row["source_role"] for row in validation_failures],
    }
    summary_path = output / "source_gate_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if strict and summary["status"] != "success":
        raise ValueError("source gate failed for " + dataset_row_id)
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dataset-row-id", required=True)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)
    try:
        summary = run_task(args.repo_root, args.output_dir, args.dataset_row_id, strict=args.strict)
    except (OSError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
