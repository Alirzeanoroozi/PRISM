#!/usr/bin/env python3
"""Reconcile all large-refinement checkpoints without dropping failures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def corrected_global_dockq(stage: dict[str, object]) -> dict[str, object]:
    """Return the bounded multimer score while preserving the interface sum.

    DockQ's ``best_dockq`` is the sum of the selected native-interface scores
    for a multimer.  It is useful provenance, but it is not a normalized
    per-model score.  The official JSON also contains ``GlobalDockQ``, which
    is the interface-count-normalized score and is therefore the primary
    result for comparisons.
    """
    raw_path = stage.get("raw_dockq_json") or stage.get("dockq_json")
    legacy_sum = stage.get("dockq_sum", stage.get("dockq", ""))
    result: dict[str, object] = {
        "dockq": "",
        "dockq_legacy": stage.get("dockq", ""),
        "dockq_sum": legacy_sum,
        "dockq_contract": "missing_raw_json",
        "dockq_json": str(raw_path or ""),
        "dockq_json_sha256": stage.get("raw_dockq_json_sha256", ""),
        "interface_count": "",
    }
    if not raw_path:
        return result

    path = Path(str(raw_path))
    if not path.is_file():
        return result
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        result["dockq_contract"] = "invalid_raw_json"
        return result
    if not isinstance(payload, dict):
        result["dockq_contract"] = "invalid_raw_json"
        return result

    if "best_dockq" in payload:
        result["dockq_sum"] = payload.get("best_dockq", "")
    best_result = payload.get("best_result")
    if isinstance(best_result, (dict, list)):
        result["interface_count"] = len(best_result)
    if "GlobalDockQ" not in payload:
        result["dockq_contract"] = "valid_unscored"
    else:
        try:
            value = float(payload["GlobalDockQ"])
        except (TypeError, ValueError):
            result["dockq_contract"] = "invalid_global_dockq"
        else:
            if 0.0 <= value <= 1.0:
                result["dockq"] = value
                result["dockq_contract"] = "GlobalDockQ"
            else:
                result["dockq_contract"] = "invalid_global_dockq"
    result["dockq_json_sha256"] = sha256(path)
    return result


def flatten(checkpoint: dict[str, object]) -> dict[str, object]:
    row: dict[str, object] = {
        "key": checkpoint.get("key", ""),
        "manifest_index": checkpoint.get("candidate", {}).get("manifest_index", ""),
        "origin_manifest_index": checkpoint.get("candidate", {}).get("origin_manifest_index", ""),
        "index": checkpoint.get("index", ""),
        "pipeline": checkpoint.get("candidate", {}).get("pipeline", ""),
        "case_id": checkpoint.get("candidate", {}).get("case_id", ""),
        "template": checkpoint.get("candidate", {}).get("template", ""),
        "orientation": checkpoint.get("candidate", {}).get("orientation", ""),
        "status": checkpoint.get("status", ""),
        "failed_stages": ";".join(checkpoint.get("failed_stages", []) or []),
        "slurm_job_id": checkpoint.get("slurm_job_id", ""),
        "slurm_array_job_id": checkpoint.get("slurm_array_job_id", ""),
        "slurm_array_task_id": checkpoint.get("slurm_array_task_id", ""),
        "started_at": checkpoint.get("started_at", ""),
        "finished_at": checkpoint.get("finished_at", ""),
    }
    for stage_name, stage in sorted((checkpoint.get("stages") or {}).items()):
        if not isinstance(stage, dict):
            continue
        row[f"{stage_name}_status"] = stage.get("status", "")
        row[f"{stage_name}_elapsed_seconds"] = stage.get("elapsed_seconds", "")
        for key in (
            "energy", "dockq", "dockq_global", "dockq_sum", "fnat", "dockq_irmsd",
            "dockq_lrmsd", "irmsd_backbone", "refined_model", "refined_model_sha256",
            "raw_dockq_json", "raw_dockq_json_sha256", "output_prefix", "error", "error_class",
            "reason", "native_pdb", "native_sha256", "native_available_chains",
            "native_missing_chains", "native_mapping_status", "scoring_scope", "mapping",
            "model_receptor_chains", "model_ligand_chains", "model_mapping_mode",
            "model_chain_blocks",
        ):
            if key in stage:
                row[f"{stage_name}_{key}"] = stage[key]
        if stage_name in {"dockq_fiberdock", "dockq_rosetta"}:
            corrected = corrected_global_dockq(stage)
            for key, value in corrected.items():
                row[f"{stage_name}_{key}"] = value
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-root", required=True, type=Path)
    parser.add_argument("--shard-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    shard_manifest = json.loads(args.shard_manifest.read_text(encoding="utf-8"))
    expected = int(shard_manifest["selected_count"])
    checkpoints = sorted((args.comparison_root / "checkpoints").glob("*.json"))
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    duplicate_keys: list[str] = []
    for path in checkpoints:
        try:
            checkpoint = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        row = flatten(checkpoint)
        key = str(row["key"])
        if not key:
            continue
        if key in seen:
            duplicate_keys.append(key)
            continue
        seen.add(key)
        rows.append(row)
    rows.sort(key=lambda row: (int(row["manifest_index"]) if str(row["manifest_index"]).isdigit() else 10**18, str(row["key"])))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = sorted(set().union(*(set(row) for row in rows))) if rows else ["key", "status"]
    table = args.output_dir / "refinement_results.csv"
    with table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    statuses = Counter(str(row.get("status", "missing")) for row in rows)

    def values(field: str) -> list[float]:
        output: list[float] = []
        for row in rows:
            try:
                value = float(row[field])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(value):
                output.append(value)
        return output

    def metric(field: str) -> dict[str, object]:
        observed = values(field)
        out_of_range = [value for value in observed if value < 0.0 or value > 1.0]
        return {
            "count": len(observed),
            "mean": statistics.mean(observed) if observed else None,
            "median": statistics.median(observed) if observed else None,
            "min": min(observed) if observed else None,
            "max": max(observed) if observed else None,
            "out_of_range_count": len(out_of_range),
            "out_of_range_mean": statistics.mean(out_of_range) if out_of_range else None,
            "out_of_range_median": statistics.median(out_of_range) if out_of_range else None,
        }

    summary = {
        "status": "complete" if len(rows) == expected and not duplicate_keys else "incomplete",
        "expected_rows": expected,
        "checkpoint_rows": len(rows),
        "missing_rows": expected - len(rows),
        "duplicate_keys": duplicate_keys,
        "status_counts": dict(sorted(statuses.items())),
        "dockq_contract_counts": {
            stage: dict(
                sorted(
                    Counter(str(row.get(f"{stage}_dockq_contract", "missing")) for row in rows).items()
                )
            )
            for stage in ("dockq_fiberdock", "dockq_rosetta")
        },
        "dockq_metrics": {
            stage: {
                "global_dockq": metric(f"{stage}_dockq"),
                "interface_sum": metric(f"{stage}_dockq_sum"),
            }
            for stage in ("dockq_fiberdock", "dockq_rosetta")
        },
        "checkpoint_count_on_disk": len(checkpoints),
        "results_csv": str(table.resolve()),
        "results_csv_sha256": sha256(table),
    }
    (args.output_dir / "aggregation_status.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0 if summary["status"] == "complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())
