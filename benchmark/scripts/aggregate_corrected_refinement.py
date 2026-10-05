#!/usr/bin/env python3
"""Aggregate common-refinement checkpoints with the corrected DockQ contract.

The active large-refinement worker predates the current compact scorer and
stores a scalar ``dockq`` alongside the raw DockQ JSON.  That scalar may be an
interface sum/best-interface value.  This aggregator reads ``GlobalDockQ`` and
the requested receptor--ligand interfaces from the raw JSON, preserves the
legacy scalar only as a diagnostic, and never converts missing scores to zero.
It is read-only with respect to the refinement tree.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


STAGES = ("dockq_fiberdock", "dockq_rosetta")
STAGE_LABELS = {
    "dockq_fiberdock": "fiberdock",
    "dockq_rosetta": "external_rosetta",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def requested_interfaces(stage: dict[str, Any], best_result: dict[str, Any]) -> list[str]:
    receptor = str(stage.get("native_receptor_chains", ""))
    ligand = str(stage.get("native_ligand_chains", ""))
    keys = [f"{r}{l}" for r in receptor for l in ligand]
    return [key for key in keys if key in best_result]


def score_stage(stage: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Return one compact score record and zero or more interface records."""

    label = STAGE_LABELS.get(str(stage.get("stage", "")), str(stage.get("stage", "")))
    prefix = f"{label}_"
    base: dict[str, Any] = {
        f"{prefix}status": str(stage.get("status", "")),
        f"{prefix}reason": str(stage.get("reason", stage.get("error", "")) or ""),
        f"{prefix}elapsed_seconds": stage.get("elapsed_seconds", ""),
        f"{prefix}raw_dockq_json": str(stage.get("raw_dockq_json", "") or ""),
        f"{prefix}raw_dockq_json_sha256": str(stage.get("raw_dockq_json_sha256", "") or ""),
        f"{prefix}global_dockq": "",
        f"{prefix}global_dockq_status": "not_scoreable",
        f"{prefix}cross_best": "",
        f"{prefix}cross_mean": "",
        f"{prefix}cross_fnat_mean": "",
        f"{prefix}cross_irmsd_mean": "",
        f"{prefix}cross_lrmsd_mean": "",
        f"{prefix}cross_interface_count": 0,
        f"{prefix}legacy_dockq_diagnostic": stage.get("dockq", ""),
        f"{prefix}legacy_dockq_sum_diagnostic": stage.get("dockq_sum", ""),
        f"{prefix}dockq_best_internal_diagnostic": "",
        f"{prefix}mapping": str(stage.get("mapping", "") or ""),
        f"{prefix}evaluator_version": str(stage.get("evaluator_version", "") or ""),
    }
    interfaces: list[dict[str, Any]] = []
    raw_value = stage.get("raw_dockq_json")
    if not raw_value:
        return base, interfaces
    raw = Path(str(raw_value))
    if not raw.is_file():
        base[f"{prefix}status"] = "valid_unscored"
        base[f"{prefix}reason"] = "raw_dockq_json_missing"
        base[f"{prefix}global_dockq_status"] = "missing_raw_json"
        return base, interfaces
    try:
        payload = json.loads(raw.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        base[f"{prefix}status"] = "score_failed"
        base[f"{prefix}reason"] = f"invalid_raw_dockq_json:{exc}"
        return base, interfaces

    raw_hash = sha256_file(raw)
    base[f"{prefix}raw_dockq_json_sha256"] = raw_hash
    global_dockq = number(payload.get("GlobalDockQ"))
    best_internal = number(payload.get("best_dockq"))
    base[f"{prefix}global_dockq"] = "" if global_dockq is None else global_dockq
    base[f"{prefix}dockq_best_internal_diagnostic"] = "" if best_internal is None else best_internal
    base[f"{prefix}global_dockq_status"] = "scored" if global_dockq is not None else "valid_unscored"
    base[f"{prefix}status"] = "scored" if global_dockq is not None else "valid_unscored"
    base[f"{prefix}mapping"] = str(payload.get("best_mapping_str") or stage.get("mapping") or "")

    best_result = payload.get("best_result") or {}
    if not isinstance(best_result, dict):
        best_result = {}
    cross_keys = requested_interfaces(stage, best_result)
    cross_rows = []
    for interface in cross_keys:
        component = best_result.get(interface)
        if not isinstance(component, dict):
            continue
        row = {
            "stage": label,
            "interface": interface,
            "requested_cross_interface": True,
            "GlobalDockQ": "" if global_dockq is None else global_dockq,
            "DockQ": component.get("DockQ", ""),
            "Fnat": component.get("fnat", component.get("Fnat", "")),
            "iRMSD": component.get("iRMSD", ""),
            "LRMSD": component.get("LRMSD", ""),
            "mapping": base[f"{prefix}mapping"],
            "raw_dockq_json": str(raw),
            "raw_dockq_json_sha256": raw_hash,
        }
        cross_rows.append(row)
        interfaces.append(row)

    dockq_values = [number(row.get("DockQ")) for row in cross_rows]
    fnat_values = [number(row.get("Fnat")) for row in cross_rows]
    irmsd_values = [number(row.get("iRMSD")) for row in cross_rows]
    lrmsd_values = [number(row.get("LRMSD")) for row in cross_rows]
    dockq_values = [value for value in dockq_values if value is not None]
    fnat_values = [value for value in fnat_values if value is not None]
    irmsd_values = [value for value in irmsd_values if value is not None]
    lrmsd_values = [value for value in lrmsd_values if value is not None]
    base[f"{prefix}cross_best"] = max(dockq_values) if dockq_values else ""
    base[f"{prefix}cross_mean"] = mean(dockq_values) if dockq_values else ""
    base[f"{prefix}cross_fnat_mean"] = mean(fnat_values) if fnat_values else ""
    base[f"{prefix}cross_irmsd_mean"] = mean(irmsd_values) if irmsd_values else ""
    base[f"{prefix}cross_lrmsd_mean"] = mean(lrmsd_values) if lrmsd_values else ""
    base[f"{prefix}cross_interface_count"] = len(cross_rows)
    return base, interfaces


def flatten(checkpoint: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    candidate = checkpoint.get("candidate") or {}
    stage_records = checkpoint.get("stages") or {}
    failed_stages = [str(name) for name, stage in stage_records.items() if isinstance(stage, dict) and stage.get("status") == "failed"]
    failed_stages.extend(str(name) for name in (checkpoint.get("failed_stages") or []) if str(name) not in failed_stages)
    failure_reasons = {
        str(name): str(stage.get("error", stage.get("reason", "")) or "")
        for name, stage in stage_records.items()
        if isinstance(stage, dict) and stage.get("status") == "failed"
    }
    row: dict[str, Any] = {
        "key": checkpoint.get("key", ""),
        "manifest_index": candidate.get("manifest_index", checkpoint.get("index", "")),
        "pipeline": candidate.get("pipeline", ""),
        "dataset": candidate.get("dataset", ""),
        "split": candidate.get("split", ""),
        "case_id": candidate.get("case_id", ""),
        "template": candidate.get("template", ""),
        "orientation": candidate.get("orientation", ""),
        "selection_rank": candidate.get("selection_rank", ""),
        "candidate_status": checkpoint.get("status", ""),
        "failed_stages": ";".join(failed_stages),
        "failure_reasons": json.dumps(failure_reasons, sort_keys=True),
        "slurm_job_id": checkpoint.get("slurm_job_id", ""),
        "started_at": checkpoint.get("started_at", ""),
        "finished_at": checkpoint.get("finished_at", ""),
    }
    interfaces: list[dict[str, Any]] = []
    for stage_name in STAGES:
        stage = dict(stage_records.get(stage_name) or {})
        stage["stage"] = stage_name
        compact, stage_interfaces = score_stage(stage)
        row.update(compact)
        for interface in stage_interfaces:
            interfaces.append({
                "key": row["key"],
                "manifest_index": row["manifest_index"],
                "pipeline": row["pipeline"],
                "case_id": row["case_id"],
                "template": row["template"],
                **interface,
            })

    transformed = number(row.get("fiberdock_global_dockq"))
    refined = number(row.get("external_rosetta_global_dockq"))
    row["paired_global_dockq_status"] = "paired" if transformed is not None and refined is not None else "unpaired"
    row["delta_external_rosetta_minus_fiberdock_global_dockq"] = "" if transformed is None or refined is None else refined - transformed
    return row, interfaces


def write_rows(path: Path, rows: list[dict[str, Any]], delimiter: str = "\t") -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row}) or ["status"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return sha256_file(path)


def aggregate(comparison_root: Path, output_dir: Path, expected_count: int | None) -> dict[str, Any]:
    checkpoints = sorted((comparison_root / "checkpoints").glob("*.json"))
    rows: list[dict[str, Any]] = []
    interfaces: list[dict[str, Any]] = []
    seen: set[str] = set()
    duplicate_keys: list[str] = []
    invalid_checkpoints = 0
    for path in checkpoints:
        try:
            checkpoint = json.loads(path.read_text(encoding="utf-8"))
            row, stage_interfaces = flatten(checkpoint)
        except (OSError, json.JSONDecodeError, TypeError, ValueError):
            invalid_checkpoints += 1
            continue
        key = str(row.get("key", ""))
        if not key:
            continue
        if key in seen:
            duplicate_keys.append(key)
            continue
        seen.add(key)
        rows.append(row)
        interfaces.extend(stage_interfaces)
    rows.sort(key=lambda row: (int(row["manifest_index"]) if str(row["manifest_index"]).isdigit() else 10**18, str(row["key"])))
    interfaces.sort(key=lambda row: (int(row["manifest_index"]) if str(row["manifest_index"]).isdigit() else 10**18, str(row["key"]), str(row["stage"]), str(row["interface"])))

    output_dir.mkdir(parents=True, exist_ok=True)
    results_path = output_dir / "refinement_comparison.tsv"
    interface_path = output_dir / "refinement_interfaces.tsv"
    results_hash = write_rows(results_path, rows)
    interface_hash = write_rows(interface_path, interfaces)
    status_counts = Counter(str(row.get("candidate_status", "missing")) for row in rows)
    score_counts = {
        label: dict(sorted(Counter(str(row.get(f"{label}_status", "missing")) for row in rows).items()))
        for label in STAGE_LABELS.values()
    }
    complete = expected_count is not None and len(rows) == expected_count and not duplicate_keys and not invalid_checkpoints
    summary = {
        "schema_version": "prism-corrected-refinement-aggregation/v1",
        "status": "complete" if complete else "incomplete",
        "expected_count": expected_count,
        "checkpoint_count_on_disk": len(checkpoints),
        "unique_rows": len(rows),
        "invalid_checkpoints": invalid_checkpoints,
        "duplicate_keys": duplicate_keys,
        "candidate_status_counts": dict(sorted(status_counts.items())),
        "dockq_stage_status_counts": score_counts,
        "interface_rows": len(interfaces),
        "results_path": str(results_path.resolve()),
        "results_sha256": results_hash,
        "interfaces_path": str(interface_path.resolve()),
        "interfaces_sha256": interface_hash,
        "cleanup_eligible": bool(complete),
        "cleanup_note": "No files were deleted by this read-only aggregator.",
    }
    (output_dir / "aggregation_status.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-count", type=int)
    args = parser.parse_args()
    summary = aggregate(args.comparison_root.resolve(), args.output_dir.resolve(), args.expected_count)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["status"] == "complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())
