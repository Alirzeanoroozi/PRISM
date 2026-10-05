#!/usr/bin/env python3
"""Aggregate the 40-candidate FiberDock/Rosetta/DockQ comparison."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
from collections import Counter
from pathlib import Path


def flatten_stage(row: dict[str, object], stage: dict[str, object], prefix: str) -> None:
    row[f"{prefix}_status"] = stage.get("status", "missing")
    row[f"{prefix}_elapsed_seconds"] = stage.get("elapsed_seconds", "")
    row[f"{prefix}_error"] = stage.get("error", "")


def corrected_global_dockq(stage: dict[str, object], label: str) -> tuple[object, str, object, object, object]:
    """Read bounded GlobalDockQ from preserved raw JSON, not best_dockq sum."""
    raw_path = stage.get("dockq_json")
    if not raw_path and stage.get("model_pdb"):
        raw_dir = Path(str(stage["model_pdb"])).resolve().parent / label
        candidates = sorted(raw_dir.glob("*.json"), key=lambda path: path.stat().st_mtime)
        raw_path = str(candidates[-1]) if candidates else ""
    if not raw_path:
        return "", "missing_raw_json", "", "", ""
    path = Path(str(raw_path))
    if not path.is_file():
        return "", "missing_raw_json", "", str(path), ""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        value = float(payload["GlobalDockQ"])
    except (KeyError, TypeError, ValueError, OSError, json.JSONDecodeError):
        return "", "valid_unscored", payload.get("best_dockq", "") if "payload" in locals() else "", str(path), ""
    if not 0.0 <= value <= 1.0:
        return "", "invalid_global_dockq", payload.get("best_dockq", ""), str(path), ""
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return value, "GlobalDockQ", payload.get("best_dockq", ""), str(path), digest


def capri_class(value: object) -> str:
    if value == "" or value is None:
        return "Unknown"
    value = float(value)
    if value < 0.23:
        return "Incorrect"
    if value < 0.49:
        return "Acceptable"
    if value < 0.80:
        return "Medium"
    return "High"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-root", required=True, type=Path)
    parser.add_argument("--selected-csv", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--job-id", default="")
    args = parser.parse_args()

    root = args.comparison_root.resolve()
    selected = list(csv.DictReader(args.selected_csv.open(newline="", encoding="utf-8")))
    checkpoints = sorted((root / "checkpoints").glob("v2-chain-normalized-*.json"))
    records = [json.loads(path.read_text(encoding="utf-8")) for path in checkpoints]
    records.sort(key=lambda item: int(item["index"]))
    indices = [int(item["index"]) for item in records]
    expected = list(range(len(selected)))
    if indices != expected:
        raise SystemExit(f"no-drop validation failed: checkpoint indices={indices}, expected={expected}")
    if len(selected) != 40 or len(records) != 40:
        raise SystemExit(f"expected exactly 40 selected/checkpoint rows, got {len(selected)}/{len(records)}")

    rows: list[dict[str, object]] = []
    for record in records:
        candidate = record["candidate"]
        stages = record["stages"]
        row: dict[str, object] = {
            "index": record["index"],
            "key": record["key"],
            "pipeline": candidate.get("pipeline", ""),
            "selection_rank": candidate.get("selection_rank", ""),
            "case_id": candidate.get("case_id", ""),
            "split": candidate.get("split", ""),
            "template": candidate.get("template", ""),
            "orientation": candidate.get("orientation", ""),
            "query_left": candidate.get("query_left", ""),
            "query_right": candidate.get("query_right", ""),
            "confidence_tm_min": candidate.get("confidence_tm_min", ""),
            "confidence_tm_mean": candidate.get("confidence_tm_mean", ""),
            "match_count_total": candidate.get("match_count_total", ""),
            "benchmark_irmsd_A": candidate.get("benchmark_irmsd_A", ""),
            "native_receptor_chains": candidate.get("native_receptor_chains", ""),
            "native_ligand_chains": candidate.get("native_ligand_chains", ""),
            "native_pdb": candidate.get("native_pdb", ""),
            "candidate_status": record.get("status", ""),
            "candidate_started_at": record.get("started_at", ""),
            "candidate_finished_at": record.get("finished_at", ""),
        }
        for prefix, stage_name in (
            ("input_normalization", "input_normalization"),
            ("fiberdock", "fiberdock"),
            ("external_rosetta", "external_rosetta"),
            ("dockq_fiberdock", "dockq_fiberdock"),
            ("dockq_rosetta", "dockq_rosetta"),
        ):
            flatten_stage(row, stages.get(stage_name, {}), prefix)
        fiber = stages.get("fiberdock", {})
        rosetta = stages.get("external_rosetta", {})
        fiber_dockq = stages.get("dockq_fiberdock", {})
        rosetta_dockq = stages.get("dockq_rosetta", {})
        fiber_global, fiber_contract, fiber_sum, fiber_raw, fiber_raw_sha = corrected_global_dockq(fiber_dockq, "fiberdock")
        rosetta_global, rosetta_contract, rosetta_sum, rosetta_raw, rosetta_raw_sha = corrected_global_dockq(rosetta_dockq, "rosetta")
        row.update(
            fiberdock_energy=fiber.get("energy", ""),
            fiberdock_model=fiber.get("refined_model", ""),
            fiberdock_model_sha256=fiber.get("refined_model_sha256", ""),
            rosetta_total_score=rosetta.get("totalscore", ""),
            rosetta_interaction_score=rosetta.get("interaction_score", ""),
            rosetta_model=rosetta.get("refined_model", ""),
            rosetta_model_sha256=rosetta.get("refined_model_sha256", ""),
            fiberdock_dockq=fiber_global,
            fiberdock_dockq_legacy=fiber_dockq.get("dockq", ""),
            fiberdock_dockq_sum=fiber_sum,
            fiberdock_dockq_contract=fiber_contract,
            fiberdock_dockq_json=fiber_raw,
            fiberdock_dockq_json_sha256=fiber_raw_sha,
            fiberdock_irmsd_backbone=fiber_dockq.get("irmsd_backbone", ""),
            fiberdock_dockq_capri=capri_class(fiber_global),
            fiberdock_dockq_mapping=fiber_dockq.get("mapping", ""),
            rosetta_dockq=rosetta_global,
            rosetta_dockq_legacy=rosetta_dockq.get("dockq", ""),
            rosetta_dockq_sum=rosetta_sum,
            rosetta_dockq_contract=rosetta_contract,
            rosetta_dockq_json=rosetta_raw,
            rosetta_dockq_json_sha256=rosetta_raw_sha,
            rosetta_irmsd_backbone=rosetta_dockq.get("irmsd_backbone", ""),
            rosetta_dockq_capri=capri_class(rosetta_global),
            rosetta_dockq_mapping=rosetta_dockq.get("mapping", ""),
        )
        row["candidate_total_elapsed_seconds"] = sum(
            float(stages.get(stage_name, {}).get("elapsed_seconds", 0.0) or 0.0)
            for stage_name in ("input_normalization", "fiberdock", "external_rosetta", "dockq_fiberdock", "dockq_rosetta")
        )
        rows.append(row)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    columns = sorted({key for row in rows for key in row})
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    def values(field: str) -> list[float]:
        output = []
        for row in rows:
            try:
                output.append(float(row[field]))
            except (TypeError, ValueError, KeyError):
                pass
        return output

    summary = {
        "status": "completed",
        "selected_rows": len(selected),
        "checkpoint_rows": len(records),
        "no_drop": indices == expected,
        "by_pipeline": dict(Counter(str(row["pipeline"]) for row in rows)),
        "candidate_status": dict(Counter(str(row["candidate_status"]) for row in rows)),
        "stage_status": {
            field: dict(Counter(str(row[f"{field}_status"]) for row in rows))
            for field in ("input_normalization", "fiberdock", "external_rosetta", "dockq_fiberdock", "dockq_rosetta")
        },
        "dockq_contract_status": {
            "fiberdock": dict(Counter(str(row["fiberdock_dockq_contract"]) for row in rows)),
            "rosetta": dict(Counter(str(row["rosetta_dockq_contract"]) for row in rows)),
        },
        "metrics": {
            field: {
                "count": len(values(field)),
                "min": min(values(field)) if values(field) else None,
                "median": statistics.median(values(field)) if values(field) else None,
                "max": max(values(field)) if values(field) else None,
                "mean": statistics.mean(values(field)) if values(field) else None,
            }
            for field in (
                "input_normalization_elapsed_seconds", "fiberdock_elapsed_seconds",
                "external_rosetta_elapsed_seconds", "dockq_fiberdock_elapsed_seconds",
                "dockq_rosetta_elapsed_seconds", "candidate_total_elapsed_seconds",
                "fiberdock_dockq", "rosetta_dockq", "fiberdock_irmsd_backbone", "rosetta_irmsd_backbone",
            )
        },
        "output_csv": str(args.output_csv.resolve()),
        "checkpoint_root": str((root / "checkpoints").resolve()),
        "run_job_id": args.job_id,
        "scientific_scope": "bounded high-TM-score probe; not a complete BM55 benchmark or promotion result",
    }
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
