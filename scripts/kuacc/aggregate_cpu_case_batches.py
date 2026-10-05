#!/usr/bin/env python3
"""Aggregate completed per-case PRISM CPU attempts without dropping artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
try:
    from valar_agent.prism_cpu_batch import ArtifactCollision, merge_artifact, sha256_file
except ModuleNotFoundError:
    sys.path.insert(0, str(SCRIPT_DIR))
    from prism_cpu_batch import ArtifactCollision, merge_artifact, sha256_file


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def link_once(destination: Path, target: Path) -> None:
    if destination.is_symlink():
        if destination.resolve() != target.resolve():
            raise ArtifactCollision(f"link collision at {destination}")
        return
    if destination.exists():
        raise ArtifactCollision(f"non-link path blocks required link {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.symlink_to(target.resolve(), target_is_directory=target.is_dir())


def copy_tree(source: Path, destination: Path) -> int:
    copied = 0
    if not source.is_dir():
        return copied
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(source)
        merge_artifact(path, destination / relative)
        copied += 1
    return copied


def canonical_inputs(entries: list[dict]) -> str:
    lines = ["Receptor,Ligand"]
    lines.extend(f"{entry['receptor']},{entry['ligand']}" for entry in entries)
    return "\n".join(lines) + "\n"


def load_task_records(manifest: dict) -> tuple[list[dict], list[dict]]:
    records = []
    incomplete = []
    for entry in manifest["entries"]:
        status_path = Path(entry["case_root"]) / "task_status.json"
        status = json.loads(status_path.read_text(encoding="utf-8")) if status_path.is_file() else {"status": "missing"}
        records.append(status)
        if status.get("status") != "completed":
            incomplete.append({"index": entry["index"], "case_id": entry["case_id"], "status": status.get("status", "missing")})
    return records, incomplete


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-manifest", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()

    manifest_path = Path(args.batch_manifest).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    output_root = Path(args.output_root).resolve()
    status_path = output_root / "aggregation_status.json"
    run_summary_path = output_root / "run_summary.json"
    if run_summary_path.is_file():
        existing = json.loads(run_summary_path.read_text(encoding="utf-8"))
        if existing.get("status") == "completed" and existing.get("batch_manifest_sha256") == sha256_file(manifest_path):
            print(json.dumps({"status": "resumed", "run_root": str(output_root)}, sort_keys=True))
            return 0

    records, incomplete = load_task_records(manifest)
    if incomplete:
        atomic_json(status_path, {"status": "blocked", "reason": "incomplete_case_tasks", "tasks": incomplete, "updated_at": now()})
        print(json.dumps({"status": "blocked", "incomplete": len(incomplete)}, sort_keys=True))
        return 2

    started = time.perf_counter()
    atomic_json(status_path, {"status": "running", "started_at": now(), "job_id": os.environ.get("SLURM_JOB_ID", "")})
    output_root.mkdir(parents=True, exist_ok=True)
    entries = manifest["entries"]
    try:
        dataset_dir = Path(manifest["dataset_dir"])
        pipeline_repo = Path(manifest["pipeline_repo"])
        template_list = Path(manifest["template_list"])
        link_once(output_root / "prism.py", pipeline_repo / "prism.py")
        link_once(output_root / "src", pipeline_repo / "src")
        link_once(output_root / "external_tools", pipeline_repo / "external_tools")
        link_once(output_root / "prism_timed_runner.py", Path(manifest["timed_runner"]))
        link_once(output_root / "usalign_prism_wrapper.py", Path(manifest["usalign_wrapper"]))
        atomic_text(output_root / "inputs.csv", canonical_inputs(entries))
        merge_artifact(dataset_dir / "dataset_manifest.json", output_root / "dataset_manifest.json")
        merge_artifact(template_list, output_root / "templates" / "calculated_templates.txt")
        merge_artifact(template_list, output_root / "templates" / "checked_templates.txt")
        for entry in entries:
            merge_artifact(Path(entry["native_source"]), output_root / "processed" / "pdbs" / entry["native_name"])

        transformation_count = 0
        alignment_count = 0
        evidence_count = 0
        timing_rows = []
        effective_aligner = "tmalign" if manifest["aligner"] == "usalign" else manifest["aligner"]
        for entry, record in zip(entries, records):
            attempt_root = Path(record["attempt_root"])
            transformation = attempt_root / "processed" / "transformation"
            for artifact in sorted(transformation.glob("*")):
                if artifact.is_file():
                    merge_artifact(artifact, output_root / "processed" / "transformation" / artifact.name)
                    transformation_count += artifact.name.endswith("_L.pdb")
            source_alignment = attempt_root / "processed" / f"alignment_{effective_aligner}"
            alignment_destination = output_root / "processed" / f"alignment_{effective_aligner}" / "batches" / f"{entry['index']:04d}"
            alignment_count += copy_tree(source_alignment, alignment_destination)
            evidence_destination = output_root / "evidence" / "batches" / f"{entry['index']:04d}"
            for name in ("run_summary.json", "command_timing.json", "surface_timing.json", "stage_status.jsonl", "pipeline.console.log"):
                source = attempt_root / name
                if source.is_file():
                    merge_artifact(source, evidence_destination / name)
                    evidence_count += 1
            controller_log = Path(entry["case_root"]) / "controller_logs" / f"attempt_{record.get('attempt', 0):04d}.console.log"
            if controller_log.is_file():
                merge_artifact(controller_log, evidence_destination / controller_log.name)
                evidence_count += 1
            summary = json.loads((attempt_root / "run_summary.json").read_text(encoding="utf-8"))
            timing_rows.append({
                "index": entry["index"],
                "case_id": entry["case_id"],
                "attempt_root": str(attempt_root),
                "wall_seconds": summary.get("wall_seconds", 0.0),
                "alignment_records": summary.get("alignment_records", 0),
                "successful_alignment_records": summary.get("successful_alignment_records", 0),
                "failed_alignment_records": summary.get("failed_alignment_records", 0),
                "transformed_model_pairs": summary.get("transformed_model_pairs", 0),
            })

        aggregate_alignment_records = sum(row["alignment_records"] for row in timing_rows)
        aggregate_success = sum(row["successful_alignment_records"] for row in timing_rows)
        aggregate_failed = sum(row["failed_alignment_records"] for row in timing_rows)
        expected_records = manifest["expected_alignment_records"]
        if aggregate_alignment_records != expected_records:
            raise RuntimeError(
                f"alignment record total {aggregate_alignment_records} does not match expected {expected_records}"
            )
        timing = {
            "stage": "full_bm55_case_aggregation",
            "case_rows": timing_rows,
            "case_wall_seconds_sum": sum(row["wall_seconds"] for row in timing_rows),
            "case_wall_seconds_max": max((row["wall_seconds"] for row in timing_rows), default=0.0),
            "aggregation_seconds": time.perf_counter() - started,
        }
        atomic_json(output_root / "timing.json", timing)
        summary = {
            "status": "completed",
            "dataset": manifest["dataset"],
            "scenario": "bm55_full_all_templates",
            "aligner": manifest["aligner"],
            "rank": False,
            "rank_method": "baseline",
            "prodigy": False,
            "refinement": False,
            "template_list": manifest["template_list"],
            "template_sha256": manifest["template_sha256"],
            "template_count": manifest["template_count"],
            "pair_count": manifest["pair_count"],
            "expected_alignment_records": expected_records,
            "alignment_records": aggregate_alignment_records,
            "successful_alignment_records": aggregate_success,
            "failed_alignment_records": aggregate_failed,
            "transformed_model_pairs": transformation_count,
            "alignment_artifact_count": alignment_count,
            "evidence_artifact_count": evidence_count,
            "case_task_count": len(entries),
            "case_task_status": "all_completed",
            "batch_manifest": str(manifest_path),
            "batch_manifest_sha256": sha256_file(manifest_path),
            "aggregation_seconds": time.perf_counter() - started,
            "job_id": os.environ.get("SLURM_JOB_ID", ""),
            "completed_at": now(),
        }
        atomic_json(run_summary_path, summary)
        atomic_json(status_path, {"status": "completed", "run_summary": str(run_summary_path), "updated_at": now()})
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
    except (ArtifactCollision, OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
        atomic_json(status_path, {"status": "failed", "error": f"{type(exc).__name__}: {exc}", "updated_at": now()})
        print(f"aggregation failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
