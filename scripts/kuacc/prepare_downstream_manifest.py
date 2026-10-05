#!/usr/bin/env python3
"""Prepare auditable per-candidate downstream batches from compact PRISM results.

This stage reads the validated candidate table and per-case manifests, checks
the actual transformed PDB inputs and native-chain contracts, and writes small
JSONL batches. It never copies or modifies source structures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def chain_order(path: Path) -> list[str]:
    seen: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                chain = line[21]
                if chain not in seen:
                    seen.append(chain)
    return seen


def manifest_for(submission: dict[str, Any], run_root: Path, pipeline: str) -> Path:
    value = (submission.get("manifests") or {}).get(pipeline)
    candidate = Path(value) if value else run_root / pipeline / "batch_manifest_corrected.json"
    if not candidate.is_file():
        raise FileNotFoundError(f"missing corrected manifest for {pipeline}: {candidate}")
    return candidate


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected object: {path}")
    return value


def candidate_key(row: dict[str, str]) -> str:
    raw = "|".join(row.get(field, "") for field in (
        "pipeline", "case_index", "case_id", "template", "orientation", "query_left", "query_right"
    ))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:24]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submission-manifest", required=True, type=Path)
    parser.add_argument("--candidate-csv", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--batch-size", type=int, default=500)
    args = parser.parse_args()
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be positive")

    submission_path = args.submission_manifest.resolve()
    candidate_path = args.candidate_csv.resolve()
    run_root = submission_path.parent.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    batches_root = output_root / "batches"
    batches_root.mkdir(parents=True, exist_ok=True)
    submission = load_json(submission_path)

    pipeline_entries: dict[str, dict[int, dict[str, Any]]] = {}
    pipeline_manifest_paths: dict[str, Path] = {}
    for pipeline in ("tmalign", "multiprot", "usalign"):
        manifest_path = manifest_for(submission, run_root, pipeline)
        pipeline_manifest_paths[pipeline] = manifest_path
        manifest = load_json(manifest_path)
        pipeline_entries[pipeline] = {int(entry["index"]): entry for entry in manifest.get("entries", [])}

    case_cache: dict[tuple[str, int], dict[str, Any]] = {}
    chain_audit: dict[tuple[str, int], dict[str, Any]] = {}
    invalid_rows: list[dict[str, Any]] = []
    valid_rows: list[dict[str, Any]] = []

    with candidate_path.open(newline="", encoding="utf-8") as handle:
        for source_row in csv.DictReader(handle):
            pipeline = source_row["pipeline"]
            case_index = int(source_row["case_index"])
            entry = pipeline_entries[pipeline][case_index]
            cache_key = (pipeline, case_index)
            case = case_cache.get(cache_key)
            if case is None:
                task_path = Path(entry["case_root"]) / "task_status.json"
                task = load_json(task_path)
                attempt_root = Path(str(task.get("attempt_root", "")))
                summary_path = attempt_root / "run_summary.json"
                dataset_manifest_path = attempt_root / "dataset_manifest.json"
                dataset_manifest = load_json(dataset_manifest_path)
                pair = next(
                    (
                        item for item in dataset_manifest.get("pairs", [])
                        if str(item.get("receptor", "")).strip() == str(entry.get("receptor", "")).strip()
                        and str(item.get("ligand", "")).strip() == str(entry.get("ligand", "")).strip()
                    ),
                    None,
                )
                if pair is None:
                    raise RuntimeError(f"dataset pair missing for {pipeline} case {case_index}")
                native_path = attempt_root / "processed" / "pdbs" / entry["native_name"]
                query_left_path = attempt_root / "processed" / "pdbs" / entry["query_source_names"][0]
                query_right_path = attempt_root / "processed" / "pdbs" / entry["query_source_names"][1]
                query_left_chains = chain_order(query_left_path) if query_left_path.is_file() else []
                query_right_chains = chain_order(query_right_path) if query_right_path.is_file() else []
                case = {
                    "entry": entry,
                    "task": task,
                    "attempt_root": attempt_root,
                    "summary_path": summary_path,
                    "native_path": native_path,
                    "query_left_chains": query_left_chains,
                    "query_right_chains": query_right_chains,
                    "transformed_shape_checked": False,
                    "transformed_shape_valid": False,
                    "pair": pair,
                }
                case_cache[cache_key] = case

                native_chains = chain_order(native_path) if native_path.is_file() else []
                expected_receptor = str(pair.get("native_receptor_chains", "")).strip()
                expected_ligand = str(pair.get("native_ligand_chains", "")).strip()
                chain_audit[cache_key] = {
                    "pipeline": pipeline,
                    "case_index": case_index,
                    "case_id": entry.get("case_id", ""),
                    "receptor": entry.get("receptor", ""),
                    "ligand": entry.get("ligand", ""),
                    "native_pdb": str(native_path),
                    "native_chain_ids": "".join(native_chains),
                    "native_chain_count": len(native_chains),
                    "expected_receptor_chains": expected_receptor,
                    "expected_ligand_chains": expected_ligand,
                    "expected_receptor_chain_count": len(expected_receptor),
                    "expected_ligand_chain_count": len(expected_ligand),
                    "query_left_chain_ids": "".join(query_left_chains),
                    "query_right_chain_ids": "".join(query_right_chains),
                    "query_left_chain_count": len(query_left_chains),
                    "query_right_chain_count": len(query_right_chains),
                }

            query_left = source_row["query_left"].strip()
            query_right = source_row["query_right"].strip()
            orientation = source_row["orientation"].strip()
            template = source_row["template"].strip()
            attempt_root = case["attempt_root"]
            left = attempt_root / "processed" / "transformation" / f"{template}_{query_left}_{query_right}_{orientation}_L.pdb"
            right = attempt_root / "processed" / "transformation" / f"{template}_{query_left}_{query_right}_{orientation}_R.pdb"
            native = case["native_path"]
            if not case["transformed_shape_checked"]:
                transformed_left_chains = chain_order(left) if left.is_file() else []
                transformed_right_chains = chain_order(right) if right.is_file() else []
                if left.is_file() and right.is_file():
                    case["transformed_shape_valid"] = (
                        len(transformed_left_chains) == len(case["query_left_chains"])
                        and len(transformed_right_chains) == len(case["query_right_chains"])
                    )
                    case["transformed_shape_checked"] = True
                    audit = chain_audit[cache_key]
                    audit["transformed_sample_left_chain_ids"] = "".join(transformed_left_chains)
                    audit["transformed_sample_right_chain_ids"] = "".join(transformed_right_chains)
                    audit["transformed_sample_left_chain_count"] = len(transformed_left_chains)
                    audit["transformed_sample_right_chain_count"] = len(transformed_right_chains)
                    audit["transformed_shape_matches_query"] = case["transformed_shape_valid"]
            left_chains = case["query_left_chains"]
            right_chains = case["query_right_chains"]
            pair = case["pair"]
            expected_left = str(pair.get("native_receptor_chains", "")).strip()
            expected_right = str(pair.get("native_ligand_chains", "")).strip()
            chain_shape_valid = (
                bool(left_chains) and bool(right_chains) and bool(native.is_file())
                and len(left_chains) == len(expected_left)
                and len(right_chains) == len(expected_right)
                and case["transformed_shape_valid"]
            )
            item = {
                "task_id": len(valid_rows),
                "key": candidate_key(source_row),
                "pipeline": pipeline,
                "dataset": "bm55_full",
                "split": entry.get("split", ""),
                "case_index": case_index,
                "case_id": entry.get("case_id", ""),
                "template": template,
                "orientation": orientation,
                "receptor": str(entry.get("receptor", "")).strip(),
                "ligand": str(entry.get("ligand", "")).strip(),
                "candidate_status": source_row.get("status", ""),
                "error_reason": source_row.get("error_reason", ""),
                "left": str(left),
                "right": str(right),
                "native_pdb": str(native),
                "model_left_chains": "".join(left_chains),
                "model_right_chains": "".join(right_chains),
                "native_receptor_chains": expected_left,
                "native_ligand_chains": expected_right,
                "benchmark_irmsd_A": pair.get("benchmark_irmsd_A"),
                "input_shape_valid": chain_shape_valid,
            }
            if not chain_shape_valid:
                invalid_rows.append({**item, "invalid_reason": "missing_input_or_chain_count_mismatch"})
                continue
            valid_rows.append(item)

            audit = chain_audit[cache_key]
            audit["model_left_chain_ids"] = item["model_left_chains"]
            audit["model_right_chain_ids"] = item["model_right_chains"]
            audit["model_left_chain_count"] = len(left_chains)
            audit["model_right_chain_count"] = len(right_chains)
            audit["chain_shape_valid"] = chain_shape_valid

    # Re-number after invalid rows are excluded, making array task IDs dense.
    for task_id, item in enumerate(valid_rows):
        item["task_id"] = task_id

    batch_rows: list[dict[str, Any]] = []
    for batch_id, start in enumerate(range(0, len(valid_rows), args.batch_size)):
        batch = valid_rows[start:start + args.batch_size]
        batch_path = batches_root / f"batch_{batch_id:05d}.jsonl"
        with batch_path.open("w", encoding="utf-8") as handle:
            for item in batch:
                handle.write(json.dumps(item, sort_keys=True) + "\n")
        batch_rows.append({
            "batch_id": batch_id,
            "path": str(batch_path),
            "count": len(batch),
            "first_task_id": batch[0]["task_id"] if batch else "",
            "last_task_id": batch[-1]["task_id"] if batch else "",
        })

    with (output_root / "chain_audit.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = sorted({key for row in chain_audit.values() for key in row})
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(sorted(chain_audit.values(), key=lambda row: (row["pipeline"], row["case_index"])))
    with (output_root / "invalid_candidates.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = sorted({key for row in invalid_rows for key in row}) if invalid_rows else ["invalid_reason"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(invalid_rows)
    with (output_root / "batch_index.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["batch_id", "path", "count", "first_task_id", "last_task_id"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(batch_rows)

    manifest = {
        "status": "ready",
        "created_at": now(),
        "run_root": str(run_root),
        "submission_manifest": str(submission_path),
        "submission_manifest_sha256": sha256(submission_path),
        "candidate_csv": str(candidate_path),
        "candidate_csv_sha256": sha256(candidate_path),
        "task_unit": "one_candidate_record_per downstream batch item",
        "batch_size": args.batch_size,
        "candidate_rows_seen": len(valid_rows) + len(invalid_rows),
        "candidate_rows_ready": len(valid_rows),
        "candidate_rows_invalid": len(invalid_rows),
        "case_chain_rows": len(chain_audit),
        "batch_count": len(batch_rows),
        "batch_index": str(output_root / "batch_index.csv"),
        "chain_audit": str(output_root / "chain_audit.csv"),
        "invalid_candidates": str(output_root / "invalid_candidates.csv"),
        "pipeline_manifests": {name: str(path) for name, path in pipeline_manifest_paths.items()},
        "resume_contract": "batch worker appends per-candidate checkpoint records and skips completed checkpoints",
        "dockq_contract": "assemble model chains, map to native receptor/ligand chains, calculate DockQ and backbone iRMSD",
        "pyrosetta_contract": "separate fail-closed preflight; no fallback to external Rosetta",
    }
    atomic_json(output_root / "downstream_manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
