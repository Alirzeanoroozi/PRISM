#!/usr/bin/env python3
"""Prepare one resumable full-template task for every BM55 input pair."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
try:
    from valar_agent.prism_cpu_batch import build_batch_plan, sha256_file
except ModuleNotFoundError:  # remote bundle staging includes this module beside the script
    sys.path.insert(0, str(SCRIPT_DIR))
    from prism_cpu_batch import build_batch_plan, sha256_file


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


def same_plan(existing: dict, plan: dict) -> bool:
    return all(existing.get(key) == plan.get(key) for key in (
        "dataset_dir",
        "template_list",
        "template_count",
        "template_sha256",
        "template_chain_slots",
        "alignment_records_per_case",
        "aligner",
        "pair_count",
    ))


def ensure_source_link(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        if destination.resolve() != source.resolve():
            raise RuntimeError(f"source link points elsewhere: {destination}")
        return
    if destination.exists():
        raise RuntimeError(f"refusing to replace existing source path: {destination}")
    destination.symlink_to(source.resolve())


def stage_case_dataset(entry: dict, original_pair: dict) -> None:
    case_dataset = Path(entry["case_dataset_dir"])
    case_dataset.mkdir(parents=True, exist_ok=True)
    source_dir = case_dataset / "sources"
    ensure_source_link(Path(entry["native_source"]), source_dir / entry["native_name"])
    for query_name, query_source in zip(entry["query_source_names"], entry["query_sources"]):
        ensure_source_link(Path(query_source), source_dir / query_name)

    atomic_text(
        case_dataset / "inputs.csv",
        "Receptor,Ligand\n" + f"{entry['receptor']},{entry['ligand']}\n",
    )
    atomic_json(
        case_dataset / "dataset_manifest.json",
        {
            "dataset": "bm55_full",
            "pair_count": 1,
            "pairs": [original_pair],
            "prepared_from_full_manifest": True,
        },
    )


def load_original_pairs(dataset_dir: Path) -> list[dict]:
    payload = json.loads((dataset_dir / "dataset_manifest.json").read_text(encoding="utf-8"))
    return list(payload["pairs"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-dir", required=True)
    parser.add_argument("--template-list", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--aligner", choices=("tmalign", "multiprot", "usalign"), required=True)
    parser.add_argument("--pipeline-repo", required=True)
    parser.add_argument("--pipeline-python", required=True)
    parser.add_argument("--helper", required=True)
    parser.add_argument("--timed-runner", required=True)
    parser.add_argument("--usalign-wrapper", required=True)
    parser.add_argument("--expected-template-count", type=int, default=19_855)
    args = parser.parse_args()

    dataset_dir = Path(args.dataset_dir).resolve()
    template_list = Path(args.template_list).resolve()
    output_root = Path(args.output_root).resolve()
    plan = build_batch_plan(
        dataset_dir=dataset_dir,
        template_list=template_list,
        output_root=output_root,
        aligner=args.aligner,
        expected_template_count=args.expected_template_count,
    )
    manifest_path = output_root / "batch_manifest.json"
    if manifest_path.is_file():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not same_plan(existing, plan):
            raise SystemExit("existing batch manifest does not match the requested immutable inputs")
        print(json.dumps({"status": "resumed", "manifest": str(manifest_path)}, sort_keys=True))
        return 0
    if output_root.exists() and any(output_root.iterdir()):
        raise SystemExit(f"refusing to initialize a non-empty output root without a matching manifest: {output_root}")

    original_pairs = load_original_pairs(dataset_dir)
    for entry, original_pair in zip(plan["entries"], original_pairs):
        stage_case_dataset(entry, original_pair)
        status_path = Path(entry["case_root"]) / "task_status.json"
        if not status_path.exists():
            atomic_json(
                status_path,
                {
                    "status": "not_started",
                    "index": entry["index"],
                    "case_id": entry["case_id"],
                    "updated_at": datetime.now(timezone.utc).isoformat(),
                },
            )

    manifest = {
        **plan,
        "status": "ready",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "job_submission_scope": "full_bm55_all_pairs_all_templates",
        "pipeline_repo": str(Path(args.pipeline_repo).resolve()),
        "pipeline_python": str(Path(args.pipeline_python).resolve()),
        "helper": str(Path(args.helper).resolve()),
        "timed_runner": str(Path(args.timed_runner).resolve()),
        "usalign_wrapper": str(Path(args.usalign_wrapper).resolve()),
        "rank": False,
        "rank_method": "baseline",
        "prodigy": False,
        "refinement": False,
        "template_file_sha256": sha256_file(template_list),
    }
    atomic_json(manifest_path, manifest)
    print(json.dumps({"status": "ready", "manifest": str(manifest_path), "task_count": manifest["pair_count"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
