#!/usr/bin/env python3
"""Score one resumable chunk of a corrected GTalign stage."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--pipeline-repo", required=True)
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--count", type=int, required=True)
    args = parser.parse_args()

    run_root = Path(args.run_root).resolve()
    pipeline_repo = Path(args.pipeline_repo).resolve()
    manifest_path = Path(args.input_manifest).resolve()
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    task_dir = output_root / "tasks"
    task_dir.mkdir(parents=True, exist_ok=True)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    entries = payload["entries"]
    dataset_manifest = json.loads((run_root / "dataset_manifest.json").read_text(encoding="utf-8"))
    pair_index = {
        (str(item["receptor"]).strip(), str(item["ligand"]).strip()): item
        for item in dataset_manifest["pairs"]
    }
    selected = entries[args.start : args.start + args.count]
    sys.path.insert(0, str(pipeline_repo))
    from src.eval.dockq import calculate_dockq, dockq_to_capri_class
    from src.eval.irmsd_backbone import calculate_irmsd_backbone

    records = []
    task_started = time.perf_counter()
    for entry in selected:
        started = time.perf_counter()
        row = {
            **entry,
            "run_root": str(run_root),
            "status": "pending",
            "dockq": None,
            "dockq_global": None,
            "dockq_sum": None,
            "dockq_capri": "",
            "dockq_fnat": None,
            "dockq_irmsd": None,
            "dockq_lrmsd": None,
            "irmsd_backbone": None,
            "dockq_json_status": "",
            "raw_dockq_json": "",
            "raw_dockq_json_sha256": "",
            "error": "",
        }
        try:
            model = Path(str(entry["model_pdb"])).resolve()
            if not model.is_file():
                raise FileNotFoundError(f"model PDB not found: {model}")
            pair = pair_index[(str(entry["receptor"]).strip(), str(entry["ligand"]).strip())]
            native = run_root / pair["native_pdb"]
            mapping = f"{entry['model_receptor_chains']}{entry['model_ligand_chains']}:{pair['native_receptor_chains']}{pair['native_ligand_chains']}"
            row.update(case_id=pair["case_id"], benchmark_split=pair["split"], benchmark_irmsd_A=pair.get("benchmark_irmsd_A"), native_pdb=str(native), dockq_mapping=mapping)
            work_dir = task_dir / f"dockq_tmp_{int(entry['index']):06d}"
            result = calculate_dockq(str(model), str(native), mapping=mapping, work_dir=str(work_dir), n_cpu=1)
            global_score = result.get("dockq_global")
            row.update(
                status="scored" if result.get("dockq_json_status") == "valid" else "valid_unscored",
                dockq=global_score,
                dockq_global=global_score,
                dockq_sum=result.get("dockq_sum"),
                dockq_capri=dockq_to_capri_class(global_score),
                dockq_fnat=result.get("fnat"),
                dockq_irmsd=result.get("irmsd"),
                dockq_lrmsd=result.get("lrmsd"),
                dockq_json_status=result.get("dockq_json_status", ""),
                raw_dockq_json=result.get("raw_dockq_json", ""),
                raw_dockq_json_sha256=result.get("raw_dockq_json_sha256", ""),
            )
            row["irmsd_backbone"] = calculate_irmsd_backbone(
                str(model), str(entry["model_receptor_chains"]), str(entry["model_ligand_chains"]),
                str(native), str(pair["native_receptor_chains"]), str(pair["native_ligand_chains"]),
            )
        except KeyError:
            row.update(status="out_of_scope", error="pair absent from dataset manifest")
        except Exception as exc:  # preserve per-candidate failure without dropping the row
            row.update(status="score_failed", error=f"{type(exc).__name__}: {exc}")
        row["elapsed_seconds"] = time.perf_counter() - started
        records.append(row)

    task_name = f"task_{args.start:06d}_{args.count:06d}.json"
    task_path = task_dir / task_name
    task_payload = {
        "schema_version": "prism-gtalign-score-task-20260920",
        "stage": payload["stage"],
        "start": args.start,
        "count": args.count,
        "selected": len(selected),
        "status_counts": {status: sum(row["status"] == status for row in records) for status in sorted({row["status"] for row in records})},
        "elapsed_seconds": time.perf_counter() - task_started,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
        "records": records,
    }
    task_path.write_text(json.dumps(task_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"stage": payload["stage"], "task": str(task_path), "selected": len(records), "elapsed_seconds": task_payload["elapsed_seconds"], "status_counts": task_payload["status_counts"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
