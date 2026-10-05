#!/usr/bin/env python3
"""Score one resumable downstream batch with DockQ and backbone iRMSD."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
import time
from pathlib import Path
from typing import Any


def now() -> str:
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat()


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
                chain = line[21].strip() or "_"
                if chain not in seen:
                    seen.append(chain)
    return seen


def native_mapping_preflight(
    native: Path,
    native_receptor: object,
    native_ligand: object,
) -> dict[str, Any]:
    """Validate native reference chains before assembling a DockQ model."""
    receptor = re.sub(r"\*+$", "", "".join(str(native_receptor or "").split()))
    ligand = re.sub(r"\*+$", "", "".join(str(native_ligand or "").split()))
    base: dict[str, Any] = {
        "native_receptor_chains": receptor,
        "native_ligand_chains": ligand,
        "native_mapping_status": "valid",
        "scoring_scope": "complete_complex",
    }
    if not native.is_file():
        return {
            **base,
            "native_mapping_status": "invalid_mapping",
            "scoring_scope": "none",
            "reason": f"native PDB is missing: {native}",
        }
    available = chain_order(native)
    missing = [chain for chain in receptor + ligand if chain not in available]
    overlap = sorted(set(receptor) & set(ligand))
    base.update({
        "native_available_chains": "".join(available),
        "native_missing_chains": "".join(missing),
    })
    if overlap:
        return {
            **base,
            "native_mapping_status": "invalid_mapping",
            "scoring_scope": "none",
            "reason": f"native receptor/ligand groups overlap: {overlap}",
        }
    if not receptor or not ligand:
        if len(available) == 1 and len(receptor + ligand) <= 1:
            return {
                **base,
                "native_mapping_status": "valid_unscored",
                "scoring_scope": "no_protein_protein_interface",
                "reason": "native reference contains only one protein chain",
            }
        return {
            **base,
            "native_mapping_status": "invalid_mapping",
            "scoring_scope": "none",
            "reason": "both native receptor and ligand chain groups are required",
        }
    if missing:
        return {
            **base,
            "native_mapping_status": "invalid_mapping",
            "scoring_scope": "none",
            "reason": f"native mapping references absent chain(s): {missing}",
        }
    return base


def assemble(left: Path, right: Path, output: Path, receptor_chains: str, ligand_chains: str) -> tuple[str, str]:
    left_raw = chain_order(left)
    right_raw = chain_order(right)
    if len(left_raw) != len(receptor_chains) or len(right_raw) != len(ligand_chains):
        raise ValueError(
            f"chain count mismatch left={left_raw!r}/{receptor_chains!r} right={right_raw!r}/{ligand_chains!r}"
        )
    if set(receptor_chains) & set(ligand_chains):
        raise ValueError(f"native receptor/ligand chain IDs overlap: {receptor_chains!r}/{ligand_chains!r}")
    left_map = dict(zip(left_raw, receptor_chains))
    right_map = dict(zip(right_raw, ligand_chains))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for source, mapping in ((left, left_map), (right, right_map)):
            with source.open(encoding="utf-8", errors="replace") as source_handle:
                for line in source_handle:
                    if line.startswith(("ATOM", "HETATM", "TER")):
                        if len(line) > 21 and line[21] in mapping:
                            line = line[:21] + mapping[line[21]] + line[22:]
                        handle.write(line)
            handle.write("TER\n")
        handle.write("END\n")
    return receptor_chains, ligand_chains


def result_row(item: dict[str, Any]) -> dict[str, Any]:
    return {
        "task_id": item.get("task_id", ""),
        "key": item.get("key", ""),
        "pipeline": item.get("pipeline", ""),
        "dataset": item.get("dataset", ""),
        "split": item.get("split", ""),
        "case_index": item.get("case_index", ""),
        "case_id": item.get("case_id", ""),
        "template": item.get("template", ""),
        "orientation": item.get("orientation", ""),
        "receptor": item.get("receptor", ""),
        "ligand": item.get("ligand", ""),
        "candidate_status": item.get("candidate_status", ""),
        "model_left_chains": item.get("model_left_chains", ""),
        "model_right_chains": item.get("model_right_chains", ""),
        "native_receptor_chains": item.get("native_receptor_chains", ""),
        "native_ligand_chains": item.get("native_ligand_chains", ""),
        "native_available_chains": "",
        "native_missing_chains": "",
        "native_mapping_status": "",
        "scoring_scope": "",
        "benchmark_irmsd_A": item.get("benchmark_irmsd_A", ""),
        "left_pdb": item.get("left", ""),
        "right_pdb": item.get("right", ""),
        "model_pdb": "",
        "assembled_model_retained": False,
        "native_pdb": item.get("native_pdb", ""),
        "dockq_mapping": "",
        "dockq": "",
        "dockq_capri": "",
        "dockq_fnat": "",
        "dockq_irmsd": "",
        "dockq_lrmsd": "",
        "irmsd_backbone": "",
        "status": "pending",
        "error": "",
        "error_class": "",
        "elapsed_seconds": "",
        "scored_at": "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-file", required=True, type=Path)
    parser.add_argument("--checkpoint-file", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--pipeline-repo", required=True, type=Path)
    parser.add_argument("--dockq-python", required=True, type=Path)
    args = parser.parse_args()
    batch_path = args.batch_file.resolve()
    checkpoint_path = args.checkpoint_file.resolve()
    output_csv = args.output_csv.resolve()
    repo = args.pipeline_repo.resolve()
    dockq_python = args.dockq_python.resolve()
    if not dockq_python.is_file():
        raise SystemExit(f"missing DockQ Python: {dockq_python}")
    sys.path.insert(0, str(repo))
    from src.eval.dockq import calculate_dockq, dockq_to_capri_class
    from src.eval.irmsd_backbone import calculate_irmsd_backbone

    entries = []
    with batch_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                entries.append(json.loads(line))
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    completed: dict[str, dict[str, Any]] = {}
    if checkpoint_path.is_file():
        with checkpoint_path.open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    item = json.loads(line)
                    if item.get("status") in {
                        "scored", "score_failed", "scorer_error", "invalid",
                        "invalid_mapping", "valid_unscored",
                    }:
                        completed[item["key"]] = item

    started = time.perf_counter()
    with checkpoint_path.open("a", encoding="utf-8") as checkpoint_handle:
        for entry in entries:
            key = entry["key"]
            if key in completed:
                continue
            row = result_row(entry)
            item_started = time.perf_counter()
            try:
                left = Path(entry["left"])
                right = Path(entry["right"])
                native = Path(entry["native_pdb"])
                if not left.is_file() or not right.is_file() or not native.is_file():
                    raise FileNotFoundError(f"missing input left={left} right={right} native={native}")
                native_info = native_mapping_preflight(
                    native, entry["native_receptor_chains"], entry["native_ligand_chains"]
                )
                row.update({key: value for key, value in native_info.items() if key != "reason"})
                if native_info["native_mapping_status"] == "valid_unscored":
                    row.update(status="valid_unscored", error=native_info.get("reason", ""))
                    row["elapsed_seconds"] = time.perf_counter() - item_started
                    row["scored_at"] = now()
                    checkpoint_handle.write(json.dumps(row, sort_keys=True) + "\n")
                    checkpoint_handle.flush()
                    os.fsync(checkpoint_handle.fileno())
                    completed[key] = row
                    continue
                if native_info["native_mapping_status"] != "valid":
                    row.update(
                        status="invalid_mapping",
                        error=native_info.get("reason", ""),
                    )
                    row["elapsed_seconds"] = time.perf_counter() - item_started
                    row["scored_at"] = now()
                    checkpoint_handle.write(json.dumps(row, sort_keys=True) + "\n")
                    checkpoint_handle.flush()
                    os.fsync(checkpoint_handle.fileno())
                    completed[key] = row
                    continue
                native_receptor = str(native_info["native_receptor_chains"])
                native_ligand = str(native_info["native_ligand_chains"])
                with tempfile.TemporaryDirectory(prefix=f"dockq_{key}_") as temp_dir:
                    model = Path(temp_dir) / f"{key}.pdb"
                    model_rec, model_lig = assemble(
                        left, right, model,
                        native_receptor,
                        native_ligand,
                    )
                    mapping = f"{model_rec}{model_lig}:{native_receptor}{native_ligand}"
                    scores = calculate_dockq(
                        str(model), str(native), mapping=mapping,
                        work_dir=str(Path(temp_dir) / "dockq_tmp"),
                    )
                    row.update(
                        model_pdb="",
                        assembled_model_retained=False,
                        left_pdb=str(left),
                        right_pdb=str(right),
                        dockq_mapping=mapping,
                        dockq=scores.get("dockq"),
                        dockq_capri=dockq_to_capri_class(scores.get("dockq")),
                        dockq_fnat=scores.get("fnat"),
                        dockq_irmsd=scores.get("irmsd"),
                        dockq_lrmsd=scores.get("lrmsd"),
                        irmsd_backbone=calculate_irmsd_backbone(
                            str(model), model_rec, model_lig, str(native),
                            native_receptor, native_ligand,
                        ),
                        status="scored",
                    )
            except Exception as exc:
                row.update(
                    status="scorer_error",
                    error=f"{type(exc).__name__}: {exc}",
                    error_class=type(exc).__name__,
                )
            row["elapsed_seconds"] = time.perf_counter() - item_started
            row["scored_at"] = now()
            checkpoint_handle.write(json.dumps(row, sort_keys=True) + "\n")
            checkpoint_handle.flush()
            os.fsync(checkpoint_handle.fileno())
            completed[key] = row

    rows = [completed.get(entry["key"], result_row(entry)) for entry in entries]
    fields = list(rows[0]) if rows else list(result_row({}))
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_csv.with_name(output_csv.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output_csv)
    atomic_json(output_csv.with_suffix(".summary.json"), {
        "status": "completed",
        "batch_file": str(batch_path),
        "checkpoint_file": str(checkpoint_path),
        "output_csv": str(output_csv),
        "rows": len(rows),
        "scored": sum(row.get("status") == "scored" for row in rows),
        "score_failed": sum(row.get("status") in {"score_failed", "scorer_error"} for row in rows),
        "invalid_mapping": sum(row.get("status") == "invalid_mapping" for row in rows),
        "valid_unscored": sum(row.get("status") == "valid_unscored" for row in rows),
        "elapsed_seconds": time.perf_counter() - started,
        "job_id": os.environ.get("SLURM_JOB_ID", ""),
        "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
    })
    print(json.dumps({
        "status": "completed",
        "rows": len(rows),
        "scored": sum(row.get("status") == "scored" for row in rows),
        "score_failed": sum(row.get("status") in {"score_failed", "scorer_error"} for row in rows),
        "invalid_mapping": sum(row.get("status") == "invalid_mapping" for row in rows),
        "valid_unscored": sum(row.get("status") == "valid_unscored" for row in rows),
        "output_csv": str(output_csv),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
