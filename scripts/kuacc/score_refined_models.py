#!/usr/bin/env python3
"""Score full refined PRISM complexes with the same DockQ/iRMSD functions."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--pipeline-repo", required=True)
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()

    run_root = Path(args.run_root).resolve()
    pipeline_repo = Path(args.pipeline_repo).resolve()
    input_manifest = Path(args.input_manifest).resolve()
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    import sys

    sys.path.insert(0, str(run_root))
    sys.path.insert(0, str(pipeline_repo))
    from src.eval.dockq import calculate_dockq, dockq_to_capri_class
    from src.eval.irmsd_backbone import calculate_irmsd_backbone

    run_summary = json.loads((run_root / "run_summary.json").read_text(encoding="utf-8")) if (run_root / "run_summary.json").is_file() else {}
    dataset_manifest = json.loads((run_root / "dataset_manifest.json").read_text(encoding="utf-8"))
    pair_index = {
        (str(row["receptor"]).strip(), str(row["ligand"]).strip()): row
        for row in dataset_manifest.get("pairs", [])
    }
    entries = json.loads(input_manifest.read_text(encoding="utf-8"))
    rows = []
    started = time.perf_counter()
    work_dir = output_root / "dockq_tmp"
    work_dir.mkdir(parents=True, exist_ok=True)
    for entry in entries:
        row = {
            "stage": entry.get("stage", "refined"),
            "dataset": run_summary.get("dataset", dataset_manifest.get("dataset", "")),
            "scenario": run_summary.get("scenario", ""),
            "aligner": run_summary.get("aligner", ""),
            "template": entry.get("template", ""),
            "receptor": entry.get("receptor", ""),
            "ligand": entry.get("ligand", ""),
            "case_id": "",
            "benchmark_split": "",
            "benchmark_irmsd_A": None,
            "orientation": entry.get("orientation", ""),
            "model_pdb": entry.get("model_pdb", ""),
            "native_pdb": "",
            "dockq": None,
            "dockq_capri": "",
            "dockq_fnat": None,
            "dockq_irmsd": None,
            "dockq_lrmsd": None,
            "irmsd_backbone": None,
            "dockq_mapping": "",
            "model_receptor_chains": entry.get("model_receptor_chains", ""),
            "model_ligand_chains": entry.get("model_ligand_chains", ""),
            "status": "pending",
            "error": "",
        }
        pair = pair_index.get((str(entry.get("receptor", "")).strip(), str(entry.get("ligand", "")).strip()))
        if pair is None:
            row.update(status="out_of_scope", error="no dataset manifest row for model pair")
            rows.append(row)
            continue
        row.update(
            dataset=dataset_manifest.get("dataset", ""),
            case_id=pair.get("case_id", ""),
            benchmark_split=pair.get("split", ""),
            benchmark_irmsd_A=pair.get("benchmark_irmsd_A"),
        )
        native = run_root / pair["native_pdb"]
        row["native_pdb"] = str(native)
        model = Path(entry["model_pdb"]).resolve()
        mapping = f"{entry['model_receptor_chains']}{entry['model_ligand_chains']}:{pair['native_receptor_chains']}{pair['native_ligand_chains']}"
        row["dockq_mapping"] = mapping
        try:
            dockq = calculate_dockq(str(model), str(native), mapping=mapping, work_dir=str(work_dir))
            row.update(
                dockq=dockq.get("dockq"),
                dockq_capri=dockq_to_capri_class(dockq.get("dockq")),
                dockq_fnat=dockq.get("fnat"),
                dockq_irmsd=dockq.get("irmsd"),
                dockq_lrmsd=dockq.get("lrmsd"),
            )
            row["irmsd_backbone"] = calculate_irmsd_backbone(
                str(model),
                entry["model_receptor_chains"],
                entry["model_ligand_chains"],
                str(native),
                str(pair["native_receptor_chains"]),
                str(pair["native_ligand_chains"]),
            )
            row["status"] = "scored"
        except Exception as exc:
            row.update(status="score_failed", error=f"{type(exc).__name__}: {exc}")
        rows.append(row)

    columns = sorted({key for row in rows for key in row})
    output_csv = output_root / "dockq_irmsd.csv"
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "status": "completed",
        "stage": entries[0].get("stage", "refined") if entries else "refined",
        "model_groups": len(entries),
        "rows": len(rows),
        "scored": sum(row.get("status") == "scored" for row in rows),
        "score_failed": sum(row.get("status") == "score_failed" for row in rows),
        "out_of_scope": sum(row.get("status") == "out_of_scope" for row in rows),
        "unscored": sum(row.get("status") != "scored" for row in rows),
        "elapsed_seconds": time.perf_counter() - started,
        "csv": str(output_csv),
        "input_manifest": str(input_manifest),
    }
    (output_root / "score_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
