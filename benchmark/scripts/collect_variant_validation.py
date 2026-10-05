#!/usr/bin/env python3
"""Collect comparable runtime, failure, and score summaries for variant arms."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def mean(rows: list[dict[str, str]], key: str) -> float | None:
    values = []
    for row in rows:
        try:
            values.append(float(row[key]))
        except (KeyError, TypeError, ValueError):
            pass
    return sum(values) / len(values) if values else None


def summarize(name: str, run_root: Path, score_csv: Path | None) -> dict[str, object]:
    exits = list(run_root.glob("batch_*/status/exit.json"))
    exit_rows = [json.loads(path.read_text()) for path in exits]
    alignments = sum(len(list(path.glob("processed/alignment/**/*.json"))) for path in run_root.glob("batch_*"))
    gt_alignments = sum(len(list(path.glob("processed/alignment_gtalign/**/*.json"))) for path in run_root.glob("batch_*"))
    transformations = sum(len(list(path.glob("processed/transformation/*.pdb"))) for path in run_root.glob("batch_*"))
    external_models = sum(len(list(path.glob("processed/rosetta_refinement/*_rosetta_0001.pdb"))) for path in run_root.glob("batch_*"))
    pyro_models = sum(len(list(path.glob("processed/pyrosetta_refinement/structures/*_rosetta.pdb"))) for path in run_root.glob("batch_*"))
    scores = []
    if score_csv and score_csv.is_file():
        with score_csv.open(newline="", encoding="utf-8") as handle:
            scores = list(csv.DictReader(handle))
    return {
        "variant": name,
        "batches": len(exits),
        "completed_batches": sum(row.get("return_code") == 0 for row in exit_rows),
        "failed_batches": sum(row.get("return_code") != 0 for row in exit_rows),
        "elapsed_seconds": sum(int(row.get("elapsed_seconds", 0)) for row in exit_rows),
        "alignment_json": alignments + gt_alignments,
        "transformation_halves": transformations,
        "external_models": external_models,
        "pyrosetta_models": pyro_models,
        "score_rows": len(scores),
        "scoreable_rows": sum(row.get("score_status") == "scored" for row in scores),
        "dockq_cross_mean": mean(scores, "dockq_cross_mean"),
        "dockq_cross_best_mean": mean(scores, "dockq_cross_best"),
        "irmsd_grouped_mean": mean(scores, "irmsd_grouped_min"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", action="append", required=True, help="NAME=RUN_ROOT[=SCORE_CSV]")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for spec in args.variant:
        parts = spec.split("=", 2)
        if len(parts) < 2:
            parser.error("--variant must be NAME=RUN_ROOT[=SCORE_CSV]")
        rows.append(summarize(parts[0], Path(parts[1]), Path(parts[2]) if len(parts) == 3 else None))
    fields = list(rows[0]) if rows else ["variant"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} variant summaries to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
