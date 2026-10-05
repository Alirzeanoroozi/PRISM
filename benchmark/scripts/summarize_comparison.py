#!/usr/bin/env python3
"""Summarize full current-vs-legacy benchmark outcomes and model scores."""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def floats(rows, key):
    return [float(row[key]) for row in rows if row.get(key, "") not in ("", None)]


def fmt(value):
    return "n/a" if value is None else f"{value:.3f}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", type=Path, required=True)
    parser.add_argument("--pairwise", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    status = read_csv(args.status)
    pairwise = read_csv(args.pairwise)
    scores = read_csv(args.scores)

    summary_rows = []
    for pipeline in ("tmalign_rosetta", "multiprot_fiberdock"):
        rows = [r for r in scores if r.get("pipeline") == pipeline]
        ready = [r for r in rows if r.get("status") == "ready"]
        dockq = floats(ready, "dockq")
        irmsd = floats(ready, "irmsd")
        summary_rows.append({
            "pipeline": pipeline,
            "model_rows": len(rows),
            "ready_rows": len(ready),
            "dockq_n": len(dockq),
            "dockq_mean": fmt(statistics.mean(dockq) if dockq else None),
            "dockq_median": fmt(statistics.median(dockq) if dockq else None),
            "dockq_best": fmt(max(dockq) if dockq else None),
            "irmsd_n": len(irmsd),
            "irmsd_mean": fmt(statistics.mean(irmsd) if irmsd else None),
            "irmsd_median": fmt(statistics.median(irmsd) if irmsd else None),
            "irmsd_best": fmt(min(irmsd) if irmsd else None),
            "score_errors": sum(bool(r.get("score_error")) for r in ready),
        })
    with (args.output_dir / "model_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader(); writer.writerows(summary_rows)

    pair_rows = []
    for pair in pairwise:
        out = dict(pair)
        for pipeline in ("tmalign_rosetta", "multiprot_fiberdock"):
            rows = [r for r in scores if r.get("pipeline") == pipeline and r.get("pair_id") == pair["pair_id"] and r.get("status") == "ready"]
            dockq = floats(rows, "dockq")
            irmsd = floats(rows, "irmsd")
            out[f"{pipeline}_models"] = len(rows)
            out[f"{pipeline}_best_dockq"] = fmt(max(dockq) if dockq else None)
            out[f"{pipeline}_best_irmsd"] = fmt(min(irmsd) if irmsd else None)
        c = out.get("tmalign_rosetta_best_dockq")
        l = out.get("multiprot_fiberdock_best_dockq")
        out["dockq_delta_current_minus_legacy"] = fmt(float(c) - float(l)) if c != "n/a" and l != "n/a" else "n/a"
        pair_rows.append(out)
    with (args.output_dir / "pairwise_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]))
        writer.writeheader(); writer.writerows(pair_rows)

    status_counts = {}
    for row in status:
        key = (row["pipeline"], row["status"])
        status_counts[key] = status_counts.get(key, 0) + 1
    classes = {}
    for row in pairwise:
        classes[row["comparison_class"]] = classes.get(row["comparison_class"], 0) + 1
    current_score_pairs = sum(float(row["tmalign_rosetta_models"]) > 0 for row in pair_rows)
    legacy_score_pairs = sum(float(row["multiprot_fiberdock_models"]) > 0 for row in pair_rows)
    shared_score_pairs = sum(float(row["tmalign_rosetta_models"]) > 0 and float(row["multiprot_fiberdock_models"]) > 0 for row in pair_rows)

    lines = ["# Full current-vs-legacy comparison", "", "## Scope", "", "- Shared manifest: 257 pairs, 26 batches per pipeline (25×10 plus a final batch of 7).", "- PDB staging: 470/472 unique structures available; `1erk` and `4zai` remained explicit input-missing cases.", "- Current pipeline: TMalign + Rosetta refinement. Legacy pipeline: MultiProt + FiberDock.", "- Scores are model-level DockQ/iRMSD values; they are summarized separately by pipeline and joined pairwise by `pair_id`.", "", "## Batch outcomes", "", "| pipeline | status | pairs |", "|---|---:|---:|"]
    for (pipeline, state), count in sorted(status_counts.items()):
        lines.append(f"| {pipeline} | {state} | {count} |")
    lines += ["", "| pairwise class | pairs |", "|---|---:|"]
    for key, count in sorted(classes.items()):
        lines.append(f"| {key} | {count} |")
    lines += ["", "## Model-level scores", "", "| pipeline | model rows | ready | DockQ n | DockQ mean | DockQ median | DockQ best | iRMSD n | iRMSD mean | iRMSD median | iRMSD best | score errors |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in summary_rows:
        lines.append("| " + " | ".join(str(value) for value in row.values()) + " |")
    lines += ["", "## Interpretation", "", "- Observation: both pipelines completed 255/257 input pairs; the same two pairs were unavailable during staging.", f"- Observation: the current run produced {sum(int(row['model_rows']) for row in summary_rows if row['pipeline'] == 'tmalign_rosetta')} Rosetta models across {current_score_pairs} pairs with at least one scoreable model.", f"- Observation: the legacy run produced {sum(int(row['model_rows']) for row in summary_rows if row['pipeline'] == 'multiprot_fiberdock')} FiberDock models for {legacy_score_pairs} pair; there were {shared_score_pairs} pairs with scoreable models from both pipelines.", "- Inference: the dominant legacy limitation in this run is candidate-generation/refinement attrition, so model-quality averages are not a fair standalone method comparison; coverage and failure mode must be reported alongside scores.", "- Inference: current-vs-legacy score deltas are only interpretable for pairs with scoreable models and validated chain mappings. Rosetta and FiberDock energy values are not compared numerically here.", "", "Generated from `final_status.csv`, `final_pairwise.csv`, and `scored_models.csv`."]
    (args.output_dir / "FINAL_REPORT.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {args.output_dir / 'FINAL_REPORT.md'}")


if __name__ == "__main__":
    main()
