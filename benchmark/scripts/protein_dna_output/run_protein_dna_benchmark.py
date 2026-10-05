#!/usr/bin/env python3
"""
Run a curated protein-DNA benchmark manifest and emit per_prediction.csv and pair_summary.csv.
"""

import argparse
import csv
import statistics
from pathlib import Path

try:
    from benchmark.scripts.protein_dna_output.score_single_protein_dna_pair import parse_chain_list, score_one
except ImportError:
    from score_single_protein_dna_pair import parse_chain_list, score_one


def to_float(value):
    if value in ("", None, "NA"):
        return None
    try:
        return float(value)
    except Exception:
        return None


def summarize(values, best_mode):
    if not values:
        return None, None, None
    if len(values) == 1:
        return float(values[0]), 0.0, float(values[0])
    mean_value = float(sum(values)) / float(len(values))
    variance_value = float(statistics.pvariance(values))
    best_value = max(values) if best_mode == "max" else min(values)
    return mean_value, variance_value, best_value


def write_csv(path, rows):
    fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_manifest(manifest_csv, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(manifest_csv, newline="") as handle:
        manifest_rows = list(csv.DictReader(handle))

    per_prediction = []
    grouped = {}
    manifest_root = Path(manifest_csv).resolve().parent

    for row in manifest_rows:
        model_path = (manifest_root / row["model_path"]).resolve()
        native_path = (manifest_root / row["native_pdb"]).resolve()
        score = score_one(
            str(model_path),
            str(native_path),
            model_protein_chains=parse_chain_list(row.get("model_protein_chains")) or None,
            model_dna_chains=parse_chain_list(row.get("model_dna_chains")) or None,
            native_protein_chains=parse_chain_list(row.get("native_protein_chains")) or None,
            native_dna_chains=parse_chain_list(row.get("native_dna_chains")) or None,
            score_json=(manifest_root / row["score_json"]).resolve() if row.get("score_json") else None,
        )
        score.update(
            {
                "pair_id": row["pair_id"],
                "template_id": row.get("template_id", ""),
                "benchmark_case": row.get("benchmark_case", ""),
            }
        )
        per_prediction.append(score)
        grouped.setdefault(row["pair_id"], []).append(score)

    summary_rows = []
    for pair_id, rows in grouped.items():
        cp_vals = [to_float(r.get("contact_precision")) for r in rows if to_float(r.get("contact_precision")) is not None]
        cr_vals = [to_float(r.get("contact_recall")) for r in rows if to_float(r.get("contact_recall")) is not None]
        pir_vals = [to_float(r.get("protein_interface_recall")) for r in rows if to_float(r.get("protein_interface_recall")) is not None]
        ncr_vals = [to_float(r.get("nucleotide_contact_recall")) for r in rows if to_float(r.get("nucleotide_contact_recall")) is not None]
        tm_vals = [to_float(r.get("alignment_tm_score")) for r in rows if to_float(r.get("alignment_tm_score")) is not None]
        rmsd_vals = [to_float(r.get("alignment_rmsd")) for r in rows if to_float(r.get("alignment_rmsd")) is not None]

        cp_mean, cp_var, cp_best = summarize(cp_vals, "max")
        cr_mean, cr_var, cr_best = summarize(cr_vals, "max")
        pir_mean, pir_var, pir_best = summarize(pir_vals, "max")
        ncr_mean, ncr_var, ncr_best = summarize(ncr_vals, "max")
        tm_mean, tm_var, tm_best = summarize(tm_vals, "max")
        rmsd_mean, rmsd_var, rmsd_best = summarize(rmsd_vals, "min")

        best_tm_row = max((r for r in rows if to_float(r.get("alignment_tm_score")) is not None), key=lambda r: r["alignment_tm_score"], default=None)
        best_contact_row = max((r for r in rows if to_float(r.get("contact_recall")) is not None), key=lambda r: r["contact_recall"], default=None)

        summary_rows.append(
            {
                "pair_id": pair_id,
                "n_predictions_total": len(rows),
                "n_scored": len(cp_vals),
                "contact_precision_mean": cp_mean,
                "contact_precision_variance": cp_var,
                "contact_precision_best": cp_best,
                "contact_recall_mean": cr_mean,
                "contact_recall_variance": cr_var,
                "contact_recall_best": cr_best,
                "protein_interface_recall_mean": pir_mean,
                "protein_interface_recall_variance": pir_var,
                "protein_interface_recall_best": pir_best,
                "nucleotide_contact_recall_mean": ncr_mean,
                "nucleotide_contact_recall_variance": ncr_var,
                "nucleotide_contact_recall_best": ncr_best,
                "alignment_tm_score_mean": tm_mean,
                "alignment_tm_score_variance": tm_var,
                "alignment_tm_score_best": tm_best,
                "alignment_rmsd_mean": rmsd_mean,
                "alignment_rmsd_variance": rmsd_var,
                "alignment_rmsd_best_min": rmsd_best,
                "best_tm_model": best_tm_row.get("model_pdb") if best_tm_row else None,
                "best_contact_model": best_contact_row.get("model_pdb") if best_contact_row else None,
            }
        )

    write_csv(output_dir / "per_prediction.csv", per_prediction)
    write_csv(output_dir / "pair_summary.csv", summary_rows)


def main():
    parser = argparse.ArgumentParser(description="Run the curated protein-DNA benchmark manifest")
    parser.add_argument("--manifest", default="benchmark/data/protein_dna_curated_manifest.csv")
    parser.add_argument("--output-dir", default="benchmark/scripts/protein_dna_output/results")
    args = parser.parse_args()
    run_manifest(args.manifest, args.output_dir)


if __name__ == "__main__":
    main()
