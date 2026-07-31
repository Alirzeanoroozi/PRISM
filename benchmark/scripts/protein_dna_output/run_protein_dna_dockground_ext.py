#!/usr/bin/env python3
"""
Run the PRISM-main protein_dna_ext branch on the Dockground subset and score the outputs.
"""

import argparse
import csv
import json
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

try:
    from benchmark.scripts.protein_dna_output.score_single_protein_dna_pair_ext import parse_chain_list, score_one_ext
except ImportError:
    from score_single_protein_dna_pair_ext import parse_chain_list, score_one_ext


def case_to_target_ids(case_id):
    base = case_id.lower()[1:]
    return f"{base}p", f"{base}d"


def write_csv(path, rows, fieldnames=None):
    fieldnames = fieldnames or (list(rows[0].keys()) if rows else [])
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize(values, best_mode):
    if not values:
        return None, None, None
    if len(values) == 1:
        return float(values[0]), 0.0, float(values[0])
    mean_value = float(sum(values)) / float(len(values))
    variance_value = float(statistics.pvariance(values))
    best_value = max(values) if best_mode == "max" else min(values)
    return mean_value, variance_value, best_value


def normalize_chain_field(value):
    if value in ("", None, "NA"):
        return []
    if isinstance(value, list):
        return value
    text = str(value).strip()
    if text.startswith("[") and text.endswith("]"):
        try:
            parsed = json.loads(text.replace("'", "\""))
            return [str(item) for item in parsed]
        except Exception:
            pass
    return parse_chain_list(text)


def stage_workspace(manifest_csv, workspace_dir):
    manifest_path = Path(manifest_csv).resolve()
    workspace_dir = Path(workspace_dir)
    if workspace_dir.exists():
        shutil.rmtree(workspace_dir)
    workspace_dir.mkdir(parents=True, exist_ok=True)
    (workspace_dir / "processed/pdbs").mkdir(parents=True, exist_ok=True)
    (workspace_dir / "templates/pdbs").mkdir(parents=True, exist_ok=True)
    (workspace_dir / "templates").mkdir(parents=True, exist_ok=True)
    (workspace_dir / "dna_run_logs").mkdir(parents=True, exist_ok=True)

    with open(manifest_path, newline="") as handle:
        rows = list(csv.DictReader(handle))

    inputs_rows = []
    checked_templates = []
    chain_hints = {}
    for row in rows:
        protein_target, dna_target = case_to_target_ids(row["case_id"])
        inputs_rows.append({"Receptor": protein_target, "Ligand": dna_target})
        shutil.copy2(row["protein_unbound_pdb"], workspace_dir / "processed/pdbs" / f"{protein_target[:4]}.pdb")
        shutil.copy2(row["dna_unbound_pdb"], workspace_dir / "processed/pdbs" / f"{dna_target[:4]}.pdb")
        template_id = row["template_id"]
        template_path = workspace_dir / "templates/pdbs" / f"{template_id[:4]}.pdb"
        shutil.copy2(row["template_complex_pdb"], template_path)
        checked_templates.append(template_id)
        chain_hints[template_id] = {
            "protein_chain_ids": parse_chain_list(row["template_protein_chain_ids"]),
            "dna_chain_ids": parse_chain_list(row["template_dna_chain_ids"]),
        }

    write_csv(workspace_dir / "inputs.csv", inputs_rows, fieldnames=["Receptor", "Ligand"])
    with open(workspace_dir / "templates/checked_templates.txt", "w") as handle:
        for template_id in checked_templates:
            handle.write(f"{template_id}\n")
    with open(workspace_dir / "templates/dna_ext_template_hints.json", "w") as handle:
        json.dump(chain_hints, handle, indent=2)
    return rows


def run_prism_main(prism_main_root, workspace_dir, profile):
    prism_main_root = Path(prism_main_root).resolve()
    workspace_dir = Path(workspace_dir).resolve()
    stdout_path = workspace_dir / "dna_run_logs" / f"stdout_{profile}.txt"
    stderr_path = workspace_dir / "dna_run_logs" / f"stderr_{profile}.txt"
    cmd = [
        sys.executable,
        str(prism_main_root / "prism.py"),
        "--generate_templates",
        "false",
        "--interaction-mode",
        "protein_dna_ext",
        "--aligner",
        "usalign",
        "--dna-ext-threshold-profile",
        profile,
    ]
    with open(stdout_path, "w") as stdout_handle, open(stderr_path, "w") as stderr_handle:
        proc = subprocess.run(cmd, cwd=str(workspace_dir), stdout=stdout_handle, stderr=stderr_handle)
    return proc.returncode, stdout_path, stderr_path


def collect_prediction_rows(workspace_dir):
    workspace_dir = Path(workspace_dir)
    refinement_csv = workspace_dir / "processed/rosetta_dna_ext_refinement/refinement_scores.csv"
    transformation_csv = workspace_dir / "processed/dna_ext_transformation/dna_ext_transformation_models.csv"
    if refinement_csv.exists():
        with open(refinement_csv, newline="") as handle:
            return list(csv.DictReader(handle)), "refinement"
    if transformation_csv.exists():
        with open(transformation_csv, newline="") as handle:
            return list(csv.DictReader(handle)), "transformation"
    return [], "missing"


def build_combined_model(workspace_dir, prediction_row):
    workspace_dir = Path(workspace_dir)
    combined_path = prediction_row.get("combined_model_path")
    if combined_path:
        path = Path(combined_path)
        if path.exists():
            return str(path)

    template_id = prediction_row["template_id"]
    protein_target = prediction_row["protein_target"]
    dna_target = prediction_row["dna_target"]
    transformed_protein = workspace_dir / "processed/dna_ext_transformation" / f"{template_id}_{protein_target}_{dna_target}_protein.pdb"
    dna_only = workspace_dir / "processed/dna_ext_features" / f"{dna_target}.dna_only.pdb"
    if not transformed_protein.exists() or not dna_only.exists():
        return None

    combined_dir = workspace_dir / "processed/dna_ext_transformation" / "combined_models"
    combined_dir.mkdir(parents=True, exist_ok=True)
    combined_path = combined_dir / f"{template_id}_{protein_target}_{dna_target}_combined.pdb"
    def non_terminal_lines(path):
        lines = []
        for line in path.read_text().splitlines():
            if line.startswith("END"):
                continue
            lines.append(line)
        return lines

    with open(combined_path, "w") as handle:
        handle.write("\n".join(non_terminal_lines(transformed_protein)))
        handle.write("\n")
        handle.write("\n".join(non_terminal_lines(dna_only)))
        handle.write("\nEND\n")
    return str(combined_path)


def score_predictions(manifest_rows, prediction_rows, workspace_dir):
    manifest_by_pair = {row["pair_id"]: row for row in manifest_rows}
    per_prediction = []
    grouped = {}
    for prediction in prediction_rows:
        pair_id = None
        for row in manifest_rows:
            protein_target, dna_target = case_to_target_ids(row["case_id"])
            if prediction.get("protein_target") == protein_target[:4] and prediction.get("dna_target") == dna_target[:4]:
                pair_id = row["pair_id"]
                break
        if pair_id is None:
            continue
        manifest_row = manifest_by_pair[pair_id]
        model_pdb = build_combined_model(workspace_dir, prediction)
        if model_pdb is None:
            continue
        score = score_one_ext(
            model_pdb,
            manifest_row["native_complex_pdb"],
            model_protein_chains=normalize_chain_field(prediction.get("protein_chain_ids")) or parse_chain_list(manifest_row["protein_chain_ids"]),
            model_dna_chains=normalize_chain_field(prediction.get("dna_chain_ids")) or parse_chain_list(manifest_row["dna_chain_ids"]),
            native_protein_chains=parse_chain_list(manifest_row["protein_chain_ids"]),
            native_dna_chains=parse_chain_list(manifest_row["dna_chain_ids"]),
            score_json=None,
        )
        score.update({"pair_id": pair_id, "template_id": manifest_row["template_id"], "case_id": manifest_row["case_id"]})
        score["status"] = prediction.get("status", score["status"])
        score["reason"] = prediction.get("reason", score["reason"])
        per_prediction.append(score)
        grouped.setdefault(pair_id, []).append(score)

    summary_rows = []
    for pair_id, rows in grouped.items():
        cp_vals = [float(r["contact_precision"]) for r in rows]
        cr_vals = [float(r["contact_recall"]) for r in rows]
        cf_vals = [float(r["contact_f1"]) for r in rows]
        drp_vals = [float(r["dna_register_contact_precision"]) for r in rows]
        drr_vals = [float(r["dna_register_contact_recall"]) for r in rows]
        drf_vals = [float(r["dna_register_contact_f1"]) for r in rows]
        pip_vals = [float(r["protein_interface_precision"]) for r in rows]
        pir_vals = [float(r["protein_interface_recall"]) for r in rows]
        pif_vals = [float(r["protein_interface_f1"]) for r in rows]
        ncp_vals = [float(r["nucleotide_contact_precision"]) for r in rows]
        ncr_vals = [float(r["nucleotide_contact_recall"]) for r in rows]
        ncf_vals = [float(r["nucleotide_contact_f1"]) for r in rows]
        tm_vals = [float(r["alignment_tm_score"]) for r in rows if r["alignment_tm_score"] is not None]
        rmsd_vals = [float(r["alignment_rmsd"]) for r in rows if r["alignment_rmsd"] is not None]
        cp_mean, cp_var, cp_best = summarize(cp_vals, "max")
        cr_mean, cr_var, cr_best = summarize(cr_vals, "max")
        cf_mean, cf_var, cf_best = summarize(cf_vals, "max")
        drp_mean, drp_var, drp_best = summarize(drp_vals, "max")
        drr_mean, drr_var, drr_best = summarize(drr_vals, "max")
        drf_mean, drf_var, drf_best = summarize(drf_vals, "max")
        pip_mean, pip_var, pip_best = summarize(pip_vals, "max")
        pir_mean, pir_var, pir_best = summarize(pir_vals, "max")
        pif_mean, pif_var, pif_best = summarize(pif_vals, "max")
        ncp_mean, ncp_var, ncp_best = summarize(ncp_vals, "max")
        ncr_mean, ncr_var, ncr_best = summarize(ncr_vals, "max")
        ncf_mean, ncf_var, ncf_best = summarize(ncf_vals, "max")
        tm_mean, tm_var, tm_best = summarize(tm_vals, "max")
        rmsd_mean, rmsd_var, rmsd_best = summarize(rmsd_vals, "min")
        summary_rows.append(
            {
                "pair_id": pair_id,
                "n_predictions_total": len(rows),
                "contact_precision_mean": cp_mean,
                "contact_precision_variance": cp_var,
                "contact_precision_best": cp_best,
                "contact_recall_mean": cr_mean,
                "contact_recall_variance": cr_var,
                "contact_recall_best": cr_best,
                "contact_f1_mean": cf_mean,
                "contact_f1_variance": cf_var,
                "contact_f1_best": cf_best,
                "dna_register_contact_precision_mean": drp_mean,
                "dna_register_contact_precision_variance": drp_var,
                "dna_register_contact_precision_best": drp_best,
                "dna_register_contact_recall_mean": drr_mean,
                "dna_register_contact_recall_variance": drr_var,
                "dna_register_contact_recall_best": drr_best,
                "dna_register_contact_f1_mean": drf_mean,
                "dna_register_contact_f1_variance": drf_var,
                "dna_register_contact_f1_best": drf_best,
                "protein_interface_precision_mean": pip_mean,
                "protein_interface_precision_variance": pip_var,
                "protein_interface_precision_best": pip_best,
                "protein_interface_recall_mean": pir_mean,
                "protein_interface_recall_variance": pir_var,
                "protein_interface_recall_best": pir_best,
                "protein_interface_f1_mean": pif_mean,
                "protein_interface_f1_variance": pif_var,
                "protein_interface_f1_best": pif_best,
                "nucleotide_contact_precision_mean": ncp_mean,
                "nucleotide_contact_precision_variance": ncp_var,
                "nucleotide_contact_precision_best": ncp_best,
                "nucleotide_contact_recall_mean": ncr_mean,
                "nucleotide_contact_recall_variance": ncr_var,
                "nucleotide_contact_recall_best": ncr_best,
                "nucleotide_contact_f1_mean": ncf_mean,
                "nucleotide_contact_f1_variance": ncf_var,
                "nucleotide_contact_f1_best": ncf_best,
                "alignment_tm_score_mean": tm_mean,
                "alignment_tm_score_variance": tm_var,
                "alignment_tm_score_best": tm_best,
                "alignment_rmsd_mean": rmsd_mean,
                "alignment_rmsd_variance": rmsd_var,
                "alignment_rmsd_best_min": rmsd_best,
            }
        )
    return per_prediction, summary_rows


def run_pipeline(manifest_csv, prism_main_root, output_root, profile):
    manifest_path = Path(manifest_csv).resolve()
    output_root = Path(output_root).resolve()
    workspace_dir = output_root / f"workspace_{profile}"
    results_dir = output_root / profile
    results_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = stage_workspace(manifest_path, workspace_dir)
    return_code, stdout_path, stderr_path = run_prism_main(prism_main_root, workspace_dir, profile)
    prediction_rows, prediction_source = collect_prediction_rows(workspace_dir)
    per_prediction, pair_summary = score_predictions(manifest_rows, prediction_rows, workspace_dir)

    write_csv(results_dir / "per_prediction.csv", per_prediction, fieldnames=[
        "pair_id", "template_id", "case_id", "model_pdb", "native_pdb", "model_protein_chains",
        "model_dna_chains", "native_protein_chains", "native_dna_chains", "contact_precision",
        "contact_recall", "contact_f1", "dna_register_contact_precision", "dna_register_contact_recall",
        "dna_register_contact_f1", "protein_interface_precision", "protein_interface_recall",
        "protein_interface_f1", "nucleotide_contact_precision", "nucleotide_contact_recall",
        "nucleotide_contact_f1", "alignment_tm_score", "alignment_rmsd", "rosetta_dna_interface_score",
        "status", "reason"
    ])
    write_csv(results_dir / "pair_summary.csv", pair_summary, fieldnames=[
        "pair_id", "n_predictions_total", "contact_precision_mean", "contact_precision_variance",
        "contact_precision_best", "contact_recall_mean", "contact_recall_variance", "contact_recall_best",
        "contact_f1_mean", "contact_f1_variance", "contact_f1_best",
        "dna_register_contact_precision_mean", "dna_register_contact_precision_variance",
        "dna_register_contact_precision_best", "dna_register_contact_recall_mean",
        "dna_register_contact_recall_variance", "dna_register_contact_recall_best",
        "dna_register_contact_f1_mean", "dna_register_contact_f1_variance",
        "dna_register_contact_f1_best",
        "protein_interface_precision_mean", "protein_interface_precision_variance", "protein_interface_precision_best",
        "protein_interface_recall_mean", "protein_interface_recall_variance", "protein_interface_recall_best",
        "protein_interface_f1_mean", "protein_interface_f1_variance", "protein_interface_f1_best",
        "nucleotide_contact_precision_mean", "nucleotide_contact_precision_variance", "nucleotide_contact_precision_best",
        "nucleotide_contact_recall_mean", "nucleotide_contact_recall_variance", "nucleotide_contact_recall_best",
        "nucleotide_contact_f1_mean", "nucleotide_contact_f1_variance", "nucleotide_contact_f1_best",
        "alignment_tm_score_mean", "alignment_tm_score_variance", "alignment_tm_score_best",
        "alignment_rmsd_mean", "alignment_rmsd_variance", "alignment_rmsd_best_min"
    ])
    summary = {
        "profile": profile,
        "return_code": return_code,
        "prediction_source": prediction_source,
        "prediction_count": len(per_prediction),
        "workspace_dir": str(workspace_dir),
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
    }
    with open(results_dir / "run_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main():
    parser = argparse.ArgumentParser(description="Run the Dockground extension benchmark workflow")
    parser.add_argument("--manifest", default="benchmark/data/protein_dna_dockground_subset_manifest.csv")
    parser.add_argument("--prism-main-root", default="/scratch/rshadi25/GitHub/PRISM-main")
    parser.add_argument("--output-root", default="benchmark/prism_processed/results/protein_dna_dockground_subset")
    parser.add_argument("--profile", choices=["strict", "dockground_calibration"], default="strict")
    args = parser.parse_args()
    summary = run_pipeline(args.manifest, args.prism_main_root, args.output_root, args.profile)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
