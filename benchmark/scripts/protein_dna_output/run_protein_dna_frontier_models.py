#!/usr/bin/env python3
"""
Run frontier protein-DNA models on the Dockground subset and score outputs.
"""

import argparse
import csv
import json
import os
import statistics
import subprocess
import sys
import shutil
import time
from pathlib import Path

try:
    from benchmark.scripts.protein_dna_output.frontier_model_adapters import (
        candidate_model_files,
        alphafold3_supported_gpu,
        load_registry,
        boltz_runtime_env,
        normalize_prediction_file,
        probe_tool,
        write_alphafold3_input,
        write_boltz_input,
        write_chai_input,
        write_csv,
        write_rosettafold_all_atom_input,
        write_rosettafoldna_inputs,
    )
    from benchmark.scripts.protein_dna_output.run_protein_dna_dockground_ext import stage_workspace
    from benchmark.scripts.protein_dna_output.score_single_protein_dna_pair_ext import parse_chain_list, score_one_ext
except ImportError:
    from frontier_model_adapters import (  # type: ignore
        candidate_model_files,
        alphafold3_supported_gpu,
        load_registry,
        boltz_runtime_env,
        normalize_prediction_file,
        probe_tool,
        write_alphafold3_input,
        write_boltz_input,
        write_chai_input,
        write_csv,
        write_rosettafold_all_atom_input,
        write_rosettafoldna_inputs,
    )
    from run_protein_dna_dockground_ext import stage_workspace  # type: ignore
    from score_single_protein_dna_pair_ext import parse_chain_list, score_one_ext  # type: ignore


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


def load_manifest(manifest_path):
    with open(manifest_path, newline="") as handle:
        return list(csv.DictReader(handle))


def build_common_inputs(tool_id, workspace_dir, manifest_rows):
    workspace_dir = Path(workspace_dir)
    staged_workspace = workspace_dir / "staged"
    staged_workspace.mkdir(parents=True, exist_ok=True)
    stage_workspace(manifest_csv=workspace_dir / "_manifest.csv", workspace_dir=staged_workspace)


def write_manifest_copy(manifest_rows, workspace_dir):
    path = Path(workspace_dir) / "_manifest.csv"
    if manifest_rows:
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0].keys()))
            writer.writeheader()
            writer.writerows(manifest_rows)
    return path


def tool_input_builder(tool_id):
    if tool_id == "chai1":
        return write_chai_input
    if tool_id == "boltz2":
        return write_boltz_input
    if tool_id == "alphafold3":
        return write_alphafold3_input
    if tool_id == "rosettafoldna":
        return write_rosettafoldna_inputs
    if tool_id == "rosettafold_all_atom":
        return write_rosettafold_all_atom_input
    raise KeyError(f"Unsupported frontier tool: {tool_id}")


def build_command(tool, artifacts, workspace_dir, output_dir, python_executable):
    tool_id = tool["tool_id"]
    if tool_id == "chai1":
        binary = artifacts["binary_path"]
        return [binary, "fold", artifacts["input_path"], str(output_dir)], workspace_dir, None
    if tool_id == "boltz2":
        binary = artifacts["binary_path"]
        return [binary, "predict", artifacts["input_path"], "--use_msa_server"], workspace_dir, boltz_runtime_env(binary)
    if tool_id == "alphafold3":
        runner = artifacts["binary_path"] or os.environ.get(tool.get("runner_env", ""), "")
        if runner and runner.endswith("/alphafold"):
            sif_runner = str(Path(runner).with_name("alphafold3.sif"))
            return [
                "singularity",
                "exec",
                "--nv",
                sif_runner,
                "python3",
                "/app/alphafold/run_alphafold.py",
                "--json_path",
                artifacts["input_path"],
                "--output_dir",
                str(output_dir),
            ], workspace_dir, None
        if runner and runner.endswith(".sif"):
            cmd = [
                "singularity",
                "exec",
                "--nv",
                runner,
                "python3",
                "/app/alphafold/run_alphafold.py",
                "--json_path",
                artifacts["input_path"],
                "--output_dir",
                str(output_dir),
            ]
            return cmd, workspace_dir, None
        runner = runner or "run_alphafold.py"
        return [runner, "--json_path", artifacts["input_path"], "--output_dir", str(output_dir)], output_dir, None
    if tool_id == "rosettafoldna":
        runner = artifacts["binary_path"] or os.environ.get(tool.get("runner_env", ""), "") or "run_RF2NA.sh"
        protein_fastas = artifacts.get("protein_fastas", [])
        dna_fastas = artifacts.get("dna_fastas", [])
        if not protein_fastas or not dna_fastas:
            return None, workspace_dir, None
        first_dna = dna_fastas[0]
        args = [runner, str(output_dir)]
        args.extend([f"P:{path}" for path in protein_fastas])
        args.append(f"D:{first_dna}")
        return args, workspace_dir, None
    if tool_id == "rosettafold_all_atom":
        return [python_executable, "-m", "rf2aa.run_inference", "-cd", artifacts["input_dir"], "--config-name", "nucleic_acid"], output_dir, None
    raise KeyError(f"Unsupported frontier tool: {tool_id}")


def run_command(cmd, cwd, stdout_path, stderr_path, env=None):
    start = time.time()
    with open(stdout_path, "w") as stdout_handle, open(stderr_path, "w") as stderr_handle:
        proc = subprocess.run(cmd, cwd=str(cwd), stdout=stdout_handle, stderr=stderr_handle, check=False, env=env)
    return proc.returncode, time.time() - start


def score_models(manifest_rows, model_rows, native_lookup, output_dir):
    per_prediction = []
    grouped = {}
    for model_row in model_rows:
        pair_id = model_row["pair_id"]
        manifest_row = native_lookup[pair_id]
        score = {
            "model_pdb": model_row.get("normalized_model_path", ""),
            "native_pdb": manifest_row["native_complex_pdb"],
            "model_protein_chains": model_row.get("model_protein_chains", ""),
            "model_dna_chains": model_row.get("model_dna_chains", ""),
            "native_protein_chains": manifest_row["protein_chain_ids"],
            "native_dna_chains": manifest_row["dna_chain_ids"],
            "contact_precision": None,
            "contact_recall": None,
            "contact_f1": None,
            "dna_register_contact_precision": None,
            "dna_register_contact_recall": None,
            "dna_register_contact_f1": None,
            "protein_interface_precision": None,
            "protein_interface_recall": None,
            "protein_interface_f1": None,
            "nucleotide_contact_precision": None,
            "nucleotide_contact_recall": None,
            "nucleotide_contact_f1": None,
            "alignment_tm_score": None,
            "alignment_rmsd": None,
            "rosetta_dna_interface_score": None,
            "runtime_seconds": model_row.get("runtime_seconds"),
            "status": model_row["status"],
            "reason": model_row["reason"],
        }
        if model_row["status"] == "passed" and model_row.get("normalized_model_path"):
            score.update(
                score_one_ext(
                    model_row["normalized_model_path"],
                    manifest_row["native_complex_pdb"],
                    model_protein_chains=parse_chain_list(model_row.get("model_protein_chains")) or parse_chain_list(manifest_row["protein_chain_ids"]),
                    model_dna_chains=parse_chain_list(model_row.get("model_dna_chains")) or parse_chain_list(manifest_row["dna_chain_ids"]),
                    native_protein_chains=parse_chain_list(manifest_row["protein_chain_ids"]),
                    native_dna_chains=parse_chain_list(manifest_row["dna_chain_ids"]),
                    score_json=None,
                )
            )
        score.update(
            {
                "pair_id": pair_id,
                "case_id": manifest_row["case_id"],
                "template_id": manifest_row["template_id"],
                "tool_id": model_row["tool_id"],
                "tool_name": model_row["tool_name"],
                "input_path": model_row["input_path"],
                "raw_model_path": model_row["raw_model_path"],
                "normalized_model_path": model_row["normalized_model_path"],
                "runtime_seconds": model_row.get("runtime_seconds"),
            }
        )
        per_prediction.append(score)
        grouped.setdefault((pair_id, model_row["tool_id"]), []).append(score)

    summary_rows = []
    for (pair_id, tool_id), rows in grouped.items():
        contact_precision_vals = [to_float(row["contact_precision"]) for row in rows if to_float(row["contact_precision"]) is not None]
        contact_recall_vals = [to_float(row["contact_recall"]) for row in rows if to_float(row["contact_recall"]) is not None]
        contact_f1_vals = [to_float(row["contact_f1"]) for row in rows if to_float(row["contact_f1"]) is not None]
        register_precision_vals = [to_float(row["dna_register_contact_precision"]) for row in rows if to_float(row["dna_register_contact_precision"]) is not None]
        register_recall_vals = [to_float(row["dna_register_contact_recall"]) for row in rows if to_float(row["dna_register_contact_recall"]) is not None]
        register_f1_vals = [to_float(row["dna_register_contact_f1"]) for row in rows if to_float(row["dna_register_contact_f1"]) is not None]
        protein_recall_vals = [to_float(row["protein_interface_recall"]) for row in rows if to_float(row["protein_interface_recall"]) is not None]
        nucleotide_recall_vals = [to_float(row["nucleotide_contact_recall"]) for row in rows if to_float(row["nucleotide_contact_recall"]) is not None]
        tm_vals = [to_float(row["alignment_tm_score"]) for row in rows if to_float(row["alignment_tm_score"]) is not None]
        rmsd_vals = [to_float(row["alignment_rmsd"]) for row in rows if to_float(row["alignment_rmsd"]) is not None]
        runtime_vals = [to_float(row["runtime_seconds"]) for row in rows if to_float(row["runtime_seconds"]) is not None]

        cp_mean, cp_var, cp_best = summarize(contact_precision_vals, "max")
        cr_mean, cr_var, cr_best = summarize(contact_recall_vals, "max")
        cf_mean, cf_var, cf_best = summarize(contact_f1_vals, "max")
        rp_mean, rp_var, rp_best = summarize(register_precision_vals, "max")
        rr_mean, rr_var, rr_best = summarize(register_recall_vals, "max")
        rf_mean, rf_var, rf_best = summarize(register_f1_vals, "max")
        pir_mean, pir_var, pir_best = summarize(protein_recall_vals, "max")
        ncr_mean, ncr_var, ncr_best = summarize(nucleotide_recall_vals, "max")
        tm_mean, tm_var, tm_best = summarize(tm_vals, "max")
        rmsd_mean, rmsd_var, rmsd_best = summarize(rmsd_vals, "min")
        runtime_mean, runtime_var, runtime_best = summarize(runtime_vals, "min")
        best_tm_row = max((row for row in rows if to_float(row["alignment_tm_score"]) is not None), key=lambda row: row["alignment_tm_score"], default=None)
        summary_rows.append(
            {
                "pair_id": pair_id,
                "tool_id": tool_id,
                "n_predictions_total": len(rows),
                "n_scored": len(contact_precision_vals),
                "contact_precision_mean": cp_mean,
                "contact_precision_variance": cp_var,
                "contact_precision_best": cp_best,
                "contact_recall_mean": cr_mean,
                "contact_recall_variance": cr_var,
                "contact_recall_best": cr_best,
                "contact_f1_mean": cf_mean,
                "contact_f1_variance": cf_var,
                "contact_f1_best": cf_best,
                "dna_register_contact_precision_mean": rp_mean,
                "dna_register_contact_precision_variance": rp_var,
                "dna_register_contact_precision_best": rp_best,
                "dna_register_contact_recall_mean": rr_mean,
                "dna_register_contact_recall_variance": rr_var,
                "dna_register_contact_recall_best": rr_best,
                "dna_register_contact_f1_mean": rf_mean,
                "dna_register_contact_f1_variance": rf_var,
                "dna_register_contact_f1_best": rf_best,
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
                "runtime_seconds_mean": runtime_mean,
                "runtime_seconds_variance": runtime_var,
                "runtime_seconds_best_min": runtime_best,
                "best_tm_model": best_tm_row.get("model_pdb") if best_tm_row else None,
            }
        )

    write_csv(output_dir / "per_prediction.csv", per_prediction)
    write_csv(output_dir / "pair_summary.csv", summary_rows)
    return per_prediction, summary_rows


def run_frontier(manifest_csv, registry_path, output_root, work_root, tools=None, execute=True, python_executable=sys.executable):
    manifest_rows = load_manifest(manifest_csv)
    registry = load_registry(registry_path)
    if tools:
        selected = [tool for tool in registry if tool["tool_id"] in tools]
    else:
        selected = registry

    output_root = Path(output_root).resolve()
    work_root = Path(work_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    native_lookup = {row["pair_id"]: row for row in manifest_rows}
    tool_capabilities = []
    tool_summaries = []

    for tool in selected:
        tool_workspace = work_root / tool["tool_id"]
        if tool_workspace.exists():
            shutil.rmtree(tool_workspace)
        tool_workspace.mkdir(parents=True, exist_ok=True)
        tool_output = output_root / tool["tool_id"]
        tool_output.mkdir(parents=True, exist_ok=True)

        manifest_copy = write_manifest_copy(manifest_rows, tool_workspace)
        stage_workspace(manifest_csv=manifest_copy, workspace_dir=tool_workspace / "staged")

        availability = probe_tool(tool, python_executable=python_executable)
        tool_capabilities.append(
            {
                "tool_id": tool["tool_id"],
                "tool_name": tool["name"],
                "integration_kind": tool["integration_kind"],
                "available": availability.available,
                "binary_path": availability.binary_path,
                "module_name": availability.module_name,
                "status": "available" if availability.available else "skipped",
                "reason": availability.reason,
                "official_url": tool.get("official_url", ""),
                "source_url": tool.get("source_url", ""),
                "notes": tool.get("notes", ""),
            }
        )

        if tool["tool_id"] == "alphafold3":
            supported, gpu_label = alphafold3_supported_gpu()
            if not supported:
                tool_capabilities[-1]["available"] = False
                tool_capabilities[-1]["status"] = "skipped"
                tool_capabilities[-1]["reason"] = gpu_label
                with open(tool_output / "run_summary.json", "w") as handle:
                    json.dump(
                        {
                            "tool_id": tool["tool_id"],
                            "tool_name": tool["name"],
                            "status": "skipped",
                            "reason": gpu_label,
                            "manifest": str(Path(manifest_csv).resolve()),
                            "registry": str(Path(registry_path).resolve()),
                        },
                        handle,
                        indent=2,
                    )
                tool_summaries.append({"tool_id": tool["tool_id"], "status": "skipped", "reason": gpu_label})
                continue

        if not availability.available or not execute:
            with open(tool_output / "run_summary.json", "w") as handle:
                json.dump(
                    {
                        "tool_id": tool["tool_id"],
                        "tool_name": tool["name"],
                        "status": "skipped" if not availability.available else "planned",
                        "reason": availability.reason if not availability.available else "dry_run",
                        "manifest": str(Path(manifest_csv).resolve()),
                        "registry": str(Path(registry_path).resolve()),
                    },
                    handle,
                    indent=2,
                )
            tool_summaries.append({"tool_id": tool["tool_id"], "status": "skipped"})
            continue

        input_builder = tool_input_builder(tool["tool_id"])
        input_artifacts = [input_builder(tool_workspace, row) for row in manifest_rows]
        model_rows = []

        for manifest_row, artifacts in zip(manifest_rows, input_artifacts):
            pair_id = manifest_row["pair_id"]
            run_dir = tool_workspace / "runs" / pair_id
            run_dir.mkdir(parents=True, exist_ok=True)
            output_dir = run_dir / "output"
            output_dir.mkdir(parents=True, exist_ok=True)
            cmd, cwd, run_env = build_command(tool, artifacts | {"binary_path": availability.binary_path}, tool_workspace, output_dir, python_executable)
            if cmd is None:
                model_rows.append(
                    {
                        "pair_id": pair_id,
                        "tool_id": tool["tool_id"],
                        "tool_name": tool["name"],
                        "input_path": artifacts.get("input_path", ""),
                        "raw_model_path": "",
                        "normalized_model_path": "",
                        "model_protein_chains": ",".join(parse_chain_list(manifest_row["protein_chain_ids"])),
                        "model_dna_chains": ",".join(parse_chain_list(manifest_row["dna_chain_ids"])),
                        "runtime_seconds": None,
                        "status": "skipped",
                        "reason": "missing_input_artifacts",
                    }
                )
                continue
            stdout_path = run_dir / "stdout.txt"
            stderr_path = run_dir / "stderr.txt"
            returncode, elapsed = run_command(cmd, cwd, stdout_path, stderr_path, env=run_env)
            if returncode != 0:
                model_rows.append(
                    {
                        "pair_id": pair_id,
                        "tool_id": tool["tool_id"],
                        "tool_name": tool["name"],
                        "input_path": artifacts.get("input_path", ""),
                        "raw_model_path": "",
                        "normalized_model_path": "",
                        "model_protein_chains": ",".join(parse_chain_list(manifest_row["protein_chain_ids"])),
                        "model_dna_chains": ",".join(parse_chain_list(manifest_row["dna_chain_ids"])),
                        "runtime_seconds": elapsed,
                        "status": "failed",
                        "reason": f"command_failed:{returncode}",
                    }
                )
                continue
            candidates = candidate_model_files(output_dir)
            if not candidates and tool["tool_id"] == "boltz2":
                boltz_workspace = tool_workspace / f"boltz_results_{pair_id}"
                candidates = candidate_model_files(boltz_workspace)
            if not candidates:
                model_rows.append(
                    {
                        "pair_id": pair_id,
                        "tool_id": tool["tool_id"],
                        "tool_name": tool["name"],
                        "input_path": artifacts.get("input_path", ""),
                        "raw_model_path": "",
                        "normalized_model_path": "",
                        "model_protein_chains": ",".join(parse_chain_list(manifest_row["protein_chain_ids"])),
                        "model_dna_chains": ",".join(parse_chain_list(manifest_row["dna_chain_ids"])),
                        "runtime_seconds": elapsed,
                        "status": "failed",
                        "reason": "no_prediction_files_found",
                    }
                )
                continue
            normalized_dir = run_dir / "normalized"
            normalized_model = normalize_prediction_file(candidates[0], normalized_dir, manifest_row)
            model_rows.append(
                {
                    "pair_id": pair_id,
                    "tool_id": tool["tool_id"],
                    "tool_name": tool["name"],
                    "input_path": artifacts.get("input_path", ""),
                    "raw_model_path": str(candidates[0].resolve()),
                    "normalized_model_path": str(Path(normalized_model).resolve()),
                    "model_protein_chains": ",".join(parse_chain_list(manifest_row["protein_chain_ids"])),
                    "model_dna_chains": ",".join(parse_chain_list(manifest_row["dna_chain_ids"])),
                    "runtime_seconds": elapsed,
                    "status": "passed",
                    "reason": "",
                }
            )

        per_prediction, summary_rows = score_models(manifest_rows, model_rows, native_lookup, tool_output)
        summary = {
            "tool_id": tool["tool_id"],
            "tool_name": tool["name"],
            "manifest": str(Path(manifest_csv).resolve()),
            "registry": str(Path(registry_path).resolve()),
            "status": "completed" if per_prediction else "skipped",
            "available": availability.available,
            "n_predictions": len(per_prediction),
            "n_summaries": len(summary_rows),
            "output_root": str(tool_output),
        }
        with open(tool_output / "run_summary.json", "w") as handle:
            json.dump(summary, handle, indent=2)
        tool_summaries.append(summary)

    write_csv(
        output_root / "tool_capability_report.csv",
        tool_capabilities,
        fieldnames=[
            "tool_id",
            "tool_name",
            "integration_kind",
            "available",
            "binary_path",
            "module_name",
            "status",
            "reason",
            "official_url",
            "source_url",
            "notes",
        ],
    )
    with open(output_root / "run_summary.json", "w") as handle:
        json.dump(
            {
                "manifest": str(Path(manifest_csv).resolve()),
                "registry": str(Path(registry_path).resolve()),
                "tool_count": len(selected),
                "available_tools": [row["tool_id"] for row in tool_capabilities if row["available"]],
                "skipped_tools": [row["tool_id"] for row in tool_capabilities if not row["available"]],
                "output_root": str(output_root),
                "work_root": str(work_root),
                "tool_summaries": tool_summaries,
            },
            handle,
            indent=2,
        )
    return tool_summaries


def main():
    parser = argparse.ArgumentParser(description="Run frontier protein-DNA models on the Dockground subset")
    parser.add_argument("--manifest", default="benchmark/data/protein_dna_dockground_subset_manifest.csv")
    parser.add_argument("--registry", default="benchmark/data/protein_dna_frontier_tools.json")
    parser.add_argument("--output-root", default="benchmark/prism_processed/results/protein_dna_frontier_models")
    parser.add_argument("--work-root", default="tmp/agent/frontier_models")
    parser.add_argument("--tools", default="", help="Comma-separated tool ids to run; default is all registry entries")
    parser.add_argument("--no-execute", action="store_true", help="Only probe and stage inputs, do not run tools")
    args = parser.parse_args()

    tool_ids = [tool.strip() for tool in args.tools.split(",") if tool.strip()] or None
    run_frontier(
        args.manifest,
        args.registry,
        args.output_root,
        args.work_root,
        tools=tool_ids,
        execute=not args.no_execute,
        python_executable=sys.executable,
    )


if __name__ == "__main__":
    main()
