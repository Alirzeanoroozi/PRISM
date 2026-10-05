#!/usr/bin/env python3
"""
Stage-by-stage validation for the frontier protein-DNA benchmark pipeline.

This script divides the frontier flow into explicit testable stages:
- manifest: validate manifest and registry shape
- workspace: stage the Dockground workspace and input files
- probe: record tool availability and AF3 GPU gating
- dry-run: exercise the shared runner without executing tools
- execute: run the shared frontier harness end to end
- af3-gpu: probe the AlphaFold 3 GPU/container path directly
"""

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

try:
    from benchmark.scripts.protein_dna_output.frontier_model_adapters import (
        alphafold3_supported_gpu,
        alphafold3_database_root,
        alphafold3_model_root,
        load_registry,
        probe_tool,
    )
    from benchmark.scripts.protein_dna_output.run_protein_dna_dockground_ext import stage_workspace
    from benchmark.scripts.protein_dna_output.run_protein_dna_frontier_models import (
        load_manifest,
        run_frontier,
        write_manifest_copy,
    )
except ImportError:
    from .frontier_model_adapters import (  # type: ignore
        alphafold3_supported_gpu,
        alphafold3_database_root,
        alphafold3_model_root,
        load_registry,
        probe_tool,
    )
    from .run_protein_dna_dockground_ext import stage_workspace  # type: ignore
    from .run_protein_dna_frontier_models import (  # type: ignore
        load_manifest,
        run_frontier,
        write_manifest_copy,
    )


REQUIRED_MANIFEST_FIELDS = [
    "pair_id",
    "case_id",
    "template_id",
    "protein_unbound_pdb",
    "dna_unbound_pdb",
    "native_complex_pdb",
    "template_complex_pdb",
    "protein_chain_ids",
    "dna_chain_ids",
    "template_protein_chain_ids",
    "template_dna_chain_ids",
]

REQUIRED_REGISTRY_FIELDS = [
    "tool_id",
    "name",
    "integration_kind",
]

STAGE_ALIASES = {
    "manifest": "manifest",
    "workspace": "workspace",
    "stage": "workspace",
    "probe": "probe",
    "dry-run": "dry-run",
    "dry_run": "dry-run",
    "execute": "execute",
    "af3-gpu": "af3-gpu",
    "af3_gpu": "af3-gpu",
}

DEFAULT_STAGES = ["manifest", "workspace", "probe", "dry-run"]


def write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)


def write_csv_rows(path, rows, fieldnames=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def normalize_stage_name(stage):
    stage = stage.strip().lower().replace("_", "-")
    if stage not in STAGE_ALIASES:
        raise KeyError(f"Unsupported stage: {stage}")
    return STAGE_ALIASES[stage]


def parse_stage_list(value):
    if not value:
        return list(DEFAULT_STAGES)
    stages = []
    for raw_stage in value.split(","):
        raw_stage = raw_stage.strip()
        if not raw_stage:
            continue
        stages.append(normalize_stage_name(raw_stage))
    ordered = []
    for stage in stages:
        if stage not in ordered:
            ordered.append(stage)
    return ordered


def validate_manifest_rows(manifest_rows, registry_rows):
    manifest_missing = sorted({field for field in REQUIRED_MANIFEST_FIELDS if any(field not in row or row[field] == "" for row in manifest_rows)})
    registry_missing = sorted({field for field in REQUIRED_REGISTRY_FIELDS if any(field not in row or row[field] == "" for row in registry_rows)})
    manifest_pair_ids = [row.get("pair_id", "") for row in manifest_rows]
    registry_tool_ids = [row.get("tool_id", "") for row in registry_rows]
    return {
        "manifest_row_count": len(manifest_rows),
        "registry_row_count": len(registry_rows),
        "manifest_pair_ids": manifest_pair_ids,
        "registry_tool_ids": registry_tool_ids,
        "missing_manifest_fields": manifest_missing,
        "missing_registry_fields": registry_missing,
        "duplicate_manifest_pair_ids": sorted({pair_id for pair_id in manifest_pair_ids if pair_id and manifest_pair_ids.count(pair_id) > 1}),
        "duplicate_registry_tool_ids": sorted({tool_id for tool_id in registry_tool_ids if tool_id and registry_tool_ids.count(tool_id) > 1}),
    }


def run_manifest_stage(manifest_csv, registry_path, output_root):
    output_dir = Path(output_root) / "stages" / "manifest"
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = load_manifest(manifest_csv)
    registry_rows = load_registry(registry_path)
    summary = validate_manifest_rows(manifest_rows, registry_rows)
    manifest_copy = output_dir / "manifest_snapshot.csv"
    registry_copy = output_dir / "registry_snapshot.json"
    write_csv_rows(manifest_copy, manifest_rows)
    write_json(registry_copy, registry_rows)
    summary.update(
        {
            "manifest": str(Path(manifest_csv).resolve()),
            "registry": str(Path(registry_path).resolve()),
            "manifest_snapshot": str(manifest_copy.resolve()),
            "registry_snapshot": str(registry_copy.resolve()),
        }
    )
    write_json(output_dir / "manifest_summary.json", summary)
    return summary


def run_workspace_stage(manifest_csv, output_root, work_root):
    output_dir = Path(output_root) / "stages" / "workspace"
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(work_root) / "workspace"
    work_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = load_manifest(manifest_csv)
    manifest_copy = write_manifest_copy(manifest_rows, work_dir)
    staged_rows = stage_workspace(manifest_csv=manifest_copy, workspace_dir=work_dir / "staged")
    summary = {
        "manifest": str(Path(manifest_csv).resolve()),
        "workspace_root": str(work_dir.resolve()),
        "manifest_copy": str(Path(manifest_copy).resolve()),
        "inputs_csv": str((work_dir / "staged" / "inputs.csv").resolve()),
        "checked_templates": str((work_dir / "staged" / "templates" / "checked_templates.txt").resolve()),
        "template_hints": str((work_dir / "staged" / "templates" / "dna_ext_template_hints.json").resolve()),
        "pair_count": len(staged_rows),
    }
    write_json(output_dir / "workspace_summary.json", summary)
    return summary


def run_probe_stage(manifest_csv, registry_path, output_root, tools=None, python_executable=sys.executable):
    output_dir = Path(output_root) / "stages" / "probe"
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = load_manifest(manifest_csv)
    registry_rows = load_registry(registry_path)
    if tools:
        registry_rows = [row for row in registry_rows if row["tool_id"] in tools]

    tool_rows = []
    for tool in registry_rows:
        availability = probe_tool(tool, python_executable=python_executable)
        row = {
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
        if tool["tool_id"] == "alphafold3":
            supported, gpu_label = alphafold3_supported_gpu()
            db_root = str(alphafold3_database_root(tool))
            model_root = str(alphafold3_model_root(tool))
            if not supported:
                row["available"] = False
                row["status"] = "skipped"
                row["reason"] = gpu_label
            row["database_root"] = db_root
            row["database_root_env"] = tool.get("database_root_env", "")
            row["database_root_default"] = tool.get("database_root_default", "")
            row["model_root"] = model_root
            row["model_root_env"] = tool.get("model_root_env", "")
            row["model_root_default"] = tool.get("model_root_default", "")
        tool_rows.append(row)

    write_csv_rows(
        output_dir / "tool_capability_report.csv",
        tool_rows,
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
            "database_root_env",
            "database_root_default",
            "database_root",
            "model_root_env",
            "model_root_default",
            "model_root",
        ],
    )
    deeppbs_rows = [row for row in tool_rows if row["tool_id"] == "deeppbs"]
    if deeppbs_rows:
        write_csv_rows(
            output_dir / "deeppbs_probe_summary.tsv",
            deeppbs_rows,
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
    summary = {
        "manifest": str(Path(manifest_csv).resolve()),
        "registry": str(Path(registry_path).resolve()),
        "pair_count": len(manifest_rows),
        "tool_count": len(tool_rows),
        "available_tools": [row["tool_id"] for row in tool_rows if row["available"]],
        "skipped_tools": [row["tool_id"] for row in tool_rows if not row["available"]],
        "tool_capability_report": str((output_dir / "tool_capability_report.csv").resolve()),
        "deeppbs_probe_summary": str((output_dir / "deeppbs_probe_summary.tsv").resolve()) if deeppbs_rows else "",
    }
    write_json(output_dir / "probe_summary.json", summary)
    return summary


def run_dry_run_stage(manifest_csv, registry_path, output_root, work_root, tools=None, python_executable=sys.executable):
    dry_output = Path(output_root) / "stages" / "dry-run"
    dry_work = Path(work_root) / "dry-run"
    dry_output.mkdir(parents=True, exist_ok=True)
    dry_work.mkdir(parents=True, exist_ok=True)
    tool_summaries = run_frontier(
        manifest_csv,
        registry_path,
        dry_output,
        dry_work,
        tools=tools,
        execute=False,
        python_executable=python_executable,
    )
    summary = {
        "manifest": str(Path(manifest_csv).resolve()),
        "registry": str(Path(registry_path).resolve()),
        "output_root": str(dry_output.resolve()),
        "work_root": str(dry_work.resolve()),
        "tool_summaries": tool_summaries,
    }
    write_json(dry_output / "stage_summary.json", summary)
    return summary


def run_execute_stage(manifest_csv, registry_path, output_root, work_root, tools=None, python_executable=sys.executable):
    exec_output = Path(output_root) / "stages" / "execute"
    exec_work = Path(work_root) / "execute"
    exec_output.mkdir(parents=True, exist_ok=True)
    exec_work.mkdir(parents=True, exist_ok=True)
    tool_summaries = run_frontier(
        manifest_csv,
        registry_path,
        exec_output,
        exec_work,
        tools=tools,
        execute=True,
        python_executable=python_executable,
    )
    summary = {
        "manifest": str(Path(manifest_csv).resolve()),
        "registry": str(Path(registry_path).resolve()),
        "output_root": str(exec_output.resolve()),
        "work_root": str(exec_work.resolve()),
        "tool_summaries": tool_summaries,
    }
    write_json(exec_output / "stage_summary.json", summary)
    return summary


def run_af3_gpu_stage(output_root):
    output_dir = Path(output_root) / "stages" / "af3-gpu"
    output_dir.mkdir(parents=True, exist_ok=True)
    supported, gpu_label = alphafold3_supported_gpu()
    report = {
        "host_supported": supported,
        "reason": gpu_label,
    }
    if supported:
        sif_path = Path("/opt/ohpc/pub/apps/alphafold/3.0.1/alphafold3.sif")
        report["container_path"] = str(sif_path)
        if sif_path.exists():
            cmd = [
                "singularity",
                "exec",
                "--nv",
                str(sif_path),
                "python3",
                "-c",
                "import jax; print('jax devices:', jax.devices()); print('jax local gpu:', jax.local_devices(backend='gpu'))",
            ]
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, universal_newlines=True)
            report["container_returncode"] = proc.returncode
            report["container_stdout"] = proc.stdout.strip()
            report["container_stderr"] = proc.stderr.strip()
        else:
            report["container_returncode"] = None
            report["container_stdout"] = ""
            report["container_stderr"] = f"missing_container:{sif_path}"
    write_json(output_dir / "af3_gpu_report.json", report)
    return report


def run_staged_validation(manifest_csv, registry_path, output_root, work_root, tools=None, stages=None, python_executable=sys.executable):
    manifest_csv = Path(manifest_csv)
    registry_path = Path(registry_path)
    output_root = Path(output_root)
    work_root = Path(work_root)
    output_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)
    stages = stages or list(DEFAULT_STAGES)
    results = {}

    if "manifest" in stages:
        results["manifest"] = run_manifest_stage(manifest_csv, registry_path, output_root)
    if "workspace" in stages:
        results["workspace"] = run_workspace_stage(manifest_csv, output_root, work_root)
    if "probe" in stages:
        results["probe"] = run_probe_stage(manifest_csv, registry_path, output_root, tools=tools, python_executable=python_executable)
    if "dry-run" in stages:
        results["dry-run"] = run_dry_run_stage(manifest_csv, registry_path, output_root, work_root, tools=tools, python_executable=python_executable)
    if "execute" in stages:
        results["execute"] = run_execute_stage(manifest_csv, registry_path, output_root, work_root, tools=tools, python_executable=python_executable)
    if "af3-gpu" in stages:
        results["af3-gpu"] = run_af3_gpu_stage(output_root)

    write_json(output_root / "staged_validation_summary.json", results)
    return results


def main():
    parser = argparse.ArgumentParser(description="Validate the frontier protein-DNA pipeline stage by stage")
    parser.add_argument("--manifest", default="benchmark/data/protein_dna_dockground_subset_manifest.csv")
    parser.add_argument("--registry", default="benchmark/data/protein_dna_frontier_tools.json")
    parser.add_argument("--output-root", default="benchmark/prism_processed/results/protein_dna_frontier_validation")
    parser.add_argument("--work-root", default="tmp/agent/frontier_validation")
    parser.add_argument("--tools", default="", help="Comma-separated tool ids to include")
    parser.add_argument(
        "--af3-db-root",
        default="",
        help="Optional AlphaFold 3 database root override; sets PRISM_AF3_DB_DIR for the validation run.",
    )
    parser.add_argument(
        "--af3-model-root",
        default="",
        help="Optional AlphaFold 3 model-parameter root override; sets PRISM_AF3_MODEL_DIR for the validation run.",
    )
    parser.add_argument(
        "--stages",
        default=",".join(DEFAULT_STAGES),
        help="Comma-separated stages to run: manifest,workspace,probe,dry-run,execute,af3-gpu",
    )
    args = parser.parse_args()

    tool_ids = [tool.strip() for tool in args.tools.split(",") if tool.strip()] or None
    if args.af3_db_root:
        os.environ["PRISM_AF3_DB_DIR"] = args.af3_db_root
    if args.af3_model_root:
        os.environ["PRISM_AF3_MODEL_DIR"] = args.af3_model_root
    stages = parse_stage_list(args.stages)
    run_staged_validation(
        args.manifest,
        args.registry,
        args.output_root,
        args.work_root,
        tools=tool_ids,
        stages=stages,
        python_executable=sys.executable,
    )


if __name__ == "__main__":
    main()
