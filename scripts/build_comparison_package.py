#!/usr/bin/env python3
"""Build a provenance-linked PRISM comparison package from retained evidence.

This script is intentionally an evidence assembler, not a benchmark runner. It
does not modify the canonical PRISM checkout, raw inputs, environments, or
existing evidence. Missing measurements remain empty and are accompanied by an
explicit status/reason field.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path


FRAMEWORK = Path(__file__).resolve().parents[2]
REPO = Path("/scratch/rshadi25/GitHub/PRISM-prescript")
EVIDENCE = FRAMEWORK / "evidence/prism-prescript"
OUT = EVIDENCE / "20260913-comparison"
STEP5 = EVIDENCE / "step5-bounded-validation-20260909"
PANEL946 = REPO / "tmp/agent/20260718-gtalign-template-panel-946-1"
PANEL20K = REPO / "tmp/agent/20260718-benchmark20k-v3"
CHECKED = REPO / "new_template/template/checked_templates.txt"
CALCULATED = REPO / "new_template/template/calculated_templates.txt"
NOTEBOOK_SOURCE = REPO / "notebooks/pipeline_orientation_comparison.ipynb"
USALIGN_WT = REPO / "tmp/agent/worktrees/run-1e6873e93a4f4082a236f7348218ecd5"
PRODIGY_WT = REPO / "tmp/agent/worktrees/run-4f0e2099ca254fde8dbd32b71a96e740"
CONTRACT_WT = Path("/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-prescript-contract-20260913")

MATRIX_FIELDS = [
    "comparison_id", "status", "execution_status", "evidence_class",
    "comparability", "pipeline", "template_panel", "template_count",
    "template_list", "template_list_sha256", "query_set", "case_id",
    "alignment_provider", "alignment_device", "alignment_version",
    "alignment_executable", "alignment_executable_sha256", "surface_method",
    "surface_version", "transformation_filter", "ranking_method", "top_k",
    "refinement_backend", "refinement_version", "evaluator",
    "evaluator_version", "source_revision", "source_branch",
    "source_worktree", "parameters", "output_path", "slurm_job_id", "return_code",
    "partition", "node", "cpu_count", "gpu_type", "initialization_seconds",
    "alignment_seconds", "transformation_seconds", "ranking_seconds",
    "refinement_seconds", "total_wall_seconds", "templates_per_second",
    "candidate_input_count", "alignment_count", "transformed_count",
    "selected_count", "refined_count", "scoreable_count", "success_count",
    "missing_count", "failed_count", "clash_rejection_count",
    "dockq_mean", "dockq_best", "irmsd_mean", "native_like_count",
    "native_like_denominator", "failure_reason", "quality_metric_scope",
    "notes", "evidence_paths",
]

PERFORMANCE_FIELDS = [
    "comparison_id", "status", "evidence_class", "comparability",
    "optimization", "template_panel", "template_count", "template_list",
    "query_set", "case_id", "alignment_provider", "alignment_device",
    "alignment_version", "surface_method", "transformation_filter",
    "ranking_method", "top_k", "refinement_backend", "evaluator",
    "source_revision", "source_worktree", "slurm_job_id", "partition",
    "node", "cpu_count", "gpu_type", "initialization_seconds",
    "alignment_seconds", "transformation_seconds", "ranking_seconds",
    "refinement_seconds", "total_wall_seconds", "baseline_seconds",
    "speedup", "templates_per_second", "candidate_throughput",
    "candidate_set_agreement", "successful_alignment_sides",
    "missing_alignment_sides", "failed_alignment_sides",
    "transformation_success", "clash_rejections", "candidate_before_ranking",
    "candidate_after_ranking", "refinement_count", "quality_scope",
    "scientific_risk", "failure_reason", "notes", "evidence_paths",
]

RANKING_FIELDS = [
    "comparison_id", "status", "evidence_class", "comparability", "case_id",
    "template_panel", "template_count", "query_set", "alignment_provider",
    "surface_method", "transformation_filter", "ranking_method", "top_k",
    "refinement_backend", "evaluator", "source_revision", "slurm_job_id",
    "candidates_entering_ranking", "candidates_forwarded", "ranking_overhead_seconds",
    "refinement_jobs_avoided", "refinement_compute_saved_seconds",
    "total_wall_seconds", "baseline_wall_seconds", "dockq_mean", "dockq_best",
    "irmsd_mean", "native_like_count", "native_like_denominator",
    "best_candidate_before_ranking", "top_ranked_quality", "ranking_regret",
    "no_contact_behavior", "failure_reason", "quality_scope", "notes",
    "evidence_paths",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def nonempty_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def selected_panel_hash(path: Path, count: int) -> str:
    lines = nonempty_lines(path)[:count]
    return hashlib.sha256(("\n".join(lines) + "\n").encode()).hexdigest()


def selected_ids_hash_from_pairs(path: Path) -> str:
    ids = []
    seen = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[1] not in seen:
            seen.add(fields[1])
            ids.append(fields[1])
    return hashlib.sha256(("\n".join(ids) + "\n").encode()).hexdigest()


def missing_side_count(entries: list[dict[str, object]]) -> int:
    return sum(len(entry.get("chains", [])) for entry in entries)


def missing_interfaces(path: Path, interface_root: Path) -> list[dict[str, object]]:
    return missing_interfaces_for_lines(nonempty_lines(path), interface_root)


def missing_interfaces_for_lines(lines: list[str], interface_root: Path) -> list[dict[str, object]]:
    missing = []
    for template_id in lines:
        chains = [
            chain for chain in template_id[4:]
            if not (interface_root / f"{template_id}_{chain}_int.pdb").is_file()
            or (interface_root / f"{template_id}_{chain}_int.pdb").stat().st_size == 0
        ]
        if chains:
            missing.append({"template_id": template_id, "chains": chains})
    return missing


def git(args: list[str]) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(REPO), *args], check=True, text=True,
            capture_output=True, timeout=20,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return "UNKNOWN"


def wt_info(path: Path) -> dict[str, object]:
    try:
        head = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"],
                              check=True, text=True, capture_output=True,
                              timeout=20).stdout.strip()
        status = subprocess.run(["git", "-C", str(path), "status", "--short"],
                                check=True, text=True, capture_output=True,
                                timeout=20).stdout.splitlines()
        return {"path": str(path), "head": head, "dirty_entries": len(status),
                "status": "VALID_WORKTREE"}
    except (OSError, subprocess.SubprocessError):
        return {"path": str(path), "head": "UNKNOWN", "dirty_entries": None,
                "status": "UNAVAILABLE"}


def blank() -> dict[str, object]:
    return {field: "" for field in MATRIX_FIELDS}


def put(row: dict[str, object], **values: object) -> dict[str, object]:
    result = blank()
    result.update(values)
    result["pipeline"] = "alignment → transformation/filtering → ranking → refinement → evaluation"
    return result


def evidence(*paths: Path) -> str:
    return ";".join(str(p) for p in paths)


def tool_version(name: str) -> str:
    return {
        "TMalign": "20220412",
        "GTalign": "0.19.00",
        "USalign": "20241108",
        "PRODIGY": "2.4.0",
        "NACCESS": "Naccess2.1 S.J.Hubbard June 1996",
        "Rosetta": "2022.42",
        "DockQ": "2.1.3",
    }.get(name, "UNKNOWN")


def append_step5(rows: list[dict[str, object]]) -> None:
    source = STEP5 / "evidence/tool-matrix.json"
    raw = json.loads(source.read_text(encoding="utf-8"))["rows"]
    for item in raw:
        cancelled = item.get("slurm_job_id") == "1657005"
        execution = "CANCELLED" if cancelled else "NOT_RUN"
        if item.get("status") == "completed":
            status = "COMPLETED"
            execution = "COMPLETED"
        elif item.get("status") == "not_comparable":
            status = "NOT_COMPARABLE"
            execution = "COMPLETED" if item.get("return_code") == 0 else "UNKNOWN"
        else:
            status = "NOT_RUN"
        aligner = item.get("aligner") or "none"
        align_device = "GPU" if aligner == "GTalign GPU" else ("CPU" if aligner == "GTalign CPU" else "UNKNOWN")
        values = dict(
            comparison_id=item.get("combination_id"), status=status,
            execution_status=execution, evidence_class=(
                "TESTED_RUNTIME_VALIDATED_BOUNDED" if item.get("status") == "completed" else "PLANNED_OR_RETAINED"
            ), comparability=("FIXTURE_ONLY" if item.get("status") == "completed" else "HISTORICAL_OR_PLANNED"),
            template_panel="UNKNOWN", template_count="UNKNOWN",
            query_set=item.get("query_id"), case_id=item.get("case_id"),
            alignment_provider=aligner, alignment_device=align_device,
            alignment_version=item.get("aligner_version") or "UNKNOWN",
            alignment_executable=("/home/rshadi25/.conda/envs/gtalign_env/bin/USalign" if aligner == "USalign" else "UNKNOWN"),
            surface_method=item.get("surface_tool") or "none",
            surface_version=item.get("surface_version") or "UNKNOWN",
            transformation_filter=item.get("evaluator") or "UNKNOWN",
            ranking_method=item.get("ranker") or "none", top_k="UNKNOWN",
            refinement_backend=item.get("refiner") or "none",
            refinement_version=item.get("refiner_version") or "UNKNOWN",
            evaluator=item.get("evaluator") or "UNKNOWN",
            evaluator_version=item.get("evaluator_version") or "UNKNOWN",
            source_revision="1026a7fc0609e19f19c48a68510e7ca39d67573d" if item.get("slurm_job_id") == "1657005" else "UNKNOWN",
            source_branch="feature/prism-cli-parity" if item.get("slurm_job_id") == "1657005" else "UNKNOWN",
            source_worktree=str(USALIGN_WT) if aligner == "USalign" else str(REPO),
            parameters=item.get("exact_command") or "UNKNOWN",
            output_path=item.get("output_root") or "UNKNOWN",
            slurm_job_id=item.get("slurm_job_id") or "",
            return_code=("0:0" if cancelled else item.get("return_code") if item.get("return_code") is not None else "UNKNOWN"),
            partition=item.get("partition") or "UNKNOWN", node=item.get("node") or "UNKNOWN",
            cpu_count=item.get("cpu_count") or "UNKNOWN", gpu_type=item.get("gpu_type") or "UNKNOWN",
            candidate_input_count="UNKNOWN", alignment_count=item.get("alignment_count") or "UNKNOWN",
            transformed_count=item.get("transformation_count") or "UNKNOWN",
            selected_count=item.get("selected_count") or "UNKNOWN",
            refined_count=item.get("refined_count") or "UNKNOWN",
            scoreable_count=item.get("evaluated_count") or "UNKNOWN",
            clash_rejection_count=item.get("clash_rejection_count") or "UNKNOWN",
            failure_reason=("CANCELLED before execution; accounting detail: CANCELLED by 1365446"
                            if cancelled else item.get("failure_reason") or ""),
            quality_metric_scope="NONE" if not item.get("evaluator") else item.get("evaluator"),
            notes=(f"chain_ids={item.get('chain_ids')}; accounting_state=CANCELLED; "
                   "accounting_detail=CANCELLED by 1365446; return_code=0:0"
                   if cancelled else f"chain_ids={item.get('chain_ids')}; return_code={item.get('return_code')}"),
            evidence_paths=evidence(source, STEP5 / "report.md"),
        )
        rows.append(put(blank(), **values))


def append_small_panel(rows: list[dict[str, object]]) -> None:
    for site, job, runtime_t, runtime_g, tps_t, tps_g, mean_diff, p90_diff in (
        ("ai", "1368365", 23.07655195798725, 95.95650402922183, 81.98798518273183, 19.717266892339318, 0.03029016384778011, 0.06787000000000001),
        ("cosbi", "1368366", 22.271294339007, 72.43766740197316, 84.95240425637327, 26.119007801574558, 0.029999487315010582, 0.068),
    ):
        common = dict(
            status="COMPLETED", execution_status="COMPLETED",
            evidence_class="RUNTIME_VALIDATED_ALIGNMENT_ONLY",
            comparability="ALIGNMENT_ONLY_SAME_PANEL",
            template_panel="small_logical_panel_946", template_count=946,
            template_list=str(PANEL946 / f"{site}/workspace/templates/calculated_templates.txt"),
            template_list_sha256=selected_panel_hash(PANEL946 / f"{site}/workspace/templates/calculated_templates.txt", 946),
            query_set="1fgnH", case_id="1fgnH_vs_946_templates",
            surface_method="unknown", surface_version="UNKNOWN",
            transformation_filter="NOT_RUN", ranking_method="none", top_k="all",
            refinement_backend="none", evaluator="alignment score comparison",
            evaluator_version="UNKNOWN", source_revision="UNKNOWN",
            source_branch="UNKNOWN", source_worktree=str(PANEL946 / site / "workspace"),
            output_path=str(PANEL946 / site / "results"), slurm_job_id=job,
            partition="ai" if site == "ai" else "cosbi", node="UNKNOWN",
            cpu_count="UNKNOWN", gpu_type="UNKNOWN", alignment_count=1892,
            candidate_input_count=1892, transformed_count=0, selected_count=0,
            refined_count=0, scoreable_count=0, success_count=1892,
            missing_count=0, failed_count=0, failure_reason="",
            quality_metric_scope="TMalign/GTalign score and match-count agreement only",
            notes=f"Both providers covered all 1,892 query/template-side pairs; mean absolute TM-score difference={mean_diff}; p90={p90_diff}; GTalign device was not recorded.",
            evidence_paths=evidence(PANEL946 / site / "results/summary.json", PANEL946 / site / "exit.json"),
        )
        rows.append(put(blank(), comparison_id=f"small946-{site}-tmalign", alignment_provider="TMalign",
                        alignment_device="CPU_OR_UNKNOWN", alignment_version="20220412",
                        alignment_executable=str(REPO / "external_tools/TMalign"),
                        alignment_executable_sha256=sha256(Path("/home/rshadi25/.conda/envs/gtalign_env/bin/TMalign")),
                        alignment_seconds=runtime_t, total_wall_seconds=runtime_t,
                        templates_per_second=946 / runtime_t, parameters="template_limit=946; score comparison",
                        **common))
        rows.append(put(blank(), comparison_id=f"small946-{site}-gtalign", alignment_provider="GTalign",
                        alignment_device="UNKNOWN", alignment_version="0.19.00",
                        alignment_executable="/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_(device_unrecorded)",
                        alignment_executable_sha256="UNKNOWN", alignment_seconds=runtime_g,
                        total_wall_seconds=runtime_g, templates_per_second=946 / runtime_g,
                        parameters="template_limit=946; score comparison", **common))


def append_large_historical(rows: list[dict[str, object]]) -> None:
    root = PANEL20K / "smoke"
    arms = [
        ("large19855-gtgpu-ext", "GTalign GPU", "GPU", "external Rosetta", "1369613", 906, 60, 20, 0.883, 0.955, 0.67, 19_561, "All 20 scored models high quality; historical report says 60 transforms."),
        ("large19855-gtgpu-pyro", "GTalign GPU", "GPU", "PyRosetta", "1368049", 906, 64, 22, 0.880, 0.955, 0.68, 19_561, "All 22 scored models high quality; historical report says 64 transforms."),
        ("large19855-tm-ext", "TMalign", "CPU", "external Rosetta", "1369054", 20_733, 1_474, 617, "", "", "", 714_780, "Exit manifest reports 20,733 seconds and 617 valid refined models; report gives 1,474 transforms."),
        ("large19855-tm-pyro", "TMalign", "CPU", "PyRosetta", "1369055", 13_257, 1_474, 462, "", "", "", 714_780, "Exit manifest reports 13,257 seconds and 462 valid refined models; report gives 1,474 transforms."),
        ("large19855-gtcpu-ext", "GTalign CPU", "CPU", "external Rosetta", "1369614", 4_145, 0, 0, "", "", "", 27, "Completed with no predictions; 27 alignment JSONs versus 19,561 for GPU."),
    ]
    for cid, aligner, device, refiner, job, wall, transformed, refined, dockq, best, irmsd, alignments, note in arms:
        evaluator = "DockQ 2.1.3" if dockq != "" else "none"
        rows.append(put(blank(), comparison_id=cid, status="COMPLETED", execution_status="COMPLETED",
                        evidence_class="HISTORICAL_RUNTIME_OR_QUALITY_EVIDENCE",
                        comparability="HISTORICAL_OBSERVATIONAL_NOT_CURRENT_PANEL",
                        template_panel="historical_large_panel_19855", template_count=19_855,
                        template_list=str(PANEL20K / "batches/shared_manifest.csv"),
                        template_list_sha256="UNKNOWN", query_set="Benchmark 5.5 rigid (257 pairs)" if aligner != "GTalign CPU" else "Benchmark 5.5 rigid",
                        case_id="20260718-benchmark20k-v3", alignment_provider=aligner,
                        alignment_device=device, alignment_version=tool_version("GTalign" if aligner.startswith("GTalign") else "TMalign"),
                        alignment_executable=("/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_gpu" if aligner == "GTalign GPU" else "/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu" if aligner == "GTalign CPU" else "/home/rshadi25/.conda/envs/gtalign_env/bin/TMalign"),
                        surface_method="unknown", surface_version="UNKNOWN", transformation_filter="TM-score threshold 0.4; historical run",
                        ranking_method="baseline", top_k="all", refinement_backend=refiner,
                        refinement_version="2022.42" if refiner == "external Rosetta" else "UNKNOWN", evaluator=evaluator,
                        evaluator_version="2.1.3" if evaluator else "UNKNOWN", source_revision="UNKNOWN",
                        source_branch="UNKNOWN", source_worktree=str(PANEL20K), parameters="Benchmark 5.5 rigid; historical run; exact full command not retained in matrix",
                        output_path=str(root), slurm_job_id=job, partition="ai", node="UNKNOWN",
                        return_code=0,
                        cpu_count="UNKNOWN", gpu_type="UNKNOWN", total_wall_seconds=wall,
                        alignment_count=alignments, transformed_count=transformed, selected_count="UNKNOWN",
                        refined_count=refined, scoreable_count=(20 if cid.endswith("gtgpu-ext") else 22 if cid.endswith("gtgpu-pyro") else "UNKNOWN"),
                        dockq_mean=dockq if dockq != "" else "", dockq_best=best if best != "" else "",
                        irmsd_mean=irmsd if irmsd != "" else "", failure_reason="" if cid != "large19855-gtcpu-ext" else "completed_no_predictions",
                        quality_metric_scope="Native DockQ spot/aggregate reported by historical artifact" if evaluator else "No quality metric",
                        notes=note, evidence_paths=evidence(PANEL20K / "validation_report.md", PANEL20K / "benchmark_and_multiprot_report.md", root)))


def append_current_exact_panel_gaps(rows: list[dict[str, object]]) -> None:
    for label, count, path, phash in (
        ("small_logical_panel_946", 946, str(PANEL946 / "ai/workspace/templates/calculated_templates.txt"), selected_panel_hash(PANEL946 / "ai/workspace/templates/calculated_templates.txt", 946)),
        ("large_checked_panel_19948", 19_948, str(CHECKED), sha256(CHECKED)),
        ("large_calculated_panel_19062", 19_062, str(CALCULATED), sha256(CALCULATED)),
    ):
        for aligner, device in (("TMalign", "CPU"), ("GTalign GPU", "GPU"), ("USalign", "CPU")):
            cid = f"{label}-{aligner.lower().replace(' ', '-')}-clean-matched"
            rows.append(put(blank(), comparison_id=cid, status="NOT_RUN", execution_status="NOT_RUN",
                            evidence_class="MISSING_CLEAN_MATCHED_RUN", comparability="NOT_COMPARABLE",
                            template_panel=label, template_count=count, template_list=path,
                            template_list_sha256=phash, query_set="Current matched panel not frozen",
                            case_id="clean-current-panel", alignment_provider=aligner,
                            alignment_device=device, alignment_version=tool_version("USalign" if aligner == "USalign" else "GTalign" if aligner.startswith("GTalign") else "TMalign"),
                            alignment_executable=("/home/rshadi25/.conda/envs/gtalign_env/bin/USalign" if aligner == "USalign" else "UNKNOWN"),
                            surface_method="NACCESS", transformation_filter="Current stable gates; exact replay not run",
                            ranking_method="baseline", top_k="all", refinement_backend="external Rosetta",
                            evaluator="DockQ 2.1.3", source_revision=git(["rev-parse", "HEAD"]),
                            source_branch=git(["branch", "--show-current"]), source_worktree=str(REPO),
                            parameters="Not submitted: exact clean same-input/same-template/same-evaluator run required",
                            failure_reason="No clean matched current-panel run; retained historical artifacts use 19,855 or logical 946 panel",
                            quality_metric_scope="UNKNOWN", notes="Explicit gap; not inferred from a sibling arm.",
                            evidence_paths=evidence(CHECKED, CALCULATED, STEP5 / "evidence/tool-matrix.json")))


def append_smoke_rows(rows: list[dict[str, object]]) -> None:
    for aligner, device, job in (("TMalign", "CPU", "1658925"), ("GTalign CPU", "CPU", "1658928"), ("GTalign GPU", "GPU", "1658929")):
        rows.append(put(blank(), comparison_id=f"current-smoke-{aligner.lower().replace(' ', '-')}", status="COMPLETED",
                        execution_status="COMPLETED", evidence_class="DETERMINISTIC_SMOKE",
                        comparability="SMOKE_ONLY", template_panel="single-template-smoke", template_count=1,
                        template_list="current smoke template fixture", query_set="current smoke fixture",
                        case_id="20260912-current-smoke", alignment_provider=aligner,
                        alignment_device=device, alignment_version=tool_version("GTalign" if aligner.startswith("GTalign") else "TMalign"),
                        source_revision=git(["rev-parse", "HEAD"]), source_branch=git(["branch", "--show-current"]), source_worktree=str(REPO),
                        parameters="bash benchmark/scripts/run_prism_pipeline_smoke.sh; --no-refine",
                        output_path="canonical recent smoke artifact", slurm_job_id=job, alignment_count=4,
                        transformed_count=0, selected_count=0, refined_count=0, return_code=0,
                        failure_reason="zero candidates is expected for this smoke; not a quality result",
                        quality_metric_scope="NONE", notes="Successful execution and stage wiring only.",
                        evidence_paths=evidence(REPO / "docs/STABLE_PIPELINE.md")))


def append_new_gtalign_runs(rows: list[dict[str, object]]) -> None:
    """Add only completed/failed manifests from the isolated live benchmark jobs."""
    roots = (
        ("checked_prefix_946", Path("/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-gtalign-20260913-small946-r1")),
        ("checked_full_19948", Path("/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-gtalign-20260913-large19948-r1")),
        ("calculated_19062", Path("/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-gtalign-20260913-large19062-r1")),
    )
    for label, root in roots:
        manifest_path = root / "run_manifest.json"
        if not manifest_path.exists():
            error_logs = sorted((root / "logs").glob("slurm-*.err"))
            if label == "checked_full_19948" and error_logs:
                error = error_logs[-1].read_text(encoding="utf-8").strip()
                rows.append(put(blank(), comparison_id="exact-checked_full_19948-gtalign-gpu-1659357",
                                status="FAILED", execution_status="FAILED", evidence_class="NEW_ALIGNMENT_ONLY_RUN",
                                comparability="EXACT_CURRENT_PANEL_STAGING_FAILED", template_panel="current_checked_full_19948",
                                template_count=19_948, template_list=str(CHECKED), template_list_sha256=sha256(CHECKED),
                                query_set="1fgnH", case_id="20260913-gtalign-checked_full_19948",
                                alignment_provider="GTalign GPU", alignment_device="GPU", alignment_version="0.19.00",
                                alignment_executable="/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_gpu",
                                surface_method="precomputed surface PDB", transformation_filter="NOT_RUN", ranking_method="none",
                                top_k="all", refinement_backend="none", evaluator="none", source_revision=git(["rev-parse", "HEAD"]),
                                source_branch=git(["branch", "--show-current"]), source_worktree=str(REPO),
                                parameters="isolated exact-panel GTalign runner; staging failed closed",
                                output_path=str(root), slurm_job_id="1659357", partition="ai", node="ai01", cpu_count=1,
                                gpu_type="requested GPU", candidate_input_count="UNKNOWN", selected_count="UNKNOWN", return_code=1,
                                missing_count=18,
                                failure_reason=error or "missing interface during staging", quality_metric_scope="NONE",
                                notes="Checked list has 18 missing interface sides across 9 templates; no alignment search was started.",
                                evidence_paths=evidence(error_logs[-1], CHECKED)))
            continue
        item = json.loads(manifest_path.read_text(encoding="utf-8"))
        selected = item.get("selected_template_count", item.get("template_limit", "UNKNOWN"))
        rows.append(put(blank(), comparison_id=f"exact-{label}-gtalign-gpu-{item.get('slurm_job_id', 'unknown')}",
                        status="COMPLETED" if item.get("status") == "COMPLETED" else "FAILED",
                        execution_status=item.get("status", "UNKNOWN"),
                        evidence_class="NEW_ALIGNMENT_ONLY_RUN",
                        comparability="EXACT_CURRENT_PANEL_ALIGNMENT_ONLY",
                        template_panel=f"current_{label}", template_count=selected,
                        template_list=item.get("template_list", "UNKNOWN"),
                        template_list_sha256=(json.loads((root / "template_list.txt.manifest.json").read_text()).get("selected_sha256", "UNKNOWN") if (root / "template_list.txt.manifest.json").exists() else "UNKNOWN"),
                        query_set="1fgnH", case_id=f"20260913-gtalign-{label}",
                        alignment_provider="GTalign GPU", alignment_device="GPU",
                        alignment_version="0.19.00", alignment_executable=item.get("gtalign_executable", "UNKNOWN"),
                        alignment_executable_sha256=sha256(Path(item["gtalign_executable"])) if Path(item.get("gtalign_executable", "")).exists() else "UNKNOWN",
                        surface_method="precomputed surface PDB", transformation_filter="NOT_RUN",
                        ranking_method="none", top_k="all", refinement_backend="none", evaluator="none",
                        source_revision=git(["rev-parse", "HEAD"]), source_branch=git(["branch", "--show-current"]),
                        source_worktree=str(REPO), parameters="isolated exact-panel GTalign runner; --nhits=1 --nalns=1 --dev-min-length=3",
                        output_path=str(root / "output"), slurm_job_id=item.get("slurm_job_id", ""), partition="ai",
                        node=item.get("hostname", "UNKNOWN"), cpu_count=1, gpu_type="requested GPU; telemetry unavailable",
                        initialization_seconds=item.get("initialization_or_staging_seconds", "UNKNOWN"),
                        alignment_seconds=item.get("alignment_search_seconds", "UNKNOWN"), total_wall_seconds=item.get("total_wall_seconds", "UNKNOWN"),
                        candidate_input_count="UNKNOWN", alignment_count="UNKNOWN", transformed_count="NOT_RUN",
                        selected_count="NOT_RUN", refined_count="NOT_RUN", scoreable_count="NOT_RUN",
                        success_count="UNKNOWN",
                        missing_count=missing_side_count(item.get("excluded_missing_interfaces", [])), failed_count="UNKNOWN",
                        return_code=item.get("return_code", "UNKNOWN"), failure_reason="" if item.get("status") == "COMPLETED" else "staging/search failed; see run manifest",
                        quality_metric_scope="NONE", notes=f"Output file count={item.get('output_count', 'UNKNOWN')}; excluded missing interfaces are recorded in the selected-list manifest.",
                        evidence_paths=evidence(manifest_path, root / "template_list.txt.manifest.json", root / "logs")))


def append_new_tmalign_runs(rows: list[dict[str, object]]) -> None:
    roots = (
        ("checked_prefix_946", Path("/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-tmalign-20260913-small946-r1")),
        ("calculated_19062", Path("/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-tmalign-20260913-large19062-r1")),
    )
    for label, root in roots:
        manifest_path = root / "run_manifest.json"
        if not manifest_path.exists():
            continue
        item = json.loads(manifest_path.read_text(encoding="utf-8"))
        counts = item.get("result_status_counts", {})
        panel_path = root / "panel_manifest.json"
        panel = json.loads(panel_path.read_text(encoding="utf-8")) if panel_path.exists() else {}
        rows.append(put(blank(), comparison_id=f"exact-{label}-tmalign-{item.get('slurm_job_id', 'unknown')}",
                        status="COMPLETED" if item.get("status") == "COMPLETED" else "FAILED",
                        execution_status=item.get("status", "UNKNOWN"), evidence_class="NEW_ALIGNMENT_ONLY_RUN",
                        comparability="EXACT_CURRENT_PANEL_ALIGNMENT_ONLY", template_panel=f"current_{label}",
                        template_count=item.get("selected_template_count", "UNKNOWN"), template_list=item.get("template_list", "UNKNOWN"),
                        template_list_sha256=panel.get("selected_template_ids_sha256", "UNKNOWN"), query_set="1fgnH",
                        case_id=f"20260913-tmalign-{label}", alignment_provider="TMalign", alignment_device="CPU",
                        alignment_version=item.get("tmalign_version", "20220412"), alignment_executable=item.get("tmalign_executable", "UNKNOWN"),
                        alignment_executable_sha256=sha256(Path(item["tmalign_executable"])) if Path(item.get("tmalign_executable", "")).exists() else "UNKNOWN",
                        surface_method="precomputed surface PDB", transformation_filter="NOT_RUN", ranking_method="none", top_k="all",
                        refinement_backend="none", evaluator="none", source_revision=git(["rev-parse", "HEAD"]),
                        source_branch=git(["branch", "--show-current"]), source_worktree=str(REPO),
                        parameters="TMalign QUERY REFERENCE -m MATRIX_PATH; 8 independent CPU workers",
                        output_path=str(root / "results"), slurm_job_id=item.get("slurm_job_id", ""), partition="ai",
                        node=item.get("hostname", "UNKNOWN"), cpu_count=item.get("worker_count", 8), gpu_type="none",
                        initialization_seconds=item.get("initialization_or_staging_seconds", "UNKNOWN"),
                        alignment_seconds=item.get("alignment_search_seconds", "UNKNOWN"), total_wall_seconds=item.get("total_wall_seconds", "UNKNOWN"),
                        candidate_input_count=item.get("selected_interface_pair_count", "UNKNOWN"), alignment_count=item.get("result_count", "UNKNOWN"),
                        transformed_count="NOT_RUN", selected_count="NOT_RUN", refined_count="NOT_RUN", scoreable_count="NOT_RUN",
                        success_count=counts.get("success", 0), missing_count=missing_side_count(item.get("excluded_missing_interfaces", [])),
                        failed_count=counts.get("failed", 0) + counts.get("unparseable", 0), return_code=item.get("return_code", "UNKNOWN"),
                        quality_metric_scope="NONE", failure_reason="" if item.get("status") == "COMPLETED" else "TMalign panel incomplete",
                        notes=f"Result status counts={json.dumps(counts, sort_keys=True)}; matrix and stdout hashes are recorded per pair.",
                        evidence_paths=evidence(manifest_path, panel_path, root / "results/alignments.jsonl", root / "logs")))


def append_new_usalign_runs(rows: list[dict[str, object]]) -> None:
    roots = (
        ("checked_prefix_946", Path("/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-usalign-20260913-small946-r1")),
        ("calculated_19062", Path("/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-usalign-20260913-large19062-r1")),
    )
    for label, root in roots:
        manifest_path = root / "run_manifest.json"
        if not manifest_path.exists():
            continue
        item = json.loads(manifest_path.read_text(encoding="utf-8"))
        counts = item.get("result_status_counts", {})
        rows.append(put(blank(), comparison_id=f"exact-{label}-usalign-{item.get('slurm_job_id', 'unknown')}",
                        status="COMPLETED" if item.get("status") == "COMPLETED" else "FAILED",
                        execution_status=item.get("status", "UNKNOWN"), evidence_class="NEW_ALIGNMENT_ONLY_RUN",
                        comparability="EXACT_CURRENT_PANEL_ALIGNMENT_ONLY", template_panel=f"current_{label}",
                        template_count=item.get("selected_template_count", "UNKNOWN"), template_list=item.get("template_list", "UNKNOWN"),
                        template_list_sha256=selected_ids_hash_from_pairs(root / "pairs.tsv"), query_set="1fgnH", case_id=f"20260913-usalign-{label}",
                        alignment_provider="USalign", alignment_device="CPU", alignment_version=item.get("usalign_version", "20241108"),
                        alignment_executable=item.get("usalign_executable", "UNKNOWN"),
                        alignment_executable_sha256=sha256(Path(item["usalign_executable"])) if Path(item.get("usalign_executable", "")).exists() else "UNKNOWN",
                        surface_method="precomputed surface PDB", transformation_filter="NOT_RUN", ranking_method="none", top_k="all",
                        refinement_backend="none", evaluator="none", source_revision=git(["rev-parse", "HEAD"]),
                        source_branch=git(["branch", "--show-current"]), source_worktree=str(REPO),
                        parameters="USalign QUERY REFERENCE -outfmt -1 -m -; 8 independent CPU workers",
                        output_path=str(root / "results"), slurm_job_id=item.get("slurm_job_id", ""), partition="ai",
                        node=item.get("hostname", "UNKNOWN"), cpu_count=item.get("worker_count", 8), gpu_type="none",
                        initialization_seconds=item.get("initialization_or_staging_seconds", "UNKNOWN"),
                        alignment_seconds=item.get("alignment_search_seconds", "UNKNOWN"), total_wall_seconds=item.get("total_wall_seconds", "UNKNOWN"),
                        candidate_input_count=item.get("selected_interface_pair_count", "UNKNOWN"),
                        alignment_count=item.get("result_count", "UNKNOWN"), transformed_count="NOT_RUN",
                        selected_count="NOT_RUN", refined_count="NOT_RUN", scoreable_count="NOT_RUN",
                        success_count=counts.get("success", 0), missing_count=missing_side_count(item.get("excluded_missing_interfaces", [])),
                        failed_count=counts.get("failed", 0) + counts.get("unparseable", 0), return_code=item.get("return_code", "UNKNOWN"),
                        quality_metric_scope="NONE", failure_reason="" if item.get("status") == "COMPLETED" else "USalign panel incomplete",
                        notes=f"Result status counts={json.dumps(counts, sort_keys=True)}; transform rows are recorded per pair.",
                        evidence_paths=evidence(manifest_path, root / "panel_manifest.json", root / "results/alignments.jsonl")))


def ranking_rows() -> list[dict[str, object]]:
    src = STEP5 / "report.md"
    rows: list[dict[str, object]] = []
    for method, topk, status, entering, forwarded, overhead, refined_saved, total, base, reason, scope, note, ev in (
        ("no ranking", "all", "NOT_RUN", "", "", "", "", "", "", "No clean same-set no-ranking arm in retained evidence", "UNKNOWN", "Required comparator", src),
        ("baseline", "all", "COMPLETED", 4, 4, "UNKNOWN", 0, 56.3, 56.3, "", "bounded current-tree load observation", "4 selected; 2 refined", src),
        ("baseline", 1, "COMPLETED", 4, 1, "UNKNOWN", 1, 62.7, 56.3, "", "bounded current-tree load observation", "Candidate load fell 75%; total wall time increased by 6.4s", src),
        ("PRODIGY", 1, "COMPLETED", 2, 1, "UNKNOWN", "UNKNOWN", "UNKNOWN", "UNKNOWN", "", "fixture-only ranking state", "Affinity -65.827; selection only; no total-time or quality claim", EVIDENCE / "step4-prodigy-20260909/report.md"),
        ("PRODIGY", 1, "COMPLETED", 2, 2, "UNKNOWN", "UNKNOWN", "UNKNOWN", "UNKNOWN", "prodigy_failed_no_contacts", "fixture-only failure preservation", "No-contact group preserved both candidates", EVIDENCE / "step4-prodigy-20260909/report.md"),
    ):
        rows.append({field: "" for field in RANKING_FIELDS})
        rows[-1].update(dict(comparison_id=f"ranking-{method.replace(' ', '-')}-{topk}-{len(rows)}", status=status,
                              evidence_class="TESTED_BOUNDED" if status == "COMPLETED" else "MISSING_CLEAN_ARM",
                              comparability="FIXTURE_OR_OBSERVATIONAL" if status == "COMPLETED" else "NOT_COMPARABLE",
                              case_id="1400311-or-prodigy-fixture", template_panel="UNKNOWN", template_count="UNKNOWN",
                              query_set="current-tree bounded example", alignment_provider="unknown",
                              surface_method="unknown", transformation_filter="unknown", ranking_method=method, top_k=topk,
                              refinement_backend="unknown", evaluator="DockQ 2.1.3" if status == "NOT_RUN" else "none",
                              source_revision="UNKNOWN", slurm_job_id="1400311" if "baseline" in method else "1656993/1656992" if method == "PRODIGY" else "",
                              candidates_entering_ranking=entering, candidates_forwarded=forwarded,
                              ranking_overhead_seconds=overhead, refinement_jobs_avoided=refined_saved,
                              total_wall_seconds=total, baseline_wall_seconds=base,
                              best_candidate_before_ranking="UNKNOWN", top_ranked_quality="UNKNOWN", ranking_regret="UNKNOWN",
                              no_contact_behavior="preserved group" if "contact" in reason else "UNKNOWN",
                              failure_reason=reason, quality_scope=scope, notes=note, evidence_paths=str(ev)))
    rows.append({field: "" for field in RANKING_FIELDS})
    rows[-1].update(dict(comparison_id="ranking-bm55-baseline-audit", status="COMPLETED",
                         evidence_class="SCORING_CONTRACT_AUDIT", comparability="RANKING_AUDIT_NOT_PRODIGY",
                         case_id="BM5.5", template_panel="historical BM5.5", template_count="UNKNOWN",
                         query_set="155 groups", alignment_provider="unknown", ranking_method="baseline/native-label audit", top_k=1,
                         candidates_entering_ranking=5828, candidates_forwarded=5827,
                         evaluator="DockQ 2.1.3", best_candidate_before_ranking="median DockQ 0.09339",
                         top_ranked_quality="median DockQ 0.02235", ranking_regret="UNKNOWN",
                         native_like_count="58 top1; oracle-positive 68", native_like_denominator=155,
                         failure_reason="", quality_scope="Ranking/scoring contract evidence; not production superiority",
                         notes="5,827 rankable candidates across 155 groups; 5,828 rows entered audit.",
                         evidence_paths=str(REPO / "tmp/agent/20260727-pipeline-comparison")))
    for method in ("no ranking", "baseline", "PRODIGY"):
        for topk in (1, 3, 5):
            if not any(r["ranking_method"] == method and str(r["top_k"]) == str(topk) for r in rows):
                row = {field: "" for field in RANKING_FIELDS}
                row.update(dict(comparison_id=f"ranking-gap-{method.replace(' ', '-')}-{topk}", status="NOT_RUN",
                                evidence_class="MISSING_CLEAN_ARM", comparability="NOT_COMPARABLE", case_id="required ranking matrix",
                                template_panel="current exact panel not frozen", template_count="UNKNOWN", query_set="same-set/native-labeled cases",
                                alignment_provider="same alignment provider as comparator", ranking_method=method, top_k=topk,
                                evaluator="DockQ 2.1.3", failure_reason="No clean same-input/same-candidate-set run",
                                quality_scope="UNKNOWN", notes="Required before any PRODIGY end-to-end conclusion.", evidence_paths=str(src)))
                rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_notebook() -> None:
    source = json.loads(NOTEBOOK_SOURCE.read_text(encoding="utf-8"))
    cells = [
        {"cell_type": "markdown", "metadata": {}, "source": [
            "# PRISM pipeline comparison (evidence-backed)\n\n",
            "This is the isolated comparison successor to `notebooks/pipeline_orientation_comparison.ipynb`. "
            "It reads the machine-readable artifacts in this package; it does not run a pipeline or mutate project data. "
            "Rows marked `NOT_RUN`, `BLOCKED`, or `NOT_COMPARABLE` are shown rather than omitted.\n"
        ]},
        {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": [
            "from pathlib import Path\nimport csv, json, os\n\n",
            "configured_root = os.environ.get('PRISM_ARTIFACT_ROOT')\n",
            "candidates = [Path(configured_root).expanduser()] if configured_root else []\n",
            "candidates.extend((Path.cwd(), Path.cwd().parent, Path.cwd().parent.parent))\n",
            "ARTIFACT_ROOT = next((candidate.resolve() for candidate in candidates if (candidate / 'PIPELINE_MATRIX.csv').exists()), None)\n",
            "if ARTIFACT_ROOT is None:\n",
            "    raise FileNotFoundError('Set PRISM_ARTIFACT_ROOT to the comparison artifact directory.')\n",
            "DATA = {}\n",
            "for name in ('PIPELINE_MATRIX.csv', 'PERFORMANCE_COMPARISON.csv', 'RANKING_COMPARISON.csv'):\n",
            "    with (ARTIFACT_ROOT / name).open(newline='', encoding='utf-8') as fh:\n",
            "        DATA[name] = list(csv.DictReader(fh))\n",
            "PROVENANCE = {}\n",
            "for name in ('manifest.json', 'TEMPLATE_PANELS.json', 'NO_DROP_SUMMARY.json', 'RECONCILIATION.json'):\n",
            "    path = ARTIFACT_ROOT / name\n",
            "    if path.exists():\n",
            "        PROVENANCE[name] = json.loads(path.read_text(encoding='utf-8'))\n",
            "print({name: len(rows) for name, rows in DATA.items()})\n",
            "print('provenance artifacts:', sorted(PROVENANCE))\n",
        ]},
        {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": [
            "def rows(name, **filters):\n",
            "    result = DATA[name]\n",
            "    for key, value in filters.items():\n",
            "        if value not in (None, '', 'ALL'):\n",
            "            if key == 'outcome':\n",
            "                if value == 'successful':\n",
            "                    result = [r for r in result if r.get('status') == 'COMPLETED' and not r.get('failure_reason')]\n",
            "                elif value == 'failed_or_missing':\n",
            "                    result = [r for r in result if r.get('status') != 'COMPLETED' or r.get('failure_reason')]\n",
            "            else:\n",
            "                result = [r for r in result if r.get(key) == str(value)]\n",
            "    return result\n\n",
            "def applicable_filters(name, filters):\n",
            "    fields = set(DATA[name][0]) if DATA[name] else set()\n",
            "    return {key: value for key, value in filters.items() if key == 'outcome' or key in fields}\n\n",
            "def funnel(record):\n",
            "    return {key: record.get(key, '') for key in ('candidate_input_count', 'alignment_count', 'transformed_count', 'selected_count', 'refined_count', 'scoreable_count')}\n\n",
            "# Parameterized filters; set values to 'ALL' to remove a filter.\n",
            "FILTERS = {'template_panel': 'ALL', 'query_set': 'ALL', 'case_id': 'ALL', 'alignment_provider': 'ALL', 'surface_method': 'ALL', 'ranking_method': 'ALL', 'top_k': 'ALL', 'refinement_backend': 'ALL', 'status': 'ALL', 'outcome': 'ALL'}\n",
            "selected = rows('PIPELINE_MATRIX.csv', **FILTERS)\n",
            "print('selected rows:', len(selected))\n",
        ]},
        {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": [
            "# Stage funnel and provenance warning.\n",
            "for record in selected[:20]:\n",
            "    print(record['comparison_id'], record['status'], funnel(record), record['comparability'])\n",
            "print('WARNING: comparisons are causal only when panel, queries, filters, evaluator, refiner, and provenance match.')\n",
        ]},
        {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": [
            "# Interactive controls backed only by the stored CSV artifacts.\n",
            "try:\n",
            "    import ipywidgets as widgets\n",
            "    from IPython.display import display\n",
            "    FILTER_FIELDS = ('template_panel', 'query_set', 'case_id', 'alignment_provider',\n",
            "                     'surface_method', 'ranking_method', 'top_k',\n",
            "                     'refinement_backend', 'status', 'outcome')\n",
            "    def widget_options(field):\n",
            "        if field == 'outcome':\n",
            "            return ['ALL', 'successful', 'failed_or_missing']\n",
            "        values = sorted({record.get(field, '') for record in DATA['PIPELINE_MATRIX.csv'] if record.get(field, '')})\n",
            "        return ['ALL'] + values\n",
            "    controls = {field: widgets.Dropdown(options=widget_options(field), value='ALL', description=field.replace('_', ' ')[:18])\n",
            "               for field in FILTER_FIELDS}\n",
            "    def show_filtered(**values):\n",
            "        selected_rows = rows('PIPELINE_MATRIX.csv', **values)\n",
            "        print('selected rows:', len(selected_rows))\n",
            "        for record in selected_rows[:20]:\n",
            "            print(record['comparison_id'], record['status'], funnel(record), record['comparability'])\n",
            "        ranking_rows = rows('RANKING_COMPARISON.csv', **applicable_filters('RANKING_COMPARISON.csv', values))\n",
            "        print('ranking rows:', len(ranking_rows))\n",
            "        performance_rows = rows('PERFORMANCE_COMPARISON.csv', **applicable_filters('PERFORMANCE_COMPARISON.csv', values))\n",
            "        print('performance rows:', len(performance_rows))\n",
            "        print('WARNING: causal comparisons require matching panel, queries, filters, evaluator, refiner, and provenance.')\n",
            "    display(widgets.VBox([widgets.HBox(list(controls.values())[i:i + 2]) for i in range(0, len(controls), 2)]))\n",
            "    display(widgets.interactive_output(show_filtered, controls))\n",
            "except ImportError:\n",
            "    print('ipywidgets unavailable; use FILTERS above for portable parameterized analysis.')\n",
        ]},
        {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": [
            "# Dependency-light comparison views. Values are parsed from the stored artifacts.\n",
            "NUMERIC_FIELDS = ('template_count', 'initialization_seconds', 'alignment_seconds',\n",
            "                  'transformation_seconds', 'ranking_seconds', 'refinement_seconds',\n",
            "                  'total_wall_seconds', 'templates_per_second', 'candidate_input_count',\n",
            "                  'alignment_count', 'transformed_count', 'selected_count', 'refined_count',\n",
            "                  'scoreable_count', 'success_count', 'missing_count', 'failed_count',\n",
            "                  'clash_rejection_count', 'dockq_mean', 'dockq_best', 'irmsd_mean',\n",
            "                  'candidates_entering_ranking', 'candidates_forwarded', 'ranking_overhead_seconds',\n",
            "                  'refinement_jobs_avoided', 'refinement_compute_saved_seconds',\n",
            "                  'total_wall_seconds', 'dockq_mean', 'dockq_best', 'irmsd_mean', 'ranking_regret')\n",
            "def number(value):\n",
            "    try:\n",
            "        return float(value) if value not in (None, '', 'UNKNOWN', 'NOT_RUN', 'NOT_MEASURED') else None\n",
            "    except (TypeError, ValueError):\n",
            "        return None\n",
            "def view(record, fields):\n",
            "    return {field: record.get(field, '') for field in fields}\n",
            "stage_fields = ('comparison_id', 'status', 'template_panel', 'alignment_provider',\n",
            "                'candidate_input_count', 'alignment_count', 'transformed_count',\n",
            "                'selected_count', 'refined_count', 'scoreable_count',\n",
            "                'initialization_seconds', 'alignment_seconds', 'transformation_seconds',\n",
            "                'ranking_seconds', 'refinement_seconds', 'total_wall_seconds',\n",
            "                'templates_per_second', 'success_count', 'missing_count', 'failed_count',\n",
            "                'clash_rejection_count', 'comparability')\n",
            "print('STAGE_FUNNEL_RUNTIME_AVAILABILITY')\n",
            "for record in selected:\n",
            "    print(view(record, stage_fields))\n",
            "print('RANKING_REDUCTION_QUALITY')\n",
            "for record in rows('RANKING_COMPARISON.csv', **applicable_filters('RANKING_COMPARISON.csv', FILTERS)):\n",
            "    print(view(record, ('comparison_id', 'status', 'ranking_method', 'top_k',\n",
            "                     'candidates_entering_ranking', 'candidates_forwarded',\n",
            "                     'ranking_overhead_seconds', 'refinement_jobs_avoided',\n",
            "                     'refinement_compute_saved_seconds', 'total_wall_seconds',\n",
            "                     'dockq_mean', 'dockq_best', 'irmsd_mean', 'ranking_regret', 'quality_scope', 'failure_reason')))\n",
            "print('PROVENANCE_WARNING: rows are directly comparable only when query, panel hash,\\n",
            "surface, thresholds, ranking, refiner, evaluator, revision, and stage manifests match.')\n",
        ]},
        {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": [
            "# Optional figures: all values come from stored CSV rows.\n",
            "try:\n",
            "    import pandas as pd\n",
            "    import matplotlib.pyplot as plt\n",
            "    from IPython.display import display\n",
            "    def coerce_numeric(frame, columns):\n",
            "        for column in columns:\n",
            "            if column not in frame:\n",
            "                frame[column] = pd.Series(index=frame.index, dtype='float64')\n",
            "            frame[column] = pd.to_numeric(frame[column], errors='coerce')\n",
            "        return frame\n",
            "    matrix = pd.DataFrame(rows('PIPELINE_MATRIX.csv', **FILTERS))\n",
            "    for col in ('template_panel', 'alignment_provider'):\n",
            "        if col not in matrix:\n",
            "            matrix[col] = ''\n",
            "    coerce_numeric(matrix, ('template_count', 'candidate_input_count', 'alignment_count', 'transformed_count', 'selected_count', 'refined_count', 'scoreable_count', 'success_count', 'missing_count', 'failed_count', 'clash_rejection_count', 'total_wall_seconds', 'dockq_mean', 'dockq_best', 'irmsd_mean'))\n",
            "    stage_cols = ['candidate_input_count', 'alignment_count', 'transformed_count', 'selected_count', 'refined_count', 'scoreable_count']\n",
            "    display(matrix.groupby(['template_panel', 'alignment_provider'], dropna=False)[stage_cols].sum(min_count=1).reset_index())\n",
            "    performance = pd.DataFrame(rows('PERFORMANCE_COMPARISON.csv', **applicable_filters('PERFORMANCE_COMPARISON.csv', FILTERS)))\n",
            "    coerce_numeric(performance, ('template_count', 'initialization_seconds', 'alignment_seconds', 'transformation_seconds', 'ranking_seconds', 'refinement_seconds', 'total_wall_seconds', 'baseline_seconds', 'speedup', 'templates_per_second'))\n",
            "    if not performance.empty:\n",
            "        performance.dropna(subset=['total_wall_seconds']).plot.barh(x='comparison_id', y='total_wall_seconds', legend=False, figsize=(10, 6), title='total wall time')\n",
            "        plt.xlabel('total wall time (s)'); plt.tight_layout(); plt.show()\n",
            "    runtime_cols = ['initialization_seconds', 'alignment_seconds', 'transformation_seconds', 'ranking_seconds', 'refinement_seconds']\n",
            "    runtime = performance.dropna(subset=['total_wall_seconds']).set_index('comparison_id')[runtime_cols]\n",
            "    if not runtime.empty:\n",
            "        runtime.plot.bar(stacked=True, figsize=(12, 5), title='runtime by stage')\n",
            "        plt.ylabel('seconds'); plt.tight_layout(); plt.show()\n",
            "    throughput = performance.dropna(subset=['templates_per_second'])\n",
            "    if not throughput.empty:\n",
            "        throughput.plot.scatter(x='template_count', y='templates_per_second', figsize=(8, 5), title='template-scale throughput')\n",
            "        plt.tight_layout(); plt.show()\n",
            "    quality = matrix.dropna(subset=['dockq_mean'])\n",
            "    if not quality.empty:\n",
            "        quality.plot.scatter(x='total_wall_seconds', y='dockq_mean', figsize=(8, 5), title='speed versus DockQ (stored quality rows only)')\n",
            "        plt.show()\n",
            "    ranking = pd.DataFrame(rows('RANKING_COMPARISON.csv', **applicable_filters('RANKING_COMPARISON.csv', FILTERS)))\n",
            "    coerce_numeric(ranking, ('template_count', 'candidates_entering_ranking', 'candidates_forwarded', 'ranking_overhead_seconds', 'refinement_jobs_avoided', 'refinement_compute_saved_seconds', 'total_wall_seconds', 'baseline_wall_seconds', 'dockq_mean', 'dockq_best', 'irmsd_mean', 'ranking_regret'))\n",
            "    ranking_numeric = ranking.dropna(subset=['candidates_entering_ranking', 'candidates_forwarded'])\n",
            "    if not ranking_numeric.empty:\n",
            "        ranking_numeric.plot.scatter(x='candidates_entering_ranking', y='candidates_forwarded', figsize=(7, 5), title='ranking candidate-load reduction')\n",
            "        plt.tight_layout(); plt.show()\n",
            "except Exception as exc:\n",
            "    print(f'Optional pandas/matplotlib figures unavailable ({type(exc).__name__}: {exc}); dependency-light tables and filters remain available.')\n",
        ]},
        {"cell_type": "markdown", "metadata": {}, "source": [
            "## Preserved orientation protocol\n\n",
            "The original orientation-comparison notebook is preserved by reference, not embedded in this report view. "
            "Its source path, SHA-256, and cell count are recorded below and in `manifest.json`; execution requires the project environment and a compute-capable notebook kernel.\n"
        ]},
    ]
    # Preserve the established protocol as a nested, reviewable appendix. This
    # keeps the successor connected to the existing notebook without allowing
    # its canonical-path execution cells to run implicitly in the report view.
    cells.append({"cell_type": "markdown", "metadata": {"source_notebook": str(NOTEBOOK_SOURCE)},
                 "source": [f"Original source: `{NOTEBOOK_SOURCE}`\n\n",
                             f"Source SHA-256: `{sha256(NOTEBOOK_SOURCE)}`\n",
                             f"Original cells retained by reference: {len(source.get('cells', []))}.\n"]})
    for index, cell in enumerate(cells, start=1):
        cell["id"] = f"prism{index:04d}"
    notebook = {"cells": cells, "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.11"},
        "prism_artifact_root": "..", "prism_source_notebook": str(NOTEBOOK_SOURCE),
    }, "nbformat": 4, "nbformat_minor": 5}
    (OUT / "notebooks").mkdir(parents=True, exist_ok=True)
    (OUT / "notebooks/prism_pipeline_comparison.ipynb").write_text(
        json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


def write_documents(rows: list[dict[str, object]], performance: list[dict[str, object]], ranking: list[dict[str, object]]) -> None:
    statuses: dict[str, int] = {}
    for row in rows:
        statuses[row["status"]] = statuses.get(row["status"], 0) + 1
    current_head = git(["rev-parse", "HEAD"])
    dirty = git(["status", "--short"]).splitlines()
    panel946 = selected_panel_hash(PANEL946 / "ai/workspace/templates/calculated_templates.txt", 946)
    new_alignment_rows = [row for row in rows if row["evidence_class"] == "NEW_ALIGNMENT_ONLY_RUN"]
    new_alignment_summary = "\n".join(
        f"| {row['comparison_id']} | {row['status']} | {row['template_count']} | {row['initialization_seconds']} | {row['alignment_seconds']} | {row['return_code']} |"
        for row in new_alignment_rows
    ) or "| none | NOT_RUN | — | — | — | — |"
    reconciliation = OUT / "RECONCILIATION.md"
    reconciliation.write_text(f"""# VALAR reconciliation checkpoint — 2026-09-13

This package preserves the two evidence namespaces while treating them as one
selected underlying project.

* PRISM-prescript job `1657005`: live accounting shows `CANCELLED by 1365446`,
  zero runtime, and no durable output. The old Step 5 `blocked/unknown`
  description is superseded for this job.
* prism-refactoring run `run-2688ab679e0f47409f2a19e182483a75`: job `1657891`
  failed with exit `127:0` and retry `1657894` failed with exit `1:0` on
  `ai22`. The direct manifest still points at retry state `PENDING`/`WAITING_FOR_JOB`,
  but checkpoints, logs, handoff, and accounting show no worker completion
  report. The discrepancy is retained, not rewritten.
* Slurm controllers were up at this checkpoint, but current AI resources were
  fragmented by unrelated user allocations. No replacement was submitted.
* New isolated runner job `1659355` failed before staging because of a shell
  initialization bug (`/etc/bashrc` with `set -u`); fixed GTalign reruns are
  separately manifested as jobs `1659356` and `1659358`. Checked-panel
  staging job `1659357` failed closed on a missing interface and remains
  preserved. The exact USalign panels completed as jobs `1659360` (946-prefix)
  and `1659361` (19,062-entry calculated panel), using the transform-producing
  `-outfmt -1 -m -` invocation and eight CPU workers.
  The matched TMalign references completed as jobs `1659448` (946-prefix) and
  `1659449` (19,062-entry calculated panel), using eight CPU workers and
  per-side JSONL/no-drop records.

Authoritative inputs: `evidence/prism-prescript/`,
`evidence/prism-refactoring/run-2688ab679e0f47409f2a19e182483a75/`, and the
live accounting output recorded during this goal.\n""", encoding="utf-8")

    (OUT / "CURRENT_STATUS.md").write_text(f"""# PRISM current status — 2026-09-13

## Bottom line

The maintained canonical pipeline is operationally implemented and its focused
software tests/smokes pass, but no clean same-input/same-template/same-evaluator
comparison has scientifically validated USalign, PRODIGY, or a current exact
19K panel. Existing high DockQ and runtime results are historical observational
evidence, not a promotion-grade causal matrix.

Canonical source: `{REPO}` at `{current_head}`, branch `{git(['branch', '--show-current'])}`, with `{len(dirty)}` dirty entries. It was not edited.

## Evidence state

| State | Meaning in this package |
|---|---|
| IMPLEMENTED | Code/path exists in canonical or an isolated candidate. |
| TESTED | Focused tests or deterministic probes passed. |
| VALIDATED | Runtime and quantitative output evidence exists within the stated scope. |
| REVIEWED | Source/evidence review completed; this does not upgrade scientific validation. |
| NOT_COMPARABLE / NOT_RUN | Inputs, evaluator, provenance, or execution are insufficient for causal comparison. |

Matrix row counts: {json.dumps(statuses, sort_keys=True)}.

## Implemented and tested

* TMalign is the maintained default; GTalign CPU/GPU are opt-in; NACCESS and
  FreeSASA are surface choices; external Rosetta is the stable refiner.
* USalign and PRODIGY changes exist only in isolated worktrees. USalign's
  transform-producing invocation is `-outfmt -1 -m -`; its bounded real probe
  produced 214 matches, TM-score 0.98359, and replay RMSD 0.7500643.
* PRODIGY 2.4.0 bounded success selected one candidate (affinity -65.827) and
  its no-contact failure preserved the full candidate group. This is state
  handling, not docking-quality validation.
* The isolated shared-contract candidate passed 28 focused contract/adapter/
  import tests and 31 tests in total. It fail-closes unavailable optional
  backends and preserves unknown/rejected alignment states; it is not wired
  into canonical provider writers and is not promoted.
* Focused canonical/isolated tests recorded previously were 22/22, 20/20,
  and 9/9 respectively. The full canonical suite recorded 344 passed and 6
  skipped in the latest retained run.

## Runtime and quality evidence

* Logical 946-template alignment panel: 1 query × 946 templates × 2 sides =
  1,892 pairs; TMalign 22.27–23.08 s and GTalign 72.44–95.96 s across two
  sites. All pairs were present; mean absolute TM-score difference was about
  0.030. GTalign device was not recorded, so this is not GPU evidence.
* Historical 19,855-template BM5.5 run: GTalign GPU wall time 906 s, 60/64
  transforms, 20/22 scored models, mean best DockQ 0.883/0.880; TMalign had
  1,474 reported transforms and much larger alignment output. GTalign CPU
  completed with no predictions (27 alignment JSONs). These are historical
  artifacts with incomplete parity and known mapping/CPU-GPU caveats.
* Current exact on-disk denominators are 19,948 checked entries (SHA-256
  `{sha256(CHECKED)}`) and 19,062 calculated entries (SHA-256
  `{sha256(CALCULATED)}`). The 946 panel is logical first-946 selection from a
  historical 19,855-entry list (selected-panel SHA-256 `{panel946}`). The new
  current checked-prefix 946 runs use a distinct selected-panel SHA-256
  `{selected_panel_hash(CHECKED, 946)}`.

## New exact-panel alignment-only runs

| Arm | Status | Selected templates | Staging seconds | Search seconds | Return code |
|---|---|---:|---:|---:|---:|
{new_alignment_summary}

These jobs use one query and precomputed surface PDBs. They validate search
execution, transform-producing output/return handling, and stage timing only;
they do not validate transformation filtering, refinement, ranking, DockQ, or
scientific superiority.

## Decision at this checkpoint

* Fastest observed completed arm: historical GTalign GPU + external Rosetta at
  906 s on the 19,855-entry BM5.5 run. It is not yet a scientifically
  defensible promoted default because current-panel parity and stage timing are
  incomplete.
* Best observed quality subset: historical GTalign GPU + external Rosetta,
  mean best DockQ 0.928 (12 scored models); GTalign GPU + PyRosetta was 0.925
  (13 scored models). The small scored subset and mapping gaps limit this claim.
* Recommended reproducible default today: retain the documented TMalign +
  NACCESS + external Rosetta recipe. Use GTalign GPU + external Rosetta as the
  provisional large-library speed candidate after exact-panel validation.
* GTalign GPU is provisionally the preferred large-library search engine from
  observed runtime/output evidence; USalign is not yet a stable alternative.
  PRODIGY has no demonstrated end-to-end compute saving or acceptable quality
  retention. A ~20K search is operationally demonstrated historically but not
  yet promotion-ready on the current exact panel.
* Highest-value optimization: instrument and reduce repeated alignment output,
  structure parsing, and staging while proving candidate-set agreement.

## Blockers and uncertainty

The maintained path lacks a uniform AlignmentResult/no-drop schema across all
providers; external Rosetta does not preserve complete per-candidate return
codes; transformation audit does not always contain actual clash counts;
MultiProt's RMSD-derived proxy is not a TM-score; exact current-panel GTalign
GPU, TMalign, and USalign evidence is alignment-only, while current exact
transformation/filter, ranking, refinement, and evaluator arms are absent; and
no PRODIGY end-to-end timing or quality retention evidence exists. Job
reconciliation is in `RECONCILIATION.md`. The isolated contract candidate
remains partial because it lacks event/parent correlation, does not wrap every
stage, and is not connected to provider writers. Existing canonical
`run_evidence`, lineage, evaluator, and completion contracts remain
authoritative; the candidate is an additive projection and must not replace
them.

## Authoritative artifacts

* `PIPELINE_MATRIX.csv/json`: all retained, planned, incomplete, and historical arms.
* `TEMPLATE_PANELS.json` and `NO_DROP_SUMMARY.json`: frozen panel definitions and retained no-drop accounting.
* `PERFORMANCE_COMPARISON.csv`: exact 946 alignment evidence, historical 19,855 results, and explicit current-panel gaps.
* `RANKING_COMPARISON.csv`: bounded baseline/PRODIGY evidence and required missing no-ranking/top-k arms.
* `notebooks/prism_pipeline_comparison.ipynb`: portable artifact-driven report view.
* `PIPELINE_COMPARISON.md`, `STABLE_PIPELINES.md`, `VALIDATION_REPORT.md`, `NEXT_ACTIONS.md`, and `CANDIDATE_CHANGES.md`.
""", encoding="utf-8")

    (OUT / "PIPELINE_COMPARISON.md").write_text("""# PRISM pipeline comparison

## Stage-by-stage interpretation

The matrix uses the common stage order `alignment → transformation/filtering →
ranking → refinement → evaluation`. Empty fields are unknown, not zero. A zero
candidate count is only a result when the stage manifest explicitly records it.

### Alignment

The retained logical 946 panel is the only same-panel alignment comparison with
complete pair counts: both TMalign and GTalign produced 1,892/1,892 pair records.
The runtime artifact does not record whether GTalign was CPU or GPU, so it cannot
answer the GPU-vs-CPU question. The historical large run suggests GTalign GPU is
operationally attractive (906 s) but it used a 19,855-entry historical list and
is not a clean comparison to current 19,948/19,062 panels. The bounded USalign
probe established the transform-producing invocation; the new exact-panel
USalign runs extend that evidence to alignment-only timing.

The new isolated alignment-only runs provide a cleaner current-panel search
comparison. On the checked-list 946-template prefix, GTalign GPU took 7.89 s
total (4.81 s search), TMalign took 12.31 s total (11.48 s search), and USalign
took 22.23 s total (21.50 s search). On the calculated panel, 19,062 source
entries yielded 19,058 materialized templates: GTalign GPU took 144.62 s total
(73.96 s search), TMalign took 157.18 s total (144.95 s search), and USalign
took 452.11 s total (439.97 s search), each CPU reference using eight workers.
These are alignment-stage throughput observations only: GTalign emits a
batched search output whereas the pairwise runners record one result per
interface, so candidate-set agreement and downstream quality are still
unmeasured.

### Transformation and filtering

The validated retained 1gte ledger contains 2,997 alignment JSONs, 560 complete
orientation attempts, 489 clash rejections, and 71 retained geometry passes.
It also records 1,163 structural missing-partner sides and 1,877 residual
unconsumed sides. The cause of 73,251 missing/unwritten potential sides is not
known. This is a historical geometry-only ledger, not a provider comparison.

### Ranking

Bounded evidence shows candidate-load reduction but not speedup. In the retained
current-tree example, baseline forwarded 4 candidates and refined 2; a top-1
ranked arm forwarded 1 and refined 1, but wall time was 62.7 s versus 56.3 s.
PRODIGY selection and no-contact preservation are tested; DockQ retention,
regret, and end-to-end savings remain unknown. MultiProt's RMSD-derived proxy is
kept out of TM-score comparisons.

### Refinement and evaluation

External Rosetta and PyRosetta have historical successful outputs, but their
per-candidate failure observability and evaluator parity are incomplete. DockQ
2.1.3 is wired in the repository-local scoring environment. Historical GTalign
GPU quality is high in the scored subset, but only 2/10 BM5.5 pairs were directly
scorable in one report because of chain mapping; this limits causal conclusions.

## Scalability and reliability

Historical TMalign generated 714,780 alignment JSONs for the large run, whereas
GTalign GPU generated 19,561 alignment hits and 60/64 transforms in the cited
arms. This indicates a likely alignment/search and filesystem/output bottleneck,
but stage-separated timing and candidate-set agreement for an exact common panel
are missing. The highest-value optimization is therefore a clean, batched,
stage-instrumented GTalign GPU versus USalign versus TMalign run at both exact
panels, with reusable preprocessing and no-drop manifests.

## Comparison rule

Only rows with matching query set, template manifest/hash, surface, thresholds,
ranking, refiner, evaluator, source revision, and complete stage manifests may
support a causal claim. All other rows remain historical, observational, or
not comparable in the CSV.

## Final decision question

The current evidence supports only a provisional answer: GTalign GPU + external
Rosetta is the fastest and best observed large-library arm, TMalign + NACCESS +
external Rosetta is the recommended reproducibility-first default, and PRODIGY
should remain opt-in pending same-set top-k DockQ/regret and end-to-end timing.
Before promotion, exact current-panel counts/hashes, matched inputs, stage
timings, USalign compatibility, GPU/CPU parity, complete evaluator mapping,
and independent review are still required.
""", encoding="utf-8")

    (OUT / "STABLE_PIPELINES.md").write_text("""# Stable PRISM pipelines

## Stable operational recipe

The only maintained recipe documented as stable is TMalign + NACCESS + external
Rosetta, using the current `gtalign_env` interpreter and defaults from
`docs/STABLE_PIPELINE.md`:

```bash
PRISM_PIPELINE_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \\
  bash benchmark/scripts/run_prism_pipeline_smoke.sh
```

The smoke is an execution/wiring check and may legitimately produce zero
candidates. It is not scientific validation.

## Bounded candidates

* USalign: isolated candidate only; use the transform-producing equivalent of
  `USalign A.pdb B.pdb -outfmt -1 -m -`. Expected output must include parsed
  match pairs, rotation matrix, translation, score fields, return code, and an
  explicit failure/empty record. It is not promoted or stable.
* PRODIGY: isolated opt-in ranking candidate only, version 2.4.0. Affinity is a
  ranking feature, not a DockQ metric. No stable command is claimed until a
  same-set top-k experiment measures overhead, savings, regret, and quality.
* GTalign CPU/GPU: current opt-in implementations pass retained software/smoke
  checks, but CPU/GPU score and candidate-set parity is unresolved. No exact
  current-panel stable benchmark command is claimed here.
* MultiProt, FiberDock, FreeSASA, and PyRosetta remain supported/experimental
  or legacy choices with explicit evidence gaps; MultiProt's proxy score must
  remain separate from TM-score.

Promotion threshold for every provider: focused regression tests, deterministic
small smoke, real Slurm run, stage/no-drop manifests, reproducibility metadata,
quantitative output checks, matched cross-arm comparison, and independent review.
""", encoding="utf-8")

    (OUT / "VALIDATION_REPORT.md").write_text("""# PRISM validation report

## IMPLEMENTED

Canonical `prism.py` supports TMalign, GTalign CPU/GPU, MultiProt, NACCESS,
FreeSASA, external Rosetta, PyRosetta, FiberDock, and optional ranking paths.
USalign and PRODIGY implementations exist in isolated candidates.

## TESTED

The latest retained evidence includes canonical focused tests (22/22), isolated
USalign tests (20/20), isolated PRODIGY tests (9/9), full canonical tests
(344 passed, 6 skipped), a real USalign transform probe, real PRODIGY success
and no-contact probes, current stable execution smokes, and the generated
comparison notebook's static/portable checks.

The notebook passed nbformat/AST validation, loaded all three CSVs plus the
four provenance JSON artifacts through `PRISM_ARTIFACT_ROOT`, and passed a
direct successful-versus-failed-or-missing filter check (15/44 of 59 matrix
rows). Kernel-backed headless execution completed with zero cell errors using
the explicit non-interactive `Agg` plotting backend; the validation record is
`NOTEBOOK_VALIDATION.json`. The default login-node Matplotlib backend is
environment-sensitive, so reproducible headless execution should set
`MPLBACKEND=Agg` and a writable `MPLCONFIGDIR`.

The isolated shared-contract candidate also passed 34 focused contract/
adapter/import/ledger tests and 37 tests in its complete clean-baseline suite.
The benchmark-side ledger is a separate, additive projection from resolved
candidate manifests; this is software evidence only and does not validate a
provider, a full pipeline, or scientific output.

## VALIDATED (bounded scope only)

USalign transform parsing/application, PRODIGY state/failure preservation,
1gte geometry-only attrition accounting, the historical 946 alignment pair
coverage, the exact current checked-prefix/calculated-panel alignment-only
timings for GTalign GPU, TMalign, and USalign, and historical 19,855-run
output/quality subsets have quantitative evidence within their stated scopes.
None establishes general scientific superiority or promotion readiness.

## REVIEWED

The canonical guidance/memory, source call paths, isolated worktrees, benchmark
manifests, notebooks, and both VALAR namespaces were reviewed. The review found
missing uniform score/status/no-drop contracts, incomplete refiner observability,
and known CPU/GPU/panel/evaluator mismatches. Independent Luna High review
confirmed that canonical `run_evidence`, `investigation_lineage`,
`investigation_contracts`, and `pipeline_completion_contract` remain the
authoritative layers; the isolated candidate is TESTED only and is not a
replacement.

## BLOCKED / UNKNOWN

The current exact panels have matched GTalign GPU/TMalign/USalign alignment-only
runs, but no clean matched full-pipeline comparison through transformation/filtering,
ranking, refinement, and evaluation; current exact TMalign timing is also
present only for alignment-stage runs. PRODIGY lacks end-to-end timing and quality/regret evidence, and no
ranking top-k matrix is complete. Job 1657005 was canceled without execution;
run-2688 had two failed worker attempts without a scientific report. These
remain explicit in the matrix.

Successful Slurm or process return codes are not being promoted to scientific
validation without stage counts, provenance, quantitative outputs, and review.
""", encoding="utf-8")

    (OUT / "NEXT_ACTIONS.md").write_text("""# Next actions, prioritized

1. Finalize the exact-panel protocol: logical 946 with its selected-entry
   manifest and calculated 19,062 with four explicit missing-interface
   exclusions are recorded; decide whether checked 19,948 is admissible only
   after resolving its 18 missing-interface sides. Record hashes, duplicate
   handling, and leakage checks.
2. Keep the existing canonical `run_evidence`, lineage, evaluator, and
   completion contracts authoritative. Use the isolated AlignmentResult/stage
   candidate only as an additive projection after exact downstream alignment
   paths and benchmark foreign keys are resolved; its adapter and
   optional-import repair are tested but are not integrated or promoted.
3. Run the missing full downstream GTalign-GPU/TMalign/USalign arms at the
   finalized panels using the highest safe live parallelism; add GTalign CPU
   only after resolving the existing CPU/GPU discrepancy. Preserve all logs
   under the run.
4. Run no-ranking, baseline, and PRODIGY top-1/3/5 on identical candidate sets
   with native labels and DockQ/iRMSD. Report overhead, avoided refinement, wall
   time, quality retention, regret, and no-contact behavior.
5. Validate exact-stage timing and candidate-set agreement before introducing
   caching, batching, query parallelism, or early filtering. The likely highest
   value is reducing repeated parsing/output overhead while preserving the
   candidate set.
6. Obtain independent review of the candidate changes and only then consider a
   promotion decision. Do not merge, push, or promote automatically.
""", encoding="utf-8")

    (OUT / "CANDIDATE_CHANGES.md").write_text(f"""# Isolated candidate changes

## USalign

Worktree: `{USALIGN_WT}`; HEAD `1026a7fc0609e19f19c48a68510e7ca39d67573d`;
feature-gated CLI branch and `src/alignment_usalign.py`; focused suite 20/20;
real probes jobs 1656964/1656966. The adapter preserves empty/failure records
and requires `-outfmt -1 -m -`. It is not on canonical `prism.py` and was not
promoted.

## PRODIGY

Worktree: `{PRODIGY_WT}`; HEAD `1026a7fc0609e19f19c48a68510e7ca39d67573d`;
feature-gated `src/prodigy_ranker.py`/candidate selector; focused suite 9/9;
real probes jobs 1656992/1656993. The candidate preserves no-contact groups
and uses explicit ranking states. It is not a quality-validated production
default and was not promoted.

## Shared stage contract

Worktree: `{CONTRACT_WT}`; pinned baseline `1026a7fc0609e19f19c48a68510e7ca39d67573d`;
isolated Luna High candidate files `src/pipeline_contract.py`,
`src/alignment_result_adapter.py`, `benchmark/scripts/
build_alignment_event_ledger.py`, the additive `prism.py` stage-event
integration, and their tests. The focused contract/adapter/import/ledger suite
passed 34 tests and the complete isolated suite passed 37 tests. The candidate
also guards absent optional backends so default TMalign CLI import remains
usable while requested unavailable backends fail explicitly. Its adapter
preserves unknown, rejection, refinement-failure, and score-gate statuses as
explicit non-success records rather than treating a score as implicit success;
it also retains optional run/attempt identifiers and provider-specific numeric
metrics such as MultiProt RMSD. The ledger preserves one explicit event per
resolved manifest row, including missing/not-run records, and rejects missing
identity or duplicate candidate IDs before output. The candidate is not
integrated into the canonical checkout and does not claim provider,
transformation, ranking, refinement, or evaluator validation.

## Review risks before promotion

The canonical path still has nonuniform provider score semantics, incomplete
per-candidate refinement return codes, possible transformation audit gaps, and
an absent `.github/copilot-instructions.md` at the requested path. The isolated
contract remains partial: it lacks complete event/parent correlation fields,
does not yet wrap every pipeline stage, and the adapter/ledger are not wired
into provider writers. The ledger is tested against synthetic resolved
manifests only and does not repair the canonical raw-output hash semantics or
foreign-key lineage by itself. The candidate diffs must be rebased/reviewed against the dirty
canonical state in a new bounded worktree; this package deliberately makes no
source merge.
""", encoding="utf-8")

    (OUT / "NOTEBOOK_VALIDATION.json").write_text(json.dumps({
        "schema_version": "prism-notebook-validation-20260913",
        "notebook": "notebooks/prism_pipeline_comparison.ipynb",
        "static_checks": {
            "status": "PASS",
            "nbformat": "4.5",
            "cell_count": 9,
            "code_cell_ast_parse": "PASS",
            "manifest_nested_file_listing": "PASS",
        },
        "portable_artifact_root_check": {
            "status": "PASS",
            "configuration": "PRISM_ARTIFACT_ROOT",
            "loaded_rows": {
                "PIPELINE_MATRIX.csv": 59,
                "PERFORMANCE_COMPARISON.csv": 60,
                "RANKING_COMPARISON.csv": 13,
            },
            "provenance_files_loaded": [
                "manifest.json", "TEMPLATE_PANELS.json", "NO_DROP_SUMMARY.json", "RECONCILIATION.json",
            ],
        },
        "filter_check": {
            "status": "PASS",
            "pipeline_matrix_successful_rows": 15,
            "pipeline_matrix_failed_or_missing_rows": 44,
            "total_rows": 59,
        },
        "headless_kernel_check": {
            "status": "PASS",
            "command": "MPLBACKEND=Agg MPLCONFIGDIR=<writable-temp-dir> PRISM_ARTIFACT_ROOT=<package> jupyter nbconvert --to notebook --execute notebooks/prism_pipeline_comparison.ipynb",
            "cell_errors": 0,
            "notes": "The default login-node Matplotlib backend emitted an environment-specific optional-figure warning; the explicit non-interactive Agg validation completed without that warning. The dependency-light tables and filters remain independent of optional plotting packages.",
        },
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    append_step5(rows)
    append_small_panel(rows)
    append_large_historical(rows)
    append_current_exact_panel_gaps(rows)
    append_smoke_rows(rows)
    append_new_gtalign_runs(rows)
    append_new_tmalign_runs(rows)
    append_new_usalign_runs(rows)
    performance = []
    for row in rows:
        if row["alignment_provider"] in {"TMalign", "GTalign", "GTalign CPU", "GTalign GPU", "USalign"}:
            perf = {field: "" for field in PERFORMANCE_FIELDS}
            alignment_only = row["evidence_class"] == "NEW_ALIGNMENT_ONLY_RUN"
            perf.update({
                "comparison_id": row["comparison_id"], "status": row["status"],
                "evidence_class": row["evidence_class"], "comparability": row["comparability"],
                "optimization": "none", "template_panel": row["template_panel"], "template_count": row["template_count"],
                "template_list": row["template_list"], "query_set": row["query_set"], "case_id": row["case_id"],
                "alignment_provider": row["alignment_provider"], "alignment_device": row["alignment_device"],
                "alignment_version": row["alignment_version"], "surface_method": row["surface_method"],
                "transformation_filter": row["transformation_filter"], "ranking_method": row["ranking_method"],
                "top_k": row["top_k"], "refinement_backend": row["refinement_backend"], "evaluator": row["evaluator"],
                "source_revision": row["source_revision"], "source_worktree": row["source_worktree"],
                "slurm_job_id": row["slurm_job_id"], "partition": row["partition"], "node": row["node"],
                "cpu_count": row["cpu_count"], "gpu_type": row["gpu_type"],
                "initialization_seconds": row["initialization_seconds"],
                "alignment_seconds": row["alignment_seconds"],
                "transformation_seconds": "NOT_RUN" if alignment_only else row["transformation_seconds"],
                "ranking_seconds": "NOT_RUN" if alignment_only else row["ranking_seconds"],
                "refinement_seconds": "NOT_RUN" if alignment_only else row["refinement_seconds"],
                "total_wall_seconds": row["total_wall_seconds"],
                "templates_per_second": row["templates_per_second"], "candidate_throughput": "NOT_MEASURED",
                "candidate_before_ranking": "NOT_MEASURED" if alignment_only else row["candidate_input_count"],
                "candidate_after_ranking": "NOT_RUN" if alignment_only else row["selected_count"],
                "refinement_count": "NOT_RUN" if alignment_only else row["refined_count"],
                "quality_scope": row["quality_metric_scope"],
                "successful_alignment_sides": row["success_count"],
                "missing_alignment_sides": row["missing_count"],
                "failed_alignment_sides": row["failed_count"],
                "transformation_success": "NOT_RUN" if alignment_only else row["transformed_count"],
                "clash_rejections": "NOT_RUN" if alignment_only else row["clash_rejection_count"],
                "failure_reason": row["failure_reason"], "notes": row["notes"], "evidence_paths": row["evidence_paths"],
            })
            performance.append(perf)
    new_performance = [row for row in performance if row["evidence_class"] == "NEW_ALIGNMENT_ONLY_RUN"]
    for panel in {row["template_panel"] for row in new_performance}:
        same_panel = [row for row in new_performance if row["template_panel"] == panel]
        total_by_provider = {}
        for row in same_panel:
            try:
                total_by_provider[row["alignment_provider"]] = float(row["total_wall_seconds"])
            except (TypeError, ValueError):
                pass
        for row in same_panel:
            try:
                total = float(row["total_wall_seconds"])
                templates = float(row["template_count"])
                search = float(row["alignment_seconds"])
                row["templates_per_second"] = f"{templates / total:.6f}"
                row["candidate_throughput"] = "NOT_MEASURED"
            except (TypeError, ValueError, ZeroDivisionError):
                row["templates_per_second"] = "UNKNOWN"
                row["candidate_throughput"] = "UNKNOWN"
            row["candidate_set_agreement"] = "NOT_MEASURED"
            if row["alignment_provider"] == "GTalign GPU" and "USalign" in total_by_provider:
                row["baseline_seconds"] = f"{total_by_provider['USalign']:.6f}"
                row["speedup"] = f"{total_by_provider['USalign'] / float(row['total_wall_seconds']):.6f}"
                row["notes"] += " Speedup is relative to USalign total wall time for the same alignment-only panel; not an end-to-end pipeline speedup."
    for label, count in (("small_logical_panel_946", 946), ("large_checked_panel_19948", 19_948), ("large_calculated_panel_19062", 19_062)):
        for opt in ("GTalign GPU batching", "safe GPU parallelism", "USalign CPU task parallelism", "template preprocessing/caching", "reusable alignment results"):
            perf = {field: "" for field in PERFORMANCE_FIELDS}
            perf.update(dict(comparison_id=f"optimization-gap-{label}-{opt.replace(' ', '-')}", status="NOT_RUN",
                             evidence_class="UNTESTED_OPTIMIZATION", comparability="NOT_COMPARABLE", optimization=opt,
                             template_panel=label, template_count=count, query_set="frozen panel required",
                             alignment_provider="GTalign GPU or USalign", failure_reason="No controlled optimization experiment run",
                             quality_scope="Candidate-set agreement required", scientific_risk="Must demonstrate no scientific candidate loss",
                             notes="Do not infer speedup from candidate-load reduction; baseline and optimized stage timings are missing.",
                             evidence_paths=str(OUT / "NEXT_ACTIONS.md")))
            performance.append(perf)
    ranking = ranking_rows()
    write_csv(OUT / "PIPELINE_MATRIX.csv", rows, MATRIX_FIELDS)
    write_csv(OUT / "PERFORMANCE_COMPARISON.csv", performance, PERFORMANCE_FIELDS)
    write_csv(OUT / "RANKING_COMPARISON.csv", ranking, RANKING_FIELDS)
    (OUT / "PIPELINE_MATRIX.json").write_text(json.dumps({
        "schema_version": "prism-pipeline-matrix-20260913", "generated_at": datetime.now(timezone.utc).isoformat(),
        "project_aliases": ["PRISM-prescript", "prism-refactoring"],
        "status_vocabulary": ["NOT_RUN", "BLOCKED", "FAILED", "NOT_COMPARABLE", "COMPLETED"],
        "rows": rows,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_documents(rows, performance, ranking)
    build_notebook()
    calculated_panel_run = Path("/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/prism-tmalign-20260913-large19062-r1/panel_manifest.json")
    calculated_panel = json.loads(calculated_panel_run.read_text(encoding="utf-8")) if calculated_panel_run.exists() else {}
    panels = {
        "schema_version": "prism-template-panel-manifest-20260913",
        "observed_at": "2026-09-13",
        "panels": [
            {
                "label": "small_logical_panel_946", "exact_count": 946,
                "selection_rule": "first 946 non-empty entries of retained historical 19,855-entry list",
                "path": str(PANEL946 / "ai/workspace/templates/calculated_templates.txt"),
                "full_source_count": len(nonempty_lines(PANEL946 / "ai/workspace/templates/calculated_templates.txt")),
                "full_source_sha256": sha256(PANEL946 / "ai/workspace/templates/calculated_templates.txt"),
                "selected_entries_sha256": selected_panel_hash(PANEL946 / "ai/workspace/templates/calculated_templates.txt", 946),
                "materializable_interface_missing": missing_interfaces(PANEL946 / "ai/workspace/templates/calculated_templates.txt", REPO / "new_template/template/interfaces"),
                "duplicate_count_selected": 0, "exclusions": "not recorded in retained artifact",
                "leakage": "historical artifact reports 0/10 benchmark target overlap; current-panel leakage not recomputed",
                "source_date": "2026-07-18", "version": "retained logical panel",
            },
            {
                "label": "current_checked_prefix_946", "exact_count": 946,
                "selection_rule": "first 946 non-empty entries of current checked_templates.txt",
                "path": str(CHECKED),
                "full_source_count": len(nonempty_lines(CHECKED)),
                "full_source_sha256": sha256(CHECKED),
                "selected_entries_sha256": selected_panel_hash(CHECKED, 946),
                "materializable_interface_missing": missing_interfaces_for_lines(nonempty_lines(CHECKED)[:946], REPO / "new_template/template/interfaces"),
                "duplicate_count_selected": len(nonempty_lines(CHECKED)[:946]) - len(set(nonempty_lines(CHECKED)[:946])),
                "exclusions": "none at list-selection stage; materialization manifest is authoritative",
                "leakage": "not recomputed for this package",
                "source_date": "filesystem observation 2026-09-13",
                "version": "current checked-list prefix used by jobs 1659356/1659360",
            },
            {
                "label": "large_checked_panel_19948", "exact_count": len(nonempty_lines(CHECKED)),
                "selection_rule": "all non-empty entries in current checked_templates.txt",
                "path": str(CHECKED), "sha256": sha256(CHECKED),
                "duplicate_count": len(nonempty_lines(CHECKED)) - len(set(nonempty_lines(CHECKED))),
                "materializable_interface_missing": missing_interfaces(CHECKED, REPO / "new_template/template/interfaces"),
                "exclusions": "not inferred; checked-list policy must be frozen before benchmark",
                "leakage": "not recomputed for this package", "source_date": "filesystem observation 2026-09-13",
                "version": "current on-disk checked list",
            },
            {
                "label": "large_calculated_panel_19062", "exact_count": len(nonempty_lines(CALCULATED)),
                "selection_rule": "all non-empty entries in current calculated_templates.txt",
                "path": str(CALCULATED), "sha256": sha256(CALCULATED),
                "duplicate_count": len(nonempty_lines(CALCULATED)) - len(set(nonempty_lines(CALCULATED))),
                "materializable_interface_missing": missing_interfaces(CALCULATED, REPO / "new_template/template/interfaces"),
                "exclusions": "calculation/preflight exclusions are not reconstructed here",
                "materialized_count": calculated_panel.get("selected_template_count", "UNKNOWN"),
                "materialized_selected_ids_sha256": calculated_panel.get("selected_template_ids_sha256", "UNKNOWN"),
                "materialized_interface_pair_count": calculated_panel.get("selected_interface_pair_count", "UNKNOWN"),
                "materialized_exclusions": calculated_panel.get("excluded_missing_interfaces", []),
                "materialized_evidence": str(calculated_panel_run),
                "leakage": "not recomputed for this package", "source_date": "filesystem observation 2026-09-13",
                "version": "current on-disk calculated list",
            },
            {
                "label": "historical_large_panel_19855", "exact_count": 19855,
                "selection_rule": "historical 2026-07-18/19 benchmark list",
                "path": str(PANEL20K), "sha256": "not available as current exact file",
                "duplicate_count": "unknown", "exclusions": "historical final_list workflow",
                "leakage": "historical report: 0/10 target overlap and 17/17 used templates present",
                "source_date": "2026-07-19", "version": "benchmark20k-v3 historical",
            },
        ],
        "warning": "historical logical 946, current checked-prefix 946, 19,062, 19,855, and 19,948 are distinct panel definitions; do not pool them silently.",
    }
    (OUT / "TEMPLATE_PANELS.json").write_text(json.dumps(panels, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    no_drop_source = STEP5 / "evidence/no-drop-ledger.json"
    no_drop = json.loads(no_drop_source.read_text(encoding="utf-8"))
    no_drop["package_scope"] = "retained 1gte ledger only; unknown stages remain null"
    no_drop["package_source"] = str(no_drop_source)
    (OUT / "NO_DROP_SUMMARY.json").write_text(json.dumps(no_drop, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "RECONCILIATION.json").write_text(json.dumps({
        "schema_version": "prism-valar-reconciliation-20260913",
        "observed_at": "2026-09-13",
        "project_aliases": ["PRISM-prescript", "prism-refactoring"],
        "prism_prescript": {"job_id": "1657005", "accounting_state": "CANCELLED", "accounting_detail": "CANCELLED by 1365446", "elapsed": "00:00:00", "return_code": "0:0", "durable_output": False, "interpretation": "never-executed validation attempt", "new_isolated_alignment_jobs": {"1659355": "FAILED runner shell initialization", "1659356": "COMPLETED GTalign GPU small checked prefix", "1659357": "FAILED GTalign GPU checked-panel staging missing interface", "1659358": "COMPLETED GTalign GPU calculated panel with recorded exclusions", "1659360": "COMPLETED USalign small checked prefix; 8 CPU workers; transform-producing invocation", "1659361": "COMPLETED USalign calculated panel; 8 CPU workers; transform-producing invocation", "1659448": "COMPLETED TMalign small checked prefix; 8 CPU workers", "1659449": "COMPLETED TMalign calculated panel; 8 CPU workers"}},
        "prism_refactoring": {"run_id": "run-2688ab679e0f47409f2a19e182483a75", "attempts": [{"job_id": "1657891", "state": "FAILED", "exit_code": "127:0", "node": "ai22"}, {"job_id": "1657894", "state": "FAILED", "exit_code": "1:0", "node": "ai22"}], "direct_manifest_state": "WAITING_FOR_JOB/PENDING (stale)", "worker_report": False, "interpretation": "no worker completion evidence"},
        "scheduler": {"controllers": "UP at checkpoint", "user_jobs": "unrelated allocations present", "new_submission": False},
        "evidence_paths": [str(FRAMEWORK / "status/prism-prescript/26-09-28-LATEST.md"), str(FRAMEWORK / "status/prism-refactoring/26-09-13-LATEST.md"), str(FRAMEWORK / "evidence/prism-refactoring/run-2688ab679e0f47409f2a19e182483a75/manifest.json"), str(FRAMEWORK / "evidence/prism-refactoring/run-2688ab679e0f47409f2a19e182483a75/checkpoints")],
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    source_inventory = {
        "schema_version": "prism-source-inventory-20260913",
        "canonical": {"path": str(REPO), "head": git(["rev-parse", "HEAD"]), "branch": git(["branch", "--show-current"]), "dirty_entries": len(git(["status", "--short"]).splitlines())},
        "guidance": [str(REPO / "AGENTS.md"), str(REPO / ".github/copilot-instructions.md"), str(REPO / "docs/STABLE_PIPELINE.md"), str(REPO / ".planning/STATE.md")],
        "missing_guidance": [str(REPO / ".github/copilot-instructions.md")],
        "memory": [str(REPO / ".agents/skills/project-memory/references/summary.md"), str(REPO / ".agents/skills/project-memory/references/decisions.md"), str(REPO / ".agents/skills/project-memory/references/open_questions.md")],
        "worktrees": [wt_info(USALIGN_WT), wt_info(PRODIGY_WT), wt_info(CONTRACT_WT)],
        "notebook_source": {"path": str(NOTEBOOK_SOURCE), "sha256": sha256(NOTEBOOK_SOURCE), "cells": len(json.loads(NOTEBOOK_SOURCE.read_text(encoding="utf-8")).get("cells", []))},
        "status_sources": [str(FRAMEWORK / "status/prism-prescript/26-09-28-LATEST.md"), str(FRAMEWORK / "status/prism-refactoring/26-09-13-LATEST.md")],
        "evidence_sources": [str(STEP5 / "evidence/tool-matrix.json"), str(EVIDENCE / "step2-transformation-attrition-20260908"), str(EVIDENCE / "step3-usalign-20260908"), str(EVIDENCE / "step4-prodigy-20260909"), str(PANEL20K / "validation_report.md")],
    }
    (OUT / "SOURCE_INVENTORY.json").write_text(json.dumps(source_inventory, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema_version": "prism-comparison-package-20260913",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "project_aliases": ["PRISM-prescript", "prism-refactoring"],
        "canonical_project": str(REPO), "canonical_head": git(["rev-parse", "HEAD"]),
        "canonical_branch": git(["branch", "--show-current"]),
        "canonical_dirty_entries": len(git(["status", "--short"]).splitlines()),
        "source_notebook": str(NOTEBOOK_SOURCE), "source_notebook_sha256": sha256(NOTEBOOK_SOURCE),
        "template_panels": {
            "small_logical_946": {"count": 946, "source": str(PANEL946 / "ai/workspace/templates/calculated_templates.txt"), "selected_hash": selected_panel_hash(PANEL946 / "ai/workspace/templates/calculated_templates.txt", 946), "full_source_count": len(nonempty_lines(PANEL946 / "ai/workspace/templates/calculated_templates.txt"))},
            "current_checked_prefix_946": {"count": 946, "source": str(CHECKED), "selected_hash": selected_panel_hash(CHECKED, 946), "full_source_count": len(nonempty_lines(CHECKED)), "full_source_sha256": sha256(CHECKED)},
            "large_checked_19948": {"count": len(nonempty_lines(CHECKED)), "source": str(CHECKED), "sha256": sha256(CHECKED)},
            "large_calculated_19062": {"count": len(nonempty_lines(CALCULATED)), "source": str(CALCULATED), "sha256": sha256(CALCULATED)},
            "historical_large_19855": {"count": 19_855, "source": str(PANEL20K), "sha256": "UNKNOWN; historical list not current exact panel"},
        },
        "executables": {name: {"path": path, "sha256": sha256(Path(path)) if Path(path).exists() else "UNKNOWN", "version": version} for name, path, version in (
            ("TMalign", "/home/rshadi25/.conda/envs/gtalign_env/bin/TMalign", "20220412"),
            ("USalign", "/home/rshadi25/.conda/envs/gtalign_env/bin/USalign", "20241108"),
            ("GTalign_CPU", "/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_cpu", "0.19.00"),
            ("GTalign_GPU", "/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign_gpu", "0.19.00"),
        )},
        "isolated_candidates": {"USalign": wt_info(USALIGN_WT), "PRODIGY": wt_info(PRODIGY_WT), "shared_stage_contract": wt_info(CONTRACT_WT)},
        "valar_status_sources": [str(FRAMEWORK / "status/prism-prescript/26-09-28-LATEST.md"), str(FRAMEWORK / "status/prism-refactoring/26-09-13-LATEST.md")],
        "reconciliation": {"prism_prescript_1657005": "CANCELLED/zero-runtime/no-output", "prism_refactoring_run_2688": "1657891 FAILED 127; 1657894 FAILED 1; direct manifest stale"},
        "files": sorted([p.name for p in OUT.iterdir() if p.is_file()] + ["notebooks/prism_pipeline_comparison.ipynb"]),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT), "matrix_rows": len(rows), "performance_rows": len(performance), "ranking_rows": len(ranking)}, indent=2))


if __name__ == "__main__":
    main()
