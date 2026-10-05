#!/usr/bin/env python3
"""Build frozen arm, experiment, task, and analysis-plan artifacts.

The manifest is a planning contract, not a result table.  It deliberately
contains blocked arms/experiments when prerequisites are unavailable rather
than silently removing them from the comparison.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark.scripts.build_investigation_source_manifest import _read_rows


ARM_FIELDS = (
    "arm_id",
    "label",
    "aligner",
    "filter_profile",
    "refiner",
    "ranking",
    "evaluator",
    "template_panel",
    "template_assets",
    "compatibility_adapter",
    "source_snapshot",
    "executable_sha256",
    "status",
    "scope",
)
EXPERIMENT_FIELDS = (
    "experiment_id",
    "objective",
    "factor_changed",
    "fixed_factors",
    "native_used_for_selection",
    "top_k",
    "cohort",
    "status",
    "required_artifacts",
)
TASK_FIELDS = (
    "task_id",
    "cohort",
    "dataset_row_id",
    "arm",
    "experiment",
    "template_shard",
    "refinement",
    "scientific_attempt",
    "scheduler_retry_id",
    "input_manifest",
    "output_root",
    "status",
)


def sha256_file(path: Path) -> str:
    if not path.is_file():
        return "missing"
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_hash(values: object) -> str:
    return hashlib.sha256(json.dumps(values, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def arm_rows(root: Path, output_dir: Path) -> list[dict[str, str]]:
    current_tmalign = root / "external_tools/TMalign"
    historical_tmalign = root / "working_version/TMalign/external_tools/TMalign"
    multiprot = root / "working_version/multiprot/external_tools/multiprot/multiprot.Linux"
    rows = [
        {
            "arm_id": "current_tmalign_rosetta",
            "label": "current root TM-align -> current filters -> Rosetta",
            "aligner": "external_tools/TMalign",
            "filter_profile": "current_root",
            "refiner": "Rosetta",
            "ranking": "native-independent current interaction/ranking score",
            "evaluator": "standardized_evaluator + frozen contracts",
            "template_panel": "matched_946",
            "template_assets": "tmp/agent/20260713-historical-current-investigation/final-preflight/current",
            "compatibility_adapter": "none",
            "source_snapshot": "prism.py;src/alignment.py;src/transformation.py;src/rosetta_refinement.py",
            "executable_sha256": sha256_file(current_tmalign),
            "status": "ready-for-smoke",
            "scope": "primary",
        },
        {
            "arm_id": "historical_tmalign_rosetta",
            "label": "working_version/TMalign -> historical filters -> Rosetta",
            "aligner": "working_version/TMalign/external_tools/TMalign",
            "filter_profile": "historical_protocol",
            "refiner": "Rosetta",
            "ranking": "native-independent historical score",
            "evaluator": "standardized_evaluator + frozen contracts",
            "template_panel": "matched_946",
            "template_assets": "tmp/agent/20260713-historical-current-investigation/final-preflight/historical",
            "compatibility_adapter": "none",
            "source_snapshot": "working_version/TMalign/prism.py;working_version/TMalign/run_files",
            "executable_sha256": sha256_file(historical_tmalign),
            "status": "ready-for-smoke",
            "scope": "primary-diagnostic",
        },
        {
            "arm_id": "historical_multiprot_fiberdock",
            "label": "working_version/multiprot -> historical filters -> FiberDock",
            "aligner": "working_version/multiprot/external_tools/multiprot/multiprot.Linux",
            "filter_profile": "historical_protocol",
            "refiner": "FiberDock",
            "ranking": "within-arm FiberDock global energy",
            "evaluator": "standardized_evaluator + paper metrics where reproducible",
            "template_panel": "matched_946",
            "template_assets": "tmp/agent/20260713-historical-current-investigation/final-preflight/historical",
            "compatibility_adapter": "benchmark/scripts/multiprot_compat_adapter.py",
            "source_snapshot": "working_version/multiprot/mainController.py;working_version/multiprot/structuralAlignment.py;working_version/multiprot/prism.ini",
            "executable_sha256": sha256_file(multiprot),
            "status": "blocked-until-compatibility-smoke",
            "scope": "primary-diagnostic",
        },
    ]
    for row in rows:
        snapshot_paths = [root / item.strip() for item in row["source_snapshot"].split(";")]
        row["source_snapshot_sha256"] = json_hash({str(path): sha256_file(path) for path in snapshot_paths})
    return rows


def experiment_rows() -> list[dict[str, str]]:
    return [
        {
            "experiment_id": "source_gate",
            "objective": "validate all four curated role files and expected chains",
            "factor_changed": "none",
            "fixed_factors": "benchmark CSVs; archive; Biopython policy",
            "native_used_for_selection": "no",
            "top_k": "",
            "cohort": "repository_bm5_bm5_5_extension",
            "status": "ready",
            "required_artifacts": "source_manifest.tsv;structure_validation.tsv;staged_sources.tsv",
        },
        {
            "experiment_id": "candidate_generation_matched_946",
            "objective": "compare current TM-align and MultiProt candidates on identical inputs",
            "factor_changed": "aligner",
            "fixed_factors": "curated inputs; matched 946 templates; no filters/refinement",
            "native_used_for_selection": "no",
            "top_k": "",
            "cohort": "repository_bm5_bm5_5_extension",
            "status": "blocked-until-smoke",
            "required_artifacts": "candidate tables;lineage.tsv;executable hashes",
        },
        {
            "experiment_id": "filter_replay",
            "objective": "replay current and historical filters over saved candidates",
            "factor_changed": "filter profile",
            "fixed_factors": "candidate streams; templates; transforms",
            "native_used_for_selection": "no",
            "top_k": "",
            "cohort": "repository_bm5_bm5_5_extension",
            "status": "blocked-until-candidate-streams",
            "required_artifacts": "filter decisions; rejection reasons; lineage.tsv",
        },
        {
            "experiment_id": "canonical_transform",
            "objective": "compare canonical heavy-atom poses before refinement",
            "factor_changed": "transformation implementation",
            "fixed_factors": "matched alignments; common filter; no refiner",
            "native_used_for_selection": "no",
            "top_k": "",
            "cohort": "repository_bm5_bm5_5_extension",
            "status": "blocked-until-candidate-streams",
            "required_artifacts": "poses.tsv;coordinate hashes;orientation records",
        },
        {
            "experiment_id": "refinement_crossover",
            "objective": "send identical canonical poses through no-refinement, Rosetta, and FiberDock",
            "factor_changed": "refiner",
            "fixed_factors": "canonical pose hashes; native-independent ranking; top-20 budget",
            "native_used_for_selection": "no",
            "top_k": "20",
            "cohort": "repository_bm5_bm5_5_extension",
            "status": "blocked-until-compatible-poses",
            "required_artifacts": "poses.tsv;tool-specific derivative hashes;energy records",
        },
        {
            "experiment_id": "ranking_budget",
            "objective": "compare native-independent within-arm ranking under a common top-20 budget",
            "factor_changed": "ranking",
            "fixed_factors": "same candidate/refined pose stream; same evaluator",
            "native_used_for_selection": "no",
            "top_k": "1,3,5,20",
            "cohort": "repository_bm5_bm5_5_extension",
            "status": "ready-after-evaluator-contract",
            "required_artifacts": "candidate ranking keys; selected poses; scores_global.tsv",
        },
        {
            "experiment_id": "evaluator_adversarial",
            "objective": "test mapping, invariance, null semantics, and multichain fixtures",
            "factor_changed": "evaluator fixture",
            "fixed_factors": "frozen mappings; raw DockQ JSON hashes",
            "native_used_for_selection": "no",
            "top_k": "",
            "cohort": "evaluator_fixtures",
            "status": "ready",
            "required_artifacts": "fixture results; mapping.tsv; null-semantics assertions",
        },
        {
            "experiment_id": "confirmatory_bm55",
            "objective": "paired unconditional best GlobalDockQ at top-20 for all 257 rows",
            "factor_changed": "arm composite (descriptive only)",
            "fixed_factors": "all source/evaluator/ranking contracts; three CSVs",
            "native_used_for_selection": "no",
            "top_k": "1,3,5,20",
            "cohort": "repository_bm5_bm5_5_extension",
            "status": "blocked-until-source-and-smoke-gates",
            "required_artifacts": "pair_summary.tsv;scores_global.tsv;scores_interfaces.tsv;claim_ledger.tsv",
        },
    ]


def task_rows(root: Path, rows: list[dict[str, object]], arms: list[dict[str, str]], experiments: list[dict[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for row in rows:
        dataset_row_id = str(row["dataset_row_id"])
        combinations = [("provenance", "source_gate")]
        combinations.extend((arm["arm_id"], experiment["experiment_id"]) for arm in arms for experiment in experiments if experiment["experiment_id"] != "source_gate")
        for arm_id, experiment_id in combinations:
            task_id = f"repository_bm5_bm5_5_extension:{dataset_row_id}:{arm_id}:{experiment_id}:attempt-0"
            output.append(
                {
                    "task_id": task_id,
                    "cohort": "repository_bm5_bm5_5_extension",
                    "dataset_row_id": dataset_row_id,
                    "arm": arm_id,
                    "experiment": experiment_id,
                    "template_shard": "matched_946" if experiment_id not in {"source_gate", "evaluator_adversarial"} else "",
                    "refinement": "none" if experiment_id in {"source_gate", "candidate_generation_matched_946", "filter_replay", "canonical_transform", "evaluator_adversarial"} else "arm-defined",
                    "scientific_attempt": "0",
                    "scheduler_retry_id": "",
                    "input_manifest": "source_manifest.tsv",
                    "output_root": f"tasks/{hashlib.sha256(task_id.encode()).hexdigest()[:16]}",
                    "status": "planned",
                }
            )
    return output


def build_manifests(repo_root: str | Path, output_dir: str | Path) -> dict[str, Path]:
    root = Path(repo_root).resolve()
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows = _read_rows(root, limit=None)
    arms = arm_rows(root, output)
    experiments = experiment_rows()
    write_tsv(output / "arm_manifest.tsv", ARM_FIELDS + ("source_snapshot_sha256",), arms)
    write_tsv(output / "experiment_manifest.tsv", EXPERIMENT_FIELDS, experiments)
    tasks = task_rows(root, rows, arms, experiments)
    write_tsv(output / "task_manifest.tsv", TASK_FIELDS, tasks)
    plan = {
        "schema_version": "prism-investigation-analysis-plan/v1",
        "cohort": "repository_bm5_bm5_5_extension",
        "dataset_source": ["benchmark/data/T_Rigid.csv", "benchmark/data/T_medium.csv", "benchmark/data/T_difficult.csv"],
        "row_count": len(rows),
        "row_counts_by_dataset": {dataset: sum(row["dataset"] == dataset for row in rows) for dataset in ("rigid", "medium", "difficult")},
        "primary_outcome": "best_GlobalDockQ_at_20_unconditional",
        "secondary_outcomes": ["top_1", "top_3", "top_5", "eligible_success_utility", "conditional_quality", "failure_stage", "runtime", "peak_rss"],
        "primary_arms": [arm["arm_id"] for arm in arms],
        "matched_template_panel": 946,
        "historical_template_inventory_secondary": 21072,
        "native_independent_ranking": True,
        "structural_metric_no_model": None,
        "no_model_unconditional_utility": 0.0,
        "cluster_units": ["normalized_input_key", "native_complex_when_available", "sequence_family_when_available"],
        "multiple_testing": {"primary": "Holm", "secondary": "Benjamini-Hochberg"},
        "paper_metrics": ["FiberDock global-energy ranking", "I-Score where available", "iRMSD", "LRMSD", "Fnat"],
        "modern_metric": "DockQ 2.1.3; secondary; never used for selection",
        "bm3_exact_paper_cohort": "blocked; not inferred from these CSVs",
        "source_hashes": {str(path): sha256_file(root / path) for path in ["benchmark/data/T_Rigid.csv", "benchmark/data/T_medium.csv", "benchmark/data/T_difficult.csv"]},
    }
    plan_path = output / "analysis_plan.json"
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    findings_fields = ("finding_id", "step", "classification", "what_checked", "validation", "supporting_reference", "difference_found", "evidence_path", "evidence_sha256", "root_cause", "uncertainty", "follow_up")
    write_tsv(output / "findings.tsv", findings_fields, [])
    claim_fields = ("claim_id", "claim", "classification", "experiment_ids", "task_ids", "source_hashes", "score_rows", "reference_sections", "status")
    write_tsv(output / "claim_ledger.tsv", claim_fields, [])
    return {
        "arm_manifest": output / "arm_manifest.tsv",
        "experiment_manifest": output / "experiment_manifest.tsv",
        "task_manifest": output / "task_manifest.tsv",
        "analysis_plan": plan_path,
        "findings": output / "findings.tsv",
        "claim_ledger": output / "claim_ledger.tsv",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    paths = build_manifests(args.repo_root, args.output_dir)
    print(f"wrote {len(paths)} investigation planning artifacts to {Path(args.output_dir).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
