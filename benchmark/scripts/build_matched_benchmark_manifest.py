#!/usr/bin/env python3
"""Build a frozen, source-gated manifest for matched PRISM comparisons.

The builder is intentionally planning-only: it never downloads structures,
rewrites benchmark inputs, or selects candidates using native scores.  It
binds every task to ``dataset_row_id`` and records the exact executable and
template-panel provenance needed by the downstream Slurm runner.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark.scripts.build_investigation_source_manifest import _read_rows
from benchmark.scripts.investigation_contracts import DEFAULT_RANKING_KEYS, validate_native_independent_ranking_keys


PILOT_ROW_IDS = (
    "rigid:000002", "rigid:000004", "rigid:000060", "rigid:000073",
    "medium:000001", "medium:000012", "medium:000014", "medium:000026",
    "difficult:000016", "difficult:000020", "difficult:000025", "difficult:000028",
)

COHORT_FIELDS = (
    "dataset_row_id", "dataset", "difficulty", "source_file", "source_row_number",
    "native_complex", "receptor", "ligand", "chain_context", "prior_model_count",
    "source_policy", "audit_only", "status",
)
ARM_FIELDS = (
    "arm_id", "label", "aligner", "refiner", "filter_profile", "ranking_keys",
    "template_workflow", "template_panel", "template_manifest", "tm_score_threshold", "status", "scope",
)
TASK_FIELDS = (
    "task_id", "dataset_row_id", "arm", "experiment", "aligner", "refiner",
    "template_workflow", "template_panel", "tm_score_threshold", "top_k", "input_manifest", "output_root", "status",
)
PIPELINE_INPUT_FIELDS = (
    "pair_id", "dataset_row_id", "dataset", "difficulty", "source_row_number",
    "native_complex", "Receptor", "Ligand", "chain_context", "source_policy",
    "audit_only",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def _chain_context(native_complex: str) -> str:
    left, right = (part.strip() for part in native_complex.split(":", 1))
    left_chains = left.split("_", 1)[1] if "_" in left else ""
    return "multichain" if len(left_chains) > 1 or len(right) > 1 else "single_chain"


def _prior_model_counts(path: Path | None) -> dict[str, int]:
    if path is None or not path.is_file():
        return {}
    counts: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("pipeline") != "tmalign_rosetta":
                continue
            pair_id = row.get("pair_id", "")
            if pair_id:
                counts[pair_id] = counts.get(pair_id, 0) + 1
    return counts


def _pair_id(row: dict[str, Any]) -> str:
    # The historical comparison manifest uses the CSV source-row number in
    # pair IDs (header is row 1), while dataset_row_id uses a zero-based
    # dataset row counter.  Preserve that distinction explicitly.
    source_number = int(row["source_row_number"])
    return f"{row['dataset']}_{source_number:04d}"


def select_rows(
    repo_root: Path,
    policy_path: Path,
    cohort: str,
    prior_models: Path | None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    decision = policy["decision"]
    excluded = set(decision["excluded_dataset_row_ids"])
    all_rows = _read_rows(repo_root, limit=None)
    by_id = {row["dataset_row_id"]: row for row in all_rows}
    if cohort == "pilot":
        wanted = list(PILOT_ROW_IDS)
    elif cohort in {"strict-clean", "strict_clean", "240"}:
        wanted = [row["dataset_row_id"] for row in all_rows if row["dataset_row_id"] not in excluded]
    elif cohort in {"all", "257"}:
        wanted = [row["dataset_row_id"] for row in all_rows]
    else:
        raise ValueError("cohort must be pilot, strict-clean, or all")
    unknown = sorted(set(wanted) - set(by_id))
    if unknown:
        raise ValueError("unknown dataset_row_id(s): " + ", ".join(unknown))
    if cohort == "pilot" and any(row_id in excluded for row_id in wanted):
        raise ValueError("pilot includes source-policy excluded rows")
    counts = _prior_model_counts(prior_models)
    rows: list[dict[str, Any]] = []
    for row_id in wanted:
        source = by_id[row_id]
        pair_id = _pair_id(source)
        rows.append({
            "dataset_row_id": row_id,
            "dataset": source["dataset"],
            "difficulty": source["difficulty"],
            "source_file": source["source_file"],
            "source_row_number": source["source_row_number"],
            "native_complex": source["native_complex"],
            "receptor": source["source_row"]["PDB ID 1"],
            "ligand": source["source_row"]["PDB ID 2"],
            "chain_context": _chain_context(source["native_complex"]),
            "prior_model_count": counts.get(pair_id, 0),
            "source_policy": "strict_clean" if row_id not in excluded else "audit_only",
            "audit_only": "true" if row_id in excluded else "false",
            "status": "eligible" if row_id not in excluded else "audit_only",
        })
    if cohort == "pilot":
        expected = {"rigid": 4, "medium": 4, "difficult": 4}
        observed = {name: sum(row["dataset"] == name for row in rows) for name in expected}
        if observed != expected:
            raise ValueError(f"pilot difficulty counts do not match {expected}: {observed}")
        if any(sum(row["chain_context"] == context and row["dataset"] == name for row in rows) != 2
               for name in expected for context in ("single_chain", "multichain")):
            raise ValueError("pilot must contain two single-chain and two multichain rows per difficulty")
    return rows, policy


def build_manifests(repo_root: str | Path, output_dir: str | Path, *, cohort: str = "pilot",
                    source_policy: str | Path | None = None, prior_models: str | Path | None = None) -> dict[str, Path]:
    root = Path(repo_root).resolve()
    output = Path(output_dir).resolve()
    policy_path = Path(source_policy) if source_policy else root / "tmp/agent/20260713-investigation-implementation/source-gate-aggregate-final/source_gate_policy.json"
    if not policy_path.is_absolute():
        policy_path = root / policy_path
    prior_path = Path(prior_models) if prior_models else root / "tmp/agent/20260712-multichain-full-comparison/scored_models.csv"
    if not prior_path.is_absolute():
        prior_path = root / prior_path
    rows, policy = select_rows(root, policy_path, cohort, prior_path)
    ranking_keys = tuple(DEFAULT_RANKING_KEYS)
    validate_native_independent_ranking_keys(ranking_keys)
    template_path = root / "new_template/template/final_list.txt"
    templates = [line.strip() for line in template_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(templates) < 10 or len(set(templates)) != len(templates):
        raise ValueError(f"expected at least 10 unique templates, found {len(templates)}")
    panel_label = f"matched_{len(templates)}"
    template_manifest = output / "template_manifest.tsv"
    write_tsv(template_manifest, ("template_id", "template_panel", "template_list_sha256"), [
        {"template_id": template, "template_panel": panel_label, "template_list_sha256": sha256_file(template_path)}
        for template in templates
    ])
    cohort_manifest = output / "cohort_manifest.tsv"
    write_tsv(cohort_manifest, COHORT_FIELDS, rows)
    pipeline_inputs = output / "pipeline_inputs.csv"
    pipeline_inputs.parent.mkdir(parents=True, exist_ok=True)
    with pipeline_inputs.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PIPELINE_INPUT_FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "pair_id": _pair_id(row),
                "dataset_row_id": row["dataset_row_id"],
                "dataset": row["dataset"],
                "difficulty": row["difficulty"],
                "source_row_number": row["source_row_number"],
                "native_complex": row["native_complex"],
                "Receptor": row["receptor"],
                "Ligand": row["ligand"],
                "chain_context": row["chain_context"],
                "source_policy": row["source_policy"],
                "audit_only": row["audit_only"],
            })
    ranking_text = ",".join(ranking_keys)
    arms = [
        {"arm_id": "tmalign_external_rosetta", "label": "TM-align -> external Rosetta", "aligner": "TM-align", "refiner": "external_rosetta", "filter_profile": "production", "ranking_keys": ranking_text, "template_workflow": "json", "template_panel": "matched_946", "template_manifest": str(template_manifest), "tm_score_threshold": "0.4", "status": "planned", "scope": "primary"},
        {"arm_id": "gtalign_external_rosetta", "label": "GTalign CPU -> external Rosetta", "aligner": "GTalign CPU", "refiner": "external_rosetta", "filter_profile": "production", "ranking_keys": ranking_text, "template_workflow": "json", "template_panel": "matched_946", "template_manifest": str(template_manifest), "tm_score_threshold": "0.4", "status": "planned", "scope": "primary_diagnostic"},
        {"arm_id": "tmalign_pyrosetta", "label": "TM-align canonical poses -> seeded PyRosetta", "aligner": "TM-align", "refiner": "pyrosetta", "filter_profile": "production", "ranking_keys": ranking_text, "template_workflow": "json", "template_panel": "matched_946", "template_manifest": str(template_manifest), "tm_score_threshold": "0.4", "status": "planned", "scope": "primary_diagnostic"},
        {"arm_id": "legacy_fiberdock", "label": "Legacy MultiProt -> FiberDock", "aligner": "MultiProt", "refiner": "FiberDock", "filter_profile": "legacy_native", "ranking_keys": ranking_text, "template_workflow": "legacy", "template_panel": "legacy_native", "template_manifest": "", "tm_score_threshold": "0.4", "status": "blocked_pending_runtime", "scope": "diagnostic"},
    ]
    arm_manifest = output / "arm_manifest.tsv"
    write_tsv(arm_manifest, ARM_FIELDS, arms)
    tasks: list[dict[str, Any]] = []
    for row in rows:
        for arm in arms:
            for experiment in ("candidate_generation", "top20_refinement", "strict_scoring"):
                task_id = f"{row['dataset_row_id']}:{arm['arm_id']}:{experiment}:attempt-0"
                tasks.append({
                    "task_id": task_id, "dataset_row_id": row["dataset_row_id"], "arm": arm["arm_id"],
                    "experiment": experiment, "aligner": arm["aligner"], "refiner": arm["refiner"],
                    "template_workflow": arm["template_workflow"], "template_panel": arm["template_panel"],
                    "tm_score_threshold": "0.4", "top_k": "20", "input_manifest": str(cohort_manifest),
                    "output_root": str(output / "tasks" / hashlib.sha256(task_id.encode()).hexdigest()[:16]),
                    "status": "planned",
                })
    task_manifest = output / "task_manifest.tsv"
    write_tsv(task_manifest, TASK_FIELDS, tasks)
    plan = {
        "schema_version": "prism-matched-benchmark/v1",
        "cohort": cohort,
        "row_count": len(rows),
        "strict_source_clean_row_count": int(policy["decision"]["strict_source_clean_row_count"]),
        "audit_only_row_count": int(policy["decision"]["audit_only_row_count"]),
        "excluded_dataset_row_ids": policy["decision"]["excluded_dataset_row_ids"],
        "template_panel": "matched_946",
        "template_count": len(templates),
        "template_list_sha256": sha256_file(template_path),
        "template_workflows": {"json": "new_template/template", "legacy": "template_old/template"},
        "tm_score_threshold": 0.4,
        "tm_score_sensitivity_thresholds": [0.30, 0.35],
        "ranking_keys": ranking_keys,
        "top_k": 20,
        "seed": "-mute all -constant_seed -jran 12345",
        "arms": [arm["arm_id"] for arm in arms],
        "primary_evaluator": "strict_no_align_after_mapping_validation",
        "diagnostic_evaluator": "alignment_enabled_only_for_failure_diagnosis",
        "source_policy_sha256": sha256_file(policy_path),
    }
    plan_path = output / "analysis_plan.json"
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"cohort_manifest": cohort_manifest, "pipeline_inputs": pipeline_inputs, "template_manifest": template_manifest, "arm_manifest": arm_manifest, "task_manifest": task_manifest, "analysis_plan": plan_path}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cohort", choices=("pilot", "strict-clean", "all"), default="pilot")
    parser.add_argument("--source-policy", type=Path)
    parser.add_argument("--prior-models", type=Path)
    args = parser.parse_args(argv)
    try:
        paths = build_manifests(args.repo_root, args.output_dir, cohort=args.cohort, source_policy=args.source_policy, prior_models=args.prior_models)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    print(json.dumps({key: str(value) for key, value in paths.items()}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
