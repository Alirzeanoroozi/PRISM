#!/usr/bin/env python3
"""Join transformed and common-refined compact PRISM comparison records.

This reducer is intentionally score-table based: it never opens a structure
or copies a run tree.  Candidate identity is the frozen case/template/
orientation/query/chain tuple.  Refinement deltas are emitted only when the
same identity has both a transformed GlobalDockQ and a primary refined
GlobalDockQ.  Missing and failed values remain explicit and are never mapped
to zero.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


IDENTITY_FIELDS = (
    "case_id", "template", "orientation", "query_left", "query_right",
    "chain_left", "chain_right",
)
REFINABLE_STATUSES = {"scored", "scored_cross_only", "valid_unscored"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def number(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def identity(row: dict[str, Any]) -> tuple[str, ...]:
    return tuple(str(row.get(field, row.get("case", "")) or "") for field in IDENTITY_FIELDS)


def identity_without_chains(row: dict[str, Any]) -> tuple[str, ...]:
    return tuple(str(row.get(field, row.get("case", "")) or "") for field in IDENTITY_FIELDS[:5])


def status(row: dict[str, Any]) -> str:
    return str(row.get("score_status", row.get("source_status", row.get("status", ""))) or "")


def global_score(row: dict[str, Any], prefix: str = "") -> float | None:
    names = (
        f"{prefix}global_dockq", f"{prefix}dockq_global", f"{prefix}GlobalDockQ",
        f"{prefix}score_dockq_global", f"{prefix}dockq",
    )
    for name in names:
        value = number(row.get(name))
        if value is not None:
            return value
    return None


def cross_score(row: dict[str, Any], prefix: str = "") -> float | None:
    for name in (f"{prefix}cross_best", f"{prefix}dockq_cross_best", f"{prefix}dockq_cross_mean"):
        value = number(row.get(name))
        if value is not None:
            return value
    return None


def primary_refined_score(row: dict[str, Any]) -> tuple[float | None, str, str]:
    """Return external-Rosetta quality, with FiberDock only as a fallback arm."""
    external = global_score(row, "external_rosetta_")
    external_status = str(row.get("external_rosetta_status", "") or "")
    if external is not None or external_status:
        return external, external_status, "external_rosetta"
    generic = global_score(row, "refined_")
    generic_status = str(row.get("refined_status", "") or "")
    if generic is not None or generic_status:
        return generic, generic_status, "refined"
    fiber = global_score(row, "fiberdock_")
    fiber_status = str(row.get("fiberdock_status", "") or "")
    return fiber, fiber_status, "fiberdock_fallback"


def quality(value: float | None) -> str:
    if value is None:
        return "unscored"
    if value < 0.23:
        return "incorrect"
    if value < 0.49:
        return "acceptable"
    if value < 0.80:
        return "medium"
    return "high"


def ranking_value(row: dict[str, Any]) -> float | None:
    for field in ("ranking_score", "confidence_score", "tm_score_mean", "confidence_tm_mean"):
        value = number(row.get(field))
        if value is not None:
            return value
    left = number(row.get("tm_score_left"))
    right = number(row.get("tm_score_right"))
    if left is not None and right is not None:
        return (left + right) / 2.0
    rank = number(row.get("selection_rank"))
    return None if rank is None else -rank


def prodigy_value(row: dict[str, Any]) -> float | None:
    """Return a retained PRODIGY affinity, where lower is better."""

    for field in (
        "prodigy_affinity_kcal_mol", "prodigy_affinity", "affinity_kcal_mol",
        "prodigy_score",
    ):
        value = number(row.get(field))
        if value is not None:
            return value
    return None


def effective_dockq(row: dict[str, Any]) -> float | None:
    """Use primary refined quality when present, otherwise transformed quality."""

    refined = number(row.get("refined_global_dockq_value"))
    if refined is not None:
        return refined
    return number(row.get("transformed_global_dockq_value"))


def selected_quality(rows: list[dict[str, Any]]) -> tuple[float | None, float | None]:
    """Return the first selected score and the best score in the selection."""

    values = [value for value in (effective_dockq(row) for row in rows) if value is not None]
    first = effective_dockq(rows[0]) if rows else None
    return first, max(values) if values else None


def merge_maps(rows: list[dict[str, str]]) -> tuple[dict[tuple[str, ...], dict[str, str]], dict[tuple[str, ...], dict[str, str]]]:
    exact: dict[tuple[str, ...], dict[str, str]] = {}
    short: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = identity(row)
        exact[key] = row
        short[identity_without_chains(row)] = row
    return exact, short


def lookup(
    row: dict[str, Any],
    exact: dict[tuple[str, ...], dict[str, str]],
    short: dict[tuple[str, ...], dict[str, str]],
) -> dict[str, str] | None:
    exact_row = exact.get(identity(row))
    if exact_row is not None:
        return exact_row
    # A short fallback is safe only when both records omit chain labels.  Do
    # not pair two candidates that agree on query/template identity but have
    # different declared chain mappings.
    if any(str(row.get(field, "") or "") for field in IDENTITY_FIELDS[-2:]):
        return None
    short_row = short.get(identity_without_chains(row))
    if short_row is not None and any(str(short_row.get(field, "") or "") for field in IDENTITY_FIELDS[-2:]):
        return None
    return short_row


def write_table(path: Path, rows: list[dict[str, Any]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row}) or ["status"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return sha256_file(path)


def confidence_interval(values: list[float]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    mean = statistics.mean(values)
    if len(values) < 2:
        return mean, mean
    half = 1.96 * statistics.stdev(values) / math.sqrt(len(values))
    return mean - half, mean + half


def aggregate(
    *,
    output_dir: Path,
    candidate_specs: list[tuple[str, Path]],
    transformed_specs: list[tuple[str, Path]],
    refined_specs: list[tuple[str, Path]],
    case_list: Path | None,
) -> dict[str, Any]:
    candidate_map = dict(candidate_specs)
    transformed_map = dict(transformed_specs)
    refined_map = dict(refined_specs)
    methods = sorted(set(candidate_map) | set(transformed_map) | set(refined_map))
    all_rows: list[dict[str, Any]] = []
    input_records: list[dict[str, Any]] = []

    for method in methods:
        candidates = read_table(candidate_map[method]) if method in candidate_map else []
        transformed = read_table(transformed_map[method]) if method in transformed_map else []
        refined = read_table(refined_map[method]) if method in refined_map else []
        transformed_exact, transformed_short = merge_maps(transformed)
        refined_exact, refined_short = merge_maps(refined)
        rows_by_key: dict[tuple[str, ...], dict[str, Any]] = {}
        for source in candidates + transformed + refined:
            key = identity(source)
            rows_by_key.setdefault(key, {"method": method, **source})
        for key, base in rows_by_key.items():
            transformed_row = lookup(base, transformed_exact, transformed_short)
            refined_row = lookup(base, refined_exact, refined_short)
            if transformed_row:
                base.update({f"transformed_{field}": value for field, value in transformed_row.items()})
            if refined_row:
                base.update({f"refined_{field}": value for field, value in refined_row.items()})
            transformed_value = global_score(transformed_row or {})
            transformed_cross = cross_score(transformed_row or {})
            refined_value, refined_status, refined_backend = primary_refined_score(refined_row or {})
            transformed_status = status(transformed_row or {})
            if transformed_row and not transformed_status:
                transformed_status = "score_present"
            delta = ""
            paired_status = "unpaired"
            if transformed_value is not None and refined_value is not None:
                delta = refined_value - transformed_value
                paired_status = "paired"
            base.update(
                {
                    "candidate_key": "|".join(key),
                    "transformed_global_dockq_value": "" if transformed_value is None else transformed_value,
                    "transformed_cross_dockq_value": "" if transformed_cross is None else transformed_cross,
                    "transformed_score_status": transformed_status or "not_available",
                    "refined_global_dockq_value": "" if refined_value is None else refined_value,
                    "refined_score_status": refined_status or "not_available",
                    "refined_backend": refined_backend if refined_row else "not_available",
                    "paired_refinement_status": paired_status,
                    "delta_dockq": delta,
                    "ranking_value": ranking_value(base),
                    "transformed_quality_category": quality(transformed_value),
                    "refined_quality_category": quality(refined_value),
                }
            )
            all_rows.append(base)
        input_records.extend(
            {"method": method, "kind": kind, "path": str(path.resolve()), "sha256": sha256_file(path)}
            for kind, path in (("candidate", candidate_map.get(method)), ("transformed", transformed_map.get(method)), ("refined", refined_map.get(method)))
            if path is not None
        )

    expected: dict[str, str] = {}
    if case_list:
        for row in read_table(case_list):
            case_id = str(row.get("case_id", row.get("case", "")) or "")
            if case_id:
                expected[case_id] = str(row.get("split", "") or "")
    for row in all_rows:
        case_id = str(row.get("case_id", row.get("case", "")) or "")
        if case_id:
            expected.setdefault(case_id, str(row.get("split", "") or ""))

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in all_rows:
        grouped[(str(row.get("method", "")), str(row.get("case_id", row.get("case", ""))))].append(row)
    case_rows: list[dict[str, Any]] = []
    ranking_rows: list[dict[str, Any]] = []
    paired_rows = [row for row in all_rows if row.get("paired_refinement_status") == "paired"]
    for method in methods:
        for case_id, split in expected.items():
            rows = grouped.get((method, case_id), [])
            transformed_values = [number(row.get("transformed_global_dockq_value")) for row in rows]
            transformed_values = [v for v in transformed_values if v is not None]
            transformed_cross_values = [number(row.get("transformed_cross_dockq_value")) for row in rows]
            transformed_cross_values = [v for v in transformed_cross_values if v is not None]
            refined_values = [number(row.get("refined_global_dockq_value")) for row in rows]
            refined_values = [v for v in refined_values if v is not None]
            refined_cross_values = [number(row.get("refined_cross_dockq_value")) for row in rows]
            refined_cross_values = [v for v in refined_cross_values if v is not None]
            transformed_irmsd = [number(row.get("transformed_cross_irmsd_mean")) for row in rows]
            transformed_irmsd = [v for v in transformed_irmsd if v is not None]
            transformed_lrmsd = [number(row.get("transformed_cross_lrmsd_mean")) for row in rows]
            transformed_lrmsd = [v for v in transformed_lrmsd if v is not None]
            refined_irmsd = [number(row.get("refined_cross_irmsd_mean")) for row in rows]
            refined_irmsd = [v for v in refined_irmsd if v is not None]
            refined_lrmsd = [number(row.get("refined_cross_lrmsd_mean")) for row in rows]
            refined_lrmsd = [v for v in refined_lrmsd if v is not None]
            deltas = [number(row.get("delta_dockq")) for row in rows]
            deltas = [v for v in deltas if v is not None]
            lo, hi = confidence_interval(deltas)
            generated = [row for row in rows if str(row.get("status", "")) == "generated"]
            candidate_count = len(generated) if generated else len(rows)
            transformed_count = sum(
                str(row.get("transformed_score_status", "not_available")) != "not_available"
                for row in rows
            )
            refined_count = sum(
                str(row.get("refined_score_status", "not_available")) != "not_available"
                for row in rows
            )
            case_rows.append(
                {
                    "method": method,
                    "case_id": case_id,
                    "split": split,
                    "candidate_count": candidate_count,
                    "aligned_candidate_count": candidate_count,
                    "transformed_count": transformed_count,
                    "transformed_scoreable_count": len(transformed_values),
                    "refined_count": refined_count,
                    "refined_scoreable_count": len(refined_values),
                    "paired_candidate_count": len(deltas),
                    "best_transformed_global_dockq": max(transformed_values) if transformed_values else "",
                    "best_transformed_cross_dockq": max(transformed_cross_values) if transformed_cross_values else "",
                    "mean_transformed_global_dockq": statistics.mean(transformed_values) if transformed_values else "",
                    "median_transformed_global_dockq": statistics.median(transformed_values) if transformed_values else "",
                    "mean_transformed_cross_dockq": statistics.mean(transformed_cross_values) if transformed_cross_values else "",
                    "best_refined_global_dockq": max(refined_values) if refined_values else "",
                    "best_refined_cross_dockq": max(refined_cross_values) if refined_cross_values else "",
                    "mean_refined_global_dockq": statistics.mean(refined_values) if refined_values else "",
                    "median_refined_global_dockq": statistics.median(refined_values) if refined_values else "",
                    "mean_refined_cross_dockq": statistics.mean(refined_cross_values) if refined_cross_values else "",
                    "mean_transformed_cross_irmsd": statistics.mean(transformed_irmsd) if transformed_irmsd else "",
                    "mean_transformed_cross_lrmsd": statistics.mean(transformed_lrmsd) if transformed_lrmsd else "",
                    "mean_refined_cross_irmsd": statistics.mean(refined_irmsd) if refined_irmsd else "",
                    "mean_refined_cross_lrmsd": statistics.mean(refined_lrmsd) if refined_lrmsd else "",
                    "best_transformed_quality_category": quality(max(transformed_values) if transformed_values else None),
                    "best_refined_quality_category": quality(max(refined_values) if refined_values else None),
                    "mean_delta_dockq": statistics.mean(deltas) if deltas else "",
                    "median_delta_dockq": statistics.median(deltas) if deltas else "",
                    "delta_ci95_low": "" if lo is None else lo,
                    "delta_ci95_high": "" if hi is None else hi,
                    "delta_fraction_improved": (sum(v > 0 for v in deltas) / len(deltas)) if deltas else "",
                    "delta_fraction_worsened": (sum(v < 0 for v in deltas) / len(deltas)) if deltas else "",
                    "coverage_status": "no_candidate" if not rows else "candidate_present",
                }
            )
            ranked = sorted(
                [row for row in rows if number(row.get("ranking_value")) is not None],
                key=lambda row: (float(row["ranking_value"]), str(row.get("candidate_key", ""))),
                reverse=True,
            )
            all_candidates = generated or rows
            prodigy_ranked = sorted(
                [row for row in rows if prodigy_value(row) is not None],
                key=lambda row: (float(prodigy_value(row)), str(row.get("candidate_key", ""))),
            )
            oracle_rows = sorted(
                [row for row in rows if effective_dockq(row) is not None],
                key=lambda row: (float(effective_dockq(row)), str(row.get("candidate_key", ""))),
                reverse=True,
            )
            oracle = effective_dockq(oracle_rows[0]) if oracle_rows else None
            strategy_specs = (
                ("all", all_candidates, "no_ranking"),
                ("deterministic", ranked[:1], "deterministic_baseline"),
                ("prodigy", prodigy_ranked[:1], "prodigy_affinity_ascending"),
                ("top1", ranked[:1], "deterministic_baseline_top1"),
                ("top3", ranked[:3], "deterministic_baseline_top3"),
                ("top5", ranked[:5], "deterministic_baseline_top5"),
                ("oracle", oracle_rows[:1], "diagnostic_oracle_upper_bound"),
            )
            for strategy, selected, interpretation in strategy_specs:
                value, best_selected = selected_quality(selected)
                selection_status = "selected" if selected else (
                    "not_available" if strategy == "prodigy" else "no_rankable_candidate"
                )
                ranking_rows.append(
                    {
                        "method": method,
                        "case_id": case_id,
                        "split": split,
                        "strategy": strategy,
                        "k": len(selected),
                        "candidate_count": candidate_count,
                        "selected_count": len(selected),
                        "selected_candidate_key": "" if not selected else selected[0].get("candidate_key", ""),
                        "selected_dockq": "" if value is None else value,
                        "best_selected_dockq": "" if best_selected is None else best_selected,
                        "oracle_best_dockq_diagnostic": "" if oracle is None else oracle,
                        "ranking_regret": "" if oracle is None or value is None else oracle - value,
                        "good_model_retained": "" if best_selected is None else best_selected >= 0.23,
                        "candidates_avoided": max(0, candidate_count - len(selected)),
                        "refinement_jobs_avoided": max(0, candidate_count - len(selected)),
                        "selection_status": selection_status,
                        "interpretation": interpretation,
                    }
                )

    method_summary_rows: list[dict[str, Any]] = []
    split_values = sorted({str(row.get("split", "") or "") for row in case_rows})
    for method in methods:
        for split in ["all", *split_values]:
            scoped = [
                row for row in case_rows
                if row["method"] == method and (split == "all" or row["split"] == split)
            ]
            candidate_counts = [int(row["candidate_count"]) for row in scoped]
            transformed_best = [number(row.get("best_transformed_global_dockq")) for row in scoped]
            transformed_best = [value for value in transformed_best if value is not None]
            refined_best = [number(row.get("best_refined_global_dockq")) for row in scoped]
            refined_best = [value for value in refined_best if value is not None]
            paired = [number(row.get("mean_delta_dockq")) for row in scoped]
            paired = [value for value in paired if value is not None]
            summary: dict[str, Any] = {
                "method": method,
                "split": split,
                "expected_cases": len(scoped),
                "cases_with_candidates": sum(row["coverage_status"] == "candidate_present" for row in scoped),
                "cases_transformed": sum(int(row["transformed_count"]) > 0 for row in scoped),
                "cases_scoreable": sum(int(row["transformed_scoreable_count"]) > 0 for row in scoped),
                "cases_refined": sum(int(row["refined_count"]) > 0 for row in scoped),
                "candidate_total": sum(candidate_counts),
                "candidate_median_per_case": statistics.median(candidate_counts) if candidate_counts else "",
                "mean_case_best_transformed_global_dockq": statistics.mean(transformed_best) if transformed_best else "",
                "median_case_best_transformed_global_dockq": statistics.median(transformed_best) if transformed_best else "",
                "mean_case_best_refined_global_dockq": statistics.mean(refined_best) if refined_best else "",
                "median_case_best_refined_global_dockq": statistics.median(refined_best) if refined_best else "",
                "paired_case_count": len(paired),
                "mean_case_delta_dockq": statistics.mean(paired) if paired else "",
                "median_case_delta_dockq": statistics.median(paired) if paired else "",
                "best_transformed_incorrect_cases": sum(row["best_transformed_quality_category"] == "incorrect" for row in scoped),
                "best_transformed_acceptable_cases": sum(row["best_transformed_quality_category"] == "acceptable" for row in scoped),
                "best_transformed_medium_cases": sum(row["best_transformed_quality_category"] == "medium" for row in scoped),
                "best_transformed_high_cases": sum(row["best_transformed_quality_category"] == "high" for row in scoped),
                "best_refined_incorrect_cases": sum(row["best_refined_quality_category"] == "incorrect" for row in scoped),
                "best_refined_acceptable_cases": sum(row["best_refined_quality_category"] == "acceptable" for row in scoped),
                "best_refined_medium_cases": sum(row["best_refined_quality_category"] == "medium" for row in scoped),
                "best_refined_high_cases": sum(row["best_refined_quality_category"] == "high" for row in scoped),
            }
            for strategy in ("all", "deterministic", "prodigy", "top1", "top3", "top5", "oracle"):
                ranked_scope = [
                    row for row in ranking_rows
                    if row["method"] == method and row["strategy"] == strategy
                    and (split == "all" or row["split"] == split)
                    and row["selection_status"] == "selected"
                ]
                selected = [number(row.get("selected_dockq")) for row in ranked_scope]
                selected = [value for value in selected if value is not None]
                retained = [bool(row.get("good_model_retained")) for row in ranked_scope if row.get("good_model_retained") != ""]
                avoided = [number(row.get("candidates_avoided")) for row in ranked_scope]
                avoided = [value for value in avoided if value is not None]
                key = strategy.replace("-", "_")
                summary[f"{key}_cases"] = len(ranked_scope)
                summary[f"{key}_mean_selected_dockq"] = statistics.mean(selected) if selected else ""
                summary[f"{key}_good_model_retention"] = (sum(retained) / len(retained)) if retained else ""
                summary[f"{key}_mean_candidates_avoided"] = statistics.mean(avoided) if avoided else ""
            method_summary_rows.append(summary)

    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_hash = write_table(output_dir / "candidate_results.tsv", all_rows)
    case_hash = write_table(output_dir / "case_results.tsv", case_rows)
    ranking_hash = write_table(output_dir / "ranking_topk.tsv", ranking_rows)
    paired_hash = write_table(output_dir / "paired_refinement.tsv", paired_rows)
    method_summary_hash = write_table(output_dir / "method_summary.tsv", method_summary_rows)
    manifest = {
        "schema_version": "prism-final-matched-comparison/v1",
        "status": "validated_compacted",
        "methods": methods,
        "expected_case_count": len(expected),
        "candidate_rows": len(all_rows),
        "case_rows": len(case_rows),
        "paired_rows": len(paired_rows),
        "ranking_rows": len(ranking_rows),
        "method_summary_rows": len(method_summary_rows),
        "candidate_results_sha256": candidate_hash,
        "case_results_sha256": case_hash,
        "ranking_topk_sha256": ranking_hash,
        "paired_refinement_sha256": paired_hash,
        "method_summary_sha256": method_summary_hash,
        "inputs": input_records,
        "primary_refinement_backend": "external_rosetta_when_available",
        "ranking_strategies": ["all", "deterministic", "prodigy", "top1", "top3", "top5", "oracle"],
        "prodigy_contract": "lower retained affinity is better; absent PRODIGY evidence remains not_available",
        "delta_contract": "exact candidate identity with transformed and primary refined GlobalDockQ",
        "confidence_interval": "normal approximation 95% CI; n<2 reported as point value",
        "missing_scores_are_empty": True,
        "oracle_is_diagnostic_only": True,
    }
    (output_dir / "aggregate_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def parse_spec(value: str) -> tuple[str, Path]:
    method, separator, raw_path = value.partition("=")
    if not separator or not method or not raw_path:
        raise argparse.ArgumentTypeError("expected METHOD=PATH")
    return method, Path(raw_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--candidate", action="append", type=parse_spec, default=[])
    parser.add_argument("--transformed", action="append", type=parse_spec, default=[])
    parser.add_argument("--refined", action="append", type=parse_spec, default=[])
    parser.add_argument("--case-list", type=Path)
    args = parser.parse_args()
    summary = aggregate(
        output_dir=args.output_dir.resolve(),
        candidate_specs=[(method, path.resolve()) for method, path in args.candidate],
        transformed_specs=[(method, path.resolve()) for method, path in args.transformed],
        refined_specs=[(method, path.resolve()) for method, path in args.refined],
        case_list=args.case_list.resolve() if args.case_list else None,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
