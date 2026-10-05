#!/usr/bin/env python3
"""Build compact case/candidate/ranking tables for the matched aligner lane.

The reducer deliberately accepts already compact candidate and score tables.
It does not inspect or copy structures, and missing/failed scores remain empty
with explicit status rather than becoming zero.  ``ranking_topk.tsv`` reports
the deterministic score-based selection when a ranking score is available and
an oracle diagnostic based on the best observed DockQ; the latter is not a
runtime ranking claim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        except csv.Error:
            dialect = csv.excel_tab if delimiter == "\t" else csv.excel
        dialect.delimiter = delimiter
        return list(csv.DictReader(handle, dialect=dialect))


def number(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def first_number(row: dict[str, Any], names: tuple[str, ...]) -> float | None:
    for name in names:
        value = number(row.get(name))
        if value is not None:
            return value
    return None


def candidate_key(row: dict[str, Any]) -> tuple[str, ...]:
    return (
        str(row.get("case_id", row.get("case", ""))),
        str(row.get("template", "")),
        str(row.get("orientation", "")),
        str(row.get("query_left", "")),
        str(row.get("query_right", "")),
    )


def score_key(row: dict[str, Any]) -> tuple[str, ...]:
    key = candidate_key(row)
    if any(key):
        return key
    return (str(row.get("candidate_index", row.get("manifest_index", ""))),)


def short_candidate_key(row: dict[str, Any]) -> tuple[str, ...]:
    return (
        str(row.get("case_id", row.get("case", ""))),
        str(row.get("template", "")),
        str(row.get("orientation", "")),
    )


def status(row: dict[str, Any]) -> str:
    return str(row.get("score_status", row.get("source_status", row.get("status", ""))) or "")


def matches_method(row: dict[str, Any], method: str) -> bool:
    value = str(row.get("method", row.get("pipeline", row.get("aligner", "")))).strip().lower()
    normalized = method.strip().lower()
    if not value:
        return True
    aliases = {
        "tmalign": {"tmalign", "tm-align", "tm_align"},
        "multiprot": {"multiprot", "multi-prot", "multi_prot"},
        "usalign": {"usalign", "u-salign", "u_salign"},
        "gtalign": {"gtalign", "gt-align", "gt_align"},
    }
    return value in aliases.get(normalized, {normalized})


def global_dockq(row: dict[str, Any]) -> float | None:
    return first_number(
        row,
        (
            "dockq_global",
            "GlobalDockQ",
            "score_dockq_global",
            "external_rosetta_global_dockq",
            "fiberdock_global_dockq",
        ),
    )


def cross_dockq(row: dict[str, Any]) -> float | None:
    return first_number(
        row,
        (
            "dockq_cross_best",
            "dockq_cross_mean",
            "score_dockq_cross_best",
            "external_rosetta_cross_best",
            "fiberdock_cross_best",
        ),
    )


def quality_category(value: float | None) -> str:
    if value is None:
        return "unscored"
    if value < 0.23:
        return "incorrect"
    if value < 0.49:
        return "acceptable"
    if value < 0.80:
        return "medium"
    return "high"


def ranking_score(row: dict[str, Any]) -> float | None:
    explicit = first_number(row, ("ranking_score", "prodigy_score", "confidence_score"))
    if explicit is not None:
        return explicit
    left = number(row.get("tm_score_left"))
    right = number(row.get("tm_score_right"))
    if left is not None and right is not None:
        return (left + right) / 2.0
    rank = number(row.get("selection_rank"))
    return None if rank is None else -rank


def merge_method_rows(
    method: str,
    candidates_path: Path | None,
    scores_path: Path | None,
) -> list[dict[str, Any]]:
    candidates = [
        row for row in (read_table(candidates_path) if candidates_path else [])
        if matches_method(row, method)
    ]
    scores = [
        row for row in (read_table(scores_path) if scores_path else [])
        if matches_method(row, method)
    ]
    score_map = {score_key(row): row for row in scores}
    short_score_map = {short_candidate_key(row): row for row in scores}
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, ...]] = set()
    for source in candidates:
        key = candidate_key(source)
        merged: dict[str, Any] = {"method": method, "candidate_origin": "candidate", **source}
        score = (
            score_map.get(key)
            or short_score_map.get(short_candidate_key(source))
            or score_map.get((str(source.get("candidate_index", "")),))
        )
        if score:
            merged.update({f"score_{k}": v for k, v in score.items() if k not in {"method"}})
            for field, value in score.items():
                merged.setdefault(field, value)
        merged["candidate_key"] = "|".join(key)
        merged["dockq_global_value"] = global_dockq(merged)
        merged["dockq_cross_value"] = cross_dockq(merged)
        merged["ranking_score_value"] = ranking_score(merged)
        rows.append(merged)
        seen.add(key)
    for source in scores:
        key = score_key(source)
        if key in seen:
            continue
        merged = {"method": method, "candidate_origin": "score_only", **source, "candidate_key": "|".join(key)}
        merged["status"] = source.get("source_status", source.get("score_status", ""))
        merged["dockq_global_value"] = global_dockq(merged)
        merged["dockq_cross_value"] = cross_dockq(merged)
        merged["ranking_score_value"] = ranking_score(merged)
        rows.append(merged)
    return rows


def write_table(path: Path, rows: list[dict[str, Any]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row}) or ["status"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return sha256_file(path)


def aggregate(
    *,
    output_dir: Path,
    candidate_specs: list[tuple[str, Path]],
    score_specs: list[tuple[str, Path]],
    case_list: Path | None,
) -> dict[str, Any]:
    candidate_map = dict(candidate_specs)
    score_map = dict(score_specs)
    methods = sorted(set(candidate_map) | set(score_map))
    all_rows: list[dict[str, Any]] = []
    for method in methods:
        all_rows.extend(merge_method_rows(method, candidate_map.get(method), score_map.get(method)))

    expected: dict[str, str] = {}
    if case_list:
        for row in read_table(case_list):
            case_id = str(row.get("case_id", row.get("case", "")))
            if case_id:
                expected[case_id] = str(row.get("split", ""))
    for row in all_rows:
        case_id = str(row.get("case_id", row.get("case", "")))
        if case_id:
            expected.setdefault(case_id, str(row.get("split", "")))

    case_rows: list[dict[str, Any]] = []
    ranking_rows: list[dict[str, Any]] = []
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in all_rows:
        grouped[(str(row["method"]), str(row.get("case_id", row.get("case", ""))))].append(row)
    for method in methods:
        for case_id, split in expected.items():
            rows = grouped.get((method, case_id), [])
            generated = [row for row in rows if str(row.get("status", "")) == "generated"]
            candidate_rows = (
                generated
                if any(row.get("candidate_origin") == "candidate" for row in rows)
                else [row for row in rows if row.get("candidate_origin") == "score_only"]
            )
            scored = [row for row in rows if status(row) in {"scored", "scored_cross_only", "valid_unscored"}]
            values = [row["dockq_global_value"] for row in rows if row.get("dockq_global_value") is not None]
            cross_values = [row["dockq_cross_value"] for row in rows if row.get("dockq_cross_value") is not None]
            best = max(values) if values else None
            case_row: dict[str, Any] = {
                "method": method,
                "case_id": case_id,
                "split": split,
                "candidate_count": len(candidate_rows),
                "audit_row_count": len(rows),
                "generated_count": len(generated),
                "scoreable_count": len(scored),
                "score_failed_count": sum(status(row) == "score_failed" for row in rows),
                "best_dockq_global": "" if best is None else best,
                "mean_dockq_global": "" if not values else statistics.mean(values),
                "median_dockq_global": "" if not values else statistics.median(values),
                "best_dockq_cross": "" if not cross_values else max(cross_values),
                "quality_category_best": quality_category(best),
                "coverage_status": "no_candidate" if not candidate_rows else "candidate_present",
            }
            case_rows.append(case_row)
            ordered = sorted(
                [row for row in rows if row.get("ranking_score_value") is not None],
                key=lambda row: (float(row["ranking_score_value"]), str(row.get("candidate_key", ""))),
                reverse=True,
            )
            for k in (1, 3, 5):
                selected = ordered[:k]
                selected_values = [row["dockq_global_value"] for row in selected if row.get("dockq_global_value") is not None]
                selected_value = selected_values[0] if selected_values else None
                ranking_rows.append(
                    {
                        "method": method,
                        "case_id": case_id,
                        "split": split,
                        "k": k,
                        "candidate_count": len(candidate_rows),
                        "selected_count": len(selected),
                        "selected_candidate_key": "" if not selected else selected[0].get("candidate_key", ""),
                        "selected_dockq_global": "" if selected_value is None else selected_value,
                        "oracle_best_dockq_global": "" if best is None else best,
                        "ranking_regret": "" if best is None or selected_value is None else best - selected_value,
                        "good_model_retained": "" if best is None or selected_value is None else selected_value >= 0.23,
                        "candidates_avoided": max(0, len(rows) - len(selected)),
                        "interpretation": "deterministic_score_selection; oracle is diagnostic upper bound",
                    }
                )

    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_hash = write_table(output_dir / "candidate_results.tsv", all_rows)
    case_hash = write_table(output_dir / "case_results.tsv", case_rows)
    ranking_hash = write_table(output_dir / "ranking_topk.tsv", ranking_rows)
    manifest = {
        "schema_version": "prism-matched-comparison-aggregate/v1",
        "status": "validated_compacted",
        "methods": methods,
        "expected_case_count": len(expected),
        "candidate_rows": len(all_rows),
        "case_rows": len(case_rows),
        "ranking_rows": len(ranking_rows),
        "candidate_results_sha256": candidate_hash,
        "case_results_sha256": case_hash,
        "ranking_topk_sha256": ranking_hash,
        "inputs": [
            {"path": str(path), "sha256": sha256_file(path), "method": method, "kind": kind}
            for kind, specs in (("candidate", candidate_specs), ("score", score_specs))
            for method, path in specs
        ],
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
    parser.add_argument("--score", action="append", type=parse_spec, default=[])
    parser.add_argument("--case-list", type=Path)
    args = parser.parse_args()
    summary = aggregate(
        output_dir=args.output_dir.resolve(),
        candidate_specs=[(method, path.resolve()) for method, path in args.candidate],
        score_specs=[(method, path.resolve()) for method, path in args.score],
        case_list=args.case_list.resolve() if args.case_list else None,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
