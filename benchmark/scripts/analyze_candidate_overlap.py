#!/usr/bin/env python3
"""Compute case-wise candidate overlap from compact generated-candidate tables."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path


def read_candidates(path: Path, method: str) -> dict[str, set[tuple[str, ...]]]:
    cases: dict[str, set[tuple[str, ...]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
        for row in csv.DictReader(handle, delimiter=delimiter):
            declared_method = (row.get("pipeline") or row.get("aligner") or "").strip().lower()
            if declared_method and declared_method not in {method.strip().lower(), "usalign" if method.lower() == "usalign" else method.strip().lower()}:
                continue
            status = row.get("status") or row.get("source_status") or row.get("score_status")
            if status not in {"generated", "scored", "score_failed", "valid_unscored"}:
                continue
            case = row.get("pair_id") or row.get("case_id")
            if not case:
                # A compact provider table without a durable case key cannot
                # support case-wise overlap; fail closed rather than pooling.
                continue
            signature = (
                row.get("template", ""),
                row.get("query_left", ""),
                row.get("query_right", ""),
                row.get("orientation", ""),
                row.get("chain_left", ""),
                row.get("chain_right", ""),
            )
            cases.setdefault(case, set()).add(signature)
    return cases


def analyze(method_paths: list[str], output_dir: Path) -> dict[str, object]:
    methods: dict[str, dict[str, set[tuple[str, ...]]]] = {}
    for spec in method_paths:
        method, raw_path = spec.split("=", 1)
        path = Path(raw_path).resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        methods[method] = read_candidates(path, method)

    cases = sorted(set().union(*(set(records) for records in methods.values())))
    pair_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for left, right in itertools.combinations(sorted(methods), 2):
        jaccards: list[float] = []
        shared_total = unique_left_total = unique_right_total = 0
        for case in cases:
            a = methods[left].get(case, set())
            b = methods[right].get(case, set())
            shared = a & b
            union = a | b
            jaccard = len(shared) / len(union) if union else None
            if jaccard is not None:
                jaccards.append(jaccard)
            shared_total += len(shared)
            unique_left_total += len(a - b)
            unique_right_total += len(b - a)
            pair_rows.append({
                "method_left": left,
                "method_right": right,
                "case_id": case,
                "left_candidates": len(a),
                "right_candidates": len(b),
                "shared_candidates": len(shared),
                "union_candidates": len(union),
                "jaccard": jaccard if jaccard is not None else "",
                "unique_left": len(a - b),
                "unique_right": len(b - a),
            })
        summary_rows.append({
            "method_left": left,
            "method_right": right,
            "case_count": len(cases),
            "shared_candidates": shared_total,
            "unique_left_candidates": unique_left_total,
            "unique_right_candidates": unique_right_total,
            "mean_case_jaccard": sum(jaccards) / len(jaccards) if jaccards else "",
            "zero_overlap_cases": sum(value == 0 for value in jaccards),
        })

    coverage_rows = []
    for method, records in sorted(methods.items()):
        counts = [len(records.get(case, set())) for case in cases]
        coverage_rows.append({
            "method": method,
            "case_count_in_union": len(cases),
            "cases_with_candidates": sum(count > 0 for count in counts),
            "candidate_total": sum(counts),
            "median_candidates_per_case": sorted(counts)[len(counts) // 2] if counts else 0,
        })

    output_dir.mkdir(parents=True, exist_ok=True)

    def write(name: str, rows: list[dict[str, object]]) -> str:
        path = output_dir / name
        fields = list(rows[0]) if rows else ["status"]
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return str(path)

    result = {
        "schema_version": "prism-candidate-overlap/v1",
        "methods": sorted(methods),
        "case_count": len(cases),
        "coverage_path": write("candidate_coverage.tsv", coverage_rows),
        "pairwise_path": write("candidate_overlap_case.tsv", pair_rows),
        "summary_path": write("candidate_overlap_summary.tsv", summary_rows),
    }
    (output_dir / "candidate_overlap_manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", action="append", required=True, help="METHOD=compact_candidate.csv")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(analyze(args.candidate, args.output_dir), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
