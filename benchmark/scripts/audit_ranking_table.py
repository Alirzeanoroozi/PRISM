#!/usr/bin/env python3
"""Audit canonical BM5.5 ranking features and labels before evaluation."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def _read(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _key(row: dict[str, str]) -> tuple[str, str]:
    return row.get("dataset_row_id", ""), row.get("source_model_sha256", "")


def audit(candidate_path: Path, score_path: Path) -> dict[str, object]:
    candidates = _read(candidate_path)
    all_scores = _read(score_path)
    scores = [
        row for row in all_scores
        if row.get("score_status") in {"scored", "scored_cross_only"}
        and row.get("source_gate_status") != "audit_only"
    ]
    explicit_not_scoreable_keys = {
        _key(row) for row in all_scores
        if row.get("score_status") == "not_scoreable"
        and row.get("source_gate_status") != "audit_only"
    }
    errors: list[str] = []

    candidate_keys = [_key(row) for row in candidates]
    score_keys = [_key(row) for row in scores]
    if len(candidate_keys) != len(set(candidate_keys)):
        errors.append("duplicate_candidate_identity")
    if len(score_keys) != len(set(score_keys)):
        errors.append("duplicate_score_identity")
    candidate_by_key = dict(zip(candidate_keys, candidates))
    rankable_candidate_keys = {
        key for key, row in candidate_by_key.items()
        if row.get("label_status") == "labeled"
    }
    score_key_set = set(score_keys)
    if rankable_candidate_keys != score_key_set:
        errors.append("candidate_score_identity_set_mismatch")

    checks = {
        "candidate_not_accepted": sum(
            row.get("status") != "refinement_accepted" for row in candidates
        ),
        "explicit_not_scoreable_retained": sum(
            _key(row) in explicit_not_scoreable_keys for row in candidates
        ),
        "label_not_attached": sum(
            row.get("label_status") != "labeled"
            and _key(row) not in explicit_not_scoreable_keys
            for row in candidates
        ),
        "wrong_label_metric": sum(
            row.get("label_status") == "labeled"
            and row.get("label_metric") != "dockq_cross_mean"
            for row in candidates
        ),
        "model_hash_mismatch": sum(
            _key(row) not in explicit_not_scoreable_keys
            and (
            not row.get("source_model_sha256")
            or row.get("observed_source_model_sha256") != row.get("source_model_sha256")
            or row.get("label_model_sha256") != row.get("source_model_sha256")
            )
            for row in candidates
        ),
        "missing_alignment_json_hash": sum(
            not row.get("alignment_left_sha256") or not row.get("alignment_right_sha256")
            for row in candidates
        ),
        "missing_alignment_raw_hash": sum(
            not row.get("alignment_left_raw_output_sha256")
            or not row.get("alignment_right_raw_output_sha256")
            for row in candidates
        ),
        "wrong_aligner": sum(
            row.get("alignment_left_aligner") != "GTalign"
            or row.get("alignment_right_aligner") != "GTalign"
            for row in candidates
        ),
        "alignment_not_success": sum(
            row.get("alignment_left_status") != "success"
            or row.get("alignment_right_status") != "success"
            for row in candidates
        ),
        "missing_features": sum(
            any(row.get(field, "") == "" for field in (
                "match_count_left", "match_count_right", "tm_score_left", "tm_score_right",
                "match_coverage_left", "match_coverage_right",
            ))
            for row in candidates
        ),
        "missing_template_provenance": sum(
            row.get("template_coverage_status") != "available"
            or not row.get("template_interface_sha256")
            or not row.get("template_size_left")
            or not row.get("template_size_right")
            for row in candidates
        ),
    }
    errors.extend(
        name for name, count in checks.items()
        if count and name != "explicit_not_scoreable_retained"
    )
    summary: dict[str, object] = {
        "audit_status": "passed" if not errors else "failed",
        "candidate_count": len(candidates),
        "eligible_score_count": len(scores),
        "rankable_candidate_count": len(rankable_candidate_keys),
        "retained_unrankable_candidate_count": len(candidates) - len(rankable_candidate_keys),
        "explicit_not_scoreable_score_count": len(explicit_not_scoreable_keys),
        "dataset_row_count": len({row.get("dataset_row_id", "") for row in candidates}),
        "batch_counts": dict(sorted(Counter(row.get("batch", "") for row in candidates).items())),
        "checks": checks,
        "errors": errors,
    }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidates", type=Path)
    parser.add_argument("scores", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = audit(args.candidates, args.scores)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if result["audit_status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
