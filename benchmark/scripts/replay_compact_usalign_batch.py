#!/usr/bin/env python3
"""Replay common transformation and compact a completed USalign batch.

The alignment provider is never rerun here.  Raw alignment JSON is consumed
once, compact success/status summaries and generated-candidate audit rows are
written, and only then are disposable alignment files removed.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
import json
from pathlib import Path
import os
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def parse_alignment_filename(name: str) -> tuple[str, str, str]:
    stem = Path(name).stem
    query, template, chain = stem.rsplit("_", 2)
    return query, template, chain


def compact_alignment_records(
    alignment_dir: Path,
    status_dir: Path,
    *,
    aligner: str,
) -> dict[str, Any]:
    status_dir.mkdir(parents=True, exist_ok=True)
    case_counts: dict[str, Counter[str]] = defaultdict(Counter)
    score_totals: dict[str, Counter[str]] = defaultdict(Counter)
    score_extrema: dict[str, dict[str, float]] = defaultdict(dict)
    raw_count = 0
    parse_failures = 0
    valid_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for path in sorted(alignment_dir.glob("*.json")):
        raw_count += 1
        query, _template, _chain = parse_alignment_filename(path.name)
        try:
            row = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            parse_failures += 1
            case_counts[query]["parse_failure"] += 1
            continue
        counts = case_counts[query]
        counts["records"] += 1
        counts["return_code_success" if row.get("return_code") == 0 else "return_code_failure"] += 1
        counts[f"status_{row.get('status', 'missing')}"] += 1
        # Explicit unavailable/truncated records remain in status accounting;
        # they are not measured zero-quality alignments and must not depress
        # numerical means.
        if row.get("status") != "success":
            continue
        for field in ("tm_score_query", "tm_score_ref", "match_count"):
            try:
                value = float(row.get(field))
            except (TypeError, ValueError):
                continue
            score_totals[query][f"{field}_sum"] += value
            valid_counts[query][f"{field}_count"] += 1
            score_extrema[query][f"{field}_min"] = min(
                score_extrema[query].get(f"{field}_min", value), value
            )
            score_extrema[query][f"{field}_max"] = max(
                score_extrema[query].get(f"{field}_max", value), value
            )

    summary_fields = [
        "query_id", "records", "return_code_success", "return_code_failure",
        "status_success", "status_mapping_truncated", "status_alignment_unavailable",
        "status_missing", "parse_failure",
        "tm_score_query_valid_count", "tm_score_ref_valid_count",
        "match_count_valid_count",
        "tm_score_query_mean", "tm_score_ref_mean", "match_count_mean",
        "tm_score_query_min", "tm_score_query_max", "tm_score_ref_min",
        "tm_score_ref_max", "match_count_min", "match_count_max",
    ]
    summary_path = status_dir / "alignment_case_summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        for query in sorted(case_counts):
            counts = case_counts[query]
            row = {
                field: query if field == "query_id" else counts.get(field, 0)
                for field in summary_fields
            }
            row.update({
                "tm_score_query_valid_count": valid_counts[query]["tm_score_query_count"],
                "tm_score_ref_valid_count": valid_counts[query]["tm_score_ref_count"],
                "match_count_valid_count": valid_counts[query]["match_count_count"],
            })
            for field in ("tm_score_query", "tm_score_ref", "match_count"):
                count = valid_counts[query][f"{field}_count"]
                row[f"{field}_mean"] = (
                    "" if count == 0 else score_totals[query][f"{field}_sum"] / count
                )
            row.update({
                field: score_extrema[query].get(field, "")
                for field in (
                    "tm_score_query_min", "tm_score_query_max",
                    "tm_score_ref_min", "tm_score_ref_max",
                    "match_count_min", "match_count_max",
                )
            })
            writer.writerow(row)

    result = {
        "schema_version": "prism-usalign-alignment-compaction/v1",
        "aligner": aligner,
        "raw_record_count": raw_count,
        "case_count": len(case_counts),
        "parse_failure_count": parse_failures,
        "summary_path": str(summary_path),
    }
    (status_dir / "alignment_compact_manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    return result


def compact_candidate_audit(
    audit_path: Path,
    status_dir: Path,
    *,
    aligner: str,
    input_rows: list[dict[str, str]] | None = None,
) -> dict[str, Any]:
    generated_path = status_dir / "candidate_generated.csv"
    rejection_path = status_dir / "candidate_rejection_summary.csv"
    input_index = {
        (row.get("Receptor", "").strip().lower(), row.get("Ligand", "").strip().lower()): row
        for row in (input_rows or [])
    }
    generated_fields = [
        "aligner", "pair_id", "benchmark_set", "source_row", "complex",
        "query_left", "query_right", "template", "chain_left",
        "chain_right", "orientation", "status", "match_count_left",
        "match_count_right", "tm_score_left", "tm_score_right",
        "tm_score_query_left", "tm_score_ref_left", "tm_score_contract_left",
        "tm_score_query_right", "tm_score_ref_right", "tm_score_contract_right",
        "match_coverage_left", "match_coverage_right", "contact_count",
        "clash_count", "error_reason", "metadata",
    ]
    rejection_counts: Counter[tuple[str, str, str, str]] = Counter()
    generated_count = 0
    total_count = 0
    with generated_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=generated_fields)
        writer.writeheader()
        if audit_path.is_file():
            with audit_path.open() as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    row = json.loads(line)
                    metadata = row.get("metadata", {}) or {}
                    total_count += 1
                    key = (
                        row.get("query_left", ""), row.get("query_right", ""),
                        row.get("orientation", ""), row.get("status", "missing"),
                    )
                    rejection_counts[key] += 1
                    if row.get("status") == "generated":
                        generated_count += 1
                        input_row = input_index.get(
                            (
                                str(row.get("query_left", "")).strip().lower(),
                                str(row.get("query_right", "")).strip().lower(),
                            ),
                            {},
                        )
                        writer.writerow({
                            field: (aligner if field == "aligner" else
                                    input_row.get(field, "")
                                    if field in {"pair_id", "benchmark_set", "source_row", "complex"} else
                                    metadata.get(field)
                                    if field in {
                                        "tm_score_query_left", "tm_score_ref_left",
                                        "tm_score_contract_left", "tm_score_query_right",
                                        "tm_score_ref_right", "tm_score_contract_right",
                                    } else
                                    json.dumps(metadata, sort_keys=True)
                                    if field == "metadata" else row.get(field))
                            for field in generated_fields
                        })

    rejection_fields = [
        "pair_id", "benchmark_set", "source_row", "complex",
        "query_left", "query_right", "orientation", "status", "count",
    ]
    with rejection_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=rejection_fields)
        writer.writeheader()
        for (query_left, query_right, orientation, status), count in sorted(rejection_counts.items()):
            input_row = input_index.get(
                (query_left.strip().lower(), query_right.strip().lower()), {}
            )
            writer.writerow({
                "pair_id": input_row.get("pair_id", ""),
                "benchmark_set": input_row.get("benchmark_set", ""),
                "source_row": input_row.get("source_row", ""),
                "complex": input_row.get("complex", ""),
                "query_left": query_left, "query_right": query_right,
                "orientation": orientation, "status": status, "count": count,
            })
    return {
        "audit_path": str(audit_path),
        "audit_record_count": total_count,
        "generated_record_count": generated_count,
        "generated_path": str(generated_path),
        "rejection_summary_path": str(rejection_path),
    }


def delete_raw_alignment_files(alignment_root: Path, run_root: Path) -> int:
    deleted = 0
    for path in sorted(alignment_root.rglob("*.json")):
        resolved = path.resolve()
        resolved.relative_to(run_root.resolve())
        path.unlink()
        deleted += 1
    for directory in sorted(alignment_root.rglob("*"), reverse=True):
        if directory.is_dir() and not directory.is_symlink():
            try:
                directory.rmdir()
            except OSError:
                pass
    return deleted


def replay_and_compact(run_root: Path) -> dict[str, Any]:
    status_dir = run_root / "status"
    exit_path = status_dir / "exit.json"
    if not exit_path.is_file():
        return {"status": "skipped", "reason": "missing_exit_json"}
    exit_record = json.loads(exit_path.read_text())
    if int(exit_record.get("return_code", 1)) != 0:
        return {"status": "skipped", "reason": "pipeline_return_code_nonzero"}

    alignment_root = run_root / "processed" / "alignment_usalign"
    alignment_dirs = [path for path in alignment_root.iterdir() if path.is_dir()] if alignment_root.is_dir() else []
    alignment_dir = next((path for path in sorted(alignment_dirs) if any(path.glob("*.json"))), None)
    if alignment_dir is None:
        return {"status": "skipped", "reason": "no_raw_alignment_json"}

    os.chdir(run_root)
    # Import after the run's environment and working directory are set so the
    # common transformer reads the same thresholds/assets as the production arm.
    from src.transformation import transformer

    templates = [line.strip() for line in Path("templates/calculated_templates.txt").read_text().splitlines() if line.strip()]
    audit_path = status_dir / "candidate_audit_replay.jsonl"
    passed_pairs = transformer(
        templates,
        alignment_dir=str(alignment_dir),
        inputs_csv="inputs.csv",
        audit_path=str(audit_path),
    )
    input_rows = []
    inputs_path = Path("inputs.csv")
    if inputs_path.is_file():
        with inputs_path.open(newline="", encoding="utf-8") as handle:
            input_rows = list(csv.DictReader(handle))
    audit_result = compact_candidate_audit(
        audit_path, status_dir, aligner="USalign", input_rows=input_rows
    )
    alignment_result = compact_alignment_records(alignment_dir, status_dir, aligner="USalign")
    deleted = delete_raw_alignment_files(alignment_root, run_root)
    result = {
        "status": "validated_compacted",
        "alignment_dir": str(alignment_dir),
        "passed_pair_count": len(passed_pairs),
        "deleted_raw_alignment_file_count": deleted,
        **alignment_result,
        **audit_result,
    }
    (status_dir / "replay_compact_status.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    args = parser.parse_args()
    result = replay_and_compact(args.run_root.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result.get("status") in {"validated_compacted", "skipped"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
