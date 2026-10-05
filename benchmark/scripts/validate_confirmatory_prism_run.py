#!/usr/bin/env python3
"""Fail-closed preflight for the PRISM confirmatory comparison.

This validator checks authorization and template-exposure provenance before a
confirmatory job is submitted.  It deliberately does not modify the frozen
source policy or infer authorization from similarity results.
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

from benchmark.scripts.build_template_source_gate import load_source_gate_policy


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_confirmatory_preflight(
    *, policy_path: Path, eligible_path: Path, template_list_path: Path,
    expected_row_count: int | None = None,
) -> dict[str, object]:
    decision = load_source_gate_policy(policy_path)
    rows = _read_rows(eligible_path)
    template_hash = sha256_file(template_list_path)
    issues: list[str] = []
    warnings: list[str] = []

    if expected_row_count is not None and len(rows) != expected_row_count:
        issues.append(f"row_count_mismatch:{len(rows)}!={expected_row_count}")

    row_ids = [row.get("dataset_row_id", "") for row in rows]
    if any(not row_id for row_id in row_ids):
        issues.append("missing_dataset_row_id")
    if len(row_ids) != len(set(row_ids)):
        issues.append("duplicate_dataset_row_id")

    template_hashes = {row.get("template_list_sha256", "") for row in rows}
    if template_hashes != {template_hash}:
        issues.append("template_list_hash_mismatch")

    unauthorized = not bool(decision.get("confirmatory_run_authorized")) or decision.get("status") != "authorized"
    if unauthorized:
        warnings.append("source_authority_not_authorized")
    else:
        excluded = {str(value) for value in decision.get("excluded_dataset_row_ids", [])}
        for row in rows:
            row_id = row.get("dataset_row_id", "")
            if row_id in excluded:
                issues.append(f"audit_only_row_present:{row_id}")
            try:
                eligible_count = int(row.get("eligible_template_count") or 0)
            except (TypeError, ValueError):
                eligible_count = 0
                issues.append(f"invalid_eligible_template_count:{row_id}")
            if eligible_count <= 0:
                issues.append(f"empty_eligible_template_list:{row_id}")
            values = [value for value in row.get("eligible_templates", "").split(",") if value]
            if len(values) != len(set(values)):
                issues.append(f"duplicate_eligible_template:{row_id}")

    status = "blocked" if unauthorized else ("invalid" if issues else "ready")
    return {
        "status": status,
        "confirmatory_run_authorized": bool(decision.get("confirmatory_run_authorized")),
        "source_gate_status": decision.get("status", ""),
        "dataset_row_count": len(rows),
        "template_list_sha256": template_hash,
        "issues": issues,
        "warnings": warnings,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--eligible-templates", type=Path, required=True)
    parser.add_argument("--template-list", type=Path, required=True)
    parser.add_argument("--expected-row-count", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = validate_confirmatory_preflight(
        policy_path=args.policy,
        eligible_path=args.eligible_templates,
        template_list_path=args.template_list,
        expected_row_count=args.expected_row_count,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))
    return 0 if report["status"] == "ready" else 2


if __name__ == "__main__":
    raise SystemExit(main())
