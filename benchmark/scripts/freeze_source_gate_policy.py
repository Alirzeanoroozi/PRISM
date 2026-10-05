#!/usr/bin/env python3
"""Freeze the source-gate inclusion policy without modifying benchmark inputs.

The aggregate source gate establishes file availability and structural contract
status, but it does not itself define whether unresolved rows may enter a
confirmatory denominator.  This script makes that decision explicit and
hash-pins the evidence used to make it.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def evidence(path: Path) -> dict[str, str]:
    return {"path": str(path.resolve()), "sha256": sha256_file(path)}


def freeze_policy(summary_path: str | Path, failures_path: str | Path, output_path: str | Path) -> dict[str, Any]:
    summary_path = Path(summary_path).resolve()
    failures_path = Path(failures_path).resolve()
    output_path = Path(output_path).resolve()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    failures = read_tsv(failures_path)

    expected = int(summary["expected_row_count"])
    failed_ids = sorted(summary["validation_failed_dataset_row_ids"])
    observed_failure_ids = sorted({row["dataset_row_id"] for row in failures})
    if failed_ids != observed_failure_ids:
        raise ValueError(
            "summary and validation failure identities disagree: "
            f"summary_only={sorted(set(failed_ids) - set(observed_failure_ids))} "
            f"failure_only={sorted(set(observed_failure_ids) - set(failed_ids))}"
        )
    if int(summary["validation_failed_dataset_row_count"]) != len(failed_ids):
        raise ValueError("summary validation failure count does not match its ID list")
    failure_keys = [(row["dataset_row_id"], row["source_role"]) for row in failures]
    if len(failure_keys) != len(set(failure_keys)):
        raise ValueError("validation failures contain duplicate dataset_row_id/source_role records")
    clean = expected - len(failed_ids)

    rows = []
    for row in failures:
        rows.append(
            {
                "dataset_row_id": row["dataset_row_id"],
                "source_role": row["source_role"],
                "native_complex": row["native_complex"],
                "raw_receptor_selector": row["raw_receptor_selector"],
                "raw_ligand_selector": row["raw_ligand_selector"],
                "archive_prefix": row["archive_prefix"],
                "parse_status": row["parse_status"],
                "expected_chain_ids": row["expected_chain_ids"],
                "polymer_chain_ids": row["polymer_chain_ids"],
                "chain_set_status": row["chain_set_status"],
                "error": row["error"],
            }
        )
    rows.sort(key=lambda row: (row["dataset_row_id"], row["source_role"]))

    evidence_files = {
        "source_gate_summary": summary_path,
        "validation_failures": failures_path,
    }
    for name in ("source_manifest.tsv", "structure_validation.tsv"):
        candidate = summary_path.parent / name
        if candidate.is_file():
            evidence_files[name.removesuffix(".tsv")] = candidate

    policy: dict[str, Any] = {
        "schema_version": "source-gate-policy/v1",
        "cohort": "repository_bm5_bm5.5_extension",
        "identity": {
            "primary_key": "dataset_row_id",
            "benchmark_row_count": expected,
            "raw_selector_authority": "benchmark CSV Complex and raw receptor/ligand selectors",
            "archive_bytes": "immutable curated archive members; no full-PDB substitution",
        },
        "decision": {
            "status": "blocked_source_authority",
            "confirmatory_run_authorized": False,
            "strict_source_clean_row_count": clean,
            "audit_only_row_count": len(failed_ids),
            "excluded_dataset_row_ids": failed_ids,
            "inclusion_rule": "A row is confirmatory-eligible only when all four curated roles pass the frozen chain contract.",
            "automatic_repair": False,
            "automatic_orientation_swap": False,
            "automatic_full_pdb_substitution": False,
        },
        "reason": {
            "confirmed": f"All {expected} curated role sets are present and hash-staged, but {len(failed_ids)} rows fail the frozen chain contract; file presence alone cannot establish receptor/ligand authority.",
            "unresolved": f"The {len(failed_ids)} rows require authoritative resolution of orientation, blank/extra chains, or selector-to-content mapping before inclusion.",
        },
        "failed_rows": rows,
        "evidence": {name: evidence(path) for name, path in sorted(evidence_files.items())},
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(policy, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return policy


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--validation-failures", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    policy = freeze_policy(args.summary, args.validation_failures, args.output)
    print(json.dumps({"output": str(args.output.resolve()), "status": policy["decision"]["status"], "eligible": policy["decision"]["strict_source_clean_row_count"], "audit_only": policy["decision"]["audit_only_row_count"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
