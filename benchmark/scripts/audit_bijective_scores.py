#!/usr/bin/env python3
"""Audit strict BM5.5 model/interface scores and provenance hashes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

SCORE_BEARING_STATUSES = {"scored", "scored_cross_only"}


def normalized_score_scope(row: dict[str, str]) -> str:
    """Normalize pre-repair complete rows while preserving explicit fallback scope."""

    if row.get("score_scope"):
        return row["score_scope"]
    if row.get("score_status") == "scored_cross_only":
        return "requested_cross_interfaces_only"
    return "complete_bijective_mapping"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit(models_path: Path, interfaces_path: Path) -> tuple[dict[str, object], list[str]]:
    models = read_tsv(models_path)
    interfaces = read_tsv(interfaces_path)
    failures: list[str] = []
    scored = [row for row in models if row.get("score_status") in SCORE_BEARING_STATUSES]
    status_counts = Counter(row.get("score_status", "") for row in models)
    if status_counts.get("score_failed"):
        failures.append(f"score_failed rows: {status_counts['score_failed']}")
    if any(row.get("source_gate_status") != "strict_clean" for row in scored):
        failures.append("a scored row is not strict_clean")
    # Pre-repair complete rows predate this explicit field; blank means the
    # grouped values were present and accepted by the original scorer.
    failed_auxiliary = sum(row.get("irmsd_status", "") not in ("", "scored") for row in scored)
    if len({row.get("source_model_path") for row in scored}) != len(scored):
        failures.append("scored source_model_path values are not unique")
    if len({row.get("source_model_sha256") for row in scored}) != len(scored):
        failures.append("scored source_model_sha256 values are not unique")

    hash_mismatches = Counter()
    for row in scored:
        for kind, path_field, hash_field in (
            ("model", "source_model_path", "source_model_sha256"),
            ("native", "native_pdb_path", "native_pdb_sha256"),
            ("raw_json", "raw_dockq_json", "raw_dockq_json_sha256"),
        ):
            path = Path(row.get(path_field, ""))
            if not path.is_file() or sha256_file(path) != row.get(hash_field):
                hash_mismatches[kind] += 1
        mapping = row.get("dockq_mapping", "")
        try:
            model_map, native_map = mapping.split(":", 1)
        except ValueError:
            failures.append(f"invalid mapping syntax: {mapping!r}")
            continue
        expected_model = row.get("model_receptor_chains", "") + row.get("model_ligand_chains", "")
        expected_native = row.get("native_receptor_chains", "") + row.get("native_ligand_chains", "")
        if len(model_map) != len(native_map) or set(model_map) != set(expected_model) or native_map != expected_native:
            failures.append(f"mapping contract mismatch: {mapping!r}")
        score_scope = normalized_score_scope(row)
        if score_scope == "requested_cross_interfaces_only":
            if row.get("dockq_global"):
                failures.append("cross-only row incorrectly reports GlobalDockQ")
            if row.get("dockq_global_status") != "unavailable_cross_only":
                failures.append("cross-only row lacks explicit GlobalDockQ-unavailable status")
            try:
                components = json.loads(row.get("dockq_component_runs", ""))
            except json.JSONDecodeError:
                components = []
            expected_component_count = len(row.get("native_receptor_chains", "")) * len(
                row.get("native_ligand_chains", "")
            )
            if len(components) != expected_component_count:
                failures.append("cross-only component provenance count mismatch")
            requested_interfaces = {
                receptor + ligand
                for receptor in row.get("native_receptor_chains", "")
                for ligand in row.get("native_ligand_chains", "")
            }
            if {str(component.get("interface", "")) for component in components} != requested_interfaces:
                failures.append("cross-only component interface set mismatch")
            for component in components:
                if component.get("status") == "scored":
                    path = Path(str(component.get("raw_json", "")))
                    if not path.is_file() or sha256_file(path) != component.get("raw_json_sha256"):
                        hash_mismatches["component_raw_json"] += 1
                elif component.get("status") != "no_native_interface":
                    failures.append("cross-only component has invalid status")
        elif score_scope != "complete_bijective_mapping":
            failures.append(f"unknown score scope: {score_scope!r}")
    if hash_mismatches:
        failures.append(f"hash mismatches: {dict(hash_mismatches)}")

    interface_by_model: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in interfaces:
        interface_by_model[row.get("staged_model_path", "")].append(row)
    missing_interfaces = 0
    for row in scored:
        records = interface_by_model.get(row.get("staged_model_path", ""), [])
        score_scope = normalized_score_scope(row)
        has_global = any(
            record.get("record_type") == "global" for record in records
        )
        if score_scope == "complete_bijective_mapping" and not has_global:
            missing_interfaces += 1
        if score_scope == "requested_cross_interfaces_only" and has_global:
            missing_interfaces += 1
        if score_scope == "requested_cross_interfaces_only":
            try:
                components = json.loads(row.get("dockq_component_runs", ""))
            except json.JSONDecodeError:
                components = []
            expected_cross = {
                str(component.get("interface", ""))
                for component in components
                if component.get("status") == "scored"
            }
        else:
            expected_cross = {
                receptor + ligand
                for receptor in row.get("native_receptor_chains", "")
                for ligand in row.get("native_ligand_chains", "")
            }
        observed_cross = {
            record.get("interface", "")
            for record in records
            if record.get("record_type") == "interface"
            and record.get("requested_cross_interface") == "True"
        }
        if score_scope == "requested_cross_interfaces_only":
            if not expected_cross or observed_cross != expected_cross:
                missing_interfaces += 1
        elif not observed_cross or not observed_cross.issubset(expected_cross):
            # DockQ omits Cartesian native chain pairs that do not form an
            # interface. Complete scoring still requires at least one actual
            # requested cross-interface record and forbids unrelated records.
            missing_interfaces += 1
    if missing_interfaces:
        failures.append(f"missing global/requested interface records: {missing_interfaces}")

    report: dict[str, object] = {
        "model_rows": len(models),
        "interface_rows": len(interfaces),
        "score_status": dict(status_counts),
        "score_scope": dict(Counter(normalized_score_scope(row) for row in scored)),
        "dockq_global_status": dict(Counter(row.get("dockq_global_status", "") for row in scored)),
        "irmsd_status": dict(Counter(row.get("irmsd_status") or "scored" for row in scored)),
        "score_bearing_rows": len(scored),
        "failed_auxiliary_irmsd_rows": failed_auxiliary,
        "source_gate_status": dict(Counter(row.get("source_gate_status", "") for row in models)),
        "scored_dataset_row_count": len({row.get("dataset_row_id") for row in scored}),
        "dockq_versions": sorted({row.get("dockq_version", "") for row in scored}),
        "mapping_validation_status": dict(Counter(row.get("mapping_validation_status", "") for row in scored)),
        "hash_mismatches": dict(hash_mismatches),
        "missing_interface_contracts": missing_interfaces,
        "audit_status": "passed" if not failures else "failed",
        "failures": failures[:100],
    }
    return report, failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", type=Path, required=True)
    parser.add_argument("--interfaces", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-model-rows", type=int)
    parser.add_argument("--expected-scored-rows", type=int)
    parser.add_argument("--expected-cross-only-rows", type=int)
    parser.add_argument("--expected-not-scoreable-rows", type=int)
    parser.add_argument("--expected-failed-auxiliary-rows", type=int)
    args = parser.parse_args()
    report, failures = audit(args.models, args.interfaces)
    if args.expected_model_rows is not None and report["model_rows"] != args.expected_model_rows:
        failures.append(f"expected {args.expected_model_rows} model rows; found {report['model_rows']}")
    observed_scored = report["score_bearing_rows"]
    if args.expected_scored_rows is not None and observed_scored != args.expected_scored_rows:
        failures.append(f"expected {args.expected_scored_rows} score-bearing rows; found {observed_scored}")
    observed_cross_only = report["score_scope"].get("requested_cross_interfaces_only", 0)
    if args.expected_cross_only_rows is not None and observed_cross_only != args.expected_cross_only_rows:
        failures.append(
            f"expected {args.expected_cross_only_rows} cross-only rows; found {observed_cross_only}"
        )
    observed_not_scoreable = report["score_status"].get("not_scoreable", 0)
    if args.expected_not_scoreable_rows is not None and observed_not_scoreable != args.expected_not_scoreable_rows:
        failures.append(
            f"expected {args.expected_not_scoreable_rows} not-scoreable rows; found {observed_not_scoreable}"
        )
    if (
        args.expected_failed_auxiliary_rows is not None
        and report["failed_auxiliary_irmsd_rows"] != args.expected_failed_auxiliary_rows
    ):
        failures.append(
            f"expected {args.expected_failed_auxiliary_rows} failed auxiliary iRMSD rows; "
            f"found {report['failed_auxiliary_irmsd_rows']}"
        )
    report["failures"] = failures[:100]
    report["audit_status"] = "passed" if not failures else "failed"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
