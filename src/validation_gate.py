#!/usr/bin/env python3
"""
Validation Gate CLI - validate subcommand.

Thin wrapper over src.provenance.run_evidence for Phase 1 quickstart scenarios.
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

from src.provenance.run_evidence import (
    ArtifactLedgerError,
    validate_before_consume,
    validate_artifact_ledger,
)


def _load_expected_inventory(path: str | Path) -> list[tuple[str, str, str]]:
    """Load the declared artifact inventory used by the consumer gate.

    JSON accepts either a list of objects or ``{"artifacts": [...]}``.
    TSV accepts the three required columns: ``dataset_row_id``,
    ``scientific_role``, and ``run_relative_path``.  The inventory is kept
    separate from the observed ledger so missing rows cannot be hidden by
    deriving expectations from whatever happened to be produced.
    """
    inventory_path = Path(path)
    if not inventory_path.is_file():
        raise ValueError(f"expected inventory does not exist: {inventory_path}")

    if inventory_path.suffix.lower() == ".json":
        try:
            payload = json.loads(inventory_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid expected-inventory JSON: {exc}") from exc
        rows = payload.get("artifacts") if isinstance(payload, dict) else payload
        if not isinstance(rows, list):
            raise ValueError("expected-inventory JSON must be a list or an artifacts list")
    else:
        with inventory_path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

    inventory: list[tuple[str, str, str]] = []
    for index, row in enumerate(rows, start=1):
        if not isinstance(row, dict):
            raise ValueError(f"expected-inventory row {index} is not an object")
        key = tuple(str(row.get(field, "")).strip() for field in (
            "dataset_row_id", "scientific_role", "run_relative_path",
        ))
        if not all(key):
            raise ValueError(
                "expected-inventory row "
                f"{index} must contain dataset_row_id, scientific_role, and run_relative_path"
            )
        inventory.append(key)

    if len(set(inventory)) != len(inventory):
        raise ValueError("expected-inventory contains duplicate artifact identities")
    return inventory


def _expected_inventory_error(args: argparse.Namespace, detail: str) -> int:
    """Emit a fail-closed result when the declared inventory is unusable."""
    result = validate_before_consume(
        args.ledger,
        [],
        run_id=args.run_id,
        staging_dir=args.run_root,
    )
    result["overall_status"] = "fail"
    result["expected_inventory_error"] = detail
    output = json.dumps(result, indent=2)
    if args.out:
        Path(args.out).write_text(output, encoding="utf-8")
    else:
        print(output)
    return 2


def validate_gate(args: argparse.Namespace) -> int:
    """Run validation gate on completed run."""
    # Load ledger and build expected inventory
    try:
        records = validate_artifact_ledger(args.ledger)
    except ArtifactLedgerError:
        # Re-enter the canonical validator so malformed or duplicate ledgers
        # produce the normal fail-closed JSON result instead of a traceback.
        result = validate_before_consume(
            args.ledger,
            [],
            run_id=args.run_id,
            staging_dir=args.run_root,
        )
        output = json.dumps(result, indent=2)
        if args.out:
            Path(args.out).write_text(output, encoding="utf-8")
        else:
            print(output)
        return 2
    
    if args.expected_inventory:
        try:
            expected_inventory = _load_expected_inventory(args.expected_inventory)
        except (OSError, ValueError) as exc:
            return _expected_inventory_error(args, str(exc))
    else:
        # Compatibility mode for legacy callers. New consumers should always
        # pass a declared inventory so missing rows are observable.
        expected_inventory = [
            (rec.dataset_row_id, rec.scientific_role, rec.run_relative_path)
            for rec in records
        ]

    # Load declared contract for configuration scanning
    configuration = {}
    if args.declared_contract_path:
        with open(args.declared_contract_path, "r", encoding="utf-8") as f:
            configuration = json.load(f)

    # Run validation
    result = validate_before_consume(
        args.ledger,
        expected_inventory,
        run_id=args.run_id,
        staging_dir=args.run_root,
        environment=None,  # Could add --env-file option later
        argv=None,         # Could add --argv option later
        configuration=configuration,
    )

    # Write output
    output = json.dumps(result, indent=2)
    if args.out:
        Path(args.out).write_text(output, encoding="utf-8")
    else:
        print(output)

    # Exit with appropriate code
    if result["overall_status"] == "fail":
        return 2
    elif result["overall_status"] == "warn":
        return 1
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="PRISM Validation Gate CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # validate subcommand
    validate = subparsers.add_parser("validate", help="Run validation gate")
    validate.add_argument("--run-root", required=True, help="Run root directory")
    validate.add_argument("--run-id", required=True, help="Run ID")
    validate.add_argument("--ledger", required=True, help="Path to artifact ledger TSV")
    validate.add_argument(
        "--expected-inventory",
        help=(
            "Declared artifact inventory (TSV or JSON). TSV requires "
            "dataset_row_id, scientific_role, and run_relative_path columns."
        ),
    )
    validate.add_argument("--declared-contract-path", help="Path to declared contract JSON (for config scanning)")
    validate.add_argument(
        "--out",
        "--output",
        dest="out",
        help="Output path for validation gate JSON",
    )
    validate.set_defaults(func=validate_gate)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
