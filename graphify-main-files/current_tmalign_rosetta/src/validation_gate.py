#!/usr/bin/env python3
"""
Validation Gate CLI - validate subcommand.

Thin wrapper over src.provenance.run_evidence for Phase 1 quickstart scenarios.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from src.provenance.run_evidence import (
    validate_before_consume,
    validate_artifact_ledger,
)


def validate_gate(args: argparse.Namespace) -> int:
    """Run validation gate on completed run."""
    # Load ledger and build expected inventory
    records = validate_artifact_ledger(args.ledger)
    
    # Build expected inventory from ledger records
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
    validate.add_argument("--declared-contract-path", help="Path to declared contract JSON (for config scanning)")
    validate.add_argument("--out", help="Output path for validation gate JSON")
    validate.set_defaults(func=validate_gate)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())