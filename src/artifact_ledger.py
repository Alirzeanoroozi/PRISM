#!/usr/bin/env python3
"""
Artifact Ledger CLI - write subcommand.

Thin wrapper over src.provenance.run_evidence for Phase 1 quickstart scenarios.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from src.provenance.run_evidence import (
    ArtifactObservation,
    append_artifact_observation,
    observe_artifact,
    validate_artifact_ledger,
    write_artifact_ledger,
)


def write_artifact(args: argparse.Namespace) -> int:
    """Write an artifact observation to the ledger."""
    # Resolve paths
    run_root = Path(args.run_root).resolve()
    abs_path = run_root / args.rel_path

    # Observe the artifact
    obs = observe_artifact(
        abs_path,
        run_id=args.run_id,
        stage=args.stage,
        dataset_row_id=args.row_id,
        scientific_role=args.role,
        run_relative_path=args.rel_path,
        produced_by=args.produced_by or "",
    )

    # Override status if provided
    if args.status:
        obs.status = args.status

    # Ledger path
    ledger_path = run_root / "artifact_ledger.tsv"

    # Append or create
    if ledger_path.is_file():
        try:
            append_artifact_observation(ledger_path, obs)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1
    else:
        write_artifact_ledger(ledger_path, [obs])

    print(f"Recorded: {obs.key} → {obs.status}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="PRISM Artifact Ledger CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # write subcommand
    write = subparsers.add_parser("write", help="Record an artifact in the ledger")
    write.add_argument("--run-id", required=True, help="Run ID")
    write.add_argument("--run-root", required=True, help="Run root directory")
    write.add_argument("--stage", required=True, help="Pipeline stage")
    write.add_argument("--row-id", required=True, help="Dataset row ID")
    write.add_argument("--role", required=True, help="Scientific role")
    write.add_argument("--rel-path", required=True, help="Relative path from run root")
    write.add_argument("--produced-by", default="", help="Tool that produced the artifact")
    write.add_argument("--status", choices=["ok", "missing", "unavailable", "failed", "partial", "not_scoreable", "error", "broken_link"], help="Artifact status (default: auto-detect)")
    write.set_defaults(func=write_artifact)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())