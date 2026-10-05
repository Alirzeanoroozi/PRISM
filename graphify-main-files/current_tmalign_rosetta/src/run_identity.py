#!/usr/bin/env python3
"""
Run Identity CLI - declare-contract and init-run subcommands.

Thin wrapper over src.provenance.run_evidence for Phase 1 quickstart scenarios.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from src.provenance.run_evidence import (
    build_declared_contract,
    build_execution_attempt,
    canonical_hash,
    canonical_json,
    _make_run_id,
)


def _load_pairs(pairs_path: str) -> list[str]:
    """Load pair IDs from file (one per line, skip empty/comments)."""
    pairs = []
    with open(pairs_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                pairs.append(line)
    return pairs


def _load_template_inventory(template_dir: str | None) -> list[dict[str, Any]]:
    """Load template inventory from directory (placeholder for future)."""
    if not template_dir:
        return []
    # TODO: Implement template inventory loading
    return []


def declare_contract(args: argparse.Namespace) -> int:
    """Declare a contract from user selectors + template inventory + parameters."""
    pairs = _load_pairs(args.pairs)
    template_inventory = _load_template_inventory(args.template_dir)

    # Build input selectors
    input_selectors = {
        "raw": pairs,
        "normalized": [],
    }

    # Build resource request
    resource_request = {
        "cpus": args.cpus,
        "memory_gb": args.memory_gb,
        "time_hours": args.time_hours,
        "partition": args.partition,
        "gpu": args.gpu,
    }

    # Build tool fingerprints (placeholder)
    tool_fingerprints = [
        {"name": "prism", "version": "1.0.0", "sha256": None},
    ]

    contract = build_declared_contract(
        pipeline_version=args.pipeline_version,
        stages_enabled=args.stages,
        aligner=args.aligner,
        refiner=args.refiner,
        input_selectors=input_selectors,
        template_inventory=template_inventory,
        parameters=args.parameters or {},
        resource_request=resource_request,
        source_inventory={
            "git_head": "a" * 40,  # placeholder
            "git_diff_hash": "b" * 64,  # placeholder
            "declared_untracked": [],
        },
        tool_fingerprints=tool_fingerprints,
    )

    # Write output
    output = canonical_json(contract)
    if args.out:
        Path(args.out).write_text(output, encoding="utf-8")
    else:
        print(output)

    return 0


def init_run(args: argparse.Namespace) -> int:
    """Initialize run identity from declared contract."""
    # Load declared contract
    with open(args.declared_contract_path, "r", encoding="utf-8") as f:
        declared_contract = json.load(f)

    # Use the contract_hash from the contract itself (computed without the contract_hash field)
    contract_hash = declared_contract.get("contract_hash", "")
    if not contract_hash:
        print("Error: declared contract missing contract_hash", file=sys.stderr)
        return 1

    # If user provided a hash, verify it matches
    if args.declared_contract_hash and args.declared_contract_hash != contract_hash:
        print(f"Error: declared_contract_hash mismatch. Expected {contract_hash}, got {args.declared_contract_hash}", file=sys.stderr)
        return 1

    # Generate run_id
    run_id = _make_run_id()

    # Build execution attempt
    attempt = build_execution_attempt(
        contract=declared_contract,
        run_id=run_id,
        parent_run_id=args.parent_run_id,
        supersedes_reason=args.supersedes_reason,
    )

    # Write output
    output = canonical_json(attempt)
    if args.out:
        Path(args.out).write_text(output, encoding="utf-8")
    else:
        print(output)

    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="PRISM Run Identity CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # declare-contract subcommand
    declare = subparsers.add_parser("declare-contract", help="Declare a pipeline contract")
    declare.add_argument("--pairs", required=True, help="Path to pairs file (one pair ID per line)")
    declare.add_argument("--template-dir", help="Path to template inventory directory")
    declare.add_argument("--aligner", required=True, choices=["tmalign", "gtalign"], help="Aligner to use")
    declare.add_argument("--refiner", required=True, choices=["pyrosetta", "external_rosetta", "none"], help="Refiner to use")
    declare.add_argument("--stages", nargs="+", required=True, help="Enabled pipeline stages")
    declare.add_argument("--pipeline-version", default="1.0.0", help="Pipeline version")
    declare.add_argument("--parameters", type=json.loads, default="{}", help="JSON parameters object")
    declare.add_argument("--cpus", type=int, default=1, help="CPUs requested")
    declare.add_argument("--memory-gb", type=int, default=4, help="Memory in GB")
    declare.add_argument("--time-hours", type=int, default=1, help="Time limit in hours")
    declare.add_argument("--partition", default="ai", help="Slurm partition")
    declare.add_argument("--gpu", action="store_true", help="Request GPU")
    declare.add_argument("--out", help="Output path for declared contract JSON")
    declare.set_defaults(func=declare_contract)

    # init-run subcommand
    init = subparsers.add_parser("init-run", help="Initialize run identity from declared contract")
    init.add_argument("--declared-contract-hash", help="SHA256 hash of declared contract (optional, verified against contract)")
    init.add_argument("--declared-contract-path", required=True, help="Path to declared contract JSON")
    init.add_argument("--parent-run-id", help="Parent run ID for retry/supersede")
    init.add_argument("--supersedes-reason", help="Reason for superseding parent run")
    init.add_argument("--out", help="Output path for run manifest JSON")
    init.set_defaults(func=init_run)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())