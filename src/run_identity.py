#!/usr/bin/env python3
"""
Run Identity CLI - declare-contract and init-run subcommands.

Thin wrapper over src.provenance.run_evidence for Phase 1 quickstart scenarios.
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from src.provenance.run_evidence import (
    build_declared_contract,
    build_execution_attempt,
    canonical_hash,
    canonical_json,
    _make_run_id,
    sha256_file,
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
    """Load a deterministic template file inventory with SHA256 hashes."""
    if not template_dir:
        return []
    root = Path(template_dir)
    paths = [root] if root.is_file() else sorted(
        path for path in root.rglob("*") if path.is_file()
    )
    base = root.parent if root.is_file() else root
    return [
        {
            "path": str(path.relative_to(base)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in paths
    ]


def _git_source_inventory() -> dict[str, Any]:
    """Record source identity without folding runtime state into contract hash."""
    def run(*command: str) -> str:
        result = subprocess.run(
            list(command), capture_output=True, text=True, check=False,
        )
        return result.stdout

    status = run("git", "status", "--short")
    tracked_diff = run("git", "diff", "HEAD", "--binary")
    index_diff = run("git", "diff", "--cached", "--binary")
    untracked = [
        line for line in status.splitlines()
        if line.startswith("?? ")
    ]
    return {
        "git_head": run("git", "rev-parse", "HEAD").strip(),
        "git_diff_hash": hashlib.sha256(tracked_diff.encode("utf-8")).hexdigest(),
        "git_index_diff_hash": hashlib.sha256(index_diff.encode("utf-8")).hexdigest(),
        "status_hash": hashlib.sha256(status.encode("utf-8")).hexdigest(),
        "declared_untracked": untracked,
    }


def _tool_fingerprints(aligner: str, refiner: str) -> list[dict[str, Any]]:
    """Resolve configured executable identities without running tools."""
    names = {
        "prism": Path("prism.py"),
        aligner: Path(os.environ.get(
            "PRISM_GTALIGN" if aligner == "gtalign" else "PRISM_TMALIGN",
            "gtalign" if aligner == "gtalign" else "external_tools/TMalign",
        )),
    }
    if refiner == "external_rosetta":
        names["rosetta_prepack"] = Path(os.environ.get(
            "PRISM_ROSETTA_PREPACK", "docking_prepack_protocol.static.linuxgccrelease"
        ))
    elif refiner == "pyrosetta":
        names["pyrosetta"] = Path("pyrosetta")
    elif refiner == "fiberdock":
        names["fiberdock"] = Path(os.environ.get(
            "PRISM_FIBERDOCK_DIR", "external_tools/fiberdock"
        )) / "FiberDock"

    fingerprints = []
    for name, configured in names.items():
        resolved = configured
        if not resolved.is_absolute():
            resolved = Path.cwd() / resolved
        if not resolved.is_file():
            found = shutil.which(str(configured))
            resolved = Path(found) if found else resolved
        fingerprints.append({
            "name": name,
            "resolved_path": str(resolved.resolve()) if resolved.is_file() else "",
            "version": None,
            "sha256": sha256_file(resolved) if resolved.is_file() else None,
        })
    return fingerprints


def declare_contract(args: argparse.Namespace) -> int:
    """Declare a contract from user selectors + template inventory + parameters."""
    pairs = _load_pairs(args.pairs)
    template_inventory = _load_template_inventory(args.template_dir)

    # Build input selectors
    input_selectors = {
        "raw": pairs,
        "normalized": sorted(set(pairs)),
    }

    # Build resource request
    resource_request = {
        "cpus": args.cpus,
        "memory_gb": args.memory_gb,
        "time_hours": args.time_hours,
        "partition": args.partition,
        "gpu": args.gpu,
    }

    tool_fingerprints = _tool_fingerprints(args.aligner, args.refiner)

    contract = build_declared_contract(
        pipeline_version=args.pipeline_version,
        stages_enabled=args.stages,
        aligner=args.aligner,
        refiner=args.refiner,
        input_selectors=input_selectors,
        template_inventory=template_inventory,
        parameters=args.parameters or {},
        resource_request=resource_request,
        source_inventory=_git_source_inventory(),
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
        command=args.command,
        output_root=args.output_root,
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
    declare.add_argument("--parameters", type=json.loads, default={}, help="JSON parameters object")
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
    init.add_argument("--output-root", default=".", help="Run output root recorded in the manifest")
    init.add_argument(
        "--command",
        action="append",
        default=[],
        help="One command argument; repeat (use --command=--flag for option-like values).",
    )
    init.add_argument("--out", help="Output path for run manifest JSON")
    init.set_defaults(func=init_run)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
