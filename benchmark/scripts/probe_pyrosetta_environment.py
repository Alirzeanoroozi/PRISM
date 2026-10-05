#!/usr/bin/env python3
"""Probe the optional PyRosetta runtime and emit a JSON environment record."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.pyrosetta_refinement import probe_environment


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="optional JSON report path")
    parser.add_argument(
        "--require-available",
        action="store_true",
        help="return nonzero when PyRosetta is unavailable",
    )
    args = parser.parse_args(argv)

    report = probe_environment()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    if args.require_available and not report["available"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
