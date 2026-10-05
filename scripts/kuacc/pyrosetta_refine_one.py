#!/usr/bin/env python3
"""Run one explicit PyRosetta refinement candidate in an isolated root."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--pipeline-repo", required=True)
    parser.add_argument("--left", required=True)
    parser.add_argument("--right", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(args.run_root).resolve()))
    sys.path.insert(0, str(Path(args.pipeline_repo).resolve()))
    os.chdir(Path(args.run_root).resolve())
    from src.pyrosetta_refinement import refine_pairs

    records = refine_pairs([(args.left, args.right)], output_root=args.output_root)
    print(json.dumps(records, indent=2, default=str))
    return 0 if records and records[0].get("status") == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
