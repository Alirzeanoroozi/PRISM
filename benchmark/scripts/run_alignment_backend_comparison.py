#!/usr/bin/env python3
"""Run a reproducible TM-align versus GTalign runtime pilot.

This measures alignment-stage execution only. It does not claim docking
quality, because no refinement or native-complex evaluator is involved.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmark.scripts.experimental_alignment.bench_alignment_backends import (
    bench_gtalign_batch,
    bench_tmalign_pairs,
    build_dataset,
)


def run(args: argparse.Namespace) -> dict:
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=False)
    queries, refs, qdir, rdir = build_dataset(output_root, args.queries, args.refs, args.residues)
    t0 = time.perf_counter()
    tm = bench_tmalign_pairs(
        args.tmalign.resolve(), queries, refs, output_root / "tmalign_work",
        limit_pairs=args.tmalign_limit_pairs,
    )
    gt = bench_gtalign_batch(
        args.gtalign.resolve(), qdir, rdir, output_root / "gtalign_out", args.refs, args.residues,
    )
    result = {
        "schema_version": "alignment-backend-comparison/v1",
        "status": "completed" if tm["ok"] and gt["ok"] else "failed",
        "dataset": {"queries": args.queries, "references": args.refs, "residues": args.residues},
        "tmalign": tm,
        "gtalign": gt,
        "total_wall_seconds": time.perf_counter() - t0,
        "quality_metrics": None,
        "interpretation": "alignment-stage runtime/output pilot only; no DockQ or iRMSD claim",
    }
    (output_root / "comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmalign", type=Path, required=True)
    parser.add_argument("--gtalign", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--queries", type=int, default=10)
    parser.add_argument("--refs", type=int, default=10)
    parser.add_argument("--residues", type=int, default=50)
    parser.add_argument("--tmalign-limit-pairs", type=int)
    args = parser.parse_args()
    result = run(args)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
