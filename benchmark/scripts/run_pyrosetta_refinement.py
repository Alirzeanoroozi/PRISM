#!/usr/bin/env python3
"""Run the opt-in PyRosetta refinement arm with explicit terminal status."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.pyrosetta_refinement import PyRosettaRefinementAdapter


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pdb", type=Path)
    parser.add_argument("--output-pdb", type=Path)
    parser.add_argument("--partners", default="A_B")
    parser.add_argument("--init-options", default="-mute all")
    parser.add_argument("--json", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--require-available", action="store_true")
    args = parser.parse_args()

    adapter = PyRosettaRefinementAdapter(init_options=args.init_options)
    if args.manifest:
        if not args.output_root:
            parser.error("--output-root is required with --manifest")
        if not args.manifest.is_file():
            parser.error(f"manifest does not exist: {args.manifest}")
        args.output_root.mkdir(parents=True, exist_ok=True)
        snapshot = args.output_root / "manifest.snapshot.csv"
        manifest_bytes = args.manifest.read_bytes()
        snapshot.write_bytes(manifest_bytes)
        manifest_sha256 = _sha256_bytes(manifest_bytes)
        results = []
        with args.manifest.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "input_pdb" not in reader.fieldnames:
                parser.error("manifest must contain an input_pdb column")
            rows = list(reader)
        seen_pose_ids = set()
        seen_outputs = set()
        for row_number, row in enumerate(rows, start=2):
                pose_id = row.get("pose_id") or row.get("dataset_row_id") or str(len(results) + 1)
                if pose_id in seen_pose_ids:
                    parser.error(f"duplicate pose_id/dataset_row_id in manifest: {pose_id}")
                seen_pose_ids.add(pose_id)
                input_path = Path(row["input_pdb"])
                if not input_path.is_absolute():
                    input_path = REPO_ROOT / input_path
                output = Path(row.get("output_pdb") or f"{pose_id}.pdb")
                if not output.is_absolute():
                    output = args.output_root / output
                output = output.resolve()
                if not output.is_relative_to(args.output_root.resolve()):
                    parser.error(f"manifest output_pdb escapes output-root: {output}")
                if str(output) in seen_outputs:
                    parser.error(f"duplicate output_pdb in manifest: {output}")
                seen_outputs.add(str(output))
                report = adapter.refine(input_path, output, partners=row.get("partners", "A_B"))
                report["pose_id"] = pose_id
                report["manifest_sha256"] = manifest_sha256
                report["manifest_row"] = {
                    "row_number": row_number,
                    "dataset_row_id": row.get("dataset_row_id"),
                    "template_id": row.get("template_id"),
                }
                results.append(report)
        (args.output_root / "results.jsonl").write_text(
            "".join(json.dumps(row, sort_keys=True) + "\n" for row in results)
        )
        aggregate = {
            "status": "success" if results and all(row.get("status") == "success" for row in results) else "blocked_or_failed",
            "rows": len(results),
            "success": sum(row.get("status") == "success" for row in results),
            "unavailable": sum(row.get("status") == "unavailable" for row in results),
            "failed": sum(row.get("status") == "failed" for row in results),
            "manifest_sha256": manifest_sha256,
        }
        (args.output_root / "summary.json").write_text(json.dumps(aggregate, indent=2, sort_keys=True) + "\n")
        print(json.dumps(aggregate, indent=2, sort_keys=True))
        return 0 if aggregate["status"] == "success" else 3 if aggregate["unavailable"] else 2

    if not args.input_pdb or not args.output_pdb or not args.json:
        parser.error("--input-pdb, --output-pdb, and --json are required without --manifest")
    report = adapter.refine(args.input_pdb, args.output_pdb, partners=args.partners)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report.get("available"):
        return 3
    return 0 if report.get("status") == "success" else 2


if __name__ == "__main__":
    raise SystemExit(main())
