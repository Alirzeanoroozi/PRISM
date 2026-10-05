#!/usr/bin/env python3
"""Create duplicate input batches for a scheduler/resource probe.

The generated batches intentionally duplicate one frozen input batch. They are
for concurrency measurements only and must not enter scientific denominators.
"""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-batch", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--count", type=int, default=8)
    args = parser.parse_args()
    source = args.source_batch.resolve()
    source_inputs = source / "inputs.csv"
    if args.count < 1:
        parser.error("--count must be positive")
    if not source_inputs.is_file():
        parser.error(f"missing source inputs.csv: {source_inputs}")
    output = args.output_root.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(f"refusing to mix probe batches into non-empty root: {output}")
    output.mkdir(parents=True, exist_ok=True)
    with source_inputs.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    with (output / "capacity_probe_manifest.tsv").open("w", encoding="utf-8") as manifest:
        manifest.write("batch\tsource_batch\trows\n")
        for index in range(1, args.count + 1):
            batch = output / f"batch_{index:04d}"
            batch.mkdir()
            with (batch / "inputs.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            pair_list = source / "pair_list"
            if pair_list.is_file():
                shutil.copy2(pair_list, batch / "pair_list")
            manifest.write(f"batch_{index:04d}\t{source}\t{len(rows)}\n")
    print(f"wrote {args.count} duplicate capacity-probe batches to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
