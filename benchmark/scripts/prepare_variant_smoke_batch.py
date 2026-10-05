#!/usr/bin/env python3
"""Select one single/single, multichain, and multi/multi pilot row."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-batch", type=Path, required=True)
    parser.add_argument("--pair-id", action="append", required=True)
    args = parser.parse_args()
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    selected = [row for row in rows if row.get("pair_id") in set(args.pair_id)]
    if len(selected) != len(set(args.pair_id)):
        found = {row.get("pair_id") for row in selected}
        missing = sorted(set(args.pair_id) - found)
        parser.error("missing pair IDs: " + ",".join(missing))
    output = args.output_batch.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(f"refusing to mix smoke inputs into non-empty directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    with (output / "inputs.csv").open("w", newline="", encoding="utf-8") as handle:
        # Retain the full pilot row so downstream scoring can recover the
        # benchmark Complex/native chain contract without a second join.
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(selected)
    with (output / "smoke_selection.tsv").open("w", encoding="utf-8") as handle:
        handle.write("pair_id\tnative_complex\tchain_context\n")
        for row in selected:
            handle.write(f"{row['pair_id']}\t{row.get('native_complex','')}\t{row.get('chain_context','')}\n")
    print(f"wrote {len(selected)} smoke rows to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
