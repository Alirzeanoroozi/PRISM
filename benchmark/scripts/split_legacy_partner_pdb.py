#!/usr/bin/env python3
"""Create an explicitly exploratory two-chain view of a legacy PDB.

Some historical PRISM outputs concatenate two partners with one chain ID and
restart residue numbering at the partner boundary. This tool does not claim
to recover the original provenance; it only makes that boundary explicit for
diagnostic scoring, and records the boundary and source hash in its log.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def split_pdb(source, destination, first_chain="A", second_chain="B"):
    source = Path(source)
    destination = Path(destination)
    raw = source.read_bytes()
    lines = raw.decode("ascii", errors="replace").splitlines(True)
    atom_indices = []
    residue_numbers = []
    for index, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 27:
            continue
        try:
            residue_number = int(line[22:26])
        except ValueError as exc:
            raise ValueError(f"invalid residue number at line {index + 1}") from exc
        atom_indices.append(index)
        residue_numbers.append(residue_number)

    resets = [
        position
        for position, (previous, current) in enumerate(
            zip(residue_numbers, residue_numbers[1:]), start=1
        )
        if current < previous
    ]
    if len(resets) != 1:
        raise ValueError(f"expected exactly one residue-number reset, found {len(resets)}")
    boundary = resets[0]
    if not first_chain or not second_chain or first_chain == second_chain:
        raise ValueError("first and second chains must be distinct non-empty IDs")

    atom_position = {line_index: position for position, line_index in enumerate(atom_indices)}
    transformed = []
    for index, line in enumerate(lines):
        if index not in atom_position:
            transformed.append(line)
            continue
        position = atom_position[index]
        chain = first_chain if position < boundary else second_chain
        transformed.append(line[:21] + chain[0] + line[22:])

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(transformed))
    return {
        "source": str(source.resolve()),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "destination": str(destination.resolve()),
        "first_chain": first_chain,
        "second_chain": second_chain,
        "atom_records": len(atom_indices),
        "first_partner_atom_records": boundary,
        "second_partner_atom_records": len(atom_indices) - boundary,
        "residue_reset_count": len(resets),
        "mode": "exploratory_boundary_split_not_proven_provenance_preserving",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("destination")
    parser.add_argument("--first-chain", default="A")
    parser.add_argument("--second-chain", default="B")
    parser.add_argument("--log", default=None)
    args = parser.parse_args()
    record = split_pdb(args.source, args.destination, args.first_chain, args.second_chain)
    if args.log:
        Path(args.log).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
