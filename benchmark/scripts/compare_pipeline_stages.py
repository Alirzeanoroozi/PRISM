#!/usr/bin/env python3
"""Build a stage-by-stage, provenance-preserving pipeline comparison ledger.

The ledger is deliberately observational: it compares retained artifacts and
does not infer causal quality differences when inputs or stage contracts are
not matched. PDB integrity checks use Bio.PDB and preserve parse failures.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from Bio.PDB import PDBParser


def count_files(root: Path, *relative_dirs: str, suffixes: tuple[str, ...] = ()) -> int:
    count = 0
    for relative in relative_dirs:
        directory = root / relative
        if not directory.is_dir():
            continue
        for path in directory.rglob("*"):
            if path.is_file() and (not suffixes or path.suffix.lower() in suffixes):
                count += 1
    return count


def read_inputs(root: Path) -> int:
    path = root / "inputs.csv"
    if not path.is_file():
        return 0
    with path.open(newline="", encoding="utf-8") as handle:
        return sum(1 for row in csv.DictReader(handle) if any(row.values()))


def alignment_counts(data_root: Path, declared_aligner: str) -> tuple[Counter[str], Counter[str]]:
    counts: Counter[str] = Counter()
    statuses: Counter[str] = Counter()
    directory = data_root / "alignment"
    if not directory.is_dir():
        return counts, statuses
    for path in directory.glob("*.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            counts["invalid_json"] += 1
            continue
        counts[str(payload.get("aligner") or declared_aligner or "unknown")] += 1
        statuses[str(payload.get("status") or "missing")] += 1
    json_names = {path.name for path in directory.glob("*.json")}
    for path in directory.iterdir():
        if path.is_file() and path.name not in json_names:
            counts[declared_aligner or "unknown"] += 1
            statuses["legacy_record"] += 1
    return counts, statuses


def pdb_integrity(data_root: Path) -> tuple[int, int, Counter[str]]:
    parser = PDBParser(QUIET=True)
    valid = invalid = 0
    chain_counts: Counter[str] = Counter()
    if not data_root.is_dir():
        return valid, invalid, chain_counts
    for path in data_root.rglob("*.pdb"):
        try:
            structure = parser.get_structure(path.stem, str(path))
            chains = tuple(sorted({chain.id for chain in structure.get_chains()}))
        except Exception:
            invalid += 1
            continue
        valid += 1
        chain_counts["".join(chains) or "<none>"] += 1
    return valid, invalid, chain_counts


def stage_signature(row: dict[str, object]) -> str:
    stages = (
        ("surface", int(row["surface_records"])),
        ("alignment", int(row["alignment_records"])),
        ("transformation", int(row["transformation_records"])),
        ("refinement", int(row["refinement_records"])),
    )
    present = [name for name, count in stages if count]
    if not present:
        return "no_derived_stage_artifact"
    return " -> ".join(present)


def collect(name: str, root: Path, kind: str, declared_aligner: str) -> dict[str, object]:
    data_root = root / "processed" if (root / "processed").is_dir() else root
    aligners, alignment_status = alignment_counts(data_root, declared_aligner)
    valid_pdb, invalid_pdb, chain_counts = pdb_integrity(data_root)
    surface = count_files(data_root, "surface_extraction", "surfaceExtract", suffixes=(".rsa", ".asa", ".pdb"))
    alignment = count_files(data_root, "alignment")
    transformation = count_files(data_root, "transformation", suffixes=(".pdb",))
    external = count_files(data_root, "rosetta_refinement", suffixes=(".pdb",))
    pyro = count_files(data_root, "pyrosetta_refinement", suffixes=(".pdb",))
    fiber = count_files(data_root, "fiberdock_refinement", "fiberdock", suffixes=(".pdb",))
    refinement = external + pyro + fiber
    row: dict[str, object] = {
        "pipeline": name,
        "kind": kind,
        "root": str(root.resolve()),
        "input_rows": read_inputs(root),
        "surface_records": surface,
        "alignment_records": alignment,
        "alignment_success": alignment_status.get("success", 0),
        "alignment_unavailable": alignment_status.get("alignment_unavailable", 0),
        "alignment_invalid_json": aligners.get("invalid_json", 0),
        "aligner_counts": json.dumps(dict(sorted(aligners.items())), sort_keys=True),
        "transformation_records": transformation,
        "external_rosetta_records": external,
        "pyrosetta_records": pyro,
        "fiberdock_records": fiber,
        "refinement_records": refinement,
        "valid_pdb_files": valid_pdb,
        "invalid_pdb_files": invalid_pdb,
        "pdb_chain_counts": json.dumps(dict(sorted(chain_counts.items())), sort_keys=True),
    }
    row["stage_signature"] = stage_signature(row)
    row["interpretation"] = (
        "observational artifact counts; not a causal method comparison"
        if kind in {"legacy", "mixed"}
        else "artifact counts from a retained current-pipeline run"
    )
    return row


def parse_spec(spec: str) -> tuple[str, Path, str, str]:
    parts = spec.split("=", 3)
    if len(parts) < 2:
        raise ValueError("spec must be NAME=ROOT[=KIND[=DECLARED_ALIGNER]]")
    name, root = parts[:2]
    kind = parts[2] if len(parts) > 2 and parts[2] else "current"
    aligner = parts[3] if len(parts) > 3 else ""
    return name, Path(root), kind, aligner


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline", action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for spec in args.pipeline:
        name, root, kind, aligner = parse_spec(spec)
        if not root.is_dir():
            raise SystemExit(f"pipeline root is missing: {root}")
        rows.append(collect(name, root, kind, aligner))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["pipeline"]
    with (args.output_dir / "pipeline_stage_ledger.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "pipeline_stage_ledger.json").write_text(
        json.dumps(rows, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    lines = [
        "# Pipeline stage comparison",
        "",
        "This report is observational. It identifies the first retained stage with materially different artifact counts, but does not claim causality unless inputs, assets, thresholds, and evaluators are matched.",
        "",
        "| Pipeline | Input rows | Surface | Alignment | Transform | Refinement | Valid PDB | Invalid PDB | Stage signature |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['pipeline']} | {row['input_rows']} | {row['surface_records']} | "
            f"{row['alignment_records']} ({row['alignment_success']} success) | "
            f"{row['transformation_records']} | {row['refinement_records']} | "
            f"{row['valid_pdb_files']} | {row['invalid_pdb_files']} | {row['stage_signature']} |"
        )
    lines += ["", "## Interpretation anchors", ""]
    for row in rows:
        lines.append(f"- `{row['pipeline']}`: {row['interpretation']}; aligners={row['aligner_counts']}; refinement split ext/PyRosetta/FiberDock={row['external_rosetta_records']}/{row['pyrosetta_records']}/{row['fiberdock_records']}.")
    (args.output_dir / "pipeline_stage_comparison.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(rows)} pipeline rows to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
