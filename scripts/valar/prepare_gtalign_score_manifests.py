#!/usr/bin/env python3
"""Prepare immutable, stage-specific manifests for corrected GTalign scoring."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def key(row: dict[str, object]) -> tuple[str, str, str, str]:
    return tuple(str(row.get(name, "")) for name in ("template", "receptor", "ligand", "orientation"))


def transformed_entries(run_root: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    source = run_root / "processed/dockq_irmsd/dockq_irmsd.csv"
    if not source.is_file():
        raise FileNotFoundError(f"historical transformed score index not found: {source}")
    entries: list[dict[str, object]] = []
    missing: list[dict[str, object]] = []
    with source.open(newline="", encoding="utf-8") as handle:
        for index, row in enumerate(csv.DictReader(handle)):
            entry = {
                "index": index,
                "stage": "transformed",
                "template": row.get("template", ""),
                "receptor": row.get("receptor", ""),
                "ligand": row.get("ligand", ""),
                "orientation": row.get("orientation", ""),
                "model_pdb": row.get("model_pdb", ""),
                "model_receptor_chains": row.get("model_receptor_chains", ""),
                "model_ligand_chains": row.get("model_ligand_chains", ""),
            }
            entries.append(entry)
            if not Path(str(entry["model_pdb"])).is_file():
                missing.append({**entry, "reason": "model_pdb_missing"})
    return entries, missing


def pyrosetta_entries(run_root: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    roots = [
        run_root / "processed/pyrosetta_refinement",
        run_root / "processed/pyrosetta_refinement_retry_20260919",
    ]
    chosen: dict[tuple[str, str, str, str], dict[str, object]] = {}
    inventory: dict[tuple[str, str, str, str], dict[str, object]] = {}
    for root in roots:
        for path in sorted(root.glob("task_*.json")):
            payload = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(payload, list):
                raise ValueError(f"task ledger must be a list: {path}")
            for record in payload:
                record_key = tuple(str(record.get(name, "")) for name in ("template", "receptor", "ligand", "orientation"))
                inventory[record_key] = {"source": str(path), **record}
                if record.get("status") not in {"success", "already_present"}:
                    continue
                output = Path(str(record.get("output_path", "")))
                if not output.is_file():
                    continue
                partners = str(record.get("partners", ""))
                if "_" not in partners:
                    continue
                left_chains, right_chains = partners.split("_", 1)
                chosen[record_key] = {
                    "stage": "pyrosetta",
                    "template": record.get("template", ""),
                    "receptor": record.get("receptor", ""),
                    "ligand": record.get("ligand", ""),
                    "orientation": record.get("orientation", ""),
                    "model_pdb": str(output.resolve()),
                    "model_receptor_chains": left_chains,
                    "model_ligand_chains": right_chains,
                }
    entries = []
    for index, item in enumerate(sorted(chosen.values(), key=key)):
        entries.append({"index": index, **item})
    missing = []
    for missing_key, record in sorted(inventory.items()):
        if missing_key not in chosen:
            missing.append({"key": list(missing_key), "reason": record.get("status", "missing"), "source": record.get("source", "")})
    return entries, missing


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--stage", choices=("transformed", "pyrosetta"), required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--chunk-size", type=int, default=400)
    args = parser.parse_args()

    run_root = Path(args.run_root).resolve()
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    if args.stage == "transformed":
        entries, missing = transformed_entries(run_root)
    else:
        entries, missing = pyrosetta_entries(run_root)
    summary_path = run_root / "run_summary.json"
    panel_path = run_root / "templates/checked_templates.txt"
    payload = {
        "schema_version": "prism-gtalign-score-manifest-20260920",
        "stage": args.stage,
        "ranking": False,
        "prodigy": False,
        "run_root": str(run_root),
        "parent_run_summary_sha256": sha256(summary_path),
        "template_count": sum(1 for _ in panel_path.open(encoding="utf-8")) if panel_path.is_file() else None,
        "template_sha256": sha256(panel_path) if panel_path.is_file() else None,
        "candidate_count": len(entries),
        "missing_count": len(missing),
        "chunk_size": args.chunk_size,
        "entries": entries,
        "missing": missing,
    }
    manifest = output_root / "input_manifest.json"
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_root / "missing_refinement.json").write_text(json.dumps(missing, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"stage": args.stage, "candidate_count": len(entries), "missing_count": len(missing), "manifest": str(manifest)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
