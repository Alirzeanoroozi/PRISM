#!/usr/bin/env python3
"""Stage current PRISM refinement outputs for canonical bijective scoring.

The PRISM-main benchmark analyzers identify template partner-chain assignments
from a historical Rosetta filename.  Current outputs keep those assignments in
their native refinement names but do not use the historical filename grammar.
This script creates *symlinks only* using that grammar, without modifying a
model PDB or its chain IDs.  Its manifest also records the query-partner chain
groups required by the strict multichain scorer.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    from benchmark.scripts.standardized_evaluator import validate_raw_pdb_chain_contract
except ModuleNotFoundError:  # Direct `python benchmark/scripts/...` invocation.
    from standardized_evaluator import validate_raw_pdb_chain_contract

TEMPLATE_TOKEN = r"[A-Za-z0-9]{4}(?:_[A-Za-z0-9]+?|[A-Za-z0-9]+?)"
CURRENT_MODEL_RE = re.compile(
    rf"^(?P<target>[A-Za-z0-9]+)_(?P<tpl1>{TEMPLATE_TOKEN})"
    rf"_(?P<tpl2>{TEMPLATE_TOKEN})_o\d+_L_.*_R_.*rosetta(?:_0001)?\.pdb$"
)
EXTERNAL_ROSETTA_MODEL_RE = re.compile(
    r"^(?P<template>[A-Za-z0-9]{6})_"
    r"(?P<target_left>[A-Za-z0-9]{4}[A-Za-z0-9]*)_"
    r"(?P<target_right>[A-Za-z0-9]{4}[A-Za-z0-9]*)_"
    r"o(?P<orientation>\d+)_L_.*_R_rosetta(?:_0001(?:_0001)?)?\.pdb$"
)
EXTERNAL_ROSETTA_SUFFIX_RE = re.compile(
    r"^(?P<pose>.*_R_rosetta)(?P<suffix>_0001_0001|_0001)?\.pdb$"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def model_chain_ids(path: Path) -> list[str]:
    chains: list[str] = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip() or "_"
                if chain not in chains:
                    chains.append(chain)
    return chains


def model_ca_counts(path: Path) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                counts[line[21].strip() or "_"] += 1
    return dict(counts)


def select_final_external_models(paths: list[Path]) -> list[Path]:
    """Select one most-refined Rosetta artifact for each generated pose."""

    suffix_rank = {"": 0, "_0001": 1, "_0001_0001": 2}
    selected: dict[tuple[Path, str], tuple[int, Path]] = {}
    for path in paths:
        match = EXTERNAL_ROSETTA_SUFFIX_RE.fullmatch(path.name)
        if match is None:
            continue
        key = (path.parent, match.group("pose"))
        rank = suffix_rank[match.group("suffix") or ""]
        previous = selected.get(key)
        if previous is None or rank > previous[0]:
            selected[key] = (rank, path)
    return sorted(path for _, path in selected.values())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_current_model_name(model_name: str, backend: str = "external_rosetta") -> dict[str, str] | None:
    """Parse a paired current refinement name from either backend grammar."""

    external_match = (
        EXTERNAL_ROSETTA_MODEL_RE.fullmatch(model_name)
        if backend == "external_rosetta"
        else None
    )
    if external_match is not None:
        fields = external_match.groupdict()
        template = fields["template"]
        return {
            "target": "",
            "target_left": fields["target_left"],
            "target_right": fields["target_right"],
            "template_1": template[:4] + template[4],
            "template_2": template[:4] + template[5],
            "orientation": fields["orientation"],
        }

    match = CURRENT_MODEL_RE.match(model_name)
    orientation = re.search(r"_o(?P<orientation>\d+)_L_", model_name)
    if match is None:
        # Current external-Rosetta outputs may put the template first and
        # preserve the raw query selectors, e.g.
        # ``1bjaAB_1DQQ_CD_3LZT__o2_L_..._R_...``.  Recover the two template
        # partner-chain tokens from the compact template prefix; the batch
        # input row remains the authority for the native complex identity.
        prefix_match = re.match(r"^(?P<template>[A-Za-z0-9]{6})_.*_o\d+_L_", model_name)
        if prefix_match is None:
            return None
        template = prefix_match.group("template")
        chains = template[4:]
        if len(chains) != 2:
            return None
        return {
            "target": "",
            "target_left": "",
            "target_right": "",
            "template_1": template[:4] + chains[0],
            "template_2": template[:4] + chains[1],
            "orientation": orientation.group("orientation") if orientation else "",
        }
    fields = match.groupdict()
    target = fields["target"]
    return {
        "target": target,
        "target_left": target,
        "target_right": target,
        "template_1": fields["tpl1"],
        "template_2": fields["tpl2"],
        "orientation": orientation.group("orientation") if orientation else "",
    }


def is_incomplete_transformation_half(model_name: str) -> bool:
    """Return whether a refinement-like file name contains only the L half."""

    return "_o" in model_name and "_L_" in model_name and "_R_" not in model_name


def normalize_target(value: str) -> str:
    value = (value or "").strip().replace(" ", "")
    value = re.sub(r"\([^)]*\)", "", value).replace("_", "")
    return value.lower()


def expected_partner_chain_count(value: str) -> int:
    """Return the number of query chains represented by a normalized target.

    PRISM target tokens are a four-character PDB identifier followed by zero
    or more chain identifiers.  An empty chain suffix denotes one selected
    chain in the BM5.5 inputs (for example ``3LZT_`` becomes ``3lzt``).
    """

    normalized = normalize_target(value)
    if len(normalized) < 4:
        raise ValueError(f"invalid target selector: {value!r}")
    return max(1, len(normalized[4:]))


def dataset_row_id(row: dict[str, str]) -> str:
    """Recover the durable BM5.5 row identity from the batch manifest."""

    existing = (row.get("dataset_row_id") or "").strip()
    if existing:
        return existing
    dataset = (row.get("benchmark_set") or row.get("dataset") or row.get("difficulty") or "").strip()
    source_row = (row.get("source_row") or row.get("source_row_number") or "").strip()
    if not dataset or not source_row:
        return ""
    try:
        row_number = int(source_row) - 1  # CSV source row 2 is dataset row 1.
    except ValueError:
        return ""
    return f"{dataset}:{row_number:06d}" if row_number > 0 else ""


def canonical_target_token(complex_id: str) -> str:
    value = (complex_id or "").strip()
    if "_" not in value or ":" not in value:
        raise ValueError(f"unexpected Complex value: {complex_id!r}")
    pdb, partners = value.split("_", 1)
    receptor, ligand = partners.split(":", 1)
    chains = "".join(ch for ch in receptor + ligand if ch.isalnum())
    if len(pdb) < 4 or not chains:
        raise ValueError(f"unexpected Complex value: {complex_id!r}")
    return pdb[:4].lower() + chains


def rows_by_batch(batch_root: Path, cohort_manifest: Path | None = None) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = {}
    cohort_rows: list[dict[str, str]] = []
    if cohort_manifest is not None:
        with cohort_manifest.open(newline="", encoding="utf-8") as handle:
            cohort_rows = list(csv.DictReader(handle, delimiter="\t"))
    for inputs in sorted(batch_root.glob("batch_*/inputs.csv")):
        rows = read_csv(inputs)
        for row in rows:
            matches = [candidate for candidate in cohort_rows if candidate.get("receptor") == row.get("Receptor") and candidate.get("ligand") == row.get("Ligand")]
            if len(matches) == 1:
                row.update(matches[0])
                row["complex"] = matches[0].get("native_complex", "").strip()
        result[inputs.parent.name] = rows
    return result


def choose_row(model_name: str, rows: list[dict[str, str]]) -> dict[str, str] | None:
    name = model_name.lower()
    compact_name = re.sub(r"[^a-z0-9]", "", name)
    hits = []
    for row in rows:
        receptor = normalize_target(row.get("Receptor", ""))
        ligand = normalize_target(row.get("Ligand", ""))
        if receptor and ligand and receptor in compact_name and ligand in compact_name:
            hits.append(row)
    return hits[0] if len(hits) == 1 else None


def stage_models(batch_root: Path, current_root: Path, output_root: Path, cohort_manifest: Path | None = None) -> list[dict[str, str]]:
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(f"refusing to mix staged models into non-empty directory: {output_root}")

    batches = rows_by_batch(batch_root, cohort_manifest=cohort_manifest)
    rows: list[dict[str, str]] = []
    seen_destinations: set[Path] = set()
    per_batch_count: defaultdict[str, int] = defaultdict(int)

    external_paths = list(current_root.glob("batch_*/processed/rosetta_refinement/*_rosetta*.pdb"))
    model_paths = set(select_final_external_models(external_paths))
    model_paths.update(path for path in external_paths if is_incomplete_transformation_half(path.name))
    model_paths.update(current_root.glob("batch_*/processed/pyrosetta_refinement/structures/*_rosetta.pdb"))
    for model in sorted(model_paths):
        batch = next((part for part in model.parts if re.fullmatch(r"batch_\d{4}", part)), "")
        record: dict[str, str] = {
            "source_model_path": str(model.resolve()),
            "source_model_sha256": sha256(model),
            "batch": batch,
            "status": "",
        }
        if is_incomplete_transformation_half(model.name):
            record.update(
                {
                    "integrity_status": "rejected",
                    "integrity_reason": "incomplete_transformation_half",
                    "status": "transformation_intermediate",
                }
            )
            rows.append(record)
            continue

        row = choose_row(model.name, batches.get(batch, []))
        if row is None:
            record["status"] = "unmatched_model"
            rows.append(record)
            continue

        set_name = row.get("benchmark_set") or row.get("dataset") or row.get("difficulty", "")
        record.update(
            {
                "pair_id": row.get("pair_id", ""),
                "dataset_row_id": dataset_row_id(row),
                "benchmark_set": set_name,
                "source_row_number": row.get("source_row") or row.get("source_row_number", ""),
                "complex": row.get("complex") or row.get("Complex") or row.get("native_complex", ""),
                "raw_receptor_selector": row.get("pdb_id_1_raw") or row.get("raw_receptor_selector", ""),
                "raw_ligand_selector": row.get("pdb_id_2_raw") or row.get("raw_ligand_selector", ""),
                "receptor": row.get("Receptor") or row.get("receptor", ""),
                "ligand": row.get("Ligand") or row.get("ligand", ""),
            }
        )

        refinement_backend = "pyrosetta" if "pyrosetta_refinement" in model.parts else "external_rosetta"
        parsed_name = parse_current_model_name(model.name, backend=refinement_backend)
        if parsed_name is None:
            record.update({"pair_id": row.get("pair_id", ""), "status": "unparseable_current_model_name"})
            rows.append(record)
            continue

        try:
            target = canonical_target_token(
                row.get("complex") or row.get("Complex") or row.get("native_complex", "")
            )
        except ValueError as exc:
            record.update({"pair_id": row.get("pair_id", ""), "status": f"invalid_complex:{exc}"})
            rows.append(record)
            continue

        per_batch_count[batch] += 1
        tpl1 = parsed_name["template_1"]
        tpl2 = parsed_name["template_2"]
        observed = model_chain_ids(model)
        ca_counts = model_ca_counts(model)
        try:
            receptor_count = expected_partner_chain_count(row.get("Receptor") or row.get("receptor", ""))
            ligand_count = expected_partner_chain_count(row.get("Ligand") or row.get("ligand", ""))
        except ValueError as exc:
            record.update(
                {
                    "pair_id": row.get("pair_id", ""),
                    "dataset_row_id": dataset_row_id(row),
                    "observed_chain_order": ",".join(observed),
                    "ca_counts": ";".join(f"{chain}:{count}" for chain, count in ca_counts.items()),
                    "integrity_status": "rejected",
                    "integrity_reason": str(exc),
                    "status": "invalid_target_chain_contract",
                }
            )
            rows.append(record)
            continue
        expected_chain_count = receptor_count + ligand_count
        model_r = "".join(observed[:receptor_count])
        model_l = "".join(observed[receptor_count:expected_chain_count])
        if len(observed) != expected_chain_count:
            reason = (
                f"unexpected_model_chain_count: observed={len(observed)} "
                f"expected={expected_chain_count} receptor={receptor_count} ligand={ligand_count}"
            )
            record.update(
                {
                    "pair_id": row.get("pair_id", ""),
                    "dataset_row_id": dataset_row_id(row),
                    "model_receptor_chains": model_r,
                    "model_ligand_chains": model_l,
                    "observed_chain_order": ",".join(observed),
                    "ca_counts": ";".join(f"{chain}:{count}" for chain, count in ca_counts.items()),
                    "chain_contract_errors": reason,
                    "integrity_status": "rejected",
                    "integrity_reason": reason,
                    "status": "invalid_model_chain_contract",
                }
            )
            rows.append(record)
            continue
        chain_contract = validate_raw_pdb_chain_contract(model, model_r, model_l)
        if not chain_contract.valid:
            record.update(
                {
                    "pair_id": row.get("pair_id", ""),
                    "model_receptor_chains": model_r,
                    "model_ligand_chains": model_l,
                    "observed_chain_order": ",".join(observed),
                    "ca_counts": ";".join(f"{chain}:{count}" for chain, count in ca_counts.items()),
                    "chain_contract_errors": "; ".join(chain_contract.errors),
                    "integrity_status": "rejected",
                    "integrity_reason": "; ".join(chain_contract.errors),
                    "status": "invalid_model_chain_contract",
                }
            )
            rows.append(record)
            continue
        destination = (
            output_root
            / set_name
            / f"joblist_{batch}"
            / f"{target}_{tpl1}_0_{tpl2}_0.rosetta_{batch}_{per_batch_count[batch]:04d}.pdb"
        )
        if destination in seen_destinations or destination.exists() or destination.is_symlink():
            raise RuntimeError(f"duplicate staged destination: {destination}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.symlink(os.path.relpath(model.resolve(), destination.parent), destination)
        seen_destinations.add(destination)
        record.update(
            {
                "pair_id": row.get("pair_id", ""),
                "dataset_row_id": dataset_row_id(row),
                "benchmark_set": set_name,
                "source_row_number": row.get("source_row") or row.get("source_row_number", ""),
                "complex": row.get("complex") or row.get("Complex") or row.get("native_complex", ""),
                "raw_receptor_selector": row.get("pdb_id_1_raw") or row.get("raw_receptor_selector", ""),
                "raw_ligand_selector": row.get("pdb_id_2_raw") or row.get("raw_ligand_selector", ""),
                "receptor": row.get("Receptor") or row.get("receptor", ""),
                "ligand": row.get("Ligand") or row.get("ligand", ""),
                "canonical_target": target,
                "template_1": tpl1,
                "template_2": tpl2,
                "target_left": parsed_name["target_left"],
                "target_right": parsed_name["target_right"],
                "orientation": parsed_name["orientation"],
                "model_receptor_chains": model_r,
                "model_ligand_chains": model_l,
                "observed_chain_order": ",".join(observed),
                "ca_counts": ";".join(f"{chain}:{count}" for chain, count in ca_counts.items()),
                "integrity_status": "valid",
                "integrity_reason": "complete_refined_pose",
                "refinement_backend": refinement_backend,
                "staged_model_path": str(destination),
                "status": "staged_symlink",
            }
        )
        rows.append(record)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", type=Path, required=True)
    parser.add_argument("--current-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cohort-manifest", type=Path, help="Optional TSV joining Receptor/Ligand to native_complex and benchmark metadata.")
    args = parser.parse_args()

    if not args.batch_root.is_dir():
        parser.error(f"batch root does not exist: {args.batch_root}")
    if not args.current_root.is_dir():
        parser.error(f"current run root does not exist: {args.current_root}")
    records = stage_models(args.batch_root, args.current_root, args.output_root, cohort_manifest=args.cohort_manifest)
    write_csv(args.manifest, records or [{"status": "no_current_models_found"}])
    staged = sum(row["status"] == "staged_symlink" for row in records)
    failed = len(records) - staged
    print(f"staged={staged} unmatched_or_invalid={failed} manifest={args.manifest}")
    return 0 if staged else 2


if __name__ == "__main__":
    raise SystemExit(main())
