#!/usr/bin/env python3
"""Assemble row-specific BM5.5 native complexes from curated bound roles.

The curated ``*_r_b.pdb`` and ``*_l_b.pdb`` archive members are the native
truth.  This adapter verifies their staged hashes and declared chain groups,
then writes a new row-specific PDB and enriches a model staging manifest with
its exact path and provenance.  Full downloaded PDB files are never used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def read_table(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_table(path: Path, rows: list[dict[str, str]], delimiter: str) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pdb_chain_ids(path: Path) -> str:
    chains: list[str] = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM"):
                chain = line[21].strip() or "_"
                if chain not in chains:
                    chains.append(chain)
    return "".join(chains)


def coordinate_lines(path: Path, expected_chains: str) -> list[str]:
    with path.open(errors="replace") as handle:
        return [
            line if line.endswith("\n") else line + "\n"
            for line in handle
            if line.startswith("ATOM") and line[21].strip() in set(expected_chains)
        ]


def unique_value(rows: list[dict[str, str]], field: str, dataset_row_id: str) -> str:
    values = {(row.get(field) or "").strip() for row in rows}
    if len(values) != 1:
        raise ValueError(f"{dataset_row_id}: inconsistent {field}: {sorted(values)}")
    return values.pop()


def assemble_native_complexes(
    stage_manifest: Path,
    source_manifest: Path,
    staged_sources: Path,
    source_gate_policy: Path,
    output_root: Path,
    assembly_manifest: Path,
    enriched_stage_manifest: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(f"refusing to mix native assemblies into non-empty directory: {output_root}")

    stages = read_table(stage_manifest, ",")
    sources = read_table(source_manifest, "\t")
    staged = read_table(staged_sources, "\t")
    policy = json.loads(source_gate_policy.read_text(encoding="utf-8"))
    excluded = set(policy["decision"]["excluded_dataset_row_ids"])

    source_by_id: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    staged_by_id_role: defaultdict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in sources:
        source_by_id[row["dataset_row_id"]].append(row)
    for row in staged:
        staged_by_id_role[(row["dataset_row_id"], row["source_role"])].append(row)

    requested_ids = sorted({row.get("dataset_row_id", "") for row in stages if row.get("status") == "staged_symlink"})
    assembly_rows: list[dict[str, str]] = []
    assembly_by_id: dict[str, dict[str, str]] = {}
    for dataset_row_id in requested_ids:
        record = {
            "dataset_row_id": dataset_row_id,
            "source_gate_status": "audit_only" if dataset_row_id in excluded else "strict_clean",
            "assembly_status": "",
            "assembly_error": "",
        }
        try:
            identity_rows = source_by_id.get(dataset_row_id, [])
            if not identity_rows:
                raise ValueError(f"{dataset_row_id}: missing source-manifest identity")
            complex_id = unique_value(identity_rows, "native_complex", dataset_row_id)
            native_r = unique_value(identity_rows, "native_receptor_chains", dataset_row_id)
            native_l = unique_value(identity_rows, "native_ligand_chains", dataset_row_id)
            if not native_r or not native_l or set(native_r).intersection(native_l):
                raise ValueError(f"{dataset_row_id}: invalid native chain groups {native_r}:{native_l}")

            role_records: dict[str, dict[str, str]] = {}
            role_paths: dict[str, Path] = {}
            for role, expected_chains in (("native_receptor", native_r), ("native_ligand", native_l)):
                candidates = [row for row in staged_by_id_role[(dataset_row_id, role)] if row.get("status") == "staged"]
                if len(candidates) != 1:
                    raise ValueError(f"{dataset_row_id}: expected one staged {role}, found {len(candidates)}")
                source_role_rows = [row for row in identity_rows if row.get("source_role") == role]
                if len(source_role_rows) != 1:
                    raise ValueError(f"{dataset_row_id}: expected one source-manifest {role}, found {len(source_role_rows)}")
                role_record = candidates[0]
                path = Path(role_record["staged_path"])
                if not path.is_file():
                    raise ValueError(f"{dataset_row_id}: missing staged {role}: {path}")
                actual_hash = sha256_file(path)
                expected_hashes = {
                    role_record.get("source_payload_sha256", ""),
                    role_record.get("staged_sha256", ""),
                    source_role_rows[0].get("sha256", ""),
                }
                if "" in expected_hashes or expected_hashes != {actual_hash}:
                    raise ValueError(f"{dataset_row_id}: hash mismatch for {role}")
                observed_chains = pdb_chain_ids(path)
                if len(observed_chains) != len(expected_chains) or set(observed_chains) != set(expected_chains):
                    raise ValueError(
                        f"{dataset_row_id}: {role} chains {observed_chains!r} do not match {expected_chains!r}"
                    )
                role_records[role] = role_record | {"actual_sha256": actual_hash}
                role_paths[role] = path

            native_path = output_root / dataset_row_id.split(":", 1)[0] / f"{dataset_row_id.replace(':', '_')}.pdb"
            native_path.parent.mkdir(parents=True, exist_ok=True)
            lines = coordinate_lines(role_paths["native_receptor"], native_r)
            lines.append("TER\n")
            lines.extend(coordinate_lines(role_paths["native_ligand"], native_l))
            lines.extend(("TER\n", "END\n"))
            native_path.write_text("".join(lines), encoding="ascii")
            assembled_chains = pdb_chain_ids(native_path)
            if len(assembled_chains) != len(native_r + native_l) or set(assembled_chains) != set(native_r + native_l):
                raise ValueError(f"{dataset_row_id}: assembled chains {assembled_chains!r} do not match {native_r}:{native_l}")
            record.update(
                {
                    "complex": complex_id,
                    "native_receptor_chains": native_r,
                    "native_ligand_chains": native_l,
                    "native_receptor_path": str(role_paths["native_receptor"]),
                    "native_receptor_archive_member": role_records["native_receptor"].get("archive_member", ""),
                    "native_receptor_sha256": role_records["native_receptor"]["actual_sha256"],
                    "native_ligand_path": str(role_paths["native_ligand"]),
                    "native_ligand_archive_member": role_records["native_ligand"].get("archive_member", ""),
                    "native_ligand_sha256": role_records["native_ligand"]["actual_sha256"],
                    "native_pdb_path": str(native_path.resolve()),
                    "native_pdb_sha256": sha256_file(native_path),
                    "assembly_status": "assembled",
                }
            )
        except Exception as exc:
            record["assembly_status"] = "assembly_failed"
            record["assembly_error"] = str(exc)
        assembly_rows.append(record)
        assembly_by_id[dataset_row_id] = record

    enriched_rows: list[dict[str, str]] = []
    for stage in stages:
        row = dict(stage)
        assembly = assembly_by_id.get(stage.get("dataset_row_id", ""), {})
        row_id = stage.get("dataset_row_id", "")
        source_gate_status = (
            "audit_only" if row_id in excluded else "strict_clean" if row_id else "unresolved"
        )
        row.update(
            {
                "source_gate_status": assembly.get("source_gate_status", source_gate_status),
                "native_assembly_status": assembly.get("assembly_status", "not_requested"),
                "native_assembly_error": assembly.get("assembly_error", ""),
                "native_pdb_path": assembly.get("native_pdb_path", ""),
                "native_pdb_sha256": assembly.get("native_pdb_sha256", ""),
                "native_receptor_archive_member": assembly.get("native_receptor_archive_member", ""),
                "native_receptor_sha256": assembly.get("native_receptor_sha256", ""),
                "native_ligand_archive_member": assembly.get("native_ligand_archive_member", ""),
                "native_ligand_sha256": assembly.get("native_ligand_sha256", ""),
            }
        )
        enriched_rows.append(row)

    write_table(assembly_manifest, assembly_rows, "\t")
    write_table(enriched_stage_manifest, enriched_rows, ",")
    return assembly_rows, enriched_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-manifest", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--staged-sources", type=Path, required=True)
    parser.add_argument("--source-gate-policy", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--assembly-manifest", type=Path, required=True)
    parser.add_argument("--enriched-stage-manifest", type=Path, required=True)
    args = parser.parse_args()
    assembly_rows, _ = assemble_native_complexes(**vars(args))
    counts: defaultdict[str, int] = defaultdict(int)
    for row in assembly_rows:
        counts[row["assembly_status"]] += 1
    print(" ".join(f"{key}={value}" for key, value in sorted(counts.items())))
    return 0 if counts["assembled"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
