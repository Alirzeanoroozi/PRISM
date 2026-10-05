#!/usr/bin/env python3
"""Build deterministic source and Biopython structure-validation manifests.

The script deliberately keeps two source contracts separate:

* raw receptor/ligand selectors are resolved only against ``benchmark/data/pdbs``;
* native receptor/ligand constituents are resolved only as exact members of
  ``benchmark/originals/benchmark5.5.tgz``.

There is no source fallback. Missing or unparsable structures remain represented
as rows, and ``--strict`` returns a non-zero status after writing both TSVs.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import re
import sys
import tarfile
import warnings
from pathlib import Path
from typing import Any, Iterable


DATASETS = (
    ("rigid", Path("benchmark/data/T_Rigid.csv")),
    ("medium", Path("benchmark/data/T_medium.csv")),
    ("difficult", Path("benchmark/data/T_difficult.csv")),
)
# The archive README documents these synthetic identifiers for repeated
# complexes.  They are candidate prefixes, not tar directory names.
ARCHIVE_PREFIX_ALIASES = {
    "1QFW": ("1QFW", "9QFW"),
    "1OYV": ("1OYV", "BOYV"),
    "3P57": ("3P57", "BP57", "CP57"),
    "3AAD": ("3AAD", "BAAD"),
}
EXPECTED_HEADER = (
    "Complex",
    "Cat.",
    "PDB ID 1",
    "Protein 1",
    "PDB ID 2",
    "Protein 2",
    "I-RMSD (Å)",
    "ΔASA(Å2)",
    "BM version introduced",
)

SOURCE_FIELDS = (
    "dataset_row_id",
    "dataset",
    "difficulty",
    "source_file",
    "source_row_number",
    "native_complex",
    "native_pdb_id",
    "native_receptor_chains_raw",
    "native_receptor_chains",
    "native_ligand_chains_raw",
    "native_ligand_chains",
    "raw_receptor_selector",
    "receptor_selector_pdb_id",
    "receptor_selector_chains_raw",
    "receptor_selector_chains",
    "raw_ligand_selector",
    "ligand_selector_pdb_id",
    "ligand_selector_chains_raw",
    "ligand_selector_chains",
    "source_role",
    "source_kind",
    "candidate_group_id",
    "candidate_count",
    "candidate_status",
    "source_path",
    "archive_member",
    "expected_archive_member",
    "archive_prefix",
    "exists",
    "size",
    "sha256",
    "resolution_status",
    "resolution_error",
    "source_row_json",
    "cohort",
    "source_scope",
    "audit_only",
    "pipeline_input",
    "truth_source",
    "archive_prefix_required",
)

VALIDATION_FIELDS = (
    "dataset_row_id",
    "dataset",
    "difficulty",
    "source_file",
    "source_row_number",
    "native_complex",
    "raw_receptor_selector",
    "raw_ligand_selector",
    "source_role",
    "source_kind",
    "cohort",
    "source_scope",
    "audit_only",
    "pipeline_input",
    "truth_source",
    "archive_prefix_required",
    "source_path",
    "archive_member",
    "archive_prefix",
    "sha256",
    "parse_status",
    "parser_name",
    "biopython_version",
    "model_count",
    "chain_count",
    "chain_ids",
    "polymer_chain_count",
    "polymer_chain_ids",
    "residue_count",
    "atom_count",
    "hetero_atom_count",
    "warning_count",
    "warning_types",
    "expected_chain_ids",
    "chain_set_status",
    "model_ids",
    "residue_ids_sha256",
    "sequence_hashes",
    "altloc_codes",
    "altloc_count",
    "disordered_atom_count",
    "disordered_residue_count",
    "duplicate_atom_id_count",
    "duplicate_residue_id_count",
    "error",
)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _clean_chain_assignment(value: str) -> str:
    """Remove benchmark residue annotations while preserving chain characters."""

    without_annotations = re.sub(r"\([^)]*\)", "", value)
    return "".join(character for character in without_annotations if character.isalnum())


def _parse_selector(raw: str) -> dict[str, str]:
    value = raw.strip()
    pdb_id, separator, chains = value.partition("_")
    return {
        "pdb_id": pdb_id.strip(),
        "chains_raw": chains if separator else "",
        "chains": _clean_chain_assignment(chains if separator else ""),
    }


def _parse_native_complex(raw: str) -> dict[str, str]:
    value = raw.strip()
    match = re.fullmatch(r"([A-Za-z0-9]{4})_([^:]*):(.*)", value)
    if match is None:
        return {
            "pdb_id": "",
            "receptor_chains_raw": "",
            "receptor_chains": "",
            "ligand_chains_raw": "",
            "ligand_chains": "",
            "error": "invalid_native_complex",
        }
    receptor_raw = match.group(2)
    ligand_raw = match.group(3).strip()
    return {
        "pdb_id": match.group(1),
        "receptor_chains_raw": receptor_raw,
        "receptor_chains": _clean_chain_assignment(receptor_raw),
        "ligand_chains_raw": ligand_raw,
        "ligand_chains": _clean_chain_assignment(ligand_raw),
        "error": "",
    }


def _read_rows(
    repo_root: Path,
    limit: int | None,
    dataset_row_ids: Iterable[str] | None = None,
) -> list[dict[str, Any]]:
    requested_ids = {str(value).strip() for value in (dataset_row_ids or ()) if str(value).strip()}
    rows: list[dict[str, Any]] = []
    for dataset, relative_path in DATASETS:
        path = repo_root / relative_path
        if not path.is_file():
            raise FileNotFoundError(f"required benchmark CSV is missing: {path}")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            try:
                header = tuple(next(reader))
            except StopIteration as exc:
                raise ValueError(f"benchmark CSV is empty: {path}") from exc
            if header != EXPECTED_HEADER:
                raise ValueError(f"unexpected header in {path}: {header!r}")
            dataset_row_number = 0
            for source_row_number, values in enumerate(reader, start=2):
                if not values or all(value == "" for value in values):
                    continue
                if len(values) != len(header):
                    raise ValueError(
                        f"row {source_row_number} in {path} has {len(values)} fields; expected {len(header)}"
                    )
                if limit is not None and not requested_ids and len(rows) >= limit:
                    return rows
                dataset_row_number += 1
                source_row = dict(zip(header, values))
                receptor = _parse_selector(source_row["PDB ID 1"])
                ligand = _parse_selector(source_row["PDB ID 2"])
                native = _parse_native_complex(source_row["Complex"])
                rows.append(
                    {
                        "dataset_row_id": f"{dataset}:{dataset_row_number:06d}",
                        "dataset": dataset,
                        "difficulty": dataset,
                        "source_file": str(relative_path),
                        "source_row_number": str(source_row_number),
                        "source_row": source_row,
                        "source_row_json": json.dumps(source_row, ensure_ascii=True, separators=(",", ":")),
                        "native_complex": source_row["Complex"],
                        "native": native,
                        "receptor": receptor,
                        "ligand": ligand,
                    }
                )
    if requested_ids:
        rows = [row for row in rows if row["dataset_row_id"] in requested_ids]
        unknown = sorted(requested_ids - {row["dataset_row_id"] for row in rows})
        if unknown:
            raise ValueError("unknown dataset_row_id(s): " + ", ".join(unknown))
    if limit is not None:
        rows = rows[:limit]
    return rows


def _local_candidates(pdb_root: Path, selector: str) -> list[Path]:
    parsed = _parse_selector(selector)
    normalized = re.sub(r"\([^)]*\)", "", selector.strip()).casefold()
    pdb_id = parsed["pdb_id"].casefold()
    candidates: list[Path] = []
    if not pdb_root.is_dir():
        return candidates
    for path in sorted(pdb_root.rglob("*.pdb"), key=lambda item: item.relative_to(pdb_root).as_posix().casefold()):
        if path.stem.casefold() == normalized:
            candidates.append(path)
            continue
        if not parsed["chains"] and path.stem.casefold() == pdb_id:
            candidates.append(path)
            continue
        if path.parent.name.casefold() == normalized and path.parent.parent.name.casefold() == "chainwise":
            candidates.append(path)
    return candidates


def _archive_index(archive_path: Path, wanted_stems: set[str]) -> tuple[dict[str, bytes], set[str], str]:
    """Return matching member bytes, observed prefixes, and an error string."""

    if not archive_path.is_file():
        return {}, set(), "archive_missing"
    members: dict[str, bytes] = {}
    prefixes: set[str] = set()
    try:
        with tarfile.open(archive_path, "r:gz") as archive:
            for member in archive.getmembers():
                name = member.name
                if not member.isfile() or "/structures/" not in name or not name.lower().endswith(".pdb"):
                    continue
                prefix, _, relative = name.partition("/")
                file_prefix = Path(relative).stem.rsplit("_", 2)[0]
                if Path(relative).name.casefold() not in wanted_stems:
                    continue
                prefixes.add(file_prefix)
                extracted = archive.extractfile(member)
                if extracted is None:
                    continue
                members[name] = extracted.read()
    except (OSError, tarfile.TarError) as exc:
        return {}, prefixes, f"archive_read_error:{type(exc).__name__}:{exc}"
    return members, prefixes, ""


def _biopython_parser_metadata(payload: bytes) -> dict[str, str]:
    try:
        import Bio
        from Bio.PDB import PDBParser
    except ImportError as exc:  # pragma: no cover - exercised only in incomplete environments
        return {
            "parse_status": "parse_error",
            "parser_name": "Bio.PDB.PDBParser",
            "biopython_version": "",
            "error": f"biopython_unavailable:{exc}",
        }

    try:
        from Bio.Data.IUPACData import protein_letters_3to1_extended
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            structure = PDBParser(QUIET=False, PERMISSIVE=True).get_structure("provenance", io.StringIO(payload.decode("utf-8", errors="replace")))
        models = list(structure.get_models())
        chains = list(structure.get_chains())
        residues = list(structure.get_residues())
        atoms = list(structure.get_atoms())
        hetero_atoms = [atom for atom in atoms if atom.get_parent().id[0] != " "]
        residue_tokens = []
        sequence_hashes: dict[str, str] = {}
        residue_id_counts: Counter[tuple[str, str, str, int, str]] = Counter()
        polymer_chain_ids: list[str] = []
        disordered_residues = 0
        disordered_atoms = 0
        one_letter = {key.upper(): value for key, value in protein_letters_3to1_extended.items()}
        for model in models:
            for chain in model.get_chains():
                if any(residue.id[0] == " " for residue in chain.get_residues()):
                    if chain.id not in polymer_chain_ids:
                        polymer_chain_ids.append(chain.id)
                sequence = []
                for residue in chain.get_residues():
                    hetflag, resseq, insertion = residue.id
                    key = (str(model.id), str(chain.id), str(hetflag), int(resseq), str(insertion or ""))
                    residue_id_counts[key] += 1
                    residue_tokens.append(
                        "|".join(
                            (
                                str(model.id),
                                str(chain.id),
                                str(hetflag),
                                str(resseq),
                                str(insertion or ""),
                                str(residue.resname).strip().upper(),
                            )
                        )
                    )
                    if residue.is_disordered():
                        disordered_residues += 1
                    if hetflag == " ":
                        sequence.append(one_letter.get(str(residue.resname).strip().upper(), "X"))
                    disordered_atoms += sum(1 for atom in residue.get_atoms() if atom.is_disordered())
                sequence_hashes[f"{model.id}:{chain.id}"] = _sha256("".join(sequence).encode("ascii"))
        raw_atom_keys: list[tuple[str, str, str, str, str, str, str]] = []
        altloc_codes: set[str] = set()
        current_model = "0"
        for line in payload.decode("utf-8", errors="replace").splitlines():
            record = line[:6].strip().upper()
            if record == "MODEL":
                current_model = line[10:14].strip() or "0"
                continue
            if record not in {"ATOM", "HETATM"} or len(line) < 27:
                continue
            altloc = line[16].strip()
            if altloc:
                altloc_codes.add(altloc)
            raw_atom_keys.append(
                (
                    current_model,
                    line[21].strip(),
                    line[22:26].strip(),
                    line[26].strip(),
                    line[12:16].strip(),
                    altloc,
                    record,
                )
            )
        warning_types = sorted({type(item.message).__name__ for item in caught})
        return {
            "parse_status": "ok",
            "parser_name": "Bio.PDB.PDBParser",
            "biopython_version": Bio.__version__,
            "model_count": str(len(models)),
            "chain_count": str(len(chains)),
            "chain_ids": ",".join(chain.id for chain in chains),
            "polymer_chain_count": str(len(polymer_chain_ids)),
            "polymer_chain_ids": ",".join(polymer_chain_ids),
            "residue_count": str(len(residues)),
            "atom_count": str(len(atoms)),
            "hetero_atom_count": str(len(hetero_atoms)),
            "warning_count": str(len(caught)),
            "warning_types": ",".join(warning_types),
            "model_ids": ",".join(str(model.id) for model in models),
            "residue_ids_sha256": _sha256("\n".join(residue_tokens).encode("utf-8")),
            "sequence_hashes": json.dumps(sequence_hashes, sort_keys=True, separators=(",", ":")),
            "altloc_codes": ",".join(sorted(altloc_codes)),
            "altloc_count": str(len(altloc_codes)),
            "disordered_atom_count": str(disordered_atoms),
            "disordered_residue_count": str(disordered_residues),
            "duplicate_atom_id_count": str(sum(count - 1 for count in Counter(raw_atom_keys).values() if count > 1)),
            "duplicate_residue_id_count": str(sum(count - 1 for count in residue_id_counts.values() if count > 1)),
            "error": "",
        }
    except Exception as exc:  # Biopython has several parser-specific exception classes.
        return {
            "parse_status": "parse_error",
            "parser_name": "Bio.PDB.PDBParser",
            "biopython_version": Bio.__version__,
            "model_ids": "",
            "residue_ids_sha256": "",
            "sequence_hashes": "",
            "altloc_codes": "",
            "altloc_count": "",
            "disordered_atom_count": "",
            "disordered_residue_count": "",
            "duplicate_atom_id_count": "",
            "duplicate_residue_id_count": "",
            "error": f"{type(exc).__name__}:{exc}",
        }


def _expected_chains(row: dict[str, Any], role: str) -> str:
    if role in {"receptor_selector", "pipeline_receptor"}:
        return row["receptor"]["chains"]
    if role in {"ligand_selector", "pipeline_ligand"}:
        return row["ligand"]["chains"]
    if role in {"native_receptor", "curated_bound_receptor"}:
        return row["native"]["receptor_chains"]
    if role in {"native_ligand", "curated_bound_ligand"}:
        return row["native"]["ligand_chains"]
    return ""


def _empty_validation(
    row: dict[str, Any],
    source: dict[str, str],
    status: str,
    error: str,
    expected_chain_ids: str = "",
) -> dict[str, str]:
    values = {
        "dataset_row_id": row["dataset_row_id"],
        "dataset": row["dataset"],
        "difficulty": row["difficulty"],
        "source_file": row["source_file"],
        "source_row_number": row["source_row_number"],
        "native_complex": row["native_complex"],
        "raw_receptor_selector": row["source_row"]["PDB ID 1"],
        "raw_ligand_selector": row["source_row"]["PDB ID 2"],
        "source_role": source["source_role"],
        "source_kind": source["source_kind"],
        "cohort": source["cohort"],
        "source_scope": source["source_scope"],
        "audit_only": source["audit_only"],
        "pipeline_input": source["pipeline_input"],
        "truth_source": source["truth_source"],
        "archive_prefix_required": source["archive_prefix_required"],
        "source_path": source["source_path"],
        "archive_member": source["archive_member"],
        "archive_prefix": source["archive_prefix"],
        "sha256": source["sha256"],
        "parse_status": status,
        "parser_name": "Bio.PDB.PDBParser",
        "biopython_version": "",
        "model_count": "",
            "chain_count": "",
            "chain_ids": "",
            "polymer_chain_count": "",
            "polymer_chain_ids": "",
        "residue_count": "",
        "atom_count": "",
        "hetero_atom_count": "",
        "warning_count": "",
        "warning_types": "",
        "expected_chain_ids": expected_chain_ids,
        "chain_set_status": "unresolved" if status != "ok" else "unknown",
        "model_ids": "",
        "residue_ids_sha256": "",
        "sequence_hashes": "",
        "altloc_codes": "",
        "altloc_count": "",
        "disordered_atom_count": "",
        "disordered_residue_count": "",
        "duplicate_atom_id_count": "",
        "duplicate_residue_id_count": "",
        "error": error,
    }
    return values


def _source_context(row: dict[str, Any], role: str) -> dict[str, str]:
    source_scope = "audit" if role in {"receptor_selector", "ligand_selector"} else "pipeline"
    pipeline_input = ""
    truth_source = ""
    if role == "pipeline_receptor":
        pipeline_input = "receptor"
    elif role == "pipeline_ligand":
        pipeline_input = "ligand"
    elif role == "native_receptor":
        truth_source = "native_receptor"
    elif role == "native_ligand":
        truth_source = "native_ligand"
    return {
        "dataset_row_id": row["dataset_row_id"],
        "dataset": row["dataset"],
        "difficulty": row["difficulty"],
        "source_file": row["source_file"],
        "source_row_number": row["source_row_number"],
        "cohort": "repository_bm5_bm5_5_extension",
        "native_complex": row["native_complex"],
        "native_pdb_id": row["native"]["pdb_id"],
        "native_receptor_chains_raw": row["native"]["receptor_chains_raw"],
        "native_receptor_chains": row["native"]["receptor_chains"],
        "native_ligand_chains_raw": row["native"]["ligand_chains_raw"],
        "native_ligand_chains": row["native"]["ligand_chains"],
        "raw_receptor_selector": row["source_row"]["PDB ID 1"],
        "receptor_selector_pdb_id": row["receptor"]["pdb_id"],
        "receptor_selector_chains_raw": row["receptor"]["chains_raw"],
        "receptor_selector_chains": row["receptor"]["chains"],
        "raw_ligand_selector": row["source_row"]["PDB ID 2"],
        "ligand_selector_pdb_id": row["ligand"]["pdb_id"],
        "ligand_selector_chains_raw": row["ligand"]["chains_raw"],
        "ligand_selector_chains": row["ligand"]["chains"],
        "source_role": role,
        "source_scope": source_scope,
        "audit_only": "true" if source_scope == "audit" else "false",
        "pipeline_input": pipeline_input,
        "truth_source": truth_source,
        "archive_prefix_required": "true" if source_scope == "pipeline" else "false",
    }


def _candidate_row(
    row: dict[str, Any],
    role: str,
    source_kind: str,
    candidate_group_id: str,
    candidate_count: int,
    candidate_status: str,
    source_path: str,
    archive_member: str,
    expected_archive_member: str,
    archive_prefix: str,
    payload: bytes | None,
    resolution_status: str,
    resolution_error: str,
) -> tuple[dict[str, str], dict[str, str]]:
    source = _source_context(row, role)
    source.update(
        {
            "source_kind": source_kind,
            "candidate_group_id": candidate_group_id,
            "candidate_count": str(candidate_count),
            "candidate_status": candidate_status,
            "source_path": source_path,
            "archive_member": archive_member,
            "expected_archive_member": expected_archive_member,
            "archive_prefix": archive_prefix,
            "exists": "true" if payload is not None else "false",
            "size": str(len(payload)) if payload is not None else "",
            "sha256": _sha256(payload) if payload is not None else "",
            "resolution_status": resolution_status,
            "resolution_error": resolution_error,
            "source_row_json": row["source_row_json"],
        }
    )
    if payload is None:
        validation = _empty_validation(row, source, "unresolved", resolution_error, _expected_chains(row, role))
    else:
        expected = _expected_chains(row, role)
        validation = _empty_validation(row, source, "parse_error", "parser_metadata_unavailable", expected)
        validation.update(_biopython_parser_metadata(payload))
        observed = set(value for value in validation.get("polymer_chain_ids", "").split(",") if value)
        requested = set(expected)
        validation["chain_set_status"] = (
            "not_declared"
            if not requested
            else "ok"
            if requested == observed
            else "superset"
            if requested.issubset(observed)
            else "mismatch"
        )
    return source, validation


def build_manifests(
    repo_root: str | Path,
    limit: int | None = None,
    dataset_row_ids: Iterable[str] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    root = Path(repo_root).resolve()
    rows = _read_rows(root, limit, dataset_row_ids)
    local_root = root / "benchmark/data/pdbs"
    archive_path = root / "benchmark/originals/benchmark5.5.tgz"
    wanted_prefixes = {
        prefix
        for row in rows
        for prefix in ARCHIVE_PREFIX_ALIASES.get(
            row["native"]["pdb_id"].upper(),
            (row["native"]["pdb_id"].upper(),),
        )
    }
    wanted_stems = {
        f"{prefix}_{side}_{unit}.pdb".casefold()
        for prefix in wanted_prefixes
        for side in ("r", "l")
        for unit in ("b", "u")
    }
    archive_bytes, observed_prefixes, archive_error = _archive_index(archive_path, wanted_stems)
    configured_prefixes = sorted(observed_prefixes or wanted_prefixes or {archive_path.name.removesuffix(".tgz")})

    def archive_chain_ids(payload: bytes) -> set[str]:
        metadata = _biopython_parser_metadata(payload)
        return {value for value in metadata.get("polymer_chain_ids", "").split(",") if value}

    archive_files: dict[tuple[str, str, str], tuple[str, bytes]] = {}
    for member_name, payload in archive_bytes.items():
        filename = Path(member_name).name
        match = re.fullmatch(r"(.+)_([rl])_([bu])\.pdb", filename, flags=re.IGNORECASE)
        if match:
            archive_files[(match.group(1).upper(), match.group(2).lower(), match.group(3).lower())] = (member_name, payload)

    row_archive_prefixes: dict[str, list[str]] = {}
    row_archive_complete_prefixes: dict[str, list[str]] = {}
    for row in rows:
        expected_receptor = set(row["native"]["receptor_chains"])
        expected_ligand = set(row["native"]["ligand_chains"])
        complete_prefixes = []
        candidates = []
        for prefix in configured_prefixes:
            receptor = archive_files.get((prefix.upper(), "r", "b"))
            ligand = archive_files.get((prefix.upper(), "l", "b"))
            unbound_receptor = archive_files.get((prefix.upper(), "r", "u"))
            unbound_ligand = archive_files.get((prefix.upper(), "l", "u"))
            if None in (receptor, ligand, unbound_receptor, unbound_ligand):
                continue
            complete_prefixes.append(prefix)
            if archive_chain_ids(receptor[1]) == expected_receptor and archive_chain_ids(ligand[1]) == expected_ligand:
                candidates.append(prefix)
        row_archive_prefixes[row["dataset_row_id"]] = sorted(candidates)
        row_archive_complete_prefixes[row["dataset_row_id"]] = sorted(complete_prefixes)

    source_rows: list[dict[str, str]] = []
    validation_rows: list[dict[str, str]] = []
    for row in rows:
        role_specs = (
            ("receptor_selector", "benchmark_data_pdbs", row["source_row"]["PDB ID 1"], ""),
            ("ligand_selector", "benchmark_data_pdbs", row["source_row"]["PDB ID 2"], ""),
            ("pipeline_receptor", "benchmark5.5_archive_member", f"{row['native']['pdb_id']}_r_u.pdb", "r_u"),
            ("pipeline_ligand", "benchmark5.5_archive_member", f"{row['native']['pdb_id']}_l_u.pdb", "l_u"),
            ("native_receptor", "benchmark5.5_archive_member", f"{row['native']['pdb_id']}_r_b.pdb", "r_b"),
            ("native_ligand", "benchmark5.5_archive_member", f"{row['native']['pdb_id']}_l_b.pdb", "l_b"),
        )
        for role, source_kind, requested_name, _side in role_specs:
            group_id = f"{row['dataset_row_id']}:{role}"
            if source_kind == "benchmark_data_pdbs":
                candidates = _local_candidates(local_root, requested_name)
                candidate_count = len(candidates)
                if candidates:
                    status = "collision" if candidate_count > 1 else "unique"
                    for candidate in candidates:
                        payload = candidate.read_bytes()
                        source, validation = _candidate_row(
                            row,
                            role,
                            source_kind,
                            group_id,
                            candidate_count,
                            status,
                            candidate.relative_to(root).as_posix(),
                            "",
                            "",
                            "",
                            payload,
                            "resolved",
                            "",
                        )
                        source_rows.append(source)
                        validation_rows.append(validation)
                else:
                    source, validation = _candidate_row(
                        row,
                        role,
                        source_kind,
                        group_id,
                        0,
                        "unresolved",
                        "",
                        "",
                        "",
                        "",
                        None,
                        "unresolved",
                        "local_selector_missing",
                    )
                    source_rows.append(source)
                    validation_rows.append(validation)
                continue

            role_side, role_unit = _side.split("_")
            matching_prefixes = row_archive_prefixes[row["dataset_row_id"]]
            complete_prefixes = row_archive_complete_prefixes[row["dataset_row_id"]]
            prefixes_to_use = matching_prefixes or complete_prefixes
            matches = []
            for prefix in prefixes_to_use:
                item = archive_files.get((prefix.upper(), role_side, role_unit))
                if item is not None:
                    matches.append((prefix, item[0], item[1]))
            if matches:
                candidate_count = len(matches)
                status = "collision" if candidate_count > 1 else "unique"
                resolution_error = "" if matching_prefixes else "archive_bound_chain_assignment_mismatch"
                for prefix, member_name, payload in matches:
                    source, validation = _candidate_row(
                        row,
                        role,
                        source_kind,
                        group_id,
                        candidate_count,
                        status,
                        archive_path.relative_to(root).as_posix(),
                        member_name,
                        member_name,
                        prefix,
                        payload,
                        "resolved",
                        resolution_error,
                    )
                    source_rows.append(source)
                    validation_rows.append(validation)
            else:
                # Keep an unresolved row for each known archive prefix. This preserves
                # the expected source identity without inventing a local substitute.
                candidate_count = 0
                error = archive_error or "archive_member_missing"
                fallback_prefixes = row_archive_prefixes[row["dataset_row_id"]] or row_archive_complete_prefixes[row["dataset_row_id"]]
                if not fallback_prefixes:
                    fallback_prefixes = list(
                        ARCHIVE_PREFIX_ALIASES.get(
                            row["native"]["pdb_id"].upper(),
                            (row["native"]["pdb_id"].upper(),),
                        )
                    )[:1]
                for prefix in fallback_prefixes:
                    member_name = (
                        f"{archive_path.name.removesuffix('.tgz')}/structures/"
                        f"{prefix}_{role_side}_{role_unit}.pdb"
                    )
                    source, validation = _candidate_row(
                        row,
                        role,
                        source_kind,
                        group_id,
                        candidate_count,
                        "unresolved",
                        archive_path.relative_to(root).as_posix() if archive_path.exists() else "",
                        "",
                        member_name,
                        prefix,
                        None,
                        "unresolved",
                        error,
                    )
                    source_rows.append(source)
                    validation_rows.append(validation)

    return source_rows, validation_rows


def _write_tsv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."), help="PRISM-prescript repository root")
    parser.add_argument("--output-dir", type=Path, required=True, help="directory for source_manifest.tsv and structure_validation.tsv")
    parser.add_argument("--strict", action="store_true", help="return 2 when any source is unresolved or any structure fails parsing")
    parser.add_argument("--limit", type=int, default=None, help="maximum number of benchmark rows to read, in dataset/file order")
    parser.add_argument(
        "--dataset-row-id",
        action="append",
        default=[],
        help="select an exact dataset_row_id; repeat for multiple rows (for example rigid:000001)",
    )
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be a positive integer")

    try:
        source_rows, validation_rows = build_manifests(
            args.repo_root,
            args.limit,
            args.dataset_row_id,
        )
    except (FileNotFoundError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    output_dir = args.output_dir
    _write_tsv(output_dir / "source_manifest.tsv", SOURCE_FIELDS, source_rows)
    _write_tsv(output_dir / "structure_validation.tsv", VALIDATION_FIELDS, validation_rows)
    unresolved = sum(
        row["source_scope"] == "pipeline"
        and (
            row["resolution_status"] != "resolved"
            or row["candidate_status"] != "unique"
            or not row["archive_prefix"]
        )
        for row in source_rows
    )
    parse_failures = sum(
        row["source_scope"] == "pipeline"
        and (row["parse_status"] != "ok" or row["chain_set_status"] not in {"ok", "not_declared"})
        for row in validation_rows
    )
    expected_curated_roles = {"pipeline_receptor", "pipeline_ligand", "native_receptor", "native_ligand"}
    curated_roles = {}
    for row in source_rows:
        if row["source_scope"] == "pipeline":
            curated_roles.setdefault(row["dataset_row_id"], set()).add(row["source_role"])
    contract_failures = sum(set(roles) != expected_curated_roles for roles in curated_roles.values())
    print(f"wrote {len(source_rows)} source rows and {len(validation_rows)} validation rows to {output_dir}")
    print(f"unresolved={unresolved} parse_failures={parse_failures} contract_failures={contract_failures}")
    return 2 if args.strict and (unresolved or parse_failures or contract_failures) else 0


if __name__ == "__main__":
    raise SystemExit(main())
