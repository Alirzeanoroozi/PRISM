#!/usr/bin/env python3
"""Stage the four hashed curated benchmark role files for each dataset row.

Only archive members recorded in ``source_manifest.tsv`` are extracted.  Local
PDB-cache candidates and downloaded/full-PDB comparators are intentionally not
eligible inputs for this staging step.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import tarfile
from pathlib import Path


ROLE_FIELDS = (
    "dataset_row_id",
    "source_role",
    "pipeline_input",
    "truth_source",
    "archive_prefix",
    "archive_member",
    "source_payload_sha256",
    "staged_path",
    "staged_sha256",
    "status",
    "error",
)
REQUIRED_ROLES = {"pipeline_receptor", "pipeline_ligand", "native_receptor", "native_ligand"}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def safe_component(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "dataset_row_id", "source_role", "source_scope", "source_kind",
        "archive_member", "expected_archive_member", "archive_prefix", "sha256",
        "resolution_status", "candidate_status",
    }
    missing = required - set(rows[0] if rows else ())
    if missing:
        raise ValueError("source manifest is missing columns: " + ", ".join(sorted(missing)))
    return rows


def stage_sources(
    manifest_path: str | Path,
    output_dir: str | Path,
    *,
    repo_root: str | Path,
    strict: bool = False,
) -> tuple[Path, int]:
    manifest = Path(manifest_path).resolve()
    output = Path(output_dir).resolve()
    root = Path(repo_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows = read_manifest(manifest)
    archive_cache: dict[Path, tarfile.TarFile] = {}
    staged_rows: list[dict[str, str]] = []
    failures = 0
    try:
        groups: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            if row.get("source_scope") == "pipeline":
                groups.setdefault(row["dataset_row_id"], []).append(row)
        for dataset_row_id, group in sorted(groups.items()):
            roles = [row["source_role"] for row in group]
            role_set = set(roles)
            if role_set != REQUIRED_ROLES or len(roles) != len(REQUIRED_ROLES):
                failures += 1
                for role in sorted(REQUIRED_ROLES - role_set):
                    staged_rows.append({
                        "dataset_row_id": dataset_row_id,
                        "source_role": role,
                        "pipeline_input": "receptor" if role == "pipeline_receptor" else "ligand" if role == "pipeline_ligand" else "",
                        "truth_source": role if role.startswith("native_") else "",
                        "archive_prefix": "",
                        "archive_member": "",
                        "source_payload_sha256": "",
                        "staged_path": "",
                        "staged_sha256": "",
                        "status": "missing_role",
                        "error": "required_curated_role_missing",
                    })
            for row in sorted(group, key=lambda item: item["source_role"]):
                record = {
                    "dataset_row_id": row["dataset_row_id"],
                    "source_role": row["source_role"],
                    "pipeline_input": row.get("pipeline_input", ""),
                    "truth_source": row.get("truth_source", ""),
                    "archive_prefix": row.get("archive_prefix", ""),
                    "archive_member": row.get("archive_member", ""),
                    "source_payload_sha256": row.get("sha256", ""),
                    "staged_path": "",
                    "staged_sha256": "",
                    "status": "",
                    "error": "",
                }
                if row.get("resolution_status") != "resolved" or row.get("candidate_status") != "unique":
                    failures += 1
                    record.update(status="source_not_unique_or_unresolved", error=row.get("resolution_error", "source_contract_failed"))
                    staged_rows.append(record)
                    continue
                if row.get("source_kind") != "benchmark5.5_archive_member":
                    failures += 1
                    record.update(status="source_kind_mismatch", error="curated_role_is_not_archive_member")
                    staged_rows.append(record)
                    continue
                archive_path = (root / row["source_path"]).resolve()
                member_name = row.get("archive_member", "")
                expected_member = row.get("expected_archive_member", "")
                expected_name = f"{row.get('archive_prefix', '')}_{'r' if row['source_role'] in {'pipeline_receptor', 'native_receptor'} else 'l'}_{'u' if row['source_role'].startswith('pipeline_') else 'b'}.pdb"
                if (
                    not archive_path.is_file()
                    or not member_name
                    or member_name != expected_member
                    or Path(member_name).name.casefold() != expected_name.casefold()
                ):
                    failures += 1
                    record.update(status="archive_member_contract_failed", error="archive_path_or_member_or_prefix_mismatch")
                    staged_rows.append(record)
                    continue
                archive = archive_cache.get(archive_path)
                if archive is None:
                    archive = tarfile.open(archive_path, "r:gz")
                    archive_cache[archive_path] = archive
                try:
                    member = archive.getmember(member_name)
                    extracted = archive.extractfile(member)
                    payload = extracted.read() if extracted is not None else None
                except (KeyError, OSError, tarfile.TarError) as exc:
                    payload = None
                    record["error"] = f"archive_extract_error:{type(exc).__name__}:{exc}"
                if payload is None:
                    failures += 1
                    record["status"] = "archive_extract_failed"
                    staged_rows.append(record)
                    continue
                observed_hash = sha256_bytes(payload)
                if observed_hash != row.get("sha256"):
                    failures += 1
                    record.update(status="hash_mismatch", staged_sha256=observed_hash, error="archive_bytes_do_not_match_manifest")
                    staged_rows.append(record)
                    continue
                row_dir = output / safe_component(dataset_row_id)
                row_dir.mkdir(parents=True, exist_ok=True)
                destination = row_dir / f"{safe_component(row['source_role'])}_{observed_hash}.pdb"
                if destination.exists() and destination.read_bytes() != payload:
                    raise FileExistsError(f"staged path contains different bytes: {destination}")
                destination.write_bytes(payload)
                record.update(staged_path=str(destination), staged_sha256=observed_hash, status="staged")
                staged_rows.append(record)
    finally:
        for archive in archive_cache.values():
            archive.close()
    staged = output / "staged_sources.tsv"
    with staged.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROLE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(staged_rows, key=lambda row: (row["dataset_row_id"], row["source_role"])))
    if strict and failures:
        raise ValueError(f"curated source staging failed for {failures} role/row checks")
    return staged, failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)
    try:
        staged, failures = stage_sources(
            args.source_manifest,
            args.output_dir,
            repo_root=args.repo_root,
            strict=args.strict,
        )
    except (OSError, ValueError, tarfile.TarError) as exc:
        print(f"error: {exc}")
        return 2
    print(f"wrote staged source manifest: {staged}")
    print(f"failures={failures}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
