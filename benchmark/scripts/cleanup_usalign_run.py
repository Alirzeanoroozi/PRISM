#!/usr/bin/env python3
"""Dry-run/apply guarded cleanup for the validated USalign run namespace."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_artifact(status_path: Path, status: dict[str, Any], path_key: str, hash_key: str, parent: Path) -> Path:
    path = Path(str(status.get(path_key, ""))).resolve()
    if not path.is_file() or path.parent != parent.resolve():
        raise RuntimeError(f"unsafe retained artifact: {path}")
    if sha256_file(path) != status.get(hash_key):
        raise RuntimeError(f"retained artifact hash mismatch: {path}")
    return path


def validate(root: Path, refinement_root: Path, final_package: Path) -> dict[str, Any]:
    if root.name != "usalign_production_19855":
        raise RuntimeError(f"unexpected USalign cleanup root: {root}")
    aggregate_dir = root / "aggregate"
    aggregate_status_path = aggregate_dir / "aggregation_status.json"
    if not aggregate_status_path.is_file():
        raise RuntimeError("USalign batch aggregation status is missing")
    aggregate_status = json.loads(aggregate_status_path.read_text(encoding="utf-8"))
    if aggregate_status.get("status") != "validated_compacted":
        raise RuntimeError("USalign batch aggregation is not validated_compacted")
    for path_key, hash_key in (
        ("candidate_generated_path", "candidate_generated_sha256"),
        ("transformed_dockq_path", "transformed_dockq_sha256"),
        ("transformed_interfaces_path", "transformed_interfaces_sha256"),
    ):
        if path_key not in aggregate_status or hash_key not in aggregate_status:
            raise RuntimeError(f"USalign aggregate is missing retained path/hash: {path_key}")
        checked_artifact(aggregate_status_path, aggregate_status, path_key, hash_key, aggregate_dir)

    handoff = root / "usalign_common_refinement_19855" / "handoff_status.json"
    if not handoff.is_file() or json.loads(handoff.read_text(encoding="utf-8")).get("status") != "validated_manifest_ready_for_resource_review":
        raise RuntimeError("USalign refinement handoff is not validated")

    refinement_status_path = refinement_root / "aggregate" / "aggregation_status.json"
    if not refinement_status_path.is_file():
        raise RuntimeError("USalign common-refinement aggregation is missing")
    refinement_status = json.loads(refinement_status_path.read_text(encoding="utf-8"))
    if refinement_status.get("status") != "complete" or not refinement_status.get("cleanup_eligible"):
        raise RuntimeError("USalign common refinement is not cleanup-eligible")
    for path_key, hash_key in (("results_path", "results_sha256"), ("interfaces_path", "interfaces_sha256")):
        checked_artifact(refinement_status_path, refinement_status, path_key, hash_key, refinement_root / "aggregate")

    final_manifest = final_package / "aggregate_manifest.json"
    if not final_manifest.is_file() or json.loads(final_manifest.read_text(encoding="utf-8")).get("status") != "validated_compacted":
        raise RuntimeError("final matched package is not validated_compacted")
    return {"batch_aggregation": aggregate_status, "refinement_aggregation": refinement_status, "final_manifest": str(final_manifest)}


def planned_paths(root: Path, refinement_root: Path) -> list[Path]:
    current = root / "current"
    paths: list[Path] = []
    for path in sorted(current.glob("batch_*")):
        if not path.is_dir():
            continue
        if path.is_symlink() or path.resolve().parent != current.resolve():
            raise RuntimeError(f"unsafe USalign batch cleanup target: {path}")
        paths.append(path)
    for relative in ("results", "adapter_inputs"):
        directory = refinement_root / relative
        if directory.is_dir():
            if directory.is_symlink() or directory.resolve().parent != refinement_root.resolve():
                raise RuntimeError(f"unsafe refinement cleanup target: {directory}")
            paths.append(directory)
    return paths


def write_record(path: Path, record: dict[str, Any]) -> None:
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def cleanup(root: Path, refinement_root: Path, final_package: Path, *, apply: bool) -> dict[str, Any]:
    root = root.resolve(); refinement_root = refinement_root.resolve(); final_package = final_package.resolve()
    evidence = validate(root, refinement_root, final_package)
    paths = planned_paths(root, refinement_root)
    record: dict[str, Any] = {
        "schema_version": "prism-usalign-cleanup/v1",
        "status": "applied" if apply else "dry_run",
        "root": str(root),
        "refinement_root": str(refinement_root),
        "final_package": str(final_package),
        "planned_paths": [str(p) for p in paths],
        "planned_count": len(paths),
        "deleted_paths": [],
        "evidence": evidence,
        "finished_at": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path = root / "cleanup_manifest.json"
    write_record(manifest_path, record)
    if apply:
        for path in paths:
            if path.is_dir(): shutil.rmtree(path)
            elif path.is_file(): path.unlink()
            record["deleted_paths"].append(str(path))
        record["finished_at"] = datetime.now(timezone.utc).isoformat()
        write_record(manifest_path, record)
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--refinement-root", type=Path, required=True)
    parser.add_argument("--final-package", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    print(json.dumps(cleanup(args.root, args.refinement_root, args.final_package, apply=args.apply), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
