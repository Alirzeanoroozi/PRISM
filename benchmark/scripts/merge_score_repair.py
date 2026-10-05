#!/usr/bin/env python3
"""Replace failed canonical score rows with an isolated repair result."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(row: dict[str, str]) -> tuple[str, str]:
    return row.get("dataset_row_id", ""), row.get("source_model_sha256", "")


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def merge(
    base_models_path: Path,
    base_interfaces_path: Path,
    repair_models_path: Path,
    repair_interfaces_path: Path,
    output_root: Path,
    *,
    allowed_base_score_statuses: set[str] | None = None,
    required_base_field: str = "",
    required_base_value: str = "",
) -> dict[str, object]:
    if output_root.exists() and any(output_root.iterdir()):
        raise ValueError(f"refusing to overwrite non-empty output root: {output_root}")
    base_models = read_tsv(base_models_path)
    base_interfaces = read_tsv(base_interfaces_path)
    repair_models = read_tsv(repair_models_path)
    repair_interfaces = read_tsv(repair_interfaces_path)

    base_keys = [identity(row) for row in base_models]
    repair_keys = [identity(row) for row in repair_models]
    if any(not all(key) for key in base_keys + repair_keys):
        raise ValueError("model row lacks durable repair identity")
    if len(base_keys) != len(set(base_keys)):
        raise ValueError("duplicate base model identity")
    if len(repair_keys) != len(set(repair_keys)):
        raise ValueError("duplicate repair model identity")
    if not set(repair_keys) <= set(base_keys):
        raise ValueError("repair identity is absent from base scores")
    base_by_key = dict(zip(base_keys, base_models))
    allowed = allowed_base_score_statuses or {"score_failed"}
    invalid = [key for key in repair_keys if base_by_key[key].get("score_status") not in allowed]
    if invalid:
        raise ValueError(f"repair attempts to replace disallowed base rows: {invalid[:5]}")
    if required_base_field:
        mismatched = [
            key for key in repair_keys
            if base_by_key[key].get(required_base_field) != required_base_value
        ]
        if mismatched:
            raise ValueError(
                f"repair base rows do not satisfy {required_base_field}={required_base_value!r}: "
                f"{mismatched[:5]}"
            )

    repair_by_key = dict(zip(repair_keys, repair_models))
    merged_models = [repair_by_key.get(key, row) for key, row in zip(base_keys, base_models)]
    if [identity(row) for row in merged_models] != base_keys:
        raise ValueError("merged model identity/order changed")

    repair_hashes = {key[1] for key in repair_keys}
    retained_interfaces = [
        row for row in base_interfaces if row.get("source_model_sha256", "") not in repair_hashes
    ]
    if any(row.get("source_model_sha256", "") not in repair_hashes for row in repair_interfaces):
        raise ValueError("repair interface does not belong to a repair model")
    merged_interfaces = retained_interfaces + repair_interfaces

    output_root.mkdir(parents=True, exist_ok=True)
    model_output = output_root / "scores_models.tsv"
    interface_output = output_root / "scores_interfaces.tsv"
    write_tsv(model_output, merged_models)
    write_tsv(interface_output, merged_interfaces)
    manifest: dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "base_models": str(base_models_path.resolve()),
        "base_models_sha256": sha256_file(base_models_path),
        "base_interfaces": str(base_interfaces_path.resolve()),
        "base_interfaces_sha256": sha256_file(base_interfaces_path),
        "repair_models": str(repair_models_path.resolve()),
        "repair_models_sha256": sha256_file(repair_models_path),
        "repair_interfaces": str(repair_interfaces_path.resolve()),
        "repair_interfaces_sha256": sha256_file(repair_interfaces_path),
        "model_rows": len(merged_models),
        "interface_rows": len(merged_interfaces),
        "replaced_rows": len(repair_models),
        "score_status": dict(Counter(row.get("score_status", "") for row in merged_models)),
    }
    (output_root / "merge_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-models", type=Path, required=True)
    parser.add_argument("--base-interfaces", type=Path, required=True)
    parser.add_argument("--repair-models", type=Path, required=True)
    parser.add_argument("--repair-interfaces", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--expected-repair-rows", type=int)
    parser.add_argument(
        "--allow-base-score-status",
        action="append",
        dest="allowed_base_score_statuses",
        help="Base score status eligible for replacement; repeatable (default: score_failed).",
    )
    parser.add_argument("--required-base-field", default="")
    parser.add_argument("--required-base-value", default="")
    args = parser.parse_args()
    if bool(args.required_base_field) != bool(args.required_base_value):
        parser.error("--required-base-field and --required-base-value must be supplied together")
    manifest = merge(
        args.base_models,
        args.base_interfaces,
        args.repair_models,
        args.repair_interfaces,
        args.output_root,
        allowed_base_score_statuses=set(args.allowed_base_score_statuses or {"score_failed"}),
        required_base_field=args.required_base_field,
        required_base_value=args.required_base_value,
    )
    if args.expected_repair_rows is not None and manifest["replaced_rows"] != args.expected_repair_rows:
        raise SystemExit(
            f"expected {args.expected_repair_rows} repair rows; found {manifest['replaced_rows']}"
        )
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
