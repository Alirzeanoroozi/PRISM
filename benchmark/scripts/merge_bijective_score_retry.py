#!/usr/bin/env python3
"""Replace retried score rows into fresh merged TSVs without altering inputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def score_key(row: dict[str, str]) -> tuple[str, str]:
    return row.get("dataset_row_id", ""), row.get("source_model_sha256", "")


def write_tsv(path: Path, rows: list[dict[str, str]], preferred_fields: list[str]) -> None:
    fields = list(preferred_fields)
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or ["record_type"], delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def merge(
    base_models_path: Path,
    base_interfaces_path: Path,
    retry_models_path: Path,
    retry_interfaces_path: Path,
    output_models_path: Path,
    output_interfaces_path: Path,
) -> tuple[int, int, int]:
    base_models = read_tsv(base_models_path)
    retry_models = read_tsv(retry_models_path)
    base_by_key = {score_key(row): row for row in base_models}
    retry_by_key = {score_key(row): row for row in retry_models}
    if len(base_by_key) != len(base_models) or any(not all(item) for item in base_by_key):
        raise ValueError("base model keys are missing or non-unique")
    if len(retry_by_key) != len(retry_models) or any(not all(item) for item in retry_by_key):
        raise ValueError("retry model keys are missing or non-unique")
    if not set(retry_by_key).issubset(base_by_key):
        raise ValueError("retry contains model keys absent from base scores")
    invalid_base = [key for key in retry_by_key if base_by_key[key].get("score_status") != "score_failed"]
    if invalid_base:
        raise ValueError(f"retry attempts to replace non-failed base rows: {invalid_base[:3]}")

    merged_models = [retry_by_key.get(score_key(row), row) for row in base_models]
    retried_paths = {row.get("staged_model_path", "") for row in retry_models}
    if "" in retried_paths:
        raise ValueError("retry model lacks staged_model_path")
    base_interfaces = read_tsv(base_interfaces_path)
    retry_interfaces = read_tsv(retry_interfaces_path)
    if any(row.get("staged_model_path", "") not in retried_paths for row in retry_interfaces):
        raise ValueError("retry interface row does not belong to a retried model")
    merged_interfaces = [
        row for row in base_interfaces if row.get("staged_model_path", "") not in retried_paths
    ] + retry_interfaces

    base_model_fields = list(base_models[0]) if base_models else []
    base_interface_fields = list(base_interfaces[0]) if base_interfaces else []
    write_tsv(output_models_path, merged_models, base_model_fields)
    write_tsv(output_interfaces_path, merged_interfaces, base_interface_fields)
    return len(merged_models), len(merged_interfaces), len(retry_models)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-models", type=Path, required=True)
    parser.add_argument("--base-interfaces", type=Path, required=True)
    parser.add_argument("--retry-models", type=Path, required=True)
    parser.add_argument("--retry-interfaces", type=Path, required=True)
    parser.add_argument("--output-models", type=Path, required=True)
    parser.add_argument("--output-interfaces", type=Path, required=True)
    args = parser.parse_args()
    models, interfaces, replaced = merge(
        args.base_models, args.base_interfaces, args.retry_models, args.retry_interfaces,
        args.output_models, args.output_interfaces,
    )
    print(f"models={models} interfaces={interfaces} replaced={replaced}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
