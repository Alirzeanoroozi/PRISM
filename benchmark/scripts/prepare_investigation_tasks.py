#!/usr/bin/env python3
"""Create an explicit one-to-ten-row task manifest for isolated KUTEM execution.

The generated CSV is intentionally small and auditable.  Each array index is
bound to exactly one dataset row, and the command is expanded before writing so
the runner never has to infer which benchmark pair a task represents.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import shlex


FIELDS = (
    "array_index",
    "task_id",
    "command",
    "config_paths",
    "input_paths",
    "output_paths",
    "scientific_retry_id",
)


def build_task_rows(
    dataset_row_ids: list[str],
    command_template: str,
    *,
    config_paths: list[str] | None = None,
    input_paths: list[str] | None = None,
    output_paths: list[str] | None = None,
    scientific_retry_id: str = "0",
) -> list[dict[str, str]]:
    if not 1 <= len(dataset_row_ids) <= 10:
        raise ValueError("between one and ten dataset_row_ids are required for array 1-N")
    cleaned = [value.strip() for value in dataset_row_ids]
    if any(not value for value in cleaned):
        raise ValueError("dataset_row_ids must be non-empty")
    if len(set(cleaned)) != len(cleaned):
        raise ValueError("dataset_row_ids must be unique")
    if not command_template.strip():
        raise ValueError("command_template must be non-empty")

    rows = []
    for index, dataset_row_id in enumerate(cleaned, start=1):
        substitutions = {
            "array_index": shlex.quote(str(index)),
            "dataset_row_id": shlex.quote(dataset_row_id),
            "task_id": shlex.quote(f"{dataset_row_id}:task-{index:04d}"),
        }
        try:
            command = command_template.format(**substitutions)
        except KeyError as exc:
            raise ValueError(f"unsupported command-template field: {exc.args[0]}") from exc
        rows.append(
            {
                "array_index": str(index),
                "task_id": substitutions["task_id"],
                "command": command,
                "config_paths": ";".join(config_paths or []),
                "input_paths": ";".join(input_paths or []),
                "output_paths": ";".join(output_paths or []),
                "scientific_retry_id": scientific_retry_id,
            }
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dataset-row-id", action="append", required=True)
    parser.add_argument("--command-template", required=True)
    parser.add_argument("--config-path", action="append", default=[])
    parser.add_argument("--input-path", action="append", default=[])
    parser.add_argument("--output-path", action="append", default=[])
    parser.add_argument("--scientific-retry-id", default="0")
    args = parser.parse_args(argv)
    try:
        rows = build_task_rows(
            args.dataset_row_id,
            args.command_template,
            config_paths=args.config_path,
            input_paths=args.input_path,
            output_paths=args.output_path,
            scientific_retry_id=args.scientific_retry_id,
        )
    except ValueError as exc:
        parser.error(str(exc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} isolated tasks to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
