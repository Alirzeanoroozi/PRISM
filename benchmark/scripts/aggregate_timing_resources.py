#!/usr/bin/env python3
"""Aggregate compact stage timing/resource records without reading structures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


def number(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def value(row: dict[str, Any], names: tuple[str, ...]) -> float | None:
    for name in names:
        parsed = number(row.get(name))
        if parsed is not None:
            return parsed
    return None


def write_table(path: Path, rows: list[dict[str, Any]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row}) or ["stage"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return sha256_file(path)


def parse_spec(raw: str) -> tuple[str, Path]:
    method, separator, path = raw.partition("=")
    if not separator or not method or not path:
        raise argparse.ArgumentTypeError("expected METHOD=PATH")
    return method, Path(path)


def timing_records(method: str, path: Path, kind: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in read_rows(path):
        if kind == "refinement":
            for prefix, label in (("fiberdock", "refinement_fiberdock"), ("external_rosetta", "refinement_external_rosetta")):
                elapsed = value(row, (f"{prefix}_elapsed_seconds",))
                if elapsed is not None:
                    records.append({"method": method, "stage": label, "wall_seconds": elapsed, "status": row.get(f"{prefix}_status", "")})
            continue
        stage = str(row.get("stage", row.get("stage_name", kind)) or kind)
        elapsed = value(row, ("wall_seconds", "elapsed_seconds", "wall_seconds_sum"))
        if elapsed is None:
            continue
        cores = value(row, ("cpus_per_task", "cpu_count", "cores", "workers"))
        gpu_hours = value(row, ("gpu_hours", "gpu_hours_used"))
        records.append({
            "method": method,
            "stage": stage,
            "wall_seconds": elapsed,
            "cores": "" if cores is None else cores,
            "gpu_hours": "" if gpu_hours is None else gpu_hours,
            "status": row.get("status", ""),
        })
    return records


def aggregate(output_dir: Path, timing_specs: list[tuple[str, Path]], refinement_specs: list[tuple[str, Path]]) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    inputs: list[dict[str, Any]] = []
    for kind, specs in (("timing", timing_specs), ("refinement", refinement_specs)):
        for method, path in specs:
            records.extend(timing_records(method, path, kind))
            inputs.append({"method": method, "kind": kind, "path": str(path.resolve()), "sha256": sha256_file(path)})

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in records:
        grouped[(str(row["method"]), str(row["stage"]))].append(row)
    summary: list[dict[str, Any]] = []
    for (method, stage), rows in sorted(grouped.items()):
        walls = [float(row["wall_seconds"]) for row in rows]
        cpu_hours = 0.0
        cpu_observed = False
        gpu_hours = 0.0
        gpu_observed = False
        for row in rows:
            cores = number(row.get("cores"))
            if cores is not None:
                cpu_hours += float(row["wall_seconds"]) * cores / 3600.0
                cpu_observed = True
            gpu = number(row.get("gpu_hours"))
            if gpu is not None:
                gpu_hours += gpu
                gpu_observed = True
        summary.append({
            "method": method,
            "stage": stage,
            "record_count": len(rows),
            "status_counts": json.dumps(dict(sorted(Counter(str(row.get("status", "")) for row in rows).items())), sort_keys=True),
            "wall_seconds_sum": sum(walls),
            "wall_seconds_mean": statistics.mean(walls),
            "wall_seconds_median": statistics.median(walls),
            "cpu_core_hours": cpu_hours if cpu_observed else "",
            "gpu_hours": gpu_hours if gpu_observed else "",
            "resource_status": "observed" if cpu_observed or gpu_observed else "timing_only",
        })

    output_dir.mkdir(parents=True, exist_ok=True)
    table_hash = write_table(output_dir / "timing_resources.tsv", summary)
    manifest = {
        "schema_version": "prism-timing-resources/v1",
        "status": "validated_compacted",
        "row_count": len(summary),
        "timing_record_count": len(records),
        "timing_resources_sha256": table_hash,
        "inputs": inputs,
        "missing_resource_fields_remain_empty": True,
    }
    (output_dir / "timing_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timing", action="append", type=parse_spec, default=[])
    parser.add_argument("--refinement", action="append", type=parse_spec, default=[])
    args = parser.parse_args()
    result = aggregate(
        args.output_dir.resolve(),
        [(m, p.resolve()) for m, p in args.timing],
        [(m, p.resolve()) for m, p in args.refinement],
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
