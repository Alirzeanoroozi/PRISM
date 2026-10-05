#!/usr/bin/env python3
"""Write the stable TSV/PDB artifacts used by the PRISM investigation.

The production pipeline can continue to emit its native files.  These helpers
provide a small adapter layer that converts those outputs into the common
investigation contract without changing candidate generation or scoring.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from benchmark.scripts.investigation_lineage import (
    LINEAGE_RECORD_FIELDS,
    LineageRecord,
    aggregate_pair_summary,
    copy_pdb_immutable,
    coordinate_sha256,
)
from benchmark.scripts.investigation_contracts import (
    SUMMARY_CONTEXT_FIELDS,
    STRUCTURAL_METRIC_FIELDS,
    summarize_top20,
)
from benchmark.scripts.standardized_evaluator import (
    DOCKQ_METRIC_FIELDS,
    evaluate_dockq_json,
)


POSE_FIELDS = (
    "pose_id",
    "pair_id",
    "template_id",
    "side",
    "alignment_id",
    "orientation",
    "filter_stage",
    "refinement_id",
    "energy_id",
    "coordinate_sha256",
    "artifact_path",
    "ranking_inputs",
    "status",
    "failure_reason",
)
PAIR_SUMMARY_FIELDS = (
    "pair_id",
    "status",
    "failure_reason",
    "primary_model_id",
    "primary_dockq",
    "primary_irmsd",
    "primary_energy",
    "dockq_mean",
    "dockq_best",
    "irmsd_mean",
    "irmsd_best",
    "energy_mean",
    "energy_best",
    "model_count",
    "scoreable_model_count",
    "secondary_model_count",
    "non_scoreable_model_count",
    "status_counts",
)
SCORE_FIELDS = (
    "record_type",
    "interface",
    "GlobalDockQ",
    *DOCKQ_METRIC_FIELDS,
    "grouped_iRMSD",
    "raw_json_sha256",
    "raw_json_path",
    "mapping_status",
)
CONTRACT_PAIR_SUMMARY_FIELDS = (
    *SUMMARY_CONTEXT_FIELDS,
    "selected_count",
    "model_valid_count",
    "scoreable_model_count",
    "best_GlobalDockQ_at_20",
    *STRUCTURAL_METRIC_FIELDS,
)


def _value(row: Mapping[str, Any] | LineageRecord, key: str, default: Any = "") -> Any:
    if isinstance(row, LineageRecord):
        return getattr(row, key, default)
    return row.get(key, default)


def _serialize(value: Any) -> Any:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def write_tsv(path: str | Path, fields: Iterable[str], rows: Iterable[Mapping[str, Any]]) -> Path:
    """Write deterministic, schema-first TSV rows."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    field_list = tuple(fields)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=field_list, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        ordered_rows = sorted(
            list(rows),
            key=lambda row: tuple(str(_serialize(row.get(field, ""))) for field in field_list),
        )
        for row in ordered_rows:
            writer.writerow({field: _serialize(row.get(field, "")) for field in field_list})
    return output


def write_lineage_tsv(path: str | Path, records: Iterable[LineageRecord | Mapping[str, Any]]) -> Path:
    rows = []
    for record in records:
        rows.append({field: _value(record, field) for field in LINEAGE_RECORD_FIELDS})
    return write_tsv(path, LINEAGE_RECORD_FIELDS, rows)


def write_poses_tsv(path: str | Path, rows: Iterable[LineageRecord | Mapping[str, Any]]) -> Path:
    output_rows = []
    for row in rows:
        output_rows.append({field: _value(row, field) for field in POSE_FIELDS})
    return write_tsv(path, POSE_FIELDS, output_rows)


def freeze_pose_artifacts(rows: Iterable[Mapping[str, Any]], destination_root: str | Path) -> list[dict[str, Any]]:
    """Copy pose PDBs into an immutable run directory and return updated rows."""

    destination = Path(destination_root)
    frozen = []
    for row in rows:
        source = Path(str(row.get("artifact_path", "")))
        pose_id = str(row.get("pose_id", source.stem))
        if not source.is_file():
            raise FileNotFoundError(source)
        if not pose_id or pose_id in {".", ".."} or Path(pose_id).name != pose_id:
            raise ValueError(f"unsafe pose_id: {pose_id!r}")
        target = destination / f"{pose_id}.pdb"
        copy_pdb_immutable(source, target)
        updated = dict(row)
        updated["artifact_path"] = str(target)
        updated["coordinate_sha256"] = coordinate_sha256(target)
        frozen.append(updated)
    return frozen


def write_dockq_tsv(
    source: str | Path | Mapping[str, Any],
    global_path: str | Path,
    interfaces_path: str | Path,
    *,
    raw_json_path: str | Path | None = None,
    mapping_status: str = "frozen",
    grouped_irmsd: float | None = None,
) -> tuple[Path, Path]:
    """Normalize one DockQ JSON document into global/interface TSV files."""

    records = evaluate_dockq_json(source, grouped_irmsd=grouped_irmsd)
    common = {
        "raw_json_path": str(raw_json_path or ""),
        "mapping_status": mapping_status,
    }
    global_rows = []
    interface_rows = []
    for record in records:
        row = dict(record)
        row.update(common)
        (global_rows if row["record_type"] == "global" else interface_rows).append(row)
    return (
        write_tsv(global_path, SCORE_FIELDS, global_rows),
        write_tsv(interfaces_path, SCORE_FIELDS, interface_rows),
    )


def write_pair_summary_tsv(path: str | Path, model_rows: Iterable[Mapping[str, Any]]) -> Path:
    summaries = aggregate_pair_summary(model_rows)
    return write_tsv(path, PAIR_SUMMARY_FIELDS, summaries)


def write_contract_pair_summary_tsv(
    path: str | Path,
    model_rows: Iterable[Mapping[str, Any]],
    *,
    foreign_keys: Mapping[str, Any],
    ranking_keys: Iterable[str] | Mapping[str, str] | None = None,
    k: int = 20,
) -> Path:
    """Write the fail-closed top-k summary used by confirmatory analysis.

    Unlike the legacy compatibility summary above, this function never turns
    absent structural metrics into zero and never allows native-derived values
    to enter ranking.  An empty candidate stream produces zero unconditional
    utility with null structural fields.
    """

    summary = summarize_top20(model_rows, foreign_keys=foreign_keys, ranking_keys=ranking_keys, k=k)
    return write_tsv(path, CONTRACT_PAIR_SUMMARY_FIELDS, [summary])
