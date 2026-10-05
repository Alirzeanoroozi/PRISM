"""Explicit evaluator and lineage contracts for controlled investigations.

The historical helpers remain available for compatibility.  This module is
the stricter adapter used when a run needs auditable foreign keys, frozen
chain mappings, and native-independent top-k selection.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from benchmark.scripts.investigation_lineage import validate_ranking_keys
from benchmark.scripts.standardized_evaluator import (
    DOCKQ_METRIC_FIELDS,
    evaluate_dockq_json,
    validate_mapping,
)


FOREIGN_KEY_FIELDS = (
    "cohort",
    "dataset_row_id",
    "arm",
    "experiment",
    "attempt",
    "pose",
    "mapping",
    "native",
    "evaluator",
)
SUMMARY_CONTEXT_FIELDS = (
    "cohort",
    "dataset_row_id",
    "arm",
    "experiment",
    "attempt",
    "native",
    "evaluator",
)

STRUCTURAL_METRIC_FIELDS = (
    "GlobalDockQ",
    *DOCKQ_METRIC_FIELDS,
    "grouped_iRMSD",
)

DEFAULT_RANKING_KEYS = (
    "tm_score_left",
    "tm_score_right",
    "match_count_left",
    "match_count_right",
)

CONTRACT_SCORE_FIELDS = (
    *FOREIGN_KEY_FIELDS,
    "record_type",
    "interface",
    *STRUCTURAL_METRIC_FIELDS,
    "raw_json_sha256",
    "raw_json_path",
    "mapping_status",
    "chain_mapping",
)


class ContractViolation(ValueError):
    """Raised when a row would violate the explicit investigation contract."""


def _missing(value: Any) -> bool:
    return value is None or (isinstance(value, str) and not value.strip())


def require_foreign_keys(
    row: Mapping[str, Any],
    *,
    fields: Sequence[str] = FOREIGN_KEY_FIELDS,
) -> dict[str, Any]:
    """Return a copy of *row* after requiring every requested FK column.

    Foreign-key values are intentionally not coerced.  A caller may use
    strings, integers, or UUID-like objects, but no row may silently carry a
    missing provenance link.
    """

    missing = [field for field in fields if field not in row or _missing(row[field])]
    if missing:
        raise ContractViolation("missing foreign-key columns/values: " + ", ".join(missing))
    return dict(row)


def add_foreign_keys(
    rows: Iterable[Mapping[str, Any]],
    foreign_keys: Mapping[str, Any],
    *,
    fields: Sequence[str] = FOREIGN_KEY_FIELDS,
) -> list[dict[str, Any]]:
    """Attach one complete FK context to rows and validate each result."""

    context = require_foreign_keys(foreign_keys, fields=fields)
    output = []
    for row in rows:
        output.append(require_foreign_keys({**dict(row), **context}, fields=fields))
    return output


@dataclass(frozen=True)
class FrozenChainMapping:
    """Canonical chain mapping validated before any score values are read."""

    mapping_id: str
    normalized_chain_mapping: tuple[tuple[str, str], ...]
    symmetric_chain_equivalence_declared: bool

    @property
    def serialized(self) -> str:
        return ";".join(f"{model}:{native}" for model, native in self.normalized_chain_mapping)

    def as_dict(self) -> dict[str, Any]:
        return {
            "mapping_id": self.mapping_id,
            "chain_mapping": dict(self.normalized_chain_mapping),
            "symmetric_chain_equivalence_declared": self.symmetric_chain_equivalence_declared,
        }


def freeze_chain_mapping(
    mapping: Mapping[str, Any] | str,
    residue_correspondence: Sequence[Mapping[str, Any]] | None = None,
    *,
    symmetric_chain_equivalence: bool = False,
) -> FrozenChainMapping:
    """Validate and canonically identify a chain mapping before scoring.

    The mapping ID is derived only from the normalized mapping and the
    explicit symmetry declaration; DockQ values and other native-derived
    measurements are not inputs to it.
    """

    validation = validate_mapping(
        mapping,
        residue_correspondence,
        symmetric_chain_equivalence=symmetric_chain_equivalence,
    )
    if not validation.valid:
        raise ContractViolation(
            "chain mapping must be valid before score evaluation: " + "; ".join(validation.errors)
        )
    normalized = tuple(validation.normalized_chain_mapping)
    payload = {
        "chain_mapping": normalized,
        "symmetric_chain_equivalence_declared": validation.symmetric_chain_equivalence_declared,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return FrozenChainMapping(
        mapping_id=f"mapping-{digest}",
        normalized_chain_mapping=normalized,
        symmetric_chain_equivalence_declared=validation.symmetric_chain_equivalence_declared,
    )


class ContractEvaluationRecords(list):
    """List-like score records carrying the frozen mapping used to create them."""

    def __init__(self, records, *, frozen_mapping: FrozenChainMapping):
        super().__init__(records)
        self.frozen_mapping = frozen_mapping


def evaluate_dockq_with_frozen_mapping(
    source: str | Path | Mapping[str, Any],
    *,
    mapping: Mapping[str, Any] | str,
    foreign_keys: Mapping[str, Any],
    residue_correspondence: Sequence[Mapping[str, Any]] | None = None,
    symmetric_chain_equivalence: bool = False,
    grouped_irmsd: float | None = None,
) -> ContractEvaluationRecords:
    """Evaluate DockQ only after freezing a complete chain-mapping context.

    ``foreign_keys`` supplies all fields except ``mapping``; the mapping FK is
    generated from the validated canonical declaration.  This prevents score
    records from being written with an unvalidated or order-dependent mapping.
    """

    frozen = freeze_chain_mapping(
        mapping,
        residue_correspondence,
        symmetric_chain_equivalence=symmetric_chain_equivalence,
    )
    context = dict(foreign_keys)
    supplied_mapping = context.pop("mapping", None)
    if supplied_mapping is not None and supplied_mapping != frozen.mapping_id:
        raise ContractViolation("supplied mapping foreign key does not match frozen mapping")
    context["mapping"] = frozen.mapping_id
    require_foreign_keys(context)

    # Keep this call after freeze_chain_mapping: malformed score input must
    # never be the first validation step for an unfrozen mapping.
    records = evaluate_dockq_json(source, grouped_irmsd=grouped_irmsd)
    enriched = []
    for record in records:
        enriched.append(
            {
                **context,
                **dict(record),
                "mapping_status": "frozen",
                "chain_mapping": frozen.serialized,
            }
        )
    return ContractEvaluationRecords(enriched, frozen_mapping=frozen)


def _native_ranking_key(key: str) -> bool:
    normalized = key.lower()
    return (
        normalized in {field.lower() for field in STRUCTURAL_METRIC_FIELDS}
        or normalized in {"dockq", "irmsd", "native_like", "native_score", "reference_score"}
        or normalized.startswith(("native_", "reference_", "primary_dockq", "primary_irmsd"))
    )


def validate_native_independent_ranking_keys(
    keys: Iterable[str] | Mapping[str, Any],
) -> tuple[str, ...]:
    """Validate ranking keys and reject native-derived measurements explicitly."""

    raw_keys = tuple(str(key) for key in (keys.keys() if isinstance(keys, Mapping) else keys))
    native = [key for key in raw_keys if _native_ranking_key(key)]
    if native:
        raise ContractViolation("native-derived fields are not valid ranking inputs: " + ", ".join(native))
    try:
        return validate_ranking_keys(raw_keys)
    except ValueError as exc:
        raise ContractViolation(str(exc)) from exc


def _ranking_spec(
    ranking_keys: Iterable[str] | Mapping[str, str] | None,
) -> tuple[tuple[str, str], ...]:
    raw = DEFAULT_RANKING_KEYS if ranking_keys is None else ranking_keys
    if isinstance(raw, Mapping):
        spec = tuple((str(key), str(direction).lower()) for key, direction in raw.items())
    else:
        spec = tuple((str(key), "desc") for key in raw)
    if not spec:
        raise ContractViolation("at least one native-independent ranking key is required")
    validate_native_independent_ranking_keys(tuple(key for key, _ in spec))
    invalid_directions = [key for key, direction in spec if direction not in {"asc", "desc"}]
    if invalid_directions:
        raise ContractViolation("ranking directions must be 'asc' or 'desc': " + ", ".join(invalid_directions))
    return spec


def _model_valid(row: Mapping[str, Any]) -> bool:
    value = row.get("model_valid")
    if value is None:
        value = row.get("valid_model")
    if value is None:
        return row.get("status") in {"generated", "refinement_accepted", "scored", "scoreable"}
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "y", "on"}:
            return True
        if normalized in {"false", "0", "no", "n", "off", "", "none", "null"}:
            return False
    return bool(value)


def _candidate_identity(row: Mapping[str, Any]) -> tuple[str, ...]:
    return tuple(str(row.get(field, "")) for field in (
        "cohort", "dataset_row_id", "arm", "experiment", "attempt", "pose", "mapping", "native", "evaluator"
    ))


def _require_context_match(row: Mapping[str, Any], context: Mapping[str, Any]) -> None:
    mismatches = [
        field
        for field in context
        if field in context and field in row and str(row[field]) != str(context[field])
    ]
    if mismatches:
        raise ContractViolation("candidate foreign keys do not match summary context: " + ", ".join(mismatches))


def _numeric_ranking_value(row: Mapping[str, Any], key: str) -> float:
    value = row.get(key)
    if isinstance(value, bool) or value is None or value == "":
        raise ContractViolation(f"candidate {row.get('pose', '<unknown>')} has no ranking value for {key}")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ContractViolation(f"ranking value {key} must be numeric") from exc
    if not math.isfinite(number):
        raise ContractViolation(f"ranking value {key} must be finite")
    return number


def select_top_k_candidates(
    rows: Iterable[Mapping[str, Any]],
    *,
    ranking_keys: Iterable[str] | Mapping[str, str] | None = None,
    k: int = 20,
) -> list[dict[str, Any]]:
    """Select valid candidates using only pre-score, native-independent data."""

    if k <= 0:
        raise ContractViolation("k must be positive")
    spec = _ranking_spec(ranking_keys)
    candidates = []
    seen = set()
    for raw in rows:
        row = require_foreign_keys(raw)
        identity = _candidate_identity(row)
        if identity in seen:
            raise ContractViolation("duplicate candidate foreign-key identity: " + "/".join(identity))
        seen.add(identity)
        if not _model_valid(row):
            continue
        values = []
        for key, direction in spec:
            number = _numeric_ranking_value(row, key)
            values.append(number if direction == "asc" else -number)
        # The tie-break uses only provenance/candidate identity, never score
        # columns.  It is stable if input rows are reordered.
        tie_break = tuple(str(row.get(field, "")) for field in (
            "dataset_row_id", "arm", "experiment", "attempt", "pose", "mapping", "evaluator"
        ))
        candidates.append((tuple(values), tie_break, row))
    candidates.sort(key=lambda item: (item[0], item[1]))
    return [dict(item[2]) for item in candidates[:k]]


def summarize_top20(
    rows: Iterable[Mapping[str, Any]],
    *,
    foreign_keys: Mapping[str, Any],
    ranking_keys: Iterable[str] | Mapping[str, str] | None = None,
    k: int = 20,
) -> dict[str, Any]:
    """Return a contract summary with explicit no-valid-model semantics.

    Structural values are taken from the highest ``GlobalDockQ`` scored model
    among the preselected candidates.  If no model is valid, every structural
    value is ``None``.  The sole unconditional metric is
    ``best_GlobalDockQ_at_20``, which is exactly ``0.0`` when no scored model
    is available.
    """

    # A pair summary is shared by many poses and mappings.  Require the
    # common pair/evaluator context here, while every candidate row still
    # carries the complete pose-level foreign-key contract.
    require_foreign_keys(foreign_keys, fields=SUMMARY_CONTEXT_FIELDS)
    context = {field: foreign_keys[field] for field in SUMMARY_CONTEXT_FIELDS}
    materialized = [require_foreign_keys(row) for row in rows]
    for row in materialized:
        _require_context_match(row, context)
    selected = select_top_k_candidates(materialized, ranking_keys=ranking_keys, k=k)
    valid_rows = [row for row in materialized if _model_valid(row)]
    scored = []
    for row in selected:
        value = row.get("GlobalDockQ")
        if value is None or value == "":
            continue
        try:
            numeric = float(value)
        except (TypeError, ValueError) as exc:
            raise ContractViolation("GlobalDockQ must be numeric when present") from exc
        if not math.isfinite(numeric):
            raise ContractViolation("GlobalDockQ must be finite when present")
        scored.append((numeric, row))
    best_scored = max(scored, key=lambda item: (item[0], _candidate_identity(item[1]))) if scored else None

    summary = {
        **context,
        "selected_count": len(selected),
        "model_valid_count": len(valid_rows),
        "scoreable_model_count": len(scored),
        "best_GlobalDockQ_at_20": 0.0 if best_scored is None else best_scored[0],
    }
    if not valid_rows or best_scored is None:
        summary.update({field: None for field in STRUCTURAL_METRIC_FIELDS})
    else:
        best_row = best_scored[1]
        summary.update({field: best_row.get(field) for field in STRUCTURAL_METRIC_FIELDS})
    return summary


def write_contract_scores_tsv(
    source: str | Path | Mapping[str, Any],
    global_path: str | Path,
    interfaces_path: str | Path,
    *,
    mapping: Mapping[str, Any] | str,
    foreign_keys: Mapping[str, Any],
    raw_json_path: str | Path | None = None,
    residue_correspondence: Sequence[Mapping[str, Any]] | None = None,
    symmetric_chain_equivalence: bool = False,
    grouped_irmsd: float | None = None,
) -> tuple[Path, Path]:
    """Write global/interface score rows with all contract FKs populated."""

    records = evaluate_dockq_with_frozen_mapping(
        source,
        mapping=mapping,
        foreign_keys=foreign_keys,
        residue_correspondence=residue_correspondence,
        symmetric_chain_equivalence=symmetric_chain_equivalence,
        grouped_irmsd=grouped_irmsd,
    )
    rows = []
    source_path = Path(source) if isinstance(source, (str, Path)) else None
    if raw_json_path is None and source_path is not None:
        raw_json_path = source_path
    if raw_json_path is not None:
        raw_path = Path(raw_json_path)
        if not raw_path.is_file():
            raise ContractViolation(f"raw DockQ JSON path does not exist: {raw_path}")
        raw_hash = hashlib.sha256(raw_path.read_bytes()).hexdigest()
        if raw_hash != records[0]["raw_json_sha256"]:
            raise ContractViolation("raw DockQ JSON hash does not match evaluator record")
    for record in records:
        rows.append({
            **dict(record),
            "raw_json_path": str(raw_json_path or ""),
        })

    def write(path: str | Path, selected: list[dict[str, Any]]) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=CONTRACT_SCORE_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in sorted(selected, key=lambda item: tuple(str(item.get(field, "")) for field in CONTRACT_SCORE_FIELDS)):
                writer.writerow({field: row.get(field, "") for field in CONTRACT_SCORE_FIELDS})
        return output

    return (
        write(global_path, [row for row in rows if row["record_type"] == "global"]),
        write(interfaces_path, [row for row in rows if row["record_type"] == "interface"]),
    )
