"""Immutable lineage and artifact helpers for investigation runs.

The module deliberately has no dependency on the production PRISM pipeline.
Lineage rows are append-only JSONL events, while PDB artifacts are copied by
creating a destination that must not already contain different bytes.  Native
scores are useful for evaluation summaries, but are explicitly excluded from
the ranking-input allowlist.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


LINEAGE_RECORD_FIELDS = (
    "schema_version",
    "record_id",
    "stage",
    "pair_id",
    "template_id",
    "side",
    "alignment_id",
    "orientation",
    "filter_stage",
    "pose_id",
    "refinement_id",
    "energy_id",
    "scoring_id",
    "artifact_path",
    "coordinate_sha256",
    "status",
    "failure_reason",
)

LINEAGE_STAGES = frozenset(
    {
        "pair",
        "template",
        "side",
        "alignment",
        "orientation",
        "filter",
        "pose",
        "refinement",
        "energy",
        "scoring",
    }
)

NON_TERMINAL_STATUSES = frozenset(
    {"pending", "started", "running", "in_progress", "alignment_pending"}
)
TERMINAL_STATUSES = frozenset(
    {
        "generated",
        "alignment_failed",
        "filter_rejected",
        "pose_created",
        "pose_written",
        "refinement_failed",
        "refinement_accepted",
        "energy_recorded",
        "energy_failed",
        "scored",
        "scoring_failed",
        "not_scoreable",
        "completed",
        "failed",
        "rejected",
        "skipped",
    }
)
_FAILURE_STATUSES = frozenset(
    {
        "alignment_failed",
        "filter_rejected",
        "refinement_failed",
        "energy_failed",
        "scoring_failed",
        "not_scoreable",
        "failed",
        "rejected",
        "skipped",
    }
)

# These are candidate-derived inputs only.  DockQ, iRMSD, native labels,
# native identifiers, and reference/native artifact fields are intentionally
# absent even though they may appear in evaluation output rows.
ALLOWED_RANKING_KEYS = frozenset(
    {
        "tm_score_left",
        "tm_score_right",
        "match_count_left",
        "match_count_right",
        "match_coverage_left",
        "match_coverage_right",
        "contact_count",
        "clash_count",
        "rosetta_interaction_score",
        "energy",
        "orientation",
        "filter_passed",
    }
)
_NATIVE_DERIVED_KEYS = frozenset(
    {
        "dockq",
        "irmsd",
        "native_like",
        "native_complex_id",
        "native_pdb",
        "native_receptor",
        "native_ligand",
        "native_interface",
        "reference_pdb",
    }
)


def _is_failure_status(status: str) -> bool:
    return status in _FAILURE_STATUSES or status.endswith(("_failed", "_rejected"))


def is_terminal_status(status: str) -> bool:
    """Return whether *status* closes one lineage event key."""
    return status in TERMINAL_STATUSES or _is_failure_status(status)


def _validate_status(status: str) -> None:
    if not isinstance(status, str) or not status:
        raise ValueError("status must be a non-empty string")
    if status not in TERMINAL_STATUSES and status not in NON_TERMINAL_STATUSES:
        raise ValueError(f"unknown lineage status: {status}")


def _stable_record_id(values: Mapping[str, Any]) -> str:
    payload = json.dumps(values, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True, slots=True)
class LineageRecord:
    """One immutable, fixed-schema event in an investigation lineage."""

    schema_version: int = 1
    record_id: str = ""
    stage: str = "pair"
    pair_id: str = ""
    template_id: str | None = None
    side: str | None = None
    alignment_id: str | None = None
    orientation: str | None = None
    filter_stage: str | None = None
    pose_id: str | None = None
    refinement_id: str | None = None
    energy_id: str | None = None
    scoring_id: str | None = None
    artifact_path: str | None = None
    coordinate_sha256: str | None = None
    status: str = "pending"
    failure_reason: str | None = None

    def __post_init__(self) -> None:
        if self.schema_version != 1:
            raise ValueError("unsupported lineage schema version")
        if not self.pair_id:
            raise ValueError("pair_id must not be empty")
        if self.stage not in LINEAGE_STAGES:
            raise ValueError(f"unknown lineage stage: {self.stage}")
        _validate_status(self.status)
        if _is_failure_status(self.status) and not self.failure_reason:
            raise ValueError(f"failure_reason is required for status {self.status}")
        if self.coordinate_sha256 is not None:
            digest = self.coordinate_sha256.lower()
            if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                raise ValueError("coordinate_sha256 must be a SHA-256 hex digest")
        if not self.record_id:
            values = {
                field: getattr(self, field)
                for field in LINEAGE_RECORD_FIELDS
                if field != "record_id"
            }
            object.__setattr__(self, "record_id", _stable_record_id(values))

    @property
    def is_terminal(self) -> bool:
        return is_terminal_status(self.status)

    def event_key(self) -> tuple[Any, ...]:
        """Key whose terminal event cannot be replaced or extended."""
        return tuple(getattr(self, field) for field in LINEAGE_RECORD_FIELDS if field not in {
            "record_id",
            "status",
            "failure_reason",
        })

    def to_dict(self) -> dict[str, Any]:
        return {field: getattr(self, field) for field in LINEAGE_RECORD_FIELDS}

    @classmethod
    def from_dict(cls, values: Mapping[str, Any]) -> "LineageRecord":
        missing = [field for field in LINEAGE_RECORD_FIELDS if field not in values]
        extra = [field for field in values if field not in LINEAGE_RECORD_FIELDS]
        if missing or extra:
            details = []
            if missing:
                details.append("missing " + ", ".join(missing))
            if extra:
                details.append("unexpected " + ", ".join(extra))
            raise ValueError("invalid lineage record schema: " + "; ".join(details))
        return cls(**{field: values[field] for field in LINEAGE_RECORD_FIELDS})


class AppendOnlyLineage:
    """Append immutable lineage events to a JSONL file.

    Re-appending the exact same record is an idempotent no-op.  A different
    record with the same id, or any event after a terminal event for the same
    event key, raises rather than rewriting history.
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._records = self._read_from_disk()
        self._record_ids = {record.record_id for record in self._records}
        self._terminal_event_keys = {record.event_key() for record in self._records if record.is_terminal}

    def _read_from_disk(self) -> tuple[LineageRecord, ...]:
        if not self.path.exists():
            return ()
        records = []
        for line_number, line in enumerate(self.path.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip():
                continue
            try:
                records.append(LineageRecord.from_dict(json.loads(line)))
            except (TypeError, ValueError, json.JSONDecodeError) as exc:
                raise ValueError(f"invalid lineage record at line {line_number}") from exc
        return tuple(records)

    def read(self) -> tuple[LineageRecord, ...]:
        return self._records

    def append(self, record: LineageRecord) -> bool:
        if not isinstance(record, LineageRecord):
            raise TypeError("append expects a LineageRecord")
        if record.record_id in self._record_ids:
            prior = next(item for item in self._records if item.record_id == record.record_id)
            if prior == record:
                return False
            raise ValueError(f"record id already contains different bytes: {record.record_id}")
        if record.event_key() in self._terminal_event_keys:
            raise ValueError("terminal lineage event cannot be extended or replaced")

        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record.to_dict(), sort_keys=True, separators=(",", ":")) + "\n")
        self._records = (*self._records, record)
        self._record_ids.add(record.record_id)
        if record.is_terminal:
            self._terminal_event_keys.add(record.event_key())
        return True

    def assert_complete(self) -> None:
        """Fail if the run still contains a non-terminal lineage event."""
        pending = [record.record_id for record in self._records if not record.is_terminal]
        if pending:
            raise ValueError("lineage contains non-terminal records: " + ", ".join(pending))


def append_lineage_record(path: str | Path, record: LineageRecord) -> bool:
    """Convenience wrapper for appending one record."""
    return AppendOnlyLineage(path).append(record)


def coordinate_sha256(path: str | Path) -> str:
    """Return the SHA-256 digest of the exact coordinate file bytes."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def create_pose_record(
    pose_path: str | Path,
    *,
    pair_id: str,
    template_id: str | None = None,
    side: str | None = None,
    pose_id: str | None = None,
    record_id: str = "",
    alignment_id: str | None = None,
    orientation: str | None = None,
    filter_stage: str | None = None,
    refinement_id: str | None = None,
    energy_id: str | None = None,
    scoring_id: str | None = None,
    status: str = "pose_created",
    failure_reason: str | None = None,
) -> LineageRecord:
    """Create a pose-stage record whose digest covers the exact PDB bytes."""
    path = Path(pose_path)
    return LineageRecord(
        record_id=record_id,
        stage="pose",
        pair_id=pair_id,
        template_id=template_id,
        side=side,
        alignment_id=alignment_id,
        orientation=orientation,
        filter_stage=filter_stage,
        pose_id=pose_id or path.stem,
        refinement_id=refinement_id,
        energy_id=energy_id,
        scoring_id=scoring_id,
        artifact_path=str(path),
        coordinate_sha256=coordinate_sha256(path),
        status=status,
        failure_reason=failure_reason,
    )


def copy_pdb_immutable(source: str | Path, destination: str | Path) -> Path:
    """Copy a PDB only when the destination is absent or byte-identical."""
    source_path = Path(source)
    destination_path = Path(destination)
    source_bytes = source_path.read_bytes()
    destination_path.parent.mkdir(parents=True, exist_ok=True)

    if destination_path.is_symlink():
        raise FileExistsError(f"refusing to overwrite symlink: {destination_path}")
    if destination_path.exists():
        if destination_path.read_bytes() == source_bytes:
            return destination_path
        raise FileExistsError(f"refusing overwrite: bytes differ at {destination_path}")

    try:
        descriptor = os.open(
            destination_path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL,
            0o644,
        )
    except FileExistsError:
        # A concurrent creator won the race.  It is still safe to accept only
        # an exact byte match; never replace what it wrote.
        if destination_path.is_file() and destination_path.read_bytes() == source_bytes:
            return destination_path
        raise FileExistsError(f"refusing overwrite: bytes differ at {destination_path}")

    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(source_bytes)
            handle.flush()
            os.fsync(handle.fileno())
    except Exception:
        try:
            destination_path.unlink()
        except FileNotFoundError:
            pass
        raise
    return destination_path


def validate_ranking_keys(keys: Iterable[str] | Mapping[str, Any]) -> tuple[str, ...]:
    """Validate that ranking inputs are allowed candidate-derived fields."""
    raw_keys = keys.keys() if isinstance(keys, Mapping) else keys
    normalized = tuple(str(key) for key in raw_keys)
    native = [
        key
        for key in normalized
        if key in _NATIVE_DERIVED_KEYS
        or key.startswith(("native_", "reference_", "primary_dockq", "primary_irmsd"))
        or key in {"dockq", "irmsd"}
    ]
    if native:
        raise ValueError("native-derived fields are not valid ranking inputs: " + ", ".join(native))
    unknown = [key for key in normalized if key not in ALLOWED_RANKING_KEYS]
    if unknown:
        raise ValueError("unsupported ranking inputs: " + ", ".join(unknown))
    return normalized


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "y", "on"}:
            return True
        if normalized in {"false", "0", "no", "n", "off", "", "none", "null"}:
            return False
    return bool(value)


def _row_value(row: Mapping[str, Any] | LineageRecord, key: str, default: Any = None) -> Any:
    if isinstance(row, LineageRecord):
        return getattr(row, key, default)
    return row.get(key, default)


def _float_values(rows: Sequence[Mapping[str, Any] | LineageRecord], key: str) -> list[float]:
    values = []
    for row in rows:
        value = _row_value(row, key)
        if value is not None and value != "":
            values.append(float(value))
    return values


def _scoreable(row: Mapping[str, Any] | LineageRecord) -> bool:
    explicit = _row_value(row, "scoreable")
    if explicit is not None:
        return _as_bool(explicit)
    if _row_value(row, "status") in {"scored", "scoreable"}:
        return True
    return any(_row_value(row, key) is not None for key in ("dockq", "irmsd"))


def _primary_sort_key(row: Mapping[str, Any] | LineageRecord, ranking_keys: tuple[str, ...]) -> tuple[Any, ...]:
    model_id = str(_row_value(row, "model_id", _row_value(row, "pose_id", "")))
    if not ranking_keys:
        return (model_id,)
    values = []
    for key in ranking_keys:
        value = _row_value(row, key)
        values.append(-float(value) if isinstance(value, (int, float)) else str(value or ""))
    return (*values, model_id)


def _mean_or_zero(values: Sequence[float]) -> float:
    return float(sum(values) / len(values)) if values else 0.0


def aggregate_pair_summary(
    rows: Iterable[Mapping[str, Any] | LineageRecord],
    *,
    pair_key: str = "pair_id",
    ranking_keys: Iterable[str] | Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Aggregate model rows into deterministic pair-level summaries.

    Primary fields describe the best scoreable model and macro statistics for
    the pair.  Secondary counts describe how many model rows were scoreable or
    not scoreable.  When none are scoreable, all numeric primary metrics are
    explicitly ``0.0``.
    """
    validated_ranking_keys = validate_ranking_keys(ranking_keys) if ranking_keys is not None else ()
    groups: dict[str, list[Mapping[str, Any] | LineageRecord]] = {}
    for row in rows:
        pair_id = _row_value(row, pair_key)
        if pair_id in (None, ""):
            raise ValueError(f"row missing {pair_key}")
        groups.setdefault(str(pair_id), []).append(row)

    summaries: list[dict[str, Any]] = []
    for pair_id in sorted(groups):
        model_rows = groups[pair_id]
        scoreable_rows = [row for row in model_rows if _scoreable(row)]
        primary = sorted(scoreable_rows, key=lambda row: _primary_sort_key(row, validated_ranking_keys))[0] if scoreable_rows else None
        dockqs = _float_values(scoreable_rows, "dockq")
        irmsds = _float_values(scoreable_rows, "irmsd")
        energies = _float_values(scoreable_rows, "energy")
        status_counts: dict[str, int] = {}
        for row in model_rows:
            status = str(_row_value(row, "status", "unknown"))
            status_counts[status] = status_counts.get(status, 0) + 1

        def primary_value(key: str) -> float:
            value = _row_value(primary, key) if primary is not None else None
            return float(value) if value is not None and value != "" else 0.0

        summaries.append(
            {
                "pair_id": pair_id,
                "status": "scoreable" if scoreable_rows else "not_scoreable",
                "failure_reason": None
                if scoreable_rows
                else next((_row_value(row, "failure_reason") for row in model_rows if _row_value(row, "failure_reason")), "no_scoreable_model"),
                "primary_model_id": None
                if primary is None
                else str(_row_value(primary, "model_id", _row_value(primary, "pose_id", ""))),
                "primary_dockq": primary_value("dockq"),
                "primary_irmsd": primary_value("irmsd"),
                "primary_energy": primary_value("energy"),
                "dockq_mean": _mean_or_zero(dockqs),
                "dockq_best": max(dockqs, default=0.0),
                "irmsd_mean": _mean_or_zero(irmsds),
                "irmsd_best": min(irmsds, default=0.0),
                "energy_mean": _mean_or_zero(energies),
                "energy_best": min(energies, default=0.0),
                "model_count": len(model_rows),
                "scoreable_model_count": len(scoreable_rows),
                "secondary_model_count": max(0, len(scoreable_rows) - 1),
                "non_scoreable_model_count": len(model_rows) - len(scoreable_rows),
                "status_counts": status_counts,
            }
        )
    return summaries


# Descriptive aliases for callers that use artifact-oriented terminology.
immutable_copy = copy_pdb_immutable
PoseRecord = LineageRecord
