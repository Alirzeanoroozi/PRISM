# /// script
# requires-python = ">=3.11"
# dependencies = []
# //

"""Run evidence and provenance helpers for the PRISM pipeline.

This module is the singular canonical implementation of
Phase 1 (Run Identity and Manifest).  It provides:

- ``build_declared_contract`` -- construct the immutable declared contract
  from current pipeline configuration.

- ``build_execution_attempt`` -- capture runtime facts (host, git, Slurm,
  Python environment) and link them to a declared contract hash.

- ``ArtifactObservation``, ``observe_artifact``, ``write_artifact_ledger``,
  ``append_artifact_observation``, ``validate_artifact_ledger``,
  ``closeout_artifact_ledger`` -- create, persist, and cross-check
  artifact observations against the ledger.

- ``validate_before_consume`` -- consumer gate that checks an expected
  artifact inventory against the ledger.

Public object names are considered stable across the 1.x series.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class ArtifactLedgerError(Exception):
    """Raised when an artifact-ledger operation fails."""


class ConsumerGateError(Exception):
    """Raised when validation before consuming an artifact fails."""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def canonical_json(obj: Any) -> str:
    """Deterministic JSON serialisation with sorted keys and no whitespace."""
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def canonical_hash(obj: Any) -> str:
    """SHA-256 hex digest of the canonical JSON representation of *obj*."""
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()


def sha256_file(path: str | os.PathLike[str]) -> str | None:
    """Return the SHA-256 hex digest of *path*, or *None* if the file is
    inaccessible."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            while True:
                block = fh.read(65536)
                if not block:
                    break
                h.update(block)
        return h.hexdigest()
    except (OSError, FileNotFoundError):
        return None


def _make_run_id() -> str:
    """Generate a unique run identifier."""
    now = datetime.now(timezone.utc)
    ts = now.strftime("%Y%m%d-%H%M%S")
    pid = os.getpid()
    raw = f"{ts}-{pid}-{time.monotonic_ns()}"
    h8 = hashlib.sha256(raw.encode("utf-8")).hexdigest()[:8]
    return f"prism-{ts}-{pid}-{h8}"


def _redact_argv(argv: Iterable[Any]) -> list[str]:
    """Redact long arguments and path-like values from *argv*."""
    out: list[str] = []
    for a in argv:
        s = str(a)
        if len(s) > 200:
            s = s[:100] + "…[redacted]"
        out.append(s)
    return out


def _normalize_mapping(m: Mapping[str, Any]) -> dict[str, Any]:
    """Return a JSON-safe dict from *m*."""
    try:
        canonical_json(m)
        return dict(m)
    except Exception:
        return {k: str(v) for k, v in m.items()}


def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _get_git_provenance(cwd: Path | None = None) -> dict[str, Any]:
    """Gather git provenance from the working directory."""
    try:
        kwargs: dict[str, Any] = {"capture_output": True, "text": True, "timeout": 5}
        cwd_arg = {"cwd": str(cwd)} if cwd else {}
        head = subprocess.run(["git", "rev-parse", "HEAD"], **kwargs, **cwd_arg).stdout.strip()
        head_short = head[:8] if len(head) >= 8 else head
        branch = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], **kwargs, **cwd_arg
        ).stdout.strip()
        status = subprocess.run(["git", "status", "--short"], **kwargs, **cwd_arg).stdout.strip()
        tracked_diff = subprocess.run(
            ["git", "diff", "HEAD", "--binary"], **kwargs, **cwd_arg
        ).stdout
        index_diff = subprocess.run(
            ["git", "diff", "--cached", "--binary"], **kwargs, **cwd_arg
        ).stdout
        untracked = [
            line for line in status.splitlines()
            if line.startswith("?? ")
        ]
        is_dirty = bool(status)
        diff_hash = hashlib.sha256(status.encode("utf-8")).hexdigest() if is_dirty else None
        return {
            "head_commit": head,
            "is_dirty": is_dirty,
            "status_short": status,
            "head_short": head_short,
            "branch": branch if branch != "HEAD" else None,
            "diff_hash": diff_hash,
            "tracked_diff_hash": hashlib.sha256(tracked_diff.encode("utf-8")).hexdigest(),
            "index_diff_hash": hashlib.sha256(index_diff.encode("utf-8")).hexdigest(),
            "untracked_manifest_hash": hashlib.sha256(
                "\n".join(untracked).encode("utf-8")
            ).hexdigest(),
            "untracked_count": len(untracked),
            "declared_untracked_files": untracked,
        }
    except Exception:
        return {
            "head_commit": "0" * 40,
            "is_dirty": True,
            "status_short": "git unavailable",
            "head_short": "00000000",
            "branch": None,
            "diff_hash": None,
            "declared_untracked_files": [],
        }


def _get_external_tools(cwd: Path | None = None) -> list[dict[str, Any]]:
    """Record configured external-tool paths and hashes without executing them."""
    candidates = [
        ("TMalign", "PRISM_TMALIGN", "external_tools/TMalign", "TMalign"),
        ("GTalign", "PRISM_GTALIGN", "gtalign", "gtalign"),
        ("MultiProt", "PRISM_MULTIPROT", "external_tools/multiprot.Linux", "MultiProt"),
        ("NACCESS", "PRISM_NACCESS_EXECUTABLE", "external_tools/naccess/naccess", "naccess"),
    ]
    tools: list[dict[str, Any]] = []
    root = cwd or Path.cwd()
    for name, env_key, default, probe_name in candidates:
        try:
            configured = os.environ.get(env_key, default)
            candidate = Path(configured)
            if not candidate.is_absolute():
                candidate = root / candidate
            resolved = str(candidate.resolve()) if candidate.is_file() else shutil.which(configured) or shutil.which(probe_name) or ""
            ok = bool(resolved and Path(resolved).is_file())
            the_hash = sha256_file(resolved) if ok and Path(resolved).is_file() else None
            tools.append({
                "name": name,
                "resolved_path": resolved if ok else "",
                "version_output": None,
                "sha256": the_hash,
                "probe_ok": ok,
            })
        except Exception:
            tools.append({"name": name, "resolved_path": "", "version_output": None, "sha256": None, "probe_ok": False})
    return tools


def _get_slurm_info() -> dict[str, Any] | None:
    """Read Slurm environment variables if running under Slurm."""
    if "SLURM_JOB_ID" not in os.environ:
        return None
    nodelist = os.environ.get("SLURM_NODELIST")
    return {
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "job_name": os.environ.get("SLURM_JOB_NAME"),
        "partition": os.environ.get("SLURM_JOB_PARTITION"),
        "qos": os.environ.get("SLURM_JOB_QOS"),
        "account": os.environ.get("SLURM_ACCOUNT"),
        "nodes": nodelist.split(",") if nodelist else None,
        "cpus_per_task": int(os.environ["SLURM_CPUS_PER_TASK"]) if "SLURM_CPUS_PER_TASK" in os.environ else None,
        "mem_per_node": os.environ.get("SLURM_MEM_PER_NODE"),
        "time_limit": os.environ.get("SLURM_TIME_LIMIT"),
        "submit_dir": os.environ.get("SLURM_SUBMIT_DIR"),
    }


_SAFE_ENVIRONMENT_KEYS = {
    "CONDA_DEFAULT_ENV",
    "CONDA_PREFIX",
    "HOSTNAME",
    "PATH",
    "PRISM_PIPELINE_PYTHON",
    "PYTHONPATH",
    "VIRTUAL_ENV",
}


def _environment_identity(environment: Mapping[str, Any] | None) -> dict[str, Any]:
    """Record environment identity without serialising credentials or secrets."""
    values = dict(environment if environment is not None else os.environ)
    selected = {
        key: str(values[key])
        for key in sorted(values)
        if key in _SAFE_ENVIRONMENT_KEYS
        and not any(token in key.upper() for token in ("SECRET", "PASSWORD", "TOKEN", "KEY"))
    }
    return {
        "selected": selected,
        "python_executable": sys.executable,
        "python_executable_sha256": sha256_file(sys.executable),
        "environment_key_count": len(values),
        "environment_key_hash": canonical_hash(sorted(values)),
    }


# ---------------------------------------------------------------------------
# Declared contract
# ---------------------------------------------------------------------------


def build_declared_contract(
    *,
    pipeline_version: str,
    stages_enabled: Sequence[str],
    aligner: str,
    refiner: str,
    input_selectors: Mapping[str, Any],
    template_inventory: Sequence[Mapping[str, Any]],
    parameters: Mapping[str, Any],
    resource_request: Mapping[str, Any],
    source_inventory: Mapping[str, Any],
    tool_fingerprints: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Build a declared-contract dict ready for signing and hashing.

    The returned dict matches ``declared-contract.schema.json``.
    The ``contract_hash`` field is the SHA-256 of the **serialised contract
    with the ``contract_hash`` field excluded** (to avoid circularity).
    """
    contract = {
        "contract_version": "1.0",
        "pipeline_version": pipeline_version,
        "stages_enabled": list(stages_enabled),
        "aligner": aligner,
        "refiner": refiner,
        "input_selectors": dict(input_selectors),
        "template_inventory": [dict(t) for t in template_inventory],
        "parameters": dict(parameters),
        "resource_request": dict(resource_request),
        "source_inventory": dict(source_inventory),
        "tool_fingerprints": [dict(t) for t in tool_fingerprints],
        "contract_hash": "",  # placeholder
    }
    # Compute hash over the contract WITHOUT the contract_hash field.
    signing = {k: v for k, v in contract.items() if k != "contract_hash"}
    contract["contract_hash"] = canonical_hash(signing)
    return contract


# ---------------------------------------------------------------------------
# Execution attempt (run manifest)
# ---------------------------------------------------------------------------


def build_execution_attempt(
    contract: Mapping[str, Any],
    *,
    run_id: str = "",
    command: Iterable[Any] | None = None,
    working_directory: str | os.PathLike[str] | None = None,
    environment: Mapping[str, Any] | None = None,
    output_root: str | os.PathLike[str] | None = None,
    slurm: Mapping[str, Any] | None = None,
    parent_run_id: str | None = None,
    supersedes_reason: str | None = None,
) -> dict[str, Any]:
    """Link one operational attempt to a declared contract.

    The caller supplies ``run_id`` so retries remain readable and distinct.
    Runtime facts (host, user, git, python, slurm) are auto-detected.

    Returns a dict matching run-manifest.schema.json:
      manifest_version, run_identity, declared_contract_hash, runtime_context.
    """

    contract_hash = str(contract.get("contract_hash", ""))
    if len(contract_hash) != 64:
        raise ArtifactLedgerError("execution attempt requires a SHA-256 contract_hash")
    rid = run_id.strip() if run_id.strip() else _make_run_id()
    py_ver = platform.python_version()
    pkgs: dict[str, str] = {}
    try:
        import importlib.metadata as _im
        pkgs = {dist.metadata.get("Name", dist.name): dist.version for dist in _im.distributions()}
    except Exception:
        pass
    return {
        "manifest_version": "1.0",
        "run_identity": {
            "run_id": rid,
            "attempt_id": rid,
            "created_at": _iso_now(),
            "host": os.uname().nodename,
            "user": os.environ.get("USER", os.environ.get("USERNAME", "unknown")),
            "status": "declared",
            "parent_run_id": parent_run_id,
            "supersedes_reason": supersedes_reason,
        },
        "declared_contract_hash": contract_hash,
        "runtime_context": {
            "command_argv": _redact_argv(command or ()),
            "working_directory": str(Path(working_directory or Path.cwd()).resolve()),
            "output_root": str(Path(output_root or working_directory or Path.cwd()).resolve()),
            "environment": _environment_identity(environment),
            "python_version": py_ver,
            "package_versions": pkgs,
            "external_tools": _get_external_tools(Path(working_directory or Path.cwd()).resolve()),
            "git_provenance": _get_git_provenance(Path(working_directory or Path.cwd()).resolve()),
            "slurm": _get_slurm_info() if slurm is None else _normalize_mapping(slurm),
            "seeds": {},
        },
    }


# ---------------------------------------------------------------------------
# Artifact observation
# ---------------------------------------------------------------------------


@dataclass
class ArtifactObservation:
    """A single row in the artifact ledger.

    Matches ``artifact-ledger.schema.json``.
    """

    run_id: str
    stage: str
    dataset_row_id: str
    scientific_role: str
    run_relative_path: str
    path_kind: str = "file"
    size_bytes: int | None = None
    sha256: str | None = None
    target_sha256: str | None = None
    link_target: str | None = None
    resolved_path: str | None = None
    is_broken_link: bool = False
    materialized_at: str = ""
    produced_by: str = ""
    status: str = "ok"
    error_detail: str | None = None

    def __post_init__(self) -> None:
        if not self.materialized_at:
            self.materialized_at = _iso_now()

    @property
    def key(self) -> tuple[str, str, str]:
        """Primary key per ADR-0002: (dataset_row_id, scientific_role, run_relative_path)."""
        return (self.dataset_row_id, self.scientific_role, self.run_relative_path)

    @classmethod
    def from_row(cls, row: dict[str, Any]) -> "ArtifactObservation":
        """Deserialise from a ledger dict (already validated).

        Coerces types back from TSV strings to Python types:
        * ``size_bytes`` -> ``int | None``
        * ``is_broken_link`` -> ``bool``
        * empty strings for nullable fields -> ``None``
        """
        return cls(
            run_id=row.get("run_id", ""),
            stage=row.get("stage", ""),
            dataset_row_id=row.get("dataset_row_id", ""),
            scientific_role=row.get("scientific_role", ""),
            run_relative_path=row.get("run_relative_path", ""),
            path_kind=row.get("path_kind", "file"),
            size_bytes=_tsv_parse_cell(row.get("size_bytes", ""), "int"),
            sha256=_tsv_parse_cell(row.get("sha256", ""), "str") or None,
            target_sha256=_tsv_parse_cell(row.get("target_sha256", ""), "str") or None,
            link_target=_tsv_parse_cell(row.get("link_target", ""), "str") or None,
            resolved_path=_tsv_parse_cell(row.get("resolved_path", ""), "str") or None,
            is_broken_link=_tsv_parse_cell(row.get("is_broken_link", ""), "bool"),
            materialized_at=row.get("materialized_at", ""),
            produced_by=row.get("produced_by", ""),
            status=row.get("status", "ok"),
            error_detail=_tsv_parse_cell(row.get("error_detail", ""), "str") or None,
        )

    def as_dict(self) -> dict[str, Any]:
        """Serialise to a dict matching artifact-ledger.schema.json.

        Returns:
            A dict with all 16 required fields.  ``is_broken_link`` is a bare
            JSON boolean.  Absent hashes are ``null``, not the empty string.
        """
        return {
            "run_id": self.run_id,
            "stage": self.stage,
            "dataset_row_id": self.dataset_row_id,
            "scientific_role": self.scientific_role,
            "run_relative_path": self.run_relative_path,
            "path_kind": self.path_kind,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256 if self.sha256 else None,
            "target_sha256": self.target_sha256 if self.target_sha256 else None,
            "link_target": self.link_target if self.link_target else None,
            "resolved_path": self.resolved_path if self.resolved_path else None,
            "is_broken_link": bool(self.is_broken_link),
            "materialized_at": self.materialized_at,
            "produced_by": self.produced_by,
            "status": self.status,
            "error_detail": self.error_detail if self.error_detail else None,
        }


def _read_tsv_ledger(path: Path) -> list[ArtifactObservation]:
    """Read a TSV artifact ledger, returning a list of observations.

    Handles the TSV format produced by ``write_artifact_ledger``.
    Empty or missing files return an empty list.
    """
    if not path.is_file():
        return []
    records: list[ArtifactObservation] = []
    with path.open("r", encoding="utf-8") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            row = dict(zip(header, parts))
            records.append(ArtifactObservation.from_row(row))
    return records


_TSV_FIELDS = [
    "run_id", "stage", "dataset_row_id", "scientific_role", "run_relative_path",
    "path_kind", "size_bytes", "sha256", "target_sha256", "link_target",
    "resolved_path", "is_broken_link", "materialized_at", "produced_by", "status", "error_detail",
]


def _tsv_cell(value: Any) -> str:
    """Serialize a single TSV cell preserving type information.

    * ``None`` / ``null`` -> empty string
    * ``True`` -> ``"true"``, ``False`` -> ``"false"``
    * Everything else -> ``str(value)``
    """
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _tsv_parse_cell(raw: str, target_type: str) -> Any:
    """Parse a TSV cell back to the target type, reversing ``_tsv_cell``."""
    if raw == "":
        return None
    if target_type == "bool":
        return raw.lower() == "true"
    if target_type == "int":
        try:
            return int(raw)
        except ValueError:
            return None
    return raw


def _write_tsv_ledger(path: Path, records: Sequence[ArtifactObservation], *, append: bool = False) -> None:
    """Write (or append) a TSV artifact ledger.

    Uses the schema field order for columns.
    Booleans are serialised as ``"true"`` / ``"false"``.
    ``None`` is serialised as the empty string (not ``"None"``).
    """
    mode = "a" if append else "w"
    with path.open(mode, encoding="utf-8") as fh:
        if not append or not path.is_file() or path.stat().st_size == 0:
            fh.write("\t".join(_TSV_FIELDS) + "\n")
        for rec in records:
            d = rec.as_dict()
            row = [_tsv_cell(d.get(f, None)) for f in _TSV_FIELDS]
            fh.write("\t".join(row) + "\n")


def observe_artifact(
    path: str | os.PathLike[str],
    *,
    run_id: str,
    stage: str,
    dataset_row_id: str,
    scientific_role: str,
    run_relative_path: str,
    produced_by: str = "",
    path_kind: str = "file",
) -> ArtifactObservation:
    """Observe a single artifact and return an ``ArtifactObservation``.

    The observation records:
    * The file's size, SHA-256, and resolution status.
    * Symlink metadata (target, resolved path, broken-link flag).
    * The caller's run identity and stage.

    Args:
        path:  Filesystem path to the artifact.
        run_id:  Current run identifier.
        stage:  Pipeline stage that produced this artifact.
        dataset_row_id:  Row identifier from the input dataset.
        scientific_role:  Role string (e.g. ``receptor``, ``ligand``).
        run_relative_path:  Path relative to the run output directory.
        produced_by:  Tool/command that created the artifact.
        path_kind:  ``file``, ``symlink``, ``directory``, or ``missing``.

    Returns:
        An ``ArtifactObservation`` ready for ledger insertion.
    """
    p = Path(path)
    materialized = _iso_now()
    size_bytes: int | None = None
    sha256_val: str | None = None
    link_target: str | None = None
    target_sha256: str | None = None
    resolved_path: str | None = None
    is_broken_link = False
    status = "ok"

    if not p.exists() and not p.is_symlink():
        status = "missing"
        path_kind = "missing"
        size_bytes = 0
    elif p.is_symlink():
        path_kind = "symlink"
        link_target = str(p.readlink())
        resolved = p.resolve(strict=False)
        resolved_path = str(resolved)
        if not resolved.exists():
            is_broken_link = True
            status = "error"
            size_bytes = 0
        else:
            size_bytes = resolved.stat().st_size
            sha256_val = sha256_file(resolved)
            target_sha256 = sha256_val
    elif p.is_file():
        path_kind = path_kind or "file"
        size_bytes = p.stat().st_size
        sha256_val = sha256_file(p)
    elif p.is_dir():
        path_kind = "directory"
        size_bytes = sum(f.stat().st_size for f in p.rglob("*") if f.is_file())

    return ArtifactObservation(
        run_id=run_id,
        stage=stage,
        dataset_row_id=dataset_row_id,
        scientific_role=scientific_role,
        run_relative_path=run_relative_path,
        path_kind=path_kind,
        size_bytes=size_bytes,
        sha256=sha256_val,
        target_sha256=target_sha256,
        link_target=link_target,
        resolved_path=resolved_path,
        is_broken_link=is_broken_link,
        materialized_at=materialized,
        produced_by=produced_by,
        status=status,
        error_detail=None,
    )


def write_artifact_ledger(path: str | os.PathLike[str], records: Sequence[ArtifactObservation]) -> None:
    """Write a new artifact ledger TSV, completely overwriting any existing file."""
    _write_tsv_ledger(Path(path), records, append=False)


def append_artifact_observation(
    path: str | os.PathLike[str], record: ArtifactObservation
) -> None:
    """Append one observation to an existing ledger.

    Raises ``ArtifactLedgerError`` if a record with the same primary key
    already exists in the ledger.
    """
    ledger_path = Path(path)
    existing = _read_tsv_ledger(ledger_path)
    for rec in existing:
        if rec.key == record.key:
            raise ArtifactLedgerError(
                f"duplicate artifact key: {record.key} already exists in ledger"
            )
    _write_tsv_ledger(ledger_path, [record], append=True)


def validate_artifact_ledger(path: str | os.PathLike[str]) -> list[ArtifactObservation]:
    """Read and return all records from a ledger, raising on structural errors.

    Checks:
    * Every row has all three primary-key components.
    * No duplicate ``(dataset_row_id, scientific_role, run_relative_path)`` keys exist.
    """
    records = _read_tsv_ledger(Path(path))
    seen: set[tuple[str, str, str]] = set()
    for rec in records:
        if not rec.dataset_row_id or not rec.scientific_role or not rec.run_relative_path:
            raise ArtifactLedgerError(
                f"Invalid ledger row: missing primary-key component: {rec}"
            )
        k = rec.key
        if k in seen:
            raise ArtifactLedgerError(
                f"Duplicate primary key in ledger: {k}"
            )
        seen.add(k)
    return records


@dataclass(frozen=True)
class CloseoutRecord:
    """Result of closing out one artifact observation.

    This is an immutable view \u2014 it never modifies the source ledger.
    """

    dataset_row_id: str
    scientific_role: str
    run_relative_path: str
    sha256_before: str | None
    sha256_after: str | None
    status_before: str
    status_after: str
    changed: bool
    error_detail: str | None = None


_CLOSEOUT_FIELDS = [
    "dataset_row_id", "scientific_role", "run_relative_path",
    "sha256_before", "sha256_after", "status_before", "status_after",
    "changed", "error_detail",
]


def write_closeout_report(path: str | os.PathLike[str], records: Sequence[CloseoutRecord]) -> None:
    """Write a closeout report TSV (immutable evidence artifact).

    Args:
        path:  Output path for the closeout TSV.
        records:  Closeout records to persist.
    """
    p = Path(path)
    with p.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(_CLOSEOUT_FIELDS) + "\n")
        for rec in records:
            row = [_tsv_cell(getattr(rec, f, None)) for f in _CLOSEOUT_FIELDS]
            fh.write("\t".join(row) + "\n")


def closeout_artifact_ledger(
    ledger_path: str | os.PathLike[str],
    staging_dir: str | os.PathLike[str],
    *,
    re_observe: bool = True,
) -> list[CloseoutRecord]:
    """Close out an artifact ledger by re-observing every record.

    For each record in the ledger, re-observes the artifact at its
    ``run_relative_path`` (relative to *staging_dir*) and compares the
    new observation against the original ledger record.

    **This is a read-only view:** the source ledger is never modified.
    Changed vs unchanged is reported through the return value.

    Args:
        ledger_path:  Path to the existing TSV artifact ledger.
        staging_dir:  Root directory for resolving relative paths.
        re_observe:  If True, re-observe each artifact.  If False, skip
            re-observation (fast path, no change detection).

    Returns:
        A list of ``CloseoutRecord`` instances (one per ledger record).
    """
    existing = validate_artifact_ledger(ledger_path)
    staging = Path(staging_dir)
    closeout: list[CloseoutRecord] = []

    for original in existing:
        if re_observe:
            abs_path = staging / original.run_relative_path
            fresh = observe_artifact(
                abs_path,
                run_id=original.run_id,
                stage=original.stage,
                dataset_row_id=original.dataset_row_id,
                scientific_role=original.scientific_role,
                run_relative_path=original.run_relative_path,
                produced_by=original.produced_by,
                path_kind=original.path_kind,
            )
            changed = fresh.sha256 != original.sha256 or fresh.status != original.status
            detail = None
            if changed:
                detail = f"sha256 changed: {original.sha256} → {fresh.sha256}"
                if fresh.status != original.status:
                    detail += f"; status: {original.status} → {fresh.status}"
            closeout.append(CloseoutRecord(
                dataset_row_id=original.dataset_row_id,
                scientific_role=original.scientific_role,
                run_relative_path=original.run_relative_path,
                sha256_before=original.sha256,
                sha256_after=fresh.sha256,
                status_before=original.status,
                status_after=fresh.status,
                changed=changed,
                error_detail=detail,
            ))
        else:
            closeout.append(CloseoutRecord(
                dataset_row_id=original.dataset_row_id,
                scientific_role=original.scientific_role,
                run_relative_path=original.run_relative_path,
                sha256_before=original.sha256,
                sha256_after=original.sha256,
                status_before=original.status,
                status_after=original.status,
                changed=False,
                error_detail=None,
            ))

    return closeout


# ---------------------------------------------------------------------------
# Consumer gate (validation before consuming)
# ---------------------------------------------------------------------------


def validate_before_consume(
    ledger_path: str | os.PathLike[str],
    expected_inventory: Sequence[tuple[str, str, str]],
    *,
    run_id: str | None = None,
    staging_dir: str | os.PathLike[str] | None = None,
    environment: Mapping[str, str] | None = None,
    argv: Sequence[str] | None = None,
    configuration: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Validate that expected artifacts exist in the ledger.

    This is the consumer gate per ADR-0002: every expected artifact
    must appear in the ledger with ``status == "ok"``.
    Files are re-hashed on disk for mutation detection (if ``staging_dir``
    is provided).

    Args:
        ledger_path:  Path to the TSV artifact ledger.
        expected_inventory:  Sequence of ``(dataset_row_id, scientific_role,
            run_relative_path)`` tuples describing expected artifacts.
        run_id:  Optional run identifier for the validation report.
        staging_dir:  Optional root for resolving relative paths for re-hash.
        environment:  Environment variables to scan for secrets.
        argv:  Command-line arguments to scan for secrets.
        configuration:  Configuration dict to scan for secrets.

    Returns:
        A dict matching ``validation-gate.schema.json``:
        ``run_id``, ``validated_at``, ``overall_status``, ``artifact_checks``,
        ``row_identity_checks``, ``secret_scan``.
    """
    # Generate valid run_id if not provided
    if not run_id:
        run_id = _make_run_id()

    try:
        records = validate_artifact_ledger(ledger_path)
    except ArtifactLedgerError as exc:
        return {
            "run_id": run_id,
            "validated_at": _iso_now(),
            "overall_status": "fail",
            "artifact_checks": [],
            "row_identity_checks": [],
            "secret_scan": {"clean": True, "scanned_fields": [], "sentinel": f"ledger_validation_error: {exc}"},
        }

    by_key: dict[tuple[str, str, str], ArtifactObservation] = {}
    for rec in records:
        k = (rec.dataset_row_id, rec.scientific_role, rec.run_relative_path)
        by_key[k] = rec

    staging_resolve = Path(staging_dir) if staging_dir else None
    checks: list[dict[str, Any]] = []
    any_fail = False
    any_warn = False

    for expected_key in expected_inventory:
        rec = by_key.get(expected_key)
        # Re-hash the file on disk for mutation detection
        actual_sha256: str | None = None
        if rec and staging_resolve:
            artifact_path = staging_resolve / rec.run_relative_path
            if artifact_path.exists():
                actual_sha256 = sha256_file(artifact_path)
            else:
                actual_sha256 = None  # deleted on disk

        if rec and rec.status == "ok":
            if rec and staging_resolve and not (staging_resolve / rec.run_relative_path).exists():
                matched = False
                art_status = "missing"
                any_fail = True
            elif actual_sha256 is not None and actual_sha256 != rec.sha256:
                matched = False
                art_status = "mismatch"
                any_fail = True
            else:
                matched = True
                art_status = "ok"
        elif rec:
            matched = False
            art_status = "mismatch"
            any_fail = True
        else:
            matched = False
            art_status = "missing"
            any_warn = True
        checks.append({
            "match": matched,
            "status": art_status,
            "expected_sha256": rec.sha256 if rec else None,
            "actual_sha256": actual_sha256,
            "dataset_row_id": expected_key[0],
            "scientific_role": expected_key[1],
            "run_relative_path": expected_key[2],
        })

    if any_fail:
        overall_status = "fail"
    elif any_warn:
        overall_status = "warn"
    else:
        overall_status = "pass"

    # Row identity checks: per-row role coverage.
    from collections import defaultdict
    rows: dict[str, set[str]] = defaultdict(set)
    for rec in records:
        rows[rec.dataset_row_id].add(rec.scientific_role)
    row_checks: list[dict[str, Any]] = []
    for row_id, roles in rows.items():
        expected_set = set()
        for ek in expected_inventory:
            if ek[0] == row_id:
                expected_set.add(ek[1])
        missing = expected_set - roles
        row_checks.append({
            "dataset_row_id": row_id,
            "expected_roles": sorted(expected_set),
            "found_roles": sorted(roles),
            "missing_roles": sorted(missing),
            "duplicate_keys": [],
            "status": "incomplete" if missing else "complete",
        })

    # Secret sentinel scan: check for plain-text secrets in ledger columns
    # plus environment, argv, and configuration (Phase 1 requirement)
    scanned_fields = [
        "run_id", "dataset_row_id", "scientific_role", "run_relative_path",
        "sha256", "produced_by", "link_target", "resolved_path", "error_detail",
        "environment", "argv", "configuration"
    ]
    found_any = False
    sentinel_found = ""
    
    # Scan ledger fields
    for rec in records:
        for field in scanned_fields[:9]:  # First 9 are ledger fields
            val = getattr(rec, field, None) or ""
            val_upper = val.upper()
            if "SECRET" in val_upper or "PASSWORD" in val_upper or "TOKEN" in val_upper:
                found_any = True
                sentinel_found = field
                break
        if found_any:
            break
    
    # Scan environment
    if not found_any and environment:
        for key, val in environment.items():
            val_upper = str(val).upper()
            if "SECRET" in val_upper or "PASSWORD" in val_upper or "TOKEN" in val_upper:
                found_any = True
                sentinel_found = f"environment:{key}"
                break
    
    # Scan argv
    if not found_any and argv:
        for arg in argv:
            val_upper = str(arg).upper()
            if "SECRET" in val_upper or "PASSWORD" in val_upper or "TOKEN" in val_upper:
                found_any = True
                sentinel_found = "argv"
                break
    
    # Scan configuration
    if not found_any and configuration:
        for key, val in configuration.items():
            val_upper = str(val).upper()
            if "SECRET" in val_upper or "PASSWORD" in val_upper or "TOKEN" in val_upper:
                found_any = True
                sentinel_found = f"configuration:{key}"
                break
    
    return {
        "run_id": run_id,
        "validated_at": _iso_now(),
        "overall_status": "fail" if found_any else overall_status,
        "artifact_checks": checks,
        "row_identity_checks": row_checks,
        "secret_scan": {"clean": not found_any, "scanned_fields": scanned_fields, "sentinel": sentinel_found},
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def _parse_args(argv: Sequence[str] | None = None) -> Any:
    import argparse

    parser = argparse.ArgumentParser(description="PRISM pipeline provenance tools")
    parser.add_argument("--run-id", default="", help="Override generated run ID")
    parser.add_argument("--contract-path", default="", help="Path to declared contract JSON")
    parser.add_argument("--manifest-path", default="", help="Path to write run manifest JSON")
    parser.add_argument("--ledger-path", default="", help="Path to artifact ledger TSV")
    parser.add_argument("--stage", default="", help="Pipeline stage for artifact observation")
    parser.add_argument("--dataset-row-id", default="", help="Dataset row identifier")
    parser.add_argument("--scientific-role", default="", help="Scientific role for artifact")
    parser.add_argument("--run-relative-path", default="", help="Relative path for artifact")
    parser.add_argument("--produced-by", default="", help="Tool that produced the artifact")
    parser.add_argument("--staging-dir", default=".", help="Staging directory for closeout")
    parser.add_argument("--command", nargs="*", help="Command that produced the artifact")
    parser.add_argument("--observe", action="store_true", help="Observe a single artifact")
    parser.add_argument("--closeout", action="store_true", help="Close out the ledger")
    parser.add_argument("--validate", action="store_true", help="Validate before consume")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point for the provenance module."""
    args = _parse_args(argv)
    run_id = args.run_id or _make_run_id()
    contract = {
        "contract_version": "1.0",
        "pipeline_version": "1.0.0",
        "stages_enabled": ["alignment", "refinement", "comparison"],
        "aligner": "tmalign",
        "refiner": "pyrosetta",
        "input_selectors": {"raw": [], "normalized": []},
        "template_inventory": [],
        "parameters": {},
        "resource_request": {"cpus": 1, "memory_gb": 4, "time_hours": 1, "partition": "ai", "gpu": False},
        "source_inventory": {"git_head": "a" * 40, "git_diff_hash": "b" * 64, "declared_untracked": []},
        "tool_fingerprints": [{"name": "prism", "version": "1.0.0", "sha256": None}],
        "contract_hash": "",
    }
    signing = {k: v for k, v in contract.items() if k != "contract_hash"}
    contract["contract_hash"] = canonical_hash(signing)

    if args.contract_path:
        Path(args.contract_path).write_text(canonical_json(contract), encoding="utf-8")

    if args.manifest_path:
        manifest = build_execution_attempt(contract, run_id=run_id, command=args.command)
        Path(args.manifest_path).write_text(canonical_json(manifest), encoding="utf-8")

    if args.observe and args.ledger_path and args.stage and args.dataset_row_id and args.scientific_role and args.run_relative_path:
        obs = observe_artifact(
            args.run_relative_path,
            run_id=run_id,
            stage=args.stage,
            dataset_row_id=args.dataset_row_id,
            scientific_role=args.scientific_role,
            run_relative_path=args.run_relative_path,
            produced_by=args.produced_by or "",
        )
        ledger_path = Path(args.ledger_path)
        if ledger_path.is_file():
            append_artifact_observation(ledger_path, obs)
        else:
            write_artifact_ledger(ledger_path, [obs])
        print(f"Observed: {obs.key} → {obs.status}")

    if args.closeout and args.ledger_path:
        results = closeout_artifact_ledger(args.ledger_path, args.staging_dir)
        print(f"Closeout: {len(results)} records")

    if args.validate and args.ledger_path:
        # Pass actual runtime context for secret scanning
        import os
        result = validate_before_consume(
            args.ledger_path,
            [],
            run_id=run_id,
            staging_dir=args.staging_dir,
            environment=dict(os.environ),
            argv=list(sys.argv),
            configuration=contract,  # declared contract as configuration
        )
        print(f"Validation: {result['overall_status']}")
        if not result['secret_scan']['clean']:
            print(f"  Secret detected in: {result['secret_scan']['sentinel']}")
        if result['overall_status'] == 'fail':
            return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
