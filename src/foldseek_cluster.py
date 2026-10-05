"""Interface-aware Foldseek multimer clustering for PRISM donor complexes.

The installed Foldseek build exposes ``easy-multimercluster`` rather than the
newer ``easy-interfacecluster`` wrapper.  This module keeps the external
workflow small and explicit: donor PDBs are indexed as two-chain, ordered
records; Foldseek's donor-level adjacency list is parsed; and each donor is
expanded back to its ordered chain pair in the shared PRISM cluster-index
contract.

No Foldseek process is started by the pure parsing/command helpers.  The CLI's
``--dry-run`` path only writes a command manifest and an empty index.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


INDEX_SCHEMA = "prism-template-cluster-index/v1"
DEFAULT_FOLDSEEK_PATH = "/home/rshadi25/.conda/envs/gtalign_env/bin/foldseek"
DEFAULT_OUTPUT_NAME = "cluster_index.json"
STRUCTURE_SUFFIXES = (".pdb", ".cif", ".mmcif", ".pdb.gz", ".cif.gz", ".mmcif.gz")


class MembershipError(ValueError):
    """Foldseek output cannot be reconciled with the input donor panel."""


@dataclass(frozen=True)
class FoldseekConfig:
    input_dir: Path
    output_dir: Path
    panel: Path | None = None
    panel_sha256: str | None = None
    foldseek_path: str = DEFAULT_FOLDSEEK_PATH
    threads: int = 1
    gpu: int = 1
    distance_threshold: float = 8.0
    extraction_mode: int = 1
    interface_lddt_threshold: float = 0.0
    chain_tm_threshold: float = 0.0
    multimer_tm_threshold: float = 0.0
    coverage_threshold: float = 0.8
    cov_mode: int = 0
    alignment_type: int = 2
    sensitivity: float = 4.0
    max_seqs: int = 300
    prefilter_mode: int = 0
    exhaustive_search: int = 0
    remove_tmp_files: int = 1
    missing_member_policy: str = "error"
    dry_run: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "input_dir", Path(self.input_dir))
        object.__setattr__(self, "output_dir", Path(self.output_dir))
        if self.panel is not None:
            object.__setattr__(self, "panel", Path(self.panel))
        if self.threads < 1:
            raise ValueError("threads must be positive")
        if self.gpu not in (0, 1):
            raise ValueError("gpu must be 0 or 1; requested GPU is never silently disabled")
        if self.extraction_mode not in (0, 1):
            raise ValueError("extraction_mode must be 0 (chains) or 1 (interfaces)")
        if self.cov_mode not in range(6):
            raise ValueError("cov_mode must be between 0 and 5")
        if self.alignment_type not in (0, 1, 2):
            raise ValueError("alignment_type must be 0, 1, or 2")
        if self.sensitivity <= 0:
            raise ValueError("sensitivity must be positive")
        if self.max_seqs < 1:
            raise ValueError("max_seqs must be positive")
        if self.prefilter_mode not in (0, 1, 2, 3):
            raise ValueError("prefilter_mode must be between 0 and 3")
        if self.exhaustive_search not in (0, 1):
            raise ValueError("exhaustive_search must be 0 or 1")
        if self.remove_tmp_files not in (0, 1):
            raise ValueError("remove_tmp_files must be 0 or 1")
        if self.missing_member_policy not in {"error", "singleton"}:
            raise ValueError("missing_member_policy must be 'error' or 'singleton'")
        if self.distance_threshold < 0:
            raise ValueError("distance_threshold must be non-negative")
        for name in (
            "interface_lddt_threshold",
            "chain_tm_threshold",
            "multimer_tm_threshold",
            "coverage_threshold",
        ):
            value = float(getattr(self, name))
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")


def _format_number(value: float | int) -> str:
    """Format a numeric CLI value without adding non-semantic zeroes."""

    return f"{float(value):g}"


def build_foldseek_command(
    config: FoldseekConfig,
    output_prefix: str | Path,
    tmp_dir: str | Path | None = None,
) -> list[str]:
    """Return the exact ``easy-multimercluster`` argv for *config*.

    ``--gpu`` is always present, including for ``gpu=0``.  Recording the
    explicit value prevents a caller from accidentally changing a requested
    GPU run into a CPU run through an omitted/default flag.
    """

    prefix = Path(output_prefix)
    temporary = Path(tmp_dir) if tmp_dir is not None else prefix.parent / ".foldseek_tmp"
    return [
        str(config.foldseek_path),
        "easy-multimercluster",
        str(config.input_dir),
        str(prefix),
        str(temporary),
        "--db-extraction-mode",
        str(config.extraction_mode),
        "--multimer-report-mode",
        "1",
        "--distance-threshold",
        _format_number(config.distance_threshold),
        "--interface-lddt-threshold",
        _format_number(config.interface_lddt_threshold),
        "--chain-tm-threshold",
        _format_number(config.chain_tm_threshold),
        "--multimer-tm-threshold",
        _format_number(config.multimer_tm_threshold),
        "-c",
        _format_number(config.coverage_threshold),
        "--cov-mode",
        str(config.cov_mode),
        "--alignment-type",
        str(config.alignment_type),
        "-s",
        _format_number(config.sensitivity),
        "--max-seqs",
        str(config.max_seqs),
        "--prefilter-mode",
        str(config.prefilter_mode),
        "--exhaustive-search",
        str(config.exhaustive_search),
        "--threads",
        str(config.threads),
        "--gpu",
        str(config.gpu),
        "--remove-tmp-files",
        str(config.remove_tmp_files),
    ]


def resolve_panel_hash(panel: str | Path | None, expected: str | None = None) -> str:
    """Hash a panel file and fail closed when an expected hash disagrees."""

    if panel is None:
        if expected:
            raise ValueError("panel_sha256 was supplied without --panel")
        return "UNKNOWN"
    source = Path(panel)
    if not source.is_file():
        raise FileNotFoundError(f"panel does not exist: {source}")
    actual = hashlib.sha256(source.read_bytes()).hexdigest()
    if expected and actual.lower() != str(expected).strip().lower():
        raise ValueError(f"panel SHA-256 mismatch: expected {expected}, observed {actual}")
    return actual


def _structure_files(input_dir: Path) -> list[Path]:
    if not input_dir.is_dir():
        raise FileNotFoundError(f"input directory does not exist: {input_dir}")
    return sorted(
        (p for p in input_dir.iterdir() if p.is_file() and p.name.lower().endswith(STRUCTURE_SUFFIXES)),
        key=lambda p: p.name,
    )


def _open_structure(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _chain_ids_from_structure(path: Path) -> list[str]:
    """Read ordered PDB chain IDs without normalizing away donor orientation."""

    if not path.name.lower().endswith(".pdb") and not path.name.lower().endswith(".pdb.gz"):
        raise ValueError(f"Foldseek donor inputs must be PDB files: {path}")
    chains: list[str] = []
    with _open_structure(path) as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain = line[21:22] if len(line) > 21 else ""
            chain = chain.strip() or "_"
            if chain not in chains:
                chains.append(chain)
    return chains


def collect_donor_inputs(input_dir: str | Path) -> list[dict[str, Any]]:
    """Collect and validate ordered, exactly-two-chain donor PDBs.

    Each record is JSON-ready and uses ``donor_id`` (the filename stem) as the
    stable Foldseek mapping key.  Chain order is the first-seen PDB order and
    is intentionally not sorted.
    """

    source_dir = Path(input_dir)
    paths = _structure_files(source_dir)
    records: list[dict[str, Any]] = []
    by_id: dict[str, Path] = {}
    for path in paths:
        donor_id = path.name
        for suffix in (".pdb.gz", ".cif.gz", ".mmcif.gz", ".pdb", ".cif", ".mmcif"):
            if donor_id.lower().endswith(suffix):
                donor_id = donor_id[: -len(suffix)]
                break
        if donor_id in by_id:
            raise MembershipError(f"duplicate donor input ID {donor_id!r}: {by_id[donor_id]} and {path}")
        chains = _chain_ids_from_structure(path)
        if len(chains) != 2:
            raise ValueError(f"donor {path} must contain exactly two chains in order; observed {chains}")
        by_id[donor_id] = path
        records.append(
            {
                "donor_id": donor_id,
                "template_id": donor_id,
                "path": str(path.resolve()),
                "chain_ids": chains,
                "orientation": f"{chains[0]}>{chains[1]}",
            }
        )
    return records


def _canonical_token(token: str) -> str:
    text = Path(str(token).strip()).name
    text = re.sub(r"\.(?:pdb|cif|mmcif)(?:\.gz)?$", "", text, flags=re.IGNORECASE)
    return text


def _resolve_donor_token(token: str, records: Sequence[Mapping[str, Any]]) -> str:
    """Map Foldseek's shortened/derived IDs back to one donor, or fail."""

    raw = _canonical_token(token)
    # Interface extraction can append ``_INT_1`` and a chain ID to headers.
    stripped = re.sub(r"_INT(?:ERFACE)?(?:_\d+)?(?:_[^_\s]+)?$", "", raw, flags=re.IGNORECASE)
    candidates: list[str] = []
    for record in records:
        donor = str(record["donor_id"])
        donor_raw = _canonical_token(donor)
        options = {donor_raw, donor_raw.lower()}
        raw_options = {raw, raw.lower(), stripped, stripped.lower()}
        if raw_options & options or stripped == donor_raw:
            candidates.append(donor)
            continue
        # MMseqs/Foldseek often uses the first PDB token (e.g. ``1ABC``) for
        # a file named ``1ABC_receptor_ligand.pdb``.  Prefix matching is only
        # accepted when it is unique; ambiguity is a hard error.
        if donor_raw.lower().startswith(raw.lower() + "_") or raw.lower().startswith(donor_raw.lower() + "_"):
            candidates.append(donor)
    if len(candidates) != 1:
        if not candidates:
            raise MembershipError(f"unaccounted Foldseek member {token!r}")
        raise MembershipError(f"ambiguous Foldseek member {token!r}: {sorted(candidates)}")
    return candidates[0]


def _read_rows(path: Path) -> list[list[str]]:
    if not path.is_file():
        return []
    # ``cluster_report`` from Foldseek 10.941cd33 includes NUL bytes between
    # selected formatted values.  Treat them as separators, not row loss.
    text = path.read_bytes().decode("utf-8", errors="replace").replace("\x00", "")
    rows: list[list[str]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            fields = line.split()
        rows.append([field.strip() for field in fields])
    return rows


def _to_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _metric_key(value: str) -> str:
    key = re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")
    aliases = {
        "interface_lddt_score": "interface_lddt",
        "interface_lddt": "interface_lddt",
        "chain_tm_score": "chain_tm",
        "chain_tm": "chain_tm",
        "multimer_tm_score": "multimer_tm",
        "multimer_tm": "multimer_tm",
        "coverage": "coverage",
        "cov": "coverage",
    }
    return aliases.get(key, key)


def _parse_report(path: Path, records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = _read_rows(path)
    if not rows:
        return []
    first = rows[0]
    first_keys = {_metric_key(first[0]), _metric_key(first[1])}
    # Foldseek's native report has no header and starts with structure IDs;
    # only treat the first row as a header when its first two fields are
    # recognisable column names.
    has_header = bool(first_keys) and first_keys <= {"query", "target", "q", "t"}
    if has_header:
        names = [_metric_key(value) for value in first[2:]]
        rows = rows[1:]
    else:
        # ``createmultimerreport`` in Foldseek 10.941cd33 emits (after query
        # and target) query/target multimer TM, query/target chain TM,
        # interface LDDT, then a rotation and translation vector.  Keep the
        # raw fields as well, but expose the first directional values under
        # the threshold names used by the command contract.
        names = [f"metric_{index}" for index in range(1, max(0, len(first) - 1))]
    parsed: list[dict[str, Any]] = []
    for row in rows:
        if len(row) < 2:
            continue
        try:
            query = _resolve_donor_token(row[0], records)
            target = _resolve_donor_token(row[1], records)
        except MembershipError:
            # Reports can include an interface-specific ID that is not in the
            # adjacency list.  It remains useful provenance, but an unknown
            # report row must not create cluster membership.
            continue
        values: dict[str, Any] = {}
        for name, raw in zip(names, row[2:]):
            number = _to_float(raw)
            values[name] = number if number is not None else raw
        if not has_header and len(row) >= 7:
            values["multimer_tm"] = _to_float(row[2])
            values["chain_tm"] = _to_float(row[4])
            values["interface_lddt"] = _to_float(row[6])
        values.update({"query": query, "target": target, "raw": row})
        parsed.append(values)
    return parsed


def _parse_representative_ids(path: Path, records: Sequence[Mapping[str, Any]]) -> set[str]:
    if not path.is_file():
        return set()
    reps: set[str] = set()
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(">"):
            continue
        token = line[1:].split()[0]
        reps.add(_resolve_donor_token(token, records))
    return reps


def _numeric(mapping: Mapping[str, Any], *keys: str, default: float | None = None) -> float | None:
    for key in keys:
        value = mapping.get(key)
        if isinstance(value, (int, float)):
            return float(value)
    return default


def select_representatives(
    members: Iterable[Mapping[str, Any]],
    *,
    coverage_threshold: float = 0.8,
    native_representative: str | None = None,
    max_representatives: int = 3,
) -> tuple[list[dict[str, Any]], str, dict[str, Any]]:
    """Select deterministic native/medoid-plus-diversity representatives.

    ``coverage`` is interpreted as the incremental covered fraction supplied
    by a candidate; this keeps the policy dependency-light and lets callers
    supply union coverage from a richer metric source.  When no metric is
    available, one candidate is considered complete coverage.  A cluster that
    still cannot reach the configured threshold after three representatives is
    explicitly marked ``SPLIT_REQUIRED``.
    """

    rows = [dict(row) for row in members]
    rows.sort(key=lambda row: str(row.get("donor_id", row.get("template_id", ""))))
    if not rows:
        return [], "EMPTY", {"coverage_threshold": coverage_threshold, "coverage_observed": 0.0}
    for row in rows:
        row.setdefault("donor_id", row.get("template_id"))
        row.setdefault("coverage", _numeric(row, "coverage_score", default=None))
        row.setdefault("medoid_score", _numeric(row, "medoid", "similarity", default=0.0) or 0.0)
        row.setdefault("diversity", _numeric(row, "diversity_score", default=0.0) or 0.0)

    native = native_representative
    if native is None:
        native = next(
            (str(row["donor_id"]) for row in rows if row.get("native_representative") is True),
            None,
        )
    if native is None or native not in {str(row["donor_id"]) for row in rows}:
        native = max(rows, key=lambda row: (float(row.get("medoid_score", 0.0)), str(row["donor_id"])))["donor_id"]

    selected: list[dict[str, Any]] = [next(row for row in rows if str(row["donor_id"]) == str(native))]
    observed = float(selected[0].get("coverage") or 0.0)
    if selected[0].get("coverage") is None:
        observed = 1.0
    while len(selected) < max_representatives and observed + 1e-12 < coverage_threshold:
        remaining = [row for row in rows if row not in selected]
        if not remaining:
            break
        # Coverage drives selection; diversity and medoid score are stable
        # tie-breakers.  The donor ID makes the result independent of input
        # filesystem/order and Python hash randomization.
        candidate = max(
            remaining,
            key=lambda row: (
                float(row.get("coverage") or 0.0),
                float(row.get("diversity") or 0.0),
                float(row.get("medoid_score") or 0.0),
                str(row["donor_id"]),
            ),
        )
        selected.append(candidate)
        observed += float(candidate.get("coverage") or 0.0)
    complete = observed + 1e-12 >= coverage_threshold
    status = "PASS" if complete else "SPLIT_REQUIRED"
    metadata = {
        "selection_method": "native_then_coverage_diversity_medoid",
        "native_representative": str(native),
        "coverage_threshold": float(coverage_threshold),
        "coverage_observed": round(observed, 12),
        "coverage_complete": complete,
        "representative_count": len(selected),
    }
    return selected, status, metadata


def _member_chain_rows(record: Mapping[str, Any]) -> list[dict[str, Any]]:
    chains = list(record["chain_ids"])
    orientation = str(record.get("orientation") or f"{chains[0]}>{chains[1]}")
    rows: list[dict[str, Any]] = []
    for index, chain in enumerate(chains):
        rows.append(
            {
                "template_id": str(record.get("template_id", record["donor_id"])),
                "chain_id": str(chain),
                "path": str(record["path"]),
                "metadata": {
                    "donor_id": str(record["donor_id"]),
                    "chain_role": "first" if index == 0 else "second",
                    "chain_order": index,
                    "orientation": orientation,
                },
            }
        )
    return rows


def _output_paths(prefix: Path) -> tuple[Path, Path, Path]:
    if prefix.is_dir():
        prefix = prefix / "foldseek"
    if prefix.name.endswith("_cluster.tsv"):
        base = prefix.with_name(prefix.name[: -len("_cluster.tsv")])
    elif prefix.name.endswith("_cluster"):
        base = prefix.with_name(prefix.name[: -len("_cluster")])
    else:
        base = prefix
    return (
        base.with_name(base.name + "_cluster.tsv"),
        base.with_name(base.name + "_cluster_report"),
        base.with_name(base.name + "_rep_seq.fasta"),
    )


def build_index_from_outputs(
    output_prefix: str | Path,
    records: Sequence[Mapping[str, Any]],
    *,
    panel_sha256: str = "UNKNOWN",
    tool_version: str = "UNKNOWN",
    parameters: Mapping[str, Any] | None = None,
    coverage_threshold: float = 0.8,
    missing_member_policy: str = "error",
) -> dict[str, Any]:
    """Parse Foldseek files and emit the shared PRISM index contract."""

    if missing_member_policy not in {"error", "singleton"}:
        raise ValueError("missing_member_policy must be 'error' or 'singleton'")

    params = dict(parameters or {})
    params.setdefault("coverage_threshold", coverage_threshold)
    params["input_state"] = "EMPTY" if not records else "NONEMPTY"
    base = Path(output_prefix)
    cluster_path, report_path, representative_path = _output_paths(base)
    if not records:
        params["input_state"] = "EMPTY"
        return {
            "schema_version": INDEX_SCHEMA,
            "tool": "foldseek",
            "tool_version": tool_version,
            "panel_sha256": panel_sha256,
            "parameters": params,
            "clusters": [],
        }
    rows = _read_rows(cluster_path)
    if not rows:
        raise MembershipError(f"missing Foldseek cluster output: {cluster_path}")
    by_id = {str(record["donor_id"]): record for record in records}
    groups: dict[str, list[str]] = {}
    seen: list[str] = []
    for row in rows:
        if len(row) < 2:
            raise ValueError(f"malformed Foldseek cluster row: {row!r}")
        root = _resolve_donor_token(row[0], records)
        member = _resolve_donor_token(row[1], records)
        groups.setdefault(root, [])
        if member in seen:
            raise MembershipError(f"duplicate Foldseek cluster member {member!r}")
        groups[root].append(member)
        seen.append(member)
    expected = set(by_id)
    observed = set(seen)
    missing = sorted(expected - observed)
    missing_set = set(missing)
    unaccounted = sorted(observed - expected)
    if missing:
        if missing_member_policy == "error":
            raise MembershipError(f"missing Foldseek cluster members: {missing}")
        # Foldseek 10.941's multimer clustering can omit short/unsupported
        # donor complexes from the adjacency output even though the input
        # directory was accepted.  Preserve exact panel accounting by adding
        # those donors as explicit conservative singletons.  They are never
        # presented as Foldseek-supported matches; routing can safely fall
        # back to the original donor for these records.
        for donor_id in missing:
            groups[donor_id] = [donor_id]
        params["missing_member_fallback_count"] = len(missing)
        params["missing_member_fallback_ids_sha256"] = hashlib.sha256(
            ("\n".join(missing) + "\n").encode("utf-8")
        ).hexdigest()
    if unaccounted:
        raise MembershipError(f"unaccounted Foldseek cluster members: {unaccounted}")
    report_rows = _parse_report(report_path, records)
    report_by_pair = {(row["query"], row["target"]): row for row in report_rows}
    fasta_reps = _parse_representative_ids(representative_path, records)

    clusters: list[dict[str, Any]] = []
    for root in sorted(groups):
        donor_ids = groups[root]
        donor_rows: list[dict[str, Any]] = []
        relevant_reports: list[dict[str, Any]] = []
        for donor_id in donor_ids:
            record = by_id[donor_id]
            metrics = dict(report_by_pair.get((root, donor_id), report_by_pair.get((donor_id, root), {})))
            if metrics:
                relevant_reports.append(metrics)
            similarity = _numeric(metrics, "multimer_tm", "chain_tm", "interface_lddt", default=0.0) or 0.0
            coverage = _numeric(metrics, "coverage", default=None)
            donor_rows.append(
                {
                    "donor_id": donor_id,
                    "path": record["path"],
                    "chain_ids": list(record["chain_ids"]),
                    "coverage": coverage,
                    "medoid_score": similarity,
                    "diversity": 1.0 - similarity,
                }
            )
        # A report without a coverage column cannot justify splitting; the
        # native representative is retained as a complete fallback.
        if all(row.get("coverage") is None for row in donor_rows):
            for row in donor_rows:
                row["coverage"] = 1.0 if row["donor_id"] == root else 0.0
        selected, selection_status, selection_metadata = select_representatives(
            donor_rows,
            coverage_threshold=coverage_threshold,
            native_representative=root if root in donor_ids else None,
        )
        member_rows = [chain for donor_id in donor_ids for chain in _member_chain_rows(by_id[donor_id])]
        representative_rows = [
            chain for donor in selected for chain in _member_chain_rows(by_id[str(donor["donor_id"])])
        ]
        status = "SINGLETON" if len(donor_ids) == 1 else selection_status
        fallback_singleton = (
            missing_member_policy == "singleton"
            and len(donor_ids) == 1
            and donor_ids[0] in missing_set
        )
        clusters.append(
            {
                "cluster_id": root,
                "status": status,
                "members": member_rows,
                "representatives": representative_rows,
                "metadata": {
                    **selection_metadata,
                    "donor_count": len(donor_ids),
                    "chain_count": len(member_rows),
                    "native_representative": root,
                    "native_representative_in_fasta": root in fasta_reps,
                    "index_membership_source": "foldseek_output" if not fallback_singleton else "singleton_fallback",
                    "fallback_reason": (
                        "foldseek_output_omitted_donor" if fallback_singleton else None
                    ),
                    "report_rows": relevant_reports,
                    "orientation_pairs": [
                        {
                            "donor_id": donor_id,
                            "orientation": by_id[donor_id]["orientation"],
                        }
                        for donor_id in donor_ids
                    ],
                },
            }
        )
    params["cluster_output"] = str(cluster_path)
    params["cluster_report"] = str(report_path)
    params["representative_output"] = str(representative_path)
    return {
        "schema_version": INDEX_SCHEMA,
        "tool": "foldseek",
        "tool_version": tool_version,
        "panel_sha256": panel_sha256,
        "parameters": params,
        "clusters": clusters,
    }


def _probe_tool(path: str) -> str:
    result = subprocess.run([path, "version"], capture_output=True, text=True, check=False)
    if result.returncode:
        message = (result.stderr or result.stdout).strip()[:1000]
        raise RuntimeError(f"Foldseek version probe failed ({result.returncode}): {message}")
    version = (result.stdout or result.stderr).strip().splitlines()
    if not version:
        raise RuntimeError("Foldseek version probe returned no version")
    return version[-1].strip()


def check_gpu_capability(path: str, *, runner=None) -> dict[str, Any]:
    """Verify that this Foldseek command exposes an explicit GPU switch."""

    if runner is None:
        runner = subprocess.run
    result = runner([path, "easy-multimercluster", "--help"], capture_output=True, text=True, check=False)
    help_text = (result.stdout or "") + "\n" + (result.stderr or "")
    if result.returncode and "usage:" not in help_text.lower():
        raise RuntimeError(f"Foldseek GPU capability probe failed ({result.returncode})")
    if not re.search(r"(?:^|\s)--gpu(?:\s|$)", help_text):
        raise RuntimeError("Foldseek build does not expose --gpu; refusing requested GPU mode")
    return {"gpu_flag": "--gpu", "supported": True}


def _gpu_failure(text: str) -> bool:
    return bool(
        re.search(
            r"(?:gpu|cuda).*(?:not available|unavailable|failed|error|fallback|falling back)|"
            r"(?:falling back|fallback).*(?:cpu|gpu)",
            text,
            flags=re.IGNORECASE | re.DOTALL,
        )
    )


def execute_foldseek(command: Sequence[str], *, gpu: int) -> subprocess.CompletedProcess[str]:
    """Execute one exact command and reject an implicit GPU-to-CPU fallback."""

    result = subprocess.run(list(command), capture_output=True, text=True, check=False)
    combined = (result.stdout or "") + "\n" + (result.stderr or "")
    if result.returncode:
        raise RuntimeError(f"Foldseek failed ({result.returncode}): {combined[-4000:]}")
    if gpu and _gpu_failure(combined):
        raise RuntimeError("Foldseek reported GPU failure/fallback for requested --gpu 1")
    return result


def _json_parameters(config: FoldseekConfig, command: Sequence[str], version: str) -> dict[str, Any]:
    return {
        "input_dir": str(config.input_dir),
        "output_dir": str(config.output_dir),
        "foldseek_path": str(config.foldseek_path),
        "threads": config.threads,
        "gpu": config.gpu,
        "distance_threshold": config.distance_threshold,
        "extraction_mode": config.extraction_mode,
        "interface_lddt_threshold": config.interface_lddt_threshold,
        "chain_tm_threshold": config.chain_tm_threshold,
        "multimer_tm_threshold": config.multimer_tm_threshold,
        "coverage_threshold": config.coverage_threshold,
        "cov_mode": config.cov_mode,
        "alignment_type": config.alignment_type,
        "sensitivity": config.sensitivity,
        "max_seqs": config.max_seqs,
        "prefilter_mode": config.prefilter_mode,
        "exhaustive_search": config.exhaustive_search,
        "remove_tmp_files": config.remove_tmp_files,
        "missing_member_policy": config.missing_member_policy,
        "command": list(command),
        "command_string": shlex.join(list(command)),
        "tool_version": version,
        "dry_run": config.dry_run,
    }


def run_build(config: FoldseekConfig) -> dict[str, Any]:
    """Run the bounded Foldseek build, or emit its no-execution dry-run index."""

    panel_hash = resolve_panel_hash(config.panel, config.panel_sha256)
    records = collect_donor_inputs(config.input_dir)
    config.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = config.output_dir / "foldseek"
    tmp_dir = config.output_dir / ".foldseek_tmp"
    command = build_foldseek_command(config, prefix, tmp_dir)
    version = "DRY_RUN" if config.dry_run else _probe_tool(config.foldseek_path)
    if not config.dry_run:
        check_gpu_capability(config.foldseek_path)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        execute_foldseek(command, gpu=config.gpu)
        index = build_index_from_outputs(
            prefix,
            records,
            panel_sha256=panel_hash,
            tool_version=version,
            parameters=_json_parameters(config, command, version),
            coverage_threshold=config.coverage_threshold,
            missing_member_policy=config.missing_member_policy,
        )
    else:
        index = build_index_from_outputs(
            prefix,
            [],
            panel_sha256=panel_hash,
            tool_version=version,
            parameters=_json_parameters(config, command, version),
            coverage_threshold=config.coverage_threshold,
            missing_member_policy=config.missing_member_policy,
        )
        index["parameters"]["input_state"] = "DRY_RUN"
    destination = config.output_dir / DEFAULT_OUTPUT_NAME
    destination.write_text(json.dumps(index, indent=2, sort_keys=True) + "\n")
    return index


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--panel-sha256")
    parser.add_argument("--foldseek-path", default=DEFAULT_FOLDSEEK_PATH)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--gpu", type=int, default=1, choices=(0, 1))
    parser.add_argument("--distance-threshold", type=float, default=8.0)
    parser.add_argument("--extraction-mode", type=int, choices=(0, 1), default=1)
    parser.add_argument(
        "--interface-lddt-threshold",
        "--interface-lddt",
        dest="interface_lddt_threshold",
        type=float,
        default=0.0,
    )
    parser.add_argument("--chain-tm-threshold", "--chain-tm", dest="chain_tm_threshold", type=float, default=0.0)
    parser.add_argument(
        "--multimer-tm-threshold",
        "--multimer-tm",
        dest="multimer_tm_threshold",
        type=float,
        default=0.0,
    )
    parser.add_argument("--coverage-threshold", "--coverage", type=float, default=0.8)
    parser.add_argument("--cov-mode", type=int, choices=range(6), default=0)
    parser.add_argument("--alignment-type", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--sensitivity", "-s", type=float, default=4.0)
    parser.add_argument("--max-seqs", type=int, default=300)
    parser.add_argument("--prefilter-mode", type=int, choices=(0, 1, 2, 3), default=0)
    parser.add_argument("--exhaustive-search", type=int, choices=(0, 1), default=0)
    parser.add_argument("--remove-tmp-files", type=int, choices=(0, 1), default=1)
    parser.add_argument(
        "--missing-member-policy",
        choices=("error", "singleton"),
        default="error",
        help="handle donors omitted by Foldseek output (default: fail closed)",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    config = FoldseekConfig(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        panel=args.panel,
        panel_sha256=args.panel_sha256,
        foldseek_path=args.foldseek_path,
        threads=args.threads,
        gpu=args.gpu,
        distance_threshold=args.distance_threshold,
        extraction_mode=args.extraction_mode,
        interface_lddt_threshold=args.interface_lddt_threshold,
        chain_tm_threshold=args.chain_tm_threshold,
        multimer_tm_threshold=args.multimer_tm_threshold,
        coverage_threshold=args.coverage_threshold,
        cov_mode=args.cov_mode,
        alignment_type=args.alignment_type,
        sensitivity=args.sensitivity,
        max_seqs=args.max_seqs,
        prefilter_mode=args.prefilter_mode,
        exhaustive_search=args.exhaustive_search,
        remove_tmp_files=args.remove_tmp_files,
        missing_member_policy=args.missing_member_policy,
        dry_run=args.dry_run,
    )
    try:
        index = run_build(config)
    except (FileNotFoundError, MembershipError, ValueError, RuntimeError) as exc:
        print(f"foldseek cluster build failed: {exc}")
        return 2
    print(json.dumps(index, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
