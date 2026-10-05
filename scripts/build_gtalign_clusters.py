#!/usr/bin/env python3
"""Build a provenance-rich, representative-aware GTalign cluster index.

The output is a deliberately small JSON contract consumed by later template
selection steps::

    {
      "schema_version": "prism-template-cluster-index/v1",
      "tool": "GTalign",
      "tool_version": "0.19.00",
      "panel_count": 123,
      "panel_sha256": "...",
      "parameters": {"threshold": 0.5, ...},
      "clusters": [{
        "cluster_id": "cluster-0001",
        "status": "COMPLETE|SINGLETON|EMPTY_DECLARED|SPLIT_REQUIRED",
        "members": [{"template_id": "...", "chain_id": "...", "path": "..."}],
        "representatives": [{"template_id": "...", "chain_id": "...", "path": "..."}],
        "metadata": {"member_metadata": {"...": {}}}
      }]
    }

``gtalign --cls`` emits ``gtalignclusters.lst`` rather than a machine-readable
format.  GTalign 0.19's post-separator lines are treated as one cluster per
line; singleton lines therefore remain observable.  For fixtures and future
versions, the parser also accepts explicit ``cluster-id: members`` lines,
tab-separated declarations, and a JSON object containing ``clusters``.

The module has no scientific dependencies.  Similarity and biological
metadata are inputs to representative selection; this code never infers a
family, species, or role from structural similarity.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import shlex
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


DEFAULT_GTALIGN_PATH = "/home/rshadi25/.conda/envs/gtalign_env/bin/gtalign"
INDEX_SCHEMA_VERSION = "prism-template-cluster-index/v1"
_STRUCTURE_SUFFIXES = {".pdb", ".cif", ".mmcif", ".gz"}


class ClusterAccountingError(ValueError):
    """Raised when GTalign output cannot be reconciled to the panel exactly."""


class PanelHashMismatch(ValueError):
    """Raised when the supplied panel digest does not match observed bytes."""


@dataclass(frozen=True)
class ClusterDeclaration:
    cluster_id: str
    members: tuple[str, ...]
    declared_status: str = "DECLARED"


@dataclass(frozen=True)
class RepresentativeSelection:
    representatives: tuple[dict[str, Any], ...]
    status: str
    metadata: dict[str, Any]


def build_command(
    input_dir: str | Path,
    raw_output: str | Path,
    cache_dir: str | Path,
    threshold: float,
    coverage: float,
    algorithm: int,
    speed: int,
    *,
    gtalign_path: str | Path = DEFAULT_GTALIGN_PATH,
    pre_score: float | None = None,
    min_length: int | None = None,
    cpu_threads_reading: int | None = None,
    sort: int | None = None,
) -> list[str]:
    """Return the exact argv used for a GTalign clustering invocation."""

    command = [
        str(gtalign_path),
        f"--cls={input_dir}",
        "-o",
        str(raw_output),
        "-c",
        str(cache_dir),
        f"--cls-threshold={threshold:g}",
        f"--cls-coverage={coverage:g}",
        f"--cls-algorithm={int(algorithm)}",
        f"--speed={int(speed)}",
    ]
    if pre_score is not None:
        command.append(f"--pre-score={float(pre_score):g}")
    if min_length is not None:
        command.append(f"--dev-min-length={int(min_length)}")
    if cpu_threads_reading is not None:
        command.append(f"--cpu-threads-reading={int(cpu_threads_reading)}")
    if sort is not None:
        command.append(f"--sort={int(sort)}")
    return command


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compute_panel_sha256(panel: str | Path) -> str:
    """Hash one panel file or a deterministic directory inventory.

    A panel manifest is hashed as its exact bytes.  A directory is hashed as
    sorted ``relative-path + NUL + file-bytes`` records, so ordering from the
    filesystem cannot change the digest.
    """

    path = Path(panel)
    if path.is_file():
        return sha256_file(path)
    if not path.is_dir():
        raise FileNotFoundError(f"panel does not exist: {path}")
    digest = hashlib.sha256()
    for child in sorted(p for p in path.rglob("*") if p.is_file()):
        relative = child.relative_to(path).as_posix().encode("utf-8")
        digest.update(relative)
        digest.update(b"\0")
        with child.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def validate_panel_hash(panel: str | Path, expected_sha256: str | None) -> str:
    """Compute and, when supplied, validate a panel SHA256 digest."""

    actual = compute_panel_sha256(panel)
    if expected_sha256 and actual.lower() != expected_sha256.lower():
        raise PanelHashMismatch(f"panel SHA256 mismatch: expected {expected_sha256}, observed {actual}")
    return actual


def _parse_json_declarations(payload: Any) -> list[ClusterDeclaration]:
    if isinstance(payload, dict) and "clusters" in payload:
        payload = payload["clusters"]
    if isinstance(payload, dict):
        payload = [{"cluster_id": key, "members": value} for key, value in payload.items()]
    if not isinstance(payload, list):
        raise ValueError("GTalign cluster JSON must be a list or an object containing clusters")
    declarations: list[ClusterDeclaration] = []
    for index, item in enumerate(payload, 1):
        if isinstance(item, dict):
            cluster_id = str(item.get("cluster_id", item.get("id", f"cluster-{index:04d}")))
            raw_members = item.get("members", item.get("member_ids", []))
            status = str(item.get("status", "DECLARED"))
        elif isinstance(item, list):
            cluster_id, raw_members, status = f"cluster-{index:04d}", item, "DECLARED"
        else:
            raise ValueError(f"invalid cluster declaration at index {index}")
        if raw_members is None:
            raw_members = []
        if not isinstance(raw_members, list):
            raise ValueError(f"members for {cluster_id} must be a list")
        members = []
        for member in raw_members:
            if isinstance(member, str):
                members.append(member)
            elif isinstance(member, dict):
                token = member.get("member_id", member.get("id", member.get("token", member.get("path"))))
                if token is None:
                    raise ValueError(f"cluster member in {cluster_id} has no id/token/path")
                members.append(Path(str(token)).stem)
            else:
                raise ValueError(f"invalid member in {cluster_id}")
        if not members:
            status = "EMPTY_DECLARED"
        elif len(members) == 1:
            status = "SINGLETON"
        declarations.append(ClusterDeclaration(cluster_id, tuple(members), status))
    return declarations


def parse_gtalign_clusters(text: str) -> list[ClusterDeclaration]:
    """Parse GTalign's ``gtalignclusters.lst`` and fixture-compatible forms."""

    stripped = text.strip()
    if not stripped:
        return []
    if stripped[:1] in "[{":
        try:
            return _parse_json_declarations(json.loads(stripped))
        except json.JSONDecodeError:
            pass

    lines = text.splitlines()
    separators = [i for i, line in enumerate(lines) if re.fullmatch(r"\s*=+\s*", line)]
    body = lines[separators[-1] + 1 :] if separators else lines
    declarations: list[ClusterDeclaration] = []
    auto_index = 1
    for raw_line in body:
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if _looks_like_gtalign_header(line):
            continue
        if re.search(r"\b0\s+structure\(s\)", line, flags=re.IGNORECASE) or re.search(
            r"no\s+clusters?", line, flags=re.IGNORECASE
        ):
            continue

        cluster_id = f"cluster-{auto_index:04d}"
        members_text = line
        declared_status = "DECLARED"
        if "\t" in raw_line:
            cluster_id, members_text = raw_line.split("\t", 1)
            cluster_id = cluster_id.strip() or f"cluster-{auto_index:04d}"
        elif re.match(r"^[^\s:]+:\s*", line):
            cluster_id, members_text = line.split(":", 1)
            cluster_id = cluster_id.strip()
        if members_text.strip() in {"", "[]", "EMPTY", "NONE", "NULL", "-"}:
            members = ()
            declared_status = "EMPTY_DECLARED"
        else:
            members = tuple(members_text.split())
            if not members:
                declared_status = "EMPTY_DECLARED"
            elif len(members) == 1:
                declared_status = "SINGLETON"
        declarations.append(ClusterDeclaration(cluster_id, members, declared_status))
        auto_index += 1

    if not declarations and re.search(r"(?:0\s+structure|empty|no\s+clusters?)", stripped, re.IGNORECASE):
        return [ClusterDeclaration("cluster-0001", (), "EMPTY_DECLARED")]
    return declarations


def _looks_like_gtalign_header(line: str) -> bool:
    prefixes = (
        "gtalign ",
        "margelevicius,",
        "command line:",
        "clustered:",
        "devices:",
        "time elapsed",
        "total residues",
    )
    if line.lower().startswith(prefixes):
        return True
    return bool(re.match(r"^\d+\s+structure\(s\)$", line, flags=re.IGNORECASE))


def _member_id(member: Mapping[str, Any]) -> str:
    for key in ("member_id", "id", "token"):
        if member.get(key) is not None:
            return str(member[key])
    if member.get("path") is not None:
        return Path(str(member["path"])).stem
    raise ValueError("panel member requires member_id, token, id, or path")


def _public_member(member: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "template_id": str(member.get("template_id", _member_id(member))),
        "chain_id": str(member.get("chain_id", "")),
        "path": str(member.get("path", "")),
    }


def _metadata_for(member: Mapping[str, Any], metadata: Mapping[str, Any] | None) -> dict[str, Any]:
    embedded = member.get("metadata")
    if not metadata:
        return copy.deepcopy(embedded) if isinstance(embedded, dict) else {}
    identifier = _member_id(member)
    value = metadata.get(identifier)
    if value is None:
        value = metadata.get(str(member.get("template_id", identifier)))
    if value is None:
        value = embedded
    return copy.deepcopy(value) if isinstance(value, dict) else {}


def _similarity(similarity: Mapping[Any, Any] | None, left: str, right: str) -> float:
    if left == right:
        return 1.0
    if not similarity:
        return 0.0
    for key in ((left, right), (right, left)):
        value = similarity.get(key)
        if isinstance(value, (int, float)):
            return float(value)
    for outer, inner in ((left, right), (right, left)):
        row = similarity.get(outer)
        if isinstance(row, Mapping) and isinstance(row.get(inner), (int, float)):
            return float(row[inner])
    return 0.0


def _category_values(member_ids: Iterable[str], metadata: Mapping[str, dict[str, Any]], field: str) -> set[str]:
    values = set()
    for member_id in member_ids:
        value = metadata.get(member_id, {}).get(field)
        if value is not None and value != "":
            values.add(json.dumps(value, sort_keys=True, default=str))
    return values


def select_representatives(
    members: Sequence[Mapping[str, Any]],
    *,
    similarity: Mapping[Any, Any] | None = None,
    metadata: Mapping[str, Any] | None = None,
    member_metadata: Mapping[str, Any] | None = None,
    coverage_threshold: float = 0.7,
    max_representatives: int = 3,
) -> RepresentativeSelection:
    """Select native/medoid then farthest/diverse representatives.

    ``similarity`` is an optional precomputed score mapping.  Missing scores
    are treated as zero only when a mapping is supplied, making coverage
    failure explicit instead of silently assuming structural equivalence.
    """

    if max_representatives < 1:
        raise ValueError("max_representatives must be positive")
    # The shared index deliberately limits adaptive selection to one through
    # three representatives; larger requests would make SPLIT_REQUIRED
    # impossible to interpret consistently downstream.
    max_representatives = min(3, int(max_representatives))
    if metadata is not None and member_metadata is not None:
        raise ValueError("pass metadata or member_metadata, not both")
    metadata = metadata if metadata is not None else member_metadata
    ordered = sorted((dict(member) for member in members), key=_member_id)
    if not ordered:
        return RepresentativeSelection((), "EMPTY_DECLARED", {"member_metadata": {}})
    member_by_id = {_member_id(member): member for member in ordered}
    if len(member_by_id) != len(ordered):
        raise ClusterAccountingError("duplicate panel member ids in representative input")
    metadata_by_id = {_member_id(member): _metadata_for(member, metadata) for member in ordered}

    native_ids = [
        identifier
        for identifier in sorted(member_by_id)
        if metadata_by_id[identifier].get("is_native") is True
        or metadata_by_id[identifier].get("native") is True
        or str(metadata_by_id[identifier].get("role", "")).lower() in {"native", "reference"}
    ]
    if native_ids:
        first_id = native_ids[0]
        first_reason = "native"
    elif similarity:
        scores = {
            identifier: sum(_similarity(similarity, identifier, other) for other in member_by_id if other != identifier)
            for identifier in member_by_id
        }
        first_id = min(scores, key=lambda identifier: (-scores[identifier], identifier))
        first_reason = "medoid"
    else:
        first_id = min(member_by_id)
        first_reason = "deterministic_first"

    selected = [first_id]
    dimensions = ("family", "species", "role")

    def unresolved_diversity() -> dict[str, list[str]]:
        result: dict[str, list[str]] = {}
        for dimension in dimensions:
            all_values = _category_values(member_by_id, metadata_by_id, dimension)
            selected_values = _category_values(selected, metadata_by_id, dimension)
            missing = sorted(all_values - selected_values)
            if missing:
                result[dimension] = missing
        return result

    def uncovered_ids() -> list[str]:
        if not similarity:
            return []
        uncovered = []
        for identifier in member_by_id:
            best = max(_similarity(similarity, identifier, representative) for representative in selected)
            if best < float(coverage_threshold):
                uncovered.append(identifier)
        return uncovered

    while len(selected) < max_representatives:
        missing_diversity = unresolved_diversity()
        uncovered = uncovered_ids()
        if not missing_diversity and not uncovered:
            break
        candidates = [identifier for identifier in member_by_id if identifier not in selected]
        if not candidates:
            break

        def candidate_rank(identifier: str) -> tuple[int, int, float, str]:
            new_categories = sum(
                json.dumps(metadata_by_id[identifier].get(dimension), sort_keys=True, default=str)
                in missing_diversity.get(dimension, [])
                for dimension in dimensions
            )
            coverage_gain = sum(
                _similarity(similarity, member_id, identifier) >= float(coverage_threshold) for member_id in uncovered
            ) if similarity else 0
            distance = 1.0 - max(_similarity(similarity, identifier, representative) for representative in selected)
            return new_categories, coverage_gain, distance, identifier

        next_id = max(candidates, key=candidate_rank)
        selected.append(next_id)

    unresolved = unresolved_diversity()
    uncovered = uncovered_ids()
    status = "SPLIT_REQUIRED" if unresolved or uncovered else ("SINGLETON" if len(ordered) == 1 else "COMPLETE")
    selection_metadata = {
        "selection_strategy": "native-first-medoid-farthest-diversity",
        "first_representative_reason": first_reason,
        "coverage_threshold": float(coverage_threshold),
        "uncovered_member_ids": uncovered,
        "uncovered_diversity": unresolved,
        "member_metadata": copy.deepcopy(metadata_by_id),
    }
    return RepresentativeSelection(tuple(member_by_id[identifier] for identifier in selected), status, selection_metadata)


def _resolve_declarations(
    declarations: Sequence[ClusterDeclaration],
    panel_members: Sequence[Mapping[str, Any]],
    *,
    missing_member_policy: str = "error",
) -> list[ClusterDeclaration]:
    if missing_member_policy not in {"error", "singleton"}:
        raise ValueError("missing_member_policy must be 'error' or 'singleton'")
    aliases: dict[str, str] = {}
    canonical: dict[str, Mapping[str, Any]] = {}
    for member in panel_members:
        identifier = _member_id(member)
        if identifier in canonical:
            raise ClusterAccountingError(f"duplicate panel member id: {identifier}")
        canonical[identifier] = member
        raw_aliases = {identifier}
        raw_aliases.update(str(alias) for alias in member.get("aliases", []) or [])
        if member.get("path"):
            source_path = Path(str(member["path"]))
            raw_aliases.update({source_path.name, str(source_path), _strip_structure_suffix(source_path.name)})
        for alias in raw_aliases:
            previous = aliases.get(alias)
            if previous is not None and previous != identifier:
                raise ClusterAccountingError(f"duplicate panel alias: {alias}")
            aliases[alias] = identifier

    resolved: list[ClusterDeclaration] = []
    seen: set[str] = set()
    cluster_ids: set[str] = set()
    for declaration in declarations:
        if declaration.cluster_id in cluster_ids:
            raise ClusterAccountingError(f"duplicate cluster id: {declaration.cluster_id}")
        cluster_ids.add(declaration.cluster_id)
        members: list[str] = []
        for token in declaration.members:
            identifier = aliases.get(token)
            if identifier is None:
                identifier = aliases.get(_strip_structure_suffix(Path(token).name))
            if identifier is None:
                raise ClusterAccountingError(f"unaccounted GTalign member: {token}")
            if identifier in seen:
                raise ClusterAccountingError(f"duplicate GTalign member: {identifier}")
            seen.add(identifier)
            members.append(identifier)
        resolved.append(ClusterDeclaration(declaration.cluster_id, tuple(members), declaration.declared_status))

    missing = sorted(set(canonical) - seen)
    if missing:
        # An explicitly declared empty cluster is valid for an empty panel,
        # but never masks omitted members from a non-empty panel.
        if missing_member_policy == "error":
            raise ClusterAccountingError(f"missing GTalign members: {', '.join(missing)}")
        # GTalign can omit structures that its classifier cannot parse or
        # that fall below its minimum usable length.  Keep the shared index
        # complete with conservative singleton fallbacks.  The special
        # declaration status is retained in cluster metadata so downstream
        # routing never confuses these with GTalign-supported clusters.
        for identifier in missing:
            cluster_id = f"fallback-{identifier}"
            suffix = 2
            while cluster_id in cluster_ids:
                cluster_id = f"fallback-{identifier}-{suffix}"
                suffix += 1
            cluster_ids.add(cluster_id)
            resolved.append(ClusterDeclaration(cluster_id, (identifier,), "UNINDEXED_FALLBACK"))
    if not declarations and canonical:
        raise ClusterAccountingError("missing GTalign cluster declarations for non-empty panel")
    return resolved


def build_index(
    *,
    panel_members: Sequence[Mapping[str, Any]],
    declarations: Sequence[ClusterDeclaration],
    panel_sha256: str,
    tool_version: str,
    parameters: Mapping[str, Any],
    provenance: Mapping[str, Any] | None = None,
    similarity_by_cluster: Mapping[str, Mapping[Any, Any]] | None = None,
    metadata_by_member: Mapping[str, Any] | None = None,
    coverage_threshold: float | None = None,
    missing_member_policy: str = "error",
) -> dict[str, Any]:
    """Resolve panel accounting and emit the shared JSON index contract."""

    resolved = _resolve_declarations(
        declarations,
        panel_members,
        missing_member_policy=missing_member_policy,
    )
    member_lookup = {_member_id(member): member for member in panel_members}
    clusters: list[dict[str, Any]] = []
    threshold = float(coverage_threshold if coverage_threshold is not None else parameters.get("coverage", 0.7))
    for declaration in resolved:
        cluster_members = [member_lookup[identifier] for identifier in declaration.members]
        selection = select_representatives(
            cluster_members,
            similarity=(similarity_by_cluster or {}).get(declaration.cluster_id),
            metadata=metadata_by_member,
            coverage_threshold=threshold,
        )
        member_objects = [_public_member(member) for member in cluster_members]
        representative_objects = [_public_member(member) for member in selection.representatives]
        metadata = copy.deepcopy(selection.metadata)
        metadata["declared_status"] = declaration.declared_status
        metadata["index_membership_source"] = (
            "singleton_fallback"
            if declaration.declared_status == "UNINDEXED_FALLBACK"
            else "gtalign_output"
        )
        if declaration.declared_status == "UNINDEXED_FALLBACK":
            metadata["fallback_reason"] = "gtalign_output_omitted_member"
        clusters.append(
            {
                "cluster_id": declaration.cluster_id,
                "status": selection.status,
                "members": member_objects,
                "representatives": representative_objects,
                "metadata": metadata,
            }
        )

    output_parameters = copy.deepcopy(dict(parameters))
    fallback_ids = sorted(
        declaration.members[0]
        for declaration in resolved
        if declaration.declared_status == "UNINDEXED_FALLBACK" and declaration.members
    )
    output_parameters["missing_member_policy"] = missing_member_policy
    output_parameters["missing_member_fallback_count"] = len(fallback_ids)
    output_parameters["missing_member_fallback_ids_sha256"] = hashlib.sha256(
        ("\n".join(fallback_ids) + "\n").encode("utf-8")
    ).hexdigest()
    return {
        "schema_version": INDEX_SCHEMA_VERSION,
        "tool": "GTalign",
        "tool_version": tool_version,
        "panel_count": len(panel_members),
        "panel_sha256": panel_sha256,
        "parameters": output_parameters,
        "clusters": clusters,
        **({"provenance": copy.deepcopy(dict(provenance))} if provenance else {}),
    }


def _strip_structure_suffix(name: str) -> str:
    stem = name
    if stem.lower().endswith(".gz"):
        stem = stem[:-3]
    for suffix in (".pdb", ".cif", ".mmcif"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem


def _chain_ids(path: Path) -> list[str]:
    chains: set[str] = set()
    try:
        with path.open("r", errors="replace") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                    chain = line[21].strip()
                    if chain:
                        chains.add(chain)
    except OSError:
        return []
    return sorted(chains)


def discover_panel_members(input_dir: str | Path) -> list[dict[str, Any]]:
    """Create one panel member per structure file and aliases for GTalign IDs."""

    root = Path(input_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"input directory does not exist: {root}")
    files = sorted(
        path for path in root.rglob("*") if path.is_file() and path.suffix.lower() in _STRUCTURE_SUFFIXES
    )
    members: list[dict[str, Any]] = []
    for path in files:
        stem = _strip_structure_suffix(path.name)
        chains = _chain_ids(path)
        suffix_match = re.match(r"^(.+)_([A-Za-z0-9])$", stem)
        # GTalign's cluster-list token identifies one chain/model even when
        # the source PDB contains several chains (for example ``1wte_A_1``).
        # Keep the first declared chain as the stable panel identity; aliases
        # below still account for every token GTalign can emit.
        chain_id = suffix_match.group(2) if suffix_match else (chains[0] if chains else "")
        template_id = suffix_match.group(1) if suffix_match else stem
        aliases = {stem}
        for chain in chains:
            aliases.add(f"{stem}_{chain}")
            aliases.add(f"{stem}_{chain}_1")
        if suffix_match:
            aliases.add(f"{stem}_1")
        members.append(
            {
                "member_id": stem,
                "template_id": template_id,
                "chain_id": chain_id,
                "path": str(path.resolve()),
                "aliases": sorted(aliases - {stem}),
            }
        )
    return members


def _load_panel_manifest(panel: Path, input_dir: Path) -> list[dict[str, Any]]:
    if panel.suffix.lower() == ".json":
        payload = json.loads(panel.read_text())
        payload = payload.get("members", payload) if isinstance(payload, dict) else payload
        if not isinstance(payload, list):
            raise ValueError("panel JSON must be a list or contain members")
        rows = payload
    else:
        rows = []
        for line in panel.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if line.lstrip().startswith("{"):
                rows.append(json.loads(line))
                continue
            fields = line.split("\t")
            if len(fields) >= 4:
                rows.append({"member_id": fields[0], "template_id": fields[1], "chain_id": fields[2], "path": fields[3]})
            elif len(fields) == 3:
                rows.append({"member_id": fields[0], "chain_id": fields[1], "path": fields[2]})
            elif len(fields) == 2:
                rows.append({"member_id": fields[0], "path": fields[1]})
            else:
                rows.append({"path": fields[0]})
    members = []
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("panel entries must be objects")
        member = dict(row)
        source = Path(str(member.get("path", "")))
        if not source.is_absolute():
            source = (panel.parent / source).resolve()
        if not source.exists():
            candidate = (input_dir / source.name).resolve()
            if candidate.exists():
                source = candidate
        if not source.exists():
            raise FileNotFoundError(f"panel member path does not exist: {source}")
        member["path"] = str(source)
        member.setdefault("member_id", _strip_structure_suffix(source.name))
        member.setdefault("template_id", member["member_id"])
        member.setdefault("chain_id", "")
        # GTalign 0.19 emits the source filename stem followed by the chain
        # identifier.  For a chain-qualified filename such as
        # ``1lyqAB_A.pdb`` this becomes ``1lyqAB_A_A``; make that repeated
        # suffix an explicit alias so JSON panel manifests remain compatible
        # with the real cluster-list output.
        aliases = {str(alias) for alias in member.get("aliases", []) or []}
        stem = _strip_structure_suffix(source.name)
        aliases.update({source.name, stem})
        chain_id = str(member.get("chain_id", "")).strip()
        if chain_id:
            aliases.update({f"{stem}_{chain_id}", f"{stem}_{chain_id}_1"})
            # GTalign may normalize the staging suffix ``_int`` away before
            # appending its repeated chain/model token.  For example, the
            # staged path ``2axtAI_I_int.pdb`` can be reported as
            # ``2axtAI_I_I``.  Keep this alias explicit and deterministic;
            # the path/member identity remains the panel's canonical key.
            normalized_stem = re.sub(r"_int$", "", stem, flags=re.IGNORECASE)
            if normalized_stem != stem:
                aliases.update(
                    {
                        f"{normalized_stem}_{chain_id}",
                        f"{normalized_stem}_{chain_id}_1",
                        f"{normalized_stem}_{chain_id}_{chain_id}",
                        f"{normalized_stem}_{chain_id}_{chain_id}_1",
                    }
                )
        member["aliases"] = sorted(aliases)
        members.append(member)
    return members


def load_panel_members(input_dir: str | Path, panel: str | Path | None) -> tuple[list[dict[str, Any]], str]:
    root = Path(input_dir)
    panel_path = Path(panel) if panel else root
    if panel and panel_path.is_file():
        return _load_panel_manifest(panel_path, root), compute_panel_sha256(panel_path)
    members = discover_panel_members(root)
    return members, compute_panel_sha256(panel_path)


def _tool_info(gtalign_path: Path) -> tuple[str, str | None, str]:
    if not gtalign_path.exists():
        return "unknown", None, ""
    executable_sha256 = sha256_file(gtalign_path)
    result = subprocess.run([str(gtalign_path), "-h"], text=True, capture_output=True, check=False)
    help_text = (result.stdout or "") + (result.stderr or "")
    match = re.search(r"gtalign\s+([0-9]+(?:\.[0-9]+)+)", help_text, flags=re.IGNORECASE)
    return (match.group(1) if match else "unknown"), executable_sha256, help_text


def _find_cluster_output(raw_output: Path) -> Path:
    preferred = [raw_output / "gtalignclusters.lst", raw_output / "gtalignclusters.txt", raw_output / "clusters.json"]
    for path in preferred:
        if path.is_file():
            return path
    candidates = sorted(path for path in raw_output.glob("*") if path.is_file() and path.suffix.lower() in {".lst", ".txt", ".json"})
    if candidates:
        return candidates[0]
    raise FileNotFoundError(f"GTalign cluster output not found in {raw_output}")


def run_clustering(
    *,
    input_dir: str | Path,
    output_dir: str | Path,
    panel: str | Path | None = None,
    panel_sha256: str | None = None,
    gtalign_path: str | Path = DEFAULT_GTALIGN_PATH,
    cache_dir: str | Path | None = None,
    threshold: float = 0.5,
    coverage: float = 0.7,
    algorithm: int = 0,
    speed: int = 13,
    pre_score: float | None = None,
    min_length: int | None = None,
    cpu_threads_reading: int | None = None,
    sort: int | None = None,
    missing_member_policy: str = "error",
    dry_run: bool = False,
) -> Path:
    """Run bounded GTalign clustering and write ``gtalign_cluster_index.json``."""

    input_path = Path(input_dir).resolve()
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    raw_output = output_path / "gtalign_raw"
    cache_path = Path(cache_dir).resolve() if cache_dir else output_path / "gtalign_cache"
    panel_members, observed_panel_hash = load_panel_members(input_path, panel)
    if panel_sha256 and observed_panel_hash.lower() != panel_sha256.lower():
        raise PanelHashMismatch(f"panel SHA256 mismatch: expected {panel_sha256}, observed {observed_panel_hash}")

    command = build_command(
        input_path,
        raw_output,
        cache_path,
        threshold,
        coverage,
        algorithm,
        speed,
        gtalign_path=gtalign_path,
        pre_score=pre_score,
        min_length=min_length,
        cpu_threads_reading=cpu_threads_reading,
        sort=sort,
    )
    started_at = _now()
    tool_version, executable_sha256, help_text = _tool_info(Path(gtalign_path))
    provenance: dict[str, Any] = {
        "status": "DRY_RUN" if dry_run else "RUNNING",
        "command": command,
        "command_text": shlex.join(command),
        "executable": str(Path(gtalign_path)),
        "executable_sha256": executable_sha256,
        "tool_help": help_text,
        "panel_count": len(panel_members),
        "panel": str(panel) if panel else None,
        "input_dir": str(input_path),
        "raw_output": str(raw_output),
        "cache_dir": str(cache_path),
        "started_at": started_at,
    }
    parameters = {
        "threshold": float(threshold),
        "coverage": float(coverage),
        "algorithm": int(algorithm),
        "speed": int(speed),
        "pre_score": pre_score,
        "min_length": min_length,
        "cpu_threads_reading": cpu_threads_reading,
        "sort": sort,
        "missing_member_policy": missing_member_policy,
        "input_dir": str(input_path),
    }
    if dry_run:
        provenance.update({"status": "DRY_RUN", "finished_at": _now()})
        index = build_index(
            panel_members=[],
            declarations=[],
            panel_sha256=observed_panel_hash,
            tool_version=tool_version,
            parameters=parameters,
            provenance=provenance,
        )
    else:
        gtalign_executable = Path(gtalign_path)
        if not gtalign_executable.is_file():
            raise FileNotFoundError(f"GTalign executable does not exist: {gtalign_executable}")
        if raw_output.exists() and any(raw_output.iterdir()):
            raise RuntimeError(f"refusing to reuse non-empty GTalign output directory: {raw_output}")
        raw_output.mkdir(parents=True, exist_ok=True)
        cache_path.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(command, text=True, capture_output=True, check=False)
        provenance["returncode"] = result.returncode
        provenance["stdout"] = result.stdout or ""
        provenance["stderr"] = result.stderr or ""
        if result.returncode != 0:
            provenance.update({"status": "FAILED", "finished_at": _now()})
            raise RuntimeError(f"GTalign failed with exit code {result.returncode}: {(result.stderr or result.stdout)[:2000]}")
        output_file = _find_cluster_output(raw_output)
        declarations = parse_gtalign_clusters(output_file.read_text(errors="replace"))
        resolved = _resolve_declarations(
            declarations,
            panel_members,
            missing_member_policy=missing_member_policy,
        )
        provenance.update({"status": "COMPLETED", "finished_at": _now(), "cluster_output": str(output_file)})
        index = build_index(
            panel_members=panel_members,
            declarations=resolved,
            panel_sha256=observed_panel_hash,
            tool_version=tool_version,
            parameters=parameters,
            provenance=provenance,
            coverage_threshold=coverage,
            missing_member_policy=missing_member_policy,
        )

    index_path = output_path / "gtalign_cluster_index.json"
    index_path.write_text(json.dumps(index, indent=2, sort_keys=True) + "\n")
    return index_path


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--panel", default=None)
    parser.add_argument("--panel-sha256", default=None)
    parser.add_argument("--gtalign-path", default=DEFAULT_GTALIGN_PATH)
    parser.add_argument("--cache-dir", default=None)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--coverage", type=float, default=0.7)
    parser.add_argument("--algorithm", type=int, choices=(0, 1), default=0)
    parser.add_argument("--speed", type=int, default=13)
    parser.add_argument("--pre-score", type=float, default=None)
    parser.add_argument("--min-length", type=int, default=None)
    parser.add_argument("--cpu-threads-reading", type=int, default=None)
    parser.add_argument("--sort", type=int, choices=range(9), default=None)
    parser.add_argument(
        "--missing-member-policy",
        choices=("error", "singleton"),
        default="error",
        help="handle GTalign-omitted panel members (default: fail closed)",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)
    try:
        output = run_clustering(**vars(args))
    except (ClusterAccountingError, PanelHashMismatch, FileNotFoundError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    else:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
