"""Auditable, read-only diagnostics for orientation comparison runs.

The helpers in this module deliberately do not change pipeline decisions.  They
turn a run directory into provenance, alignment, gate, clash, refinement,
ranking, and US-align diagnostic tables.  Missing evidence is represented as
``None`` or an explicit status instead of being treated as a pass.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable, Mapping

from .pdb_download import normalize_target_id
from .template_filtering import evaluate_protocol_candidate, load_filter_assets


GATE_NAMES = (
    "availability",
    "score_threshold",
    "match_count",
    "coverage",
    "hotspots",
    "complementary_contacts",
    "transform_materialization",
    "clash_filter",
)
DEFAULT_CLASH_DISTANCE_GRID = (2.5, 3.0, 3.5)
DEFAULT_CLASH_EVENT_GRID = (0, 1, 5, 10, 25, 50, 100)


def sha256_file(path: str | Path) -> str | None:
    """Return a file hash, or ``None`` when the artifact is unavailable."""

    path = Path(path)
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_jsonl(path: str | Path) -> list[dict[str, Any]]:
    """Load JSONL while preserving malformed-line evidence."""

    path = Path(path)
    if not path.is_file():
        return []
    records: list[dict[str, Any]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            value = json.loads(line)
            records.append(value if isinstance(value, dict) else {
                "_line": line_number,
                "_parse_error": "JSON value is not an object",
            })
        except json.JSONDecodeError as exc:
            records.append({
                "_line": line_number,
                "_parse_error": f"{type(exc).__name__}: {exc}",
            })
    return records


def read_input_pairs(path: str | Path) -> list[dict[str, str]]:
    """Read receptor/ligand selectors without silently dropping malformed rows."""

    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = {str(name).casefold(): name for name in (reader.fieldnames or [])}
        receptor = fields.get("receptor")
        ligand = fields.get("ligand")
        if not receptor or not ligand:
            raise ValueError("input CSV must contain Receptor and Ligand columns")
        rows = []
        for row_number, row in enumerate(reader, 2):
            left, right = str(row.get(receptor, "")).strip(), str(row.get(ligand, "")).strip()
            rows.append({
                "row_number": row_number,
                "query_left": left,
                "query_right": right,
                "valid": bool(left and right),
                "error": None if left and right else "missing receptor or ligand",
            })
        return rows


def orientation_requirements(
    query_left: str, query_right: str, template: str, orientation: str
) -> dict[str, str]:
    """Return the query/template-chain assignment for one orientation."""

    if len(template) != 6:
        raise ValueError(f"expected six-character template ID, got {template!r}")
    chain1, chain2 = template[4], template[5]
    if orientation == "o1":
        return {
            "query_left": query_left,
            "query_right": query_right,
            "chain_left": chain1,
            "chain_right": chain2,
        }
    if orientation == "o2":
        return {
            "query_left": query_left,
            "query_right": query_right,
            "chain_left": chain2,
            "chain_right": chain1,
        }
    raise ValueError(f"orientation must be o1 or o2, got {orientation!r}")


def _query_names(query: str) -> set[str]:
    names = {str(query)}
    try:
        names.add(normalize_target_id(query))
    except (TypeError, ValueError):
        pass
    return {name.casefold() for name in names}


def find_alignment_file(
    run_root: str | Path, query: str, template: str, chain: str
) -> Path | None:
    """Find a PRISM alignment JSON using raw and canonical query spellings."""

    root = Path(run_root)
    suffix = f"_{template}_{chain}.json".casefold()
    names = _query_names(query)
    for path in sorted((root / "processed").rglob("*.json")) if (root / "processed").exists() else []:
        if not path.name.casefold().endswith(suffix):
            continue
        prefix = path.name[:-len(suffix)].casefold()
        if prefix in names:
            return path
    return None


def _load_json(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.is_file():
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def _ca_residue_count(path: Path | None) -> int | None:
    if path is None or not path.is_file():
        return None
    residues: set[tuple[str, str, str]] = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            residues.add((line[21:22].strip(), line[22:26].strip(), line[26:27].strip()))
    return len(residues)


def template_interface_size(run_root: str | Path, template: str, chain: str) -> int | None:
    return _ca_residue_count(Path(run_root) / "templates" / "interfaces" / f"{template}_{chain}_int.pdb")


def _coverage(match_count: Any, size: int | None) -> float | None:
    try:
        if size is None or size <= 0:
            return None
        return float(match_count or 0) / float(size) * 100.0
    except (TypeError, ValueError):
        return None


def _alignment_contract_status(payload: Mapping[str, Any] | None) -> str:
    """Classify whether an alignment record is safe for gate replay."""

    if not isinstance(payload, Mapping):
        return "missing"
    if payload.get("return_code") is None:
        return "return_code_missing"
    if payload.get("return_code") != 0:
        return "return_code_nonzero"
    if "status" in payload and payload.get("status") != "success":
        return f"status_{payload.get('status')}"
    if not payload.get("raw_output_sha256"):
        return "raw_output_hash_missing"
    required = {"match_count", "match_dict", "aligner"}
    aligner = str(payload.get("aligner", "")).casefold()
    required.add("rmsd" if aligner == "multiprot" else "tm_score")
    missing = sorted(field for field in required if field not in payload)
    return "fields_missing:" + ",".join(missing) if missing else "valid"


def alignment_inventory(
    run_root: str | Path,
    inputs_csv: str | Path,
    templates: Iterable[str],
    *,
    arm: str,
    orientations: Iterable[str] = ("o1", "o2"),
    audit_path: str | Path | None = None,
) -> list[dict[str, Any]]:
    """Build one row per input/template/orientation alignment pair."""

    root = Path(run_root)
    audit_index = {}
    for record in load_jsonl(audit_path) if audit_path else []:
        key = (
            record.get("query_left"), record.get("query_right"),
            record.get("template"), record.get("orientation"),
        )
        audit_index[key] = record

    rows: list[dict[str, Any]] = []
    for pair in read_input_pairs(inputs_csv):
        for template in templates:
            for orientation in orientations:
                assignment = orientation_requirements(
                    pair["query_left"], pair["query_right"], template, orientation
                )
                left_path = find_alignment_file(
                    root, pair["query_left"], template, assignment["chain_left"]
                )
                right_path = find_alignment_file(
                    root, pair["query_right"], template, assignment["chain_right"]
                )
                left = _load_json(left_path)
                right = _load_json(right_path)
                left_size = template_interface_size(root, template, assignment["chain_left"])
                right_size = template_interface_size(root, template, assignment["chain_right"])
                audit = audit_index.get((
                    pair["query_left"], pair["query_right"], template, orientation
                ), {})
                rows.append({
                    "arm": arm,
                    "orientation": orientation,
                    "query_left": pair["query_left"],
                    "query_right": pair["query_right"],
                    "template": template,
                    "chain_left": assignment["chain_left"],
                    "chain_right": assignment["chain_right"],
                    "left_alignment_path": str(left_path) if left_path else None,
                    "right_alignment_path": str(right_path) if right_path else None,
                    "partner_left_available": left_path is not None and left is not None,
                    "partner_right_available": right_path is not None and right is not None,
                    "partner_available": left_path is not None and right_path is not None and left is not None and right is not None,
                    "alignment_contract_status_left": _alignment_contract_status(left),
                    "alignment_contract_status_right": _alignment_contract_status(right),
                    "alignment_evidence_complete": (
                        _alignment_contract_status(left) == "valid"
                        and _alignment_contract_status(right) == "valid"
                    ),
                    "aligner_left": left.get("aligner") if left else None,
                    "aligner_right": right.get("aligner") if right else None,
                    "score_gate_contract_left": (left or {}).get("score_gate_contract", (left or {}).get("tm_score_contract", "tm_score")),
                    "score_gate_contract_right": (right or {}).get("score_gate_contract", (right or {}).get("tm_score_contract", "tm_score")),
                    "return_code_left": (left or {}).get("return_code"),
                    "return_code_right": (right or {}).get("return_code"),
                    "match_count_left": (left or {}).get("match_count"),
                    "match_count_right": (right or {}).get("match_count"),
                    "coverage_left": _coverage((left or {}).get("match_count"), left_size),
                    "coverage_right": _coverage((right or {}).get("match_count"), right_size),
                    "template_size_left": left_size,
                    "template_size_right": right_size,
                    "tm_score_left": (left or {}).get("tm_score"),
                    "tm_score_right": (right or {}).get("tm_score"),
                    "native_multiprot_score_left": (left or {}).get("rmsd"),
                    "native_multiprot_score_right": (right or {}).get("rmsd"),
                    "raw_output_sha256_left": (left or {}).get("raw_output_sha256"),
                    "raw_output_sha256_right": (right or {}).get("raw_output_sha256"),
                    "alignment_json_sha256_left": sha256_file(left_path) if left_path else None,
                    "alignment_json_sha256_right": sha256_file(right_path) if right_path else None,
                    "raw_output_hash_complete": bool(
                        (left or {}).get("raw_output_sha256")
                        and (right or {}).get("raw_output_sha256")
                    ),
                    "audit_status": audit.get("status", "not_recorded"),
                    "audit_error_reason": audit.get("error_reason"),
                    "left_payload": left,
                    "right_payload": right,
                })
    return rows


def _hash_inventory(
    paths: Iterable[Path], root: Path, *, origin_by_path: Mapping[Path, str] | None = None,
    include_resolved_path: bool = False,
) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(set(paths)):
        try:
            relative = str(path.relative_to(root))
        except ValueError:
            relative = str(path)
        rows.append({
            "path": relative,
            "exists": path.is_file(),
            "bytes": path.stat().st_size if path.is_file() else None,
            "sha256": sha256_file(path),
            "origin": (origin_by_path or {}).get(path),
            **({"resolved_path": str(path.resolve())} if include_resolved_path and path.exists() else {}),
        })
    return rows


def _reference_hash_index(roots: Iterable[Path]) -> dict[str, set[str]]:
    """Index reference bytes so copied or renamed assets do not evade provenance."""

    index: dict[str, set[str]] = {}
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if path.is_file():
                digest = sha256_file(path)
                if digest:
                    index.setdefault(digest, set()).add(str(root.resolve()))
    return index


def _asset_origins(
    paths: Iterable[Path], *, current_roots: Iterable[Path], legacy_roots: Iterable[Path]
) -> dict[Path, str]:
    current_roots, legacy_roots = tuple(current_roots), tuple(legacy_roots)
    references = _reference_hash_index((*current_roots, *legacy_roots))
    origins: dict[Path, str] = {}
    for path in paths:
        digest = sha256_file(path)
        resolved = path.resolve() if path.exists() else path
        in_current = any(resolved == root.resolve() or root.resolve() in resolved.parents for root in current_roots if root.exists())
        in_legacy = any(resolved == root.resolve() or root.resolve() in resolved.parents for root in legacy_roots if root.exists())
        known_roots = references.get(digest or "", set())
        known_current = any(
            root in known_roots
            for root in (str(r.resolve()) for r in current_roots if r.exists())
        )
        known_legacy = any(
            root in known_roots
            for root in (str(r.resolve()) for r in legacy_roots if r.exists())
        )
        if in_current or (known_current and not known_legacy):
            origins[path] = "current"
        elif in_legacy or (known_legacy and not known_current):
            origins[path] = "legacy"
        elif known_current and known_legacy:
            origins[path] = "equivalent_known"
        elif digest:
            origins[path] = "unknown_unmatched"
        else:
            origins[path] = "missing"
    return origins


def asset_provenance_manifest(
    run_root: str | Path,
    inputs_csv: str | Path,
    templates: Iterable[str],
    *,
    source_root: str | Path,
    surface_backend: str,
    filter_mode: str,
    filter_asset_root: str | Path | None = None,
) -> dict[str, Any]:
    """Capture input, source, interface, target, surface, and mix evidence."""

    root, source = Path(run_root), Path(source_root)
    pairs = read_input_pairs(inputs_csv)
    queries = [q for pair in pairs for q in (pair["query_left"], pair["query_right"])]
    interface_paths = [
        root / "templates" / "interfaces" / f"{template}_{chain}_int.pdb"
        for template in templates for chain in template[4:]
    ]
    interface_list_paths = [
        root / "templates" / "interfaces_lists" / f"{template}.json"
        for template in templates
    ]
    target_paths = []
    for query in queries:
        try:
            canonical = normalize_target_id(query)
        except (TypeError, ValueError):
            canonical = query
        target_paths.extend([
            root / "processed" / "pdbs" / f"{canonical}.pdb",
            root / "processed" / "surface_extraction" / f"{canonical}.asa.pdb",
        ])
    source_paths = [source / "prism.py", source / "src" / "transformation.py", source / "src" / "alignment.py", source / "src" / "alignment_gtalign.py", source / "src" / "alignment_multiprot.py"]
    source_paths.extend([
        source / "src" / "surface_extract.py",
        source / "src" / "transformation_config.py",
        source / "src" / "stepwise_analysis.py",
    ])
    surface_paths = list((root / "processed" / "surface_extraction").glob("*.asa.pdb")) if (root / "processed" / "surface_extraction").exists() else []
    alignment_paths = list((root / "processed").glob("alignment_*/*/*.json")) if (root / "processed").exists() else []
    current_root = root / "templates"
    legacy_root = root / "template"
    reference_current_root = source / "templates"
    reference_legacy_root = source / "working_version" / "Multiprot-new" / "prism-fiberdock-cli" / "template"
    configured_filter = Path(filter_asset_root).resolve() if filter_asset_root else None
    filter_paths = []
    if filter_mode == "published_protocol" and configured_filter:
        filter_paths = [
            configured_filter / "hotspots" / f"{template}.json" for template in templates
        ] + [
            configured_filter / "contacts" / f"{template}.json" for template in templates
        ]
    consumed_paths = list(interface_paths) + interface_list_paths + filter_paths
    origins = _asset_origins(
        consumed_paths,
        current_roots=(reference_current_root,),
        legacy_roots=(reference_legacy_root,),
    )
    consumed_origins = {origins[path] for path in consumed_paths if path.is_file()}
    source_origin = {
        path: "current_source" for path in source_paths
    }
    missing_consumed = any(not path.is_file() for path in consumed_paths)
    if "current" in consumed_origins and "legacy" in consumed_origins:
        mix_status = "mixed"
    elif ("unknown_unmatched" in consumed_origins or missing_consumed) and consumed_origins & {"current", "legacy"}:
        mix_status = "unknown_unresolved"
    elif not missing_consumed and (consumed_origins == {"current"} or consumed_origins == {"equivalent_known"}):
        mix_status = "current_only"
    elif not missing_consumed and consumed_origins == {"legacy"}:
        mix_status = "legacy_only"
    elif not consumed_origins:
        mix_status = "unknown"
    else:
        mix_status = "unknown_unresolved"
    return {
        "run_root": str(root.resolve()),
        "input_csv": _hash_inventory([Path(inputs_csv)], source),
        "source_assets": _hash_inventory(source_paths, source, origin_by_path=source_origin, include_resolved_path=True),
        "template_interface_assets": _hash_inventory(interface_paths, root, origin_by_path=origins, include_resolved_path=True),
        "template_interface_list_assets": _hash_inventory(interface_list_paths, root, origin_by_path=origins, include_resolved_path=True),
        "filter_assets": _hash_inventory(filter_paths, root, origin_by_path=origins, include_resolved_path=True),
        "target_and_surface_inputs": _hash_inventory(target_paths, root),
        "surface_outputs": _hash_inventory(surface_paths, root),
        "alignment_outputs": _hash_inventory(alignment_paths, root),
        "surface_backend": surface_backend,
        "filter_mode": filter_mode,
        "filter_asset_root": str(configured_filter) if configured_filter else None,
        "asset_mix_status": mix_status,
        "current_asset_root_present": current_root.exists(),
        "legacy_asset_root_present": legacy_root.exists(),
        "consumed_asset_origins": sorted(consumed_origins),
        "input_rows": pairs,
    }


def _atom_records(path: str | Path | None, *, ca_only: bool = True) -> list[dict[str, Any]]:
    if path is None or not Path(path).is_file():
        return []
    records = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if ca_only and line[12:16].strip() != "CA":
            continue
        try:
            records.append({
                "chain": line[21:22].strip() or "_",
                "residue": f"{line[17:20].strip()}{line[22:26].strip()}{line[26:27].strip()}",
                "atom": line[12:16].strip(),
                "xyz": (
                    float(line[30:38]), float(line[38:46]), float(line[46:54])
                ),
            })
        except (ValueError, IndexError):
            continue
    return records


def _distance(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def clash_diagnostics(
    left_path: str | Path,
    right_path: str | Path,
    *,
    clash_distance: float = 3.0,
    distance_grid: Iterable[float] = DEFAULT_CLASH_DISTANCE_GRID,
    event_grid: Iterable[int] = DEFAULT_CLASH_EVENT_GRID,
) -> dict[str, Any]:
    """Calculate CA clash counts, locations, distance distribution, and grid replay."""

    left, right = _atom_records(left_path), _atom_records(right_path)
    distances = []
    events = []
    max_distance = max([clash_distance, *distance_grid])
    for left_atom in left:
        for right_atom in right:
            distance = _distance(left_atom["xyz"], right_atom["xyz"])
            if distance < max_distance:
                distances.append(distance)
            if distance < clash_distance:
                events.append({
                    "left_chain": left_atom["chain"],
                    "left_residue": left_atom["residue"],
                    "right_chain": right_atom["chain"],
                    "right_residue": right_atom["residue"],
                    "distance": distance,
                })
    sorted_grid = sorted(float(value) for value in distance_grid)
    grid_rows = []
    for distance_limit in sorted_grid:
        observed = sum(value < distance_limit for value in distances)
        for maximum in event_grid:
            grid_rows.append({
                "distance": distance_limit,
                "max_events": int(maximum),
                "observed_events": observed,
                "passes_production_style_gate": observed < int(maximum),
            })
    left_residues = {event["left_residue"] for event in events}
    right_residues = {event["right_residue"] for event in events}
    left_total = len({atom["residue"] for atom in left})
    right_total = len({atom["residue"] for atom in right})
    left_fraction = len(left_residues) / left_total if left_total else None
    right_fraction = len(right_residues) / right_total if right_total else None
    if not events:
        localization = "no_clashes"
    elif left_fraction is not None and right_fraction is not None and max(left_fraction, right_fraction) <= 0.25:
        localization = "localized_heuristic"
    else:
        localization = "distributed_heuristic"
    return {
        "left_path": str(left_path),
        "right_path": str(right_path),
        "left_ca_count": len(left),
        "right_ca_count": len(right),
        "clash_distance": float(clash_distance),
        "clash_count": len(events),
        "minimum_cross_partner_distance": min(distances) if distances else None,
        "left_clashing_residue_count": len(left_residues),
        "right_clashing_residue_count": len(right_residues),
        "left_clashing_residue_fraction": left_fraction,
        "right_clashing_residue_fraction": right_fraction,
        "localization_heuristic": localization,
        "event_residues": events,
        "grid": grid_rows,
        "all_atom_vdw_status": "not_available_without_radii",
    }


def _transformed_pair_paths(run_root: str | Path, row: Mapping[str, Any]) -> tuple[Path | None, Path | None]:
    directory = Path(run_root) / "processed" / "transformation"
    if not directory.exists():
        return None, None
    prefix = f"{row['template']}_{row['query_left']}_{row['query_right']}_{row['orientation']}"
    left = sorted(directory.glob(f"{prefix}*_L.pdb"))
    right = sorted(directory.glob(f"{prefix}*_R.pdb"))
    if not left or not right:
        return None, None
    right_by_stem = {path.name[:-6]: path for path in right}
    for left_path in left:
        match = right_by_stem.get(left_path.name[:-6])
        if match:
            return left_path, match
    return None, None


def _protocol_contacts(filter_assets: Mapping[str, Any], left_chain: str, right_chain: str) -> list[Any]:
    oriented = []
    for contact in filter_assets.get("contacts", []):
        if isinstance(contact, Mapping):
            left, right = contact.get("left", contact.get("a")), contact.get("right", contact.get("b"))
        elif isinstance(contact, (list, tuple)) and len(contact) >= 2:
            left, right = contact[0], contact[1]
        else:
            continue
        left_text, right_text = str(left), str(right)
        if left_text.startswith(f"{left_chain}.") and right_text.startswith(f"{right_chain}."):
            oriented.append([left, right])
        elif left_text.startswith(f"{right_chain}.") and right_text.startswith(f"{left_chain}."):
            oriented.append([right, left])
    return oriented


def _audit_index(path: str | Path | None) -> dict[tuple[Any, ...], dict[str, Any]]:
    records = load_jsonl(path) if path else []
    return {
        (record.get("query_left"), record.get("query_right"), record.get("template"), record.get("orientation")): record
        for record in records
    }


def gate_ledger(
    run_root: str | Path,
    alignment_rows: Iterable[Mapping[str, Any]],
    *,
    thresholds: Mapping[str, Any],
    filter_mode: str = "geometry_only_experimental",
    filter_asset_root: str | Path | None = None,
    audit_path: str | Path | None = None,
    compute_clash_diagnostics: bool = True,
) -> list[dict[str, Any]]:
    """Evaluate independent gates and cumulative replay without changing outputs."""

    audits = _audit_index(audit_path)
    rows = []
    for source in alignment_rows:
        row = dict(source)
        left, right = row.get("left_payload"), row.get("right_payload")
        available = bool(row.get("partner_available"))
        gates: dict[str, bool | None] = {name: None for name in GATE_NAMES}
        gate_status: dict[str, str] = {}
        gates["availability"] = available
        if available:
            multiprot = str(row.get("aligner_left", "")).casefold() == "multiprot" or str(row.get("aligner_right", "")).casefold() == "multiprot"
            comparable = str(thresholds.get("alignment_gate_mode", "native")).strip().lower() == "common_match_coverage"
            if comparable:
                # TM-score remains a diagnostic field, but no provider-specific
                # score gate is used in the common match/coverage contract.
                gates["score_threshold"] = True
                gate_status["score_threshold"] = "not_applied_common_match_coverage"
                minimum_matches = int(thresholds["minimum_residue_match_count"])
                minimum_coverage = float(thresholds["minimum_residue_match_percentage"])
            elif multiprot:
                gates["score_threshold"] = None
                gate_status["score_threshold"] = "not_applicable_native_multiprot_contract"
                minimum_matches = int(thresholds["multiprot_minimum_residue_match_count"])
                minimum_coverage = float(thresholds["multiprot_minimum_residue_match_percentage"])
            else:
                gates["score_threshold"] = all(float(payload.get("tm_score", 0.0) or 0.0) >= float(thresholds["tm_score_threshold"]) for payload in (left, right))
                gate_status["score_threshold"] = "tm_score"
                minimum_matches = int(thresholds["minimum_residue_match_count"])
                minimum_coverage = float(thresholds["minimum_residue_match_percentage"])
            gates["match_count"] = all(int(payload.get("match_count", 0) or 0) >= minimum_matches for payload in (left, right))
            coverage_values = [row.get("coverage_left"), row.get("coverage_right")]
            if all(value is not None for value in coverage_values):
                limits = []
                for size in (row.get("template_size_left"), row.get("template_size_right")):
                    limit = minimum_coverage - float(thresholds["diff_percentage"]) if size and size > float(thresholds["template_residue_count"]) and (comparable or not multiprot) else minimum_coverage
                    limits.append(limit)
                comparator = (lambda value, limit: value >= limit) if (multiprot or comparable) else (lambda value, limit: value > limit)
                gates["coverage"] = all(comparator(value, limit) for value, limit in zip(coverage_values, limits))
            else:
                gates["coverage"] = None
                gate_status["coverage"] = "template_interface_size_missing"

            if filter_mode != "published_protocol":
                gate_status["hotspots"] = "not_applied"
                gate_status["complementary_contacts"] = "not_applied"
            elif filter_asset_root:
                try:
                    assets = load_filter_assets(row["template"], filter_asset_root)
                    by_chain = assets.get("hotspots_by_chain", {})
                    decision = evaluate_protocol_candidate(
                        left.get("match_dict", {}), right.get("match_dict", {}),
                        by_chain.get(row["chain_left"], []), by_chain.get(row["chain_right"], []),
                        _protocol_contacts(assets, row["chain_left"], row["chain_right"]),
                        minimum_contacts=int(thresholds["contact_count_threshold"]),
                        minimum_hotspots=int(thresholds["minimum_hotspot_match_number"]),
                    )
                    gates["hotspots"] = decision.reason == "published_protocol_passed"
                    gates["complementary_contacts"] = decision.complementary_contacts >= int(thresholds["contact_count_threshold"])
                    row["hotspot_matches"] = decision.hotspot_matches
                    row["complementary_contact_count"] = decision.complementary_contacts
                except (FileNotFoundError, ValueError) as exc:
                    gate_status["hotspots"] = f"assets_unavailable:{type(exc).__name__}"
                    gate_status["complementary_contacts"] = gate_status["hotspots"]
            else:
                gate_status["hotspots"] = "asset_root_missing"
                gate_status["complementary_contacts"] = "asset_root_missing"

        # A parseable JSON file is not sufficient evidence for replay.  A
        # nonzero/missing return code or missing raw-output hash makes every
        # alignment-dependent gate unknown, even if numerical fields happen
        # to be present in the record.
        if available and row.get("alignment_evidence_complete") is not True:
            reason = (
                "alignment_contract_incomplete:"
                f"{row.get('alignment_contract_status_left')}/"
                f"{row.get('alignment_contract_status_right')}"
            )
            for name in (
                "score_threshold", "match_count", "coverage", "hotspots",
                "complementary_contacts",
            ):
                gates[name] = None
                gate_status[name] = reason

        left_transform, right_transform = _transformed_pair_paths(run_root, row)
        gates["transform_materialization"] = bool(left_transform and right_transform)
        clash = None
        if left_transform and right_transform:
            if compute_clash_diagnostics:
                clash = clash_diagnostics(
                    left_transform, right_transform,
                    clash_distance=float(thresholds["clashing_distance"]),
                )
                gates["clash_filter"] = clash["clash_count"] < int(thresholds["max_clashing_count"])
            else:
                gate_status["clash_filter"] = "deferred_set_RUN_CLASH_DIAGNOSTICS_true"
        else:
            gate_status["clash_filter"] = "transform_not_materialized"

        failures = [name for name, value in gates.items() if value is False]
        unknown = [name for name, value in gates.items() if value is None]
        cumulative = False if failures else ("unknown" if unknown else True)
        audit = audits.get((row.get("query_left"), row.get("query_right"), row.get("template"), row.get("orientation")), {})
        row.update({
            "audit_status": audit.get("status", row.get("audit_status", "not_recorded")),
            "audit_error_reason": audit.get("error_reason", row.get("audit_error_reason")),
            "independent_failures": ",".join(failures),
            "independent_unknown_gates": ",".join(unknown),
            "cumulative_pass": cumulative,
            "transformed_left_path": str(left_transform) if left_transform else None,
            "transformed_right_path": str(right_transform) if right_transform else None,
            "clash_diagnostics": clash,
        })
        for name, value in gates.items():
            row[f"gate_{name}"] = value
            row[f"gate_{name}_status"] = gate_status.get(name)
        rows.append(row)
    return rows


def leave_one_gate_out(rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Report which candidates would be rescued when one gate is omitted."""

    output = []
    for row in rows:
        gates = {name: row.get(f"gate_{name}") for name in GATE_NAMES}
        for omitted in GATE_NAMES:
            remaining = [value for name, value in gates.items() if name != omitted]
            output.append({
                "candidate": candidate_key(row),
                "arm": row.get("arm"),
                "orientation": row.get("orientation"),
                "omitted_gate": omitted,
                "rescued_if_omitted": all(value is True for value in remaining),
                "unknown_remaining_gate": any(value is None for value in remaining),
            })
    return output


def _chain_ids(path: Path | None) -> list[str]:
    if path is None or not path.is_file():
        return []
    return sorted({line[21:22].strip() for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.startswith(("ATOM", "HETATM")) and len(line) > 21})


def refinement_inventory(run_root: str | Path, gate_rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Validate refinement outputs without equating file count with quality."""

    root = Path(run_root)
    events = load_jsonl(root / "status" / "stages.jsonl")
    refinement_events = [event for event in events if event.get("stage") == "refinement"]
    stage_event = refinement_events[-1] if refinement_events else None
    output_roots = [
        root / "processed" / "rosetta_refinement",
        root / "processed" / "pyrosetta_refinement",
        root / "processed" / "fiberdock_refinement",
    ]
    rows = []
    for gate in gate_rows:
        if not gate.get("gate_transform_materialization"):
            continue
        left = Path(gate["transformed_left_path"])
        right = Path(gate["transformed_right_path"])
        stem_token = left.stem
        candidates = [path for output_root in output_roots if output_root.exists() for path in output_root.rglob("*.pdb") if stem_token in path.stem]
        output_path = candidates[0] if candidates else None
        chain_ids = _chain_ids(output_path)
        rows.append({
            "template": gate.get("template"),
            "orientation": gate.get("orientation"),
            "query_left": gate.get("query_left"),
            "query_right": gate.get("query_right"),
            "score_gate_decision": gate.get("cumulative_pass"),
            "refinement_stage_status": stage_event.get("event") if stage_event else "not_run",
            "stage_return_code": stage_event.get("return_code") if stage_event else None,
            "per_candidate_return_code": None,
            "per_candidate_return_code_status": "missing_in_external_refiner_contract",
            "output_exists": output_path is not None,
            "output_path": str(output_path) if output_path else None,
            "output_chain_ids": ",".join(chain_ids),
            "output_chain_identity_status": "observed" if chain_ids else "not_observed",
            "failure_reason": None if output_path else "refinement_output_not_found",
            "dockq_status": "not_computed",
            "dockq": None,
            "irmsd": None,
            "fnat": None,
        })
    return rows


def candidate_key(record: Mapping[str, Any]) -> str:
    return "|".join(str(record.get(field, "")) for field in ("query_left", "query_right", "template", "orientation"))


def _fingerprint(keys: Iterable[str]) -> str:
    payload = "\n".join(sorted(set(keys))).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def ranking_comparison(
    candidate_sets: Mapping[str, Iterable[str]],
    *,
    pre_ranking_sets: Mapping[str, Iterable[str]] | None = None,
    native_labels: Mapping[str, bool] | None = None,
    native_label_source: str | None = None,
    native_label_source_sha256: str | None = None,
    top_k: int = 1,
) -> list[dict[str, Any]]:
    """Compare candidate sets; affinity is never used as a native-like label."""

    native_labels = native_labels or {}
    native_label_mapping_sha256 = hashlib.sha256(
        json.dumps(sorted((str(key), bool(value)) for key, value in native_labels.items()), separators=(",", ":")).encode("utf-8")
    ).hexdigest() if native_labels else None
    rows = []
    baseline = set(
        (pre_ranking_sets or {}).get("none", candidate_sets.get("none", []))
    )
    baseline_hash = _fingerprint(baseline)
    for method, values in candidate_sets.items():
        keys = list(dict.fromkeys(values))
        selected = set(keys)
        labelled = [key for key in keys[:top_k] if key in native_labels]
        if method == "none":
            same_pre_ranking_set = _fingerprint(keys) == baseline_hash
        elif pre_ranking_sets is None or method not in pre_ranking_sets:
            same_pre_ranking_set = None
        else:
            same_pre_ranking_set = baseline_hash == _fingerprint(pre_ranking_sets[method])
        native_provenance_complete = bool(native_labels) and bool(native_label_source) and bool(native_label_source_sha256) and bool(native_label_mapping_sha256)
        recovery_allowed = same_pre_ranking_set is True and native_provenance_complete
        rows.append({
            "method": method,
            "candidate_count": len(keys),
            "candidate_set_sha256": _fingerprint(keys),
            "same_pre_ranking_set": same_pre_ranking_set,
            "forwarded_fraction": len(selected) / len(baseline) if baseline else None,
            "native_labels_available": bool(native_labels),
            "native_label_source": native_label_source,
            "native_label_source_sha256": native_label_source_sha256,
            "native_label_mapping_sha256": native_label_mapping_sha256,
            "native_label_provenance_complete": native_provenance_complete,
            "native_label_recovery_status": "validated" if recovery_allowed else "deferred_panel_or_label_provenance",
            "top_k_native_like_recovery": sum(bool(native_labels[key]) for key in labelled) / len(labelled) if recovery_allowed and labelled else None,
            "affinity_used_as_label": False,
        })
    return rows


def usalign_preflight(executable: str | Path = "USalign", *, probe: bool = False) -> dict[str, Any]:
    """Check US-align availability; never claim PRISM compatibility from presence alone."""

    raw = str(executable)
    resolved = str(Path(raw).resolve()) if Path(raw).is_file() else shutil.which(raw)
    result: dict[str, Any] = {
        "requested": raw,
        "resolved": resolved,
        "available": bool(resolved),
        "probe_requested": probe,
        "return_code": None,
        "version_or_help": None,
        "prism_json_contract": "not_validated",
    }
    if resolved and probe:
        completed = subprocess.run([resolved, "-h"], capture_output=True, text=True, check=False, timeout=30)
        result["return_code"] = completed.returncode
        result["version_or_help"] = (completed.stdout or completed.stderr)[-2000:]
    return result


def validate_usalign_contract(payload: Mapping[str, Any]) -> dict[str, Any]:
    required = {
        "match_dict",
        "rotation_mat",
        "translation",
        "tm_score",
        "tm_score_query",
        "tm_score_ref",
        "tm_score_contract",
        "aligner",
        "return_code",
        "raw_output_sha256",
    }
    missing = sorted(required - set(payload))
    return {
        "valid_prism_contract": not missing,
        "missing_fields": missing,
        "score_contract": payload.get("score_gate_contract", payload.get("tm_score_contract")),
    }
