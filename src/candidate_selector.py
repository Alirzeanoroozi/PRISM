"""Bridge between candidate audit trail and ranking — selects top-K candidates per pair."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Mapping

from .candidate_ranker import biological_baseline_score, rank_candidates


def _parse_output_pdb_path(output_pdb: str) -> tuple[str | None, str | None, str | None, str | None]:
    """Parse a pair PDB path into (receptor, ligand, template, orientation).

    The path convention from ``transformation.py`` is::

        processed/transformation/{template}_{left_query}_{right_query}_{orientation}_{L|R}.pdb

    Returns (receptor_id, ligand_id, template_id, orientation) or (None, None, None, None).
    """
    basename = os.path.basename(output_pdb)
    for suffix in ("_L.pdb", "_R.pdb", "_l.pdb", "_r.pdb"):
        if basename.endswith(suffix):
            basename = basename[: -len(suffix)]
            break
    parts = basename.rsplit("_", 1)
    if len(parts) != 2:
        return None, None, None, None
    body = parts[0]
    body_parts = body.split("_", 2)
    if len(body_parts) != 3:
        return None, None, None, None
    template, receptor, ligand = body_parts
    orientation = parts[1]
    return receptor, ligand, template, orientation


def _load_audit_records(audit_path: str) -> list[dict[str, Any]]:
    """Load all JSONL audit records from the given path."""
    records = []
    with open(audit_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return records


def _build_audit_index(records: list[dict[str, Any]]) -> dict[tuple[str, str, str, str], dict[str, Any]]:
    """Index audit records by (receptor, ligand, template, orientation)."""
    index = {}
    for rec in records:
        key = (rec.get("query_left", ""), rec.get("query_right", ""), rec.get("template", ""), rec.get("orientation", ""))
        if all(key):
            index[key] = rec
    return index


def select_top_candidates(
    passed_pairs: list[tuple[str, str]],
    audit_path: str,
    top_k: int = 5,
    min_score: float | None = None,
    rank_method: str = "baseline",
    prodigy_executable: str = "prodigy",
    prodigy_output_dir: str = "processed/ranking/prodigy",
    prodigy_distance_cutoff: float = 5.5,
    prodigy_acc_threshold: float = 0.05,
    prodigy_temperature: float = 25.0,
    prodigy_timeout: float = 120.0,
) -> list[tuple[str, str]]:
    """Select top-K candidates per receptor-ligand pair from audit trail.

    Args:
        passed_pairs: List of (left_pdb, right_pdb) tuples from transformer()
        audit_path: Path to JSONL audit file written by CandidateAudit
        top_k: Maximum candidates to keep per (receptor, ligand) pair
        min_score: Optional minimum baseline score threshold
        rank_method: ``baseline`` or the opt-in external ``prodigy`` scorer

    Returns:
        Filtered list of (left_pdb, right_pdb) tuples for the selected candidates.
    """
    if not passed_pairs:
        return []
    if top_k < 1:
        raise ValueError("top_k must be positive")
    if rank_method not in {"baseline", "prodigy"}:
        raise ValueError(f"unknown ranking method: {rank_method}")

    if rank_method == "prodigy":
        from .prodigy_ranker import select_top_candidates_with_prodigy

        return select_top_candidates_with_prodigy(
            passed_pairs,
            top_k=top_k,
            executable=prodigy_executable,
            output_dir=prodigy_output_dir,
            distance_cutoff=prodigy_distance_cutoff,
            acc_threshold=prodigy_acc_threshold,
            temperature=prodigy_temperature,
            timeout=prodigy_timeout,
        )

    if not os.path.exists(audit_path):
        print(f"  WARNING: Audit file not found at {audit_path}; skipping ranking")
        return passed_pairs

    audit_records = _load_audit_records(audit_path)
    if not audit_records:
        print(f"  WARNING: Audit file empty; skipping ranking")
        return passed_pairs

    audit_index = _build_audit_index(audit_records)

    # Group passed_pairs by (receptor, ligand)
    from collections import defaultdict
    pairs_by_receptor_ligand: dict[tuple[str, str], list[tuple[str, str, str, str]]] = defaultdict(list)

    unparseable_pairs = []
    for left_pdb, right_pdb in passed_pairs:
        receptor, ligand, template, orientation = _parse_output_pdb_path(left_pdb)
        if all(v is not None for v in (receptor, ligand, template, orientation)):
            pairs_by_receptor_ligand[(receptor, ligand)].append((left_pdb, right_pdb, template, orientation))
        else:
            # Preserve future/legacy filename conventions that this selector
            # cannot identify safely. Ranking must never make them disappear.
            unparseable_pairs.append((left_pdb, right_pdb))

    selected_pairs = list(unparseable_pairs)

    for (receptor, ligand), candidates in pairs_by_receptor_ligand.items():
        # Build rankable rows from audit data
        rankable_rows = []
        incomplete_audit = False
        for left_pdb, right_pdb, template, orientation in candidates:
            audit_key = (receptor, ligand, template, orientation)
            audit_rec = audit_index.get(audit_key)
            if audit_rec is None:
                incomplete_audit = True
                continue
            # Only rank candidates that were successfully generated
            if audit_rec.get("status") not in {"generated", "refinement_accepted"}:
                incomplete_audit = True
                continue
            score = biological_baseline_score(audit_rec)
            if score is None:
                incomplete_audit = True
                continue
            if min_score is not None and score < min_score:
                continue
            rankable_rows.append({
                "left_pdb": left_pdb,
                "right_pdb": right_pdb,
                "template": template,
                "orientation": orientation,
                "baseline_score": score,
                **audit_rec,
            })

        if incomplete_audit or not rankable_rows:
            # Audit records are append-only and may not cover a newly produced
            # candidate. Ranking is
            # optional, so preserve the generated pair rather than turning a
            # partial audit mismatch into an empty refinement run.
            selected_pairs.extend((left_pdb, right_pdb) for left_pdb, right_pdb, _, _ in candidates)
            continue

        # Rank and select top-K
        ranked = rank_candidates(rankable_rows)
        top = ranked[:top_k]
        for row in top:
            selected_pairs.append((row["left_pdb"], row["right_pdb"]))

    print(f"  Ranking: {len(passed_pairs)} input pairs -> {len(selected_pairs)} selected (top {top_k} per pair)")
    return selected_pairs


def select_top_candidates_from_audit(
    audit_path: str,
    top_k: int = 5,
    min_score: float | None = None,
) -> list[tuple[str, str]]:
    """Alternative entry point: select top-K directly from audit without passed_pairs.

    Reconstructs PDB paths from audit records. Useful when audit is the source of truth.
    """
    if not os.path.exists(audit_path):
        return []

    audit_records = _load_audit_records(audit_path)
    if not audit_records:
        return []

    # Filter to generated/refinement_accepted and group before ranking. A
    # global top-K across unrelated receptor-ligand pairs is not meaningful.
    valid = [r for r in audit_records if r.get("status") in {"generated", "refinement_accepted"}]
    if not valid:
        return []

    if top_k < 1:
        raise ValueError("top_k must be positive")
    from collections import defaultdict
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in valid:
        grouped[(row.get("query_left", ""), row.get("query_right", ""))].append(row)

    # Reconstruct PDB paths from audit records
    selected = []
    for pair_rows in grouped.values():
        ranked = rank_candidates(pair_rows)
        if min_score is not None:
            ranked = [row for row in ranked if row["baseline_score"] >= min_score]
        for row in ranked[:top_k]:
            template = row.get("template", "")
            receptor = row.get("query_left", "")
            ligand = row.get("query_right", "")
            orientation = row.get("orientation", "")
            if all(v for v in (template, receptor, ligand, orientation)):
                left_pdb = f"processed/transformation/{template}_{receptor}_{ligand}_{orientation}_L.pdb"
                right_pdb = f"processed/transformation/{template}_{receptor}_{ligand}_{orientation}_R.pdb"
                if os.path.exists(left_pdb) and os.path.exists(right_pdb):
                    selected.append((left_pdb, right_pdb))

    return selected
