#!/usr/bin/env python3
"""Deterministic template self-hit and homology gate for PRISM benchmarks."""

from __future__ import annotations

import json
from pathlib import Path

from Bio.Align import PairwiseAligner


SENSITIVITY_THRESHOLDS = (30, 40, 50, 70, 100)


def load_source_gate_policy(path: Path) -> dict[str, object]:
    """Load and validate the immutable source-authority decision record."""

    payload = json.loads(path.read_text(encoding="utf-8"))
    decision = payload.get("decision") if isinstance(payload, dict) else None
    if not isinstance(decision, dict):
        raise ValueError("source-gate policy has no decision object")
    required = {"status", "confirmatory_run_authorized", "excluded_dataset_row_ids"}
    missing = sorted(required - set(decision))
    if missing:
        raise ValueError("source-gate policy missing fields: " + ", ".join(missing))
    return decision


def confirmatory_row_eligibility(dataset_row_id: str, decision: dict[str, object]) -> tuple[bool, str]:
    """Apply the frozen source gate before template eligibility is considered."""

    excluded = {str(value) for value in decision.get("excluded_dataset_row_ids", [])}
    if dataset_row_id in excluded:
        return False, "source_gate_audit_only"
    if not bool(decision.get("confirmatory_run_authorized")):
        return False, "source_gate_not_authorized"
    if decision.get("status") != "authorized":
        return False, f"source_gate_status:{decision.get('status')}"
    return True, "source_gate_eligible"


def template_pdb_chains(template_id: str) -> tuple[str, str]:
    token = (template_id or "").strip()
    if len(token) < 5:
        raise ValueError(f"invalid template ID: {template_id!r}")
    return token[:4].lower(), token[4:].replace("_", "")


def sequence_similarity(query: str, template: str, *, aligner: PairwiseAligner | None = None) -> dict[str, float | int]:
    """Return a deterministic ungapped-match identity and bilateral coverage."""

    query, template = query.strip().upper(), template.strip().upper()
    if not query or not template:
        raise ValueError("missing sequence")
    alignment_engine = aligner or PairwiseAligner()
    alignment_engine.mode = "global"
    alignment = alignment_engine.align(query, template)[0]
    query_aligned, template_aligned = str(alignment[0]), str(alignment[1])
    paired = [(left, right) for left, right in zip(query_aligned, template_aligned) if left != "-" and right != "-"]
    matched = sum(left == right for left, right in paired)
    denominator = min(len(query), len(template))
    return {
        "identical_residues": matched,
        "identity_percent": 100.0 * matched / denominator,
        "query_coverage_percent": 100.0 * matched / len(query),
        "template_coverage_percent": 100.0 * matched / len(template),
        "shorter_coverage_percent": 100.0 * matched / denominator,
    }


def classify_template(
    *,
    target_pdb: str,
    target_chains: str,
    template_id: str,
    target_sequence: str,
    template_sequence: str,
    aligner: PairwiseAligner | None = None,
    similarity: dict[str, float | int] | None = None,
) -> dict[str, object]:
    """Return a transparent confirmatory eligibility decision for one partner."""

    template_pdb, template_chains = template_pdb_chains(template_id)
    similarity = similarity or sequence_similarity(target_sequence, template_sequence, aligner=aligner)
    self_pdb_chain = target_pdb.lower() == template_pdb and bool(set(target_chains) & set(template_chains))
    sequence_self_hit = similarity["identity_percent"] == 100.0 and min(
        similarity["query_coverage_percent"], similarity["template_coverage_percent"]
    ) >= 95.0
    if self_pdb_chain or sequence_self_hit:
        reason = "exact_self_hit"
    elif similarity["identity_percent"] > 50.0 and similarity["shorter_coverage_percent"] >= 70.0:
        reason = "homology_excluded_primary"
    else:
        reason = "eligible"
    return {
        "template_id": template_id,
        "template_pdb": template_pdb,
        "template_chains": template_chains,
        **similarity,
        "exclusion_reason": reason,
        "confirmatory_eligible": reason == "eligible",
        **{
            f"excluded_identity_gt_{threshold}": similarity["identity_percent"] > threshold
            and similarity["shorter_coverage_percent"] >= 70.0
            for threshold in SENSITIVITY_THRESHOLDS
        },
    }
