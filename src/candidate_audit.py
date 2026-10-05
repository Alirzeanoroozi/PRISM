"""Leakage-safe candidate audit records for PRISM ranking experiments.

The audit is deliberately independent of ranking or acceptance decisions. It
records candidates and failures in append-only JSONL so later models can use
the same table without silently dropping unsuccessful alignments.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import os
from typing import Any, Mapping, Optional


TERMINAL_STATUSES = {
    "generated",
    "alignment_missing",
    "alignment_failed",
    "protocol_rejected",
    "alignment_threshold_rejected",
    "transformation_failed",
    "clash_rejected",
    "refinement_failed",
    "refinement_accepted",
}


@dataclass
class CandidateRecord:
    query_left: str
    query_right: str
    template: str
    chain_left: str
    chain_right: str
    orientation: str
    source_pipeline: str = "tmalign_rosetta"
    status: str = "generated"
    match_count_left: int = 0
    match_count_right: int = 0
    tm_score_left: float = 0.0
    tm_score_right: float = 0.0
    match_coverage_left: Optional[float] = None
    match_coverage_right: Optional[float] = None
    contact_count: Optional[int] = None
    clash_count: Optional[int] = None
    rosetta_interaction_score: Optional[float] = None
    error_reason: Optional[str] = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.status not in TERMINAL_STATUSES:
            raise ValueError(f"unknown candidate status: {self.status}")
        if not self.orientation:
            raise ValueError("orientation must not be empty")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def alignment_features(
    alignment: Mapping[str, Any],
    template_residue_count: Optional[int] = None,
) -> dict[str, Any]:
    """Extract stable, numeric features from one alignment JSON object."""
    match_count = int(alignment.get("match_count", 0) or 0)
    tm_score = float(alignment.get("tm_score", 0.0) or 0.0)
    coverage = None
    if template_residue_count and template_residue_count > 0:
        coverage = 100.0 * match_count / template_residue_count
    return {
        "match_count": match_count,
        "tm_score": tm_score,
        "match_coverage": coverage,
        "mapping_count": len(alignment.get("match_dict", {}) or {}),
    }


class CandidateAudit:
    """Append candidate records without changing pipeline acceptance logic."""

    def __init__(self, path: str) -> None:
        self.path = path
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)

    def write(self, record: CandidateRecord) -> None:
        with open(self.path, "a", encoding="utf-8") as handle:
            handle.write(json.dumps(record.to_dict(), sort_keys=True) + "\n")


def record_alignment_pair(
    audit: CandidateAudit,
    *,
    query_left: str,
    query_right: str,
    template: str,
    chain_left: str,
    chain_right: str,
    orientation: str,
    left_alignment: Mapping[str, Any],
    right_alignment: Mapping[str, Any],
    template_size_left: Optional[int] = None,
    template_size_right: Optional[int] = None,
    metadata: Optional[Mapping[str, Any]] = None,
    status: str = "generated",
) -> CandidateRecord:
    left = alignment_features(left_alignment, template_size_left)
    right = alignment_features(right_alignment, template_size_right)
    record = CandidateRecord(
        query_left=query_left,
        query_right=query_right,
        template=template,
        chain_left=chain_left,
        chain_right=chain_right,
        orientation=orientation,
        status=status,
        error_reason=(
            None if status in {"generated", "refinement_accepted"} else status
        ),
        match_count_left=left["match_count"],
        match_count_right=right["match_count"],
        tm_score_left=left["tm_score"],
        tm_score_right=right["tm_score"],
        match_coverage_left=left["match_coverage"],
        match_coverage_right=right["match_coverage"],
        metadata={
            "mapping_count_left": left["mapping_count"],
            "mapping_count_right": right["mapping_count"],
            "tm_score_query_left": left_alignment.get("tm_score_query"),
            "tm_score_ref_left": left_alignment.get("tm_score_ref"),
            "tm_score_contract_left": left_alignment.get("tm_score_contract"),
            "tm_score_query_right": right_alignment.get("tm_score_query"),
            "tm_score_ref_right": right_alignment.get("tm_score_ref"),
            "tm_score_contract_right": right_alignment.get("tm_score_contract"),
            **dict(metadata or {}),
        },
    )
    audit.write(record)
    return record
