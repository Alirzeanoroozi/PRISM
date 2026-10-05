"""Deterministic biological candidate baseline used before ML reranking."""

from __future__ import annotations

from math import isfinite
from typing import Mapping


BASELINE_SCORE_VERSION = "biological-baseline/v2-real-coverage"


def _bounded(value: object, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if isfinite(number) else default


def _optional_bounded(value: object) -> float | None:
    if value in (None, ""):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if isfinite(number) else None


def biological_baseline_score(row: Mapping[str, object]) -> float | None:
    """Score a candidate without allowing failed rows into the ranking.

    TM-score and matched-residue coverage are averaged across both chains.
    Clash count is an optional penalty; missing values do not become zero
    evidence and are recorded by callers as missing features.
    """
    if row.get("status") not in {"generated", "refinement_accepted"}:
        return None
    tm = (_bounded(row.get("tm_score_left")) + _bounded(row.get("tm_score_right"))) / 2.0
    coverage_left = _optional_bounded(row.get("match_coverage_left"))
    coverage_right = _optional_bounded(row.get("match_coverage_right"))
    # A fixed match-count denominator is not biological coverage: template
    # partner lengths vary substantially. Use coverage only when both real
    # denominators were recorded; otherwise rank from the available TM-score.
    if coverage_left is None or coverage_right is None:
        evidence_score = tm
    else:
        coverage = min((coverage_left + coverage_right) / 200.0, 1.0)
        evidence_score = 0.6 * tm + 0.4 * coverage
    clash_penalty = min(_bounded(row.get("clash_count")) / 10.0, 1.0) if row.get("clash_count") is not None else 0.0
    return evidence_score - 0.2 * clash_penalty


def rank_candidates(rows: list[Mapping[str, object]]) -> list[dict]:
    """Return complete candidates sorted by baseline score, best first."""
    ranked = []
    for row in rows:
        score = biological_baseline_score(row)
        if score is not None:
            ranked.append({
                **dict(row),
                "baseline_score": score,
                "baseline_score_version": BASELINE_SCORE_VERSION,
            })
    return sorted(ranked, key=lambda row: (-row["baseline_score"], row.get("template", ""), row.get("orientation", "")))
