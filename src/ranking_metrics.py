"""Dependency-free ranking metrics for PRISM candidate experiments."""

from __future__ import annotations

from math import isnan
from typing import Mapping, Sequence


def _label(row: Mapping[str, object]) -> int:
    value = row.get("native_like")
    if value in (1, "1", True):
        return 1
    if value in (0, "0", False):
        return 0
    raise ValueError("labeled ranking metrics require native_like")


def top_k_success(rows: Sequence[Mapping[str, object]], fraction: float) -> float:
    if not 0.0 < fraction <= 1.0:
        raise ValueError("fraction must be in (0, 1]")
    ranked = sorted(rows, key=lambda row: float(row["baseline_score"]), reverse=True)
    count = max(1, round(len(ranked) * fraction))
    return sum(_label(row) for row in ranked[:count]) / count if ranked else 0.0


def enrichment_factor(rows: Sequence[Mapping[str, object]], fraction: float) -> float:
    if not rows:
        return 0.0
    prevalence = sum(_label(row) for row in rows) / len(rows)
    return top_k_success(rows, fraction) / prevalence if prevalence else 0.0


def spearman_score(rows: Sequence[Mapping[str, object]]) -> float:
    """Spearman correlation between score rank and binary biological label."""
    if len(rows) < 2:
        return 0.0
    ranked = sorted(rows, key=lambda row: float(row["baseline_score"]), reverse=True)
    score_ranks = {id(row): i + 1 for i, row in enumerate(ranked)}
    labels = [_label(row) for row in rows]
    label_order = sorted(range(len(rows)), key=lambda i: labels[i], reverse=True)
    label_ranks = {index: rank + 1 for rank, index in enumerate(label_order)}
    differences = [score_ranks[id(row)] - label_ranks[i] for i, row in enumerate(rows)]
    return 1.0 - (6.0 * sum(diff * diff for diff in differences)) / (len(rows) * (len(rows) ** 2 - 1))
