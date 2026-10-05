"""Validation and leakage-safe splitting for candidate ranking tables."""

from __future__ import annotations

import hashlib
from typing import Iterable, Mapping, Sequence


REQUIRED_FEATURES = (
    "tm_score_left",
    "tm_score_right",
    "match_count_left",
    "match_count_right",
)


def native_like_label(dockq: float, threshold: float = 0.23) -> int:
    """Return the documented CAPRI native-like label."""
    if not 0.0 <= dockq <= 1.0:
        raise ValueError("DockQ must be between 0 and 1")
    return int(dockq >= threshold)


def validate_training_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    """Validate and copy labeled rows; unlabeled rows are rejected explicitly."""
    validated: list[dict[str, object]] = []
    for index, row in enumerate(rows):
        missing = [name for name in REQUIRED_FEATURES if name not in row]
        if missing:
            raise ValueError(f"row {index} missing features: {', '.join(missing)}")
        if row.get("native_complex_id") in (None, ""):
            raise ValueError(f"row {index} missing native_complex_id")
        if row.get("dockq") is None:
            raise ValueError(f"row {index} missing dockq label")
        copied = dict(row)
        copied["native_like"] = native_like_label(float(copied["dockq"]))
        validated.append(copied)
    if not validated:
        raise ValueError("no labeled candidate rows")
    return validated


def grouped_split(
    rows: Sequence[Mapping[str, object]],
    test_fraction: float = 0.2,
    seed: int = 0,
) -> tuple[list[int], list[int]]:
    """Split row indices by native complex, never by individual decoy rows."""
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("test_fraction must be between 0 and 1")
    groups = sorted({str(row["native_complex_id"]) for row in rows})
    ranked = sorted(
        groups,
        key=lambda group: hashlib.sha256(f"{seed}:{group}".encode()).hexdigest(),
    )
    test_count = max(1, round(len(groups) * test_fraction))
    test_groups = set(ranked[:test_count])
    test = [i for i, row in enumerate(rows) if str(row["native_complex_id"]) in test_groups]
    train = [i for i, row in enumerate(rows) if str(row["native_complex_id"]) not in test_groups]
    if not train or not test:
        raise ValueError("grouped split requires at least two native complexes")
    return train, test
