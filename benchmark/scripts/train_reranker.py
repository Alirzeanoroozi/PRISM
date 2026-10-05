#!/usr/bin/env python3
"""Train the optional Stage 1 candidate reranker.

The production pipeline does not import this module. Training requires an
environment with scikit-learn and a native-labeled candidate CSV.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import pickle
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.candidate_ranker import biological_baseline_score
from src.ranking_data import REQUIRED_FEATURES, grouped_split, validate_training_rows


def _ranking_quality(rows: list[dict[str, object]], scores: list[float]) -> dict[str, float]:
    """Average ranking quality across independent native-complex groups."""
    if not rows:
        return {
            "top1_native_like": 0.0,
            "top1_dockq": 0.0,
            "top10pct_native_like": 0.0,
            "top10pct_dockq": 0.0,
            "group_count": 0.0,
        }
    if len(rows) != len(scores):
        raise ValueError("rows and scores must have equal length")
    groups: dict[str, list[tuple[float, dict[str, object]]]] = {}
    for score, row in zip(scores, rows):
        group = str(row.get("dataset_row_id") or row.get("native_complex_id") or "")
        if not group:
            raise ValueError("ranking evaluation requires dataset_row_id or native_complex_id")
        groups.setdefault(group, []).append((score, row))

    per_group = []
    for group_rows in groups.values():
        ranked = sorted(group_rows, key=lambda item: item[0], reverse=True)
        top_count = max(1, round(len(ranked) * 0.10))
        top = [row for _, row in ranked[:top_count]]
        per_group.append({
            "top1_native_like": float(ranked[0][1]["native_like"]),
            "top1_dockq": float(ranked[0][1]["dockq"]),
            "top10pct_native_like": sum(float(row["native_like"]) for row in top) / len(top),
            "top10pct_dockq": sum(float(row["dockq"]) for row in top) / len(top),
        })
    return {
        metric: sum(group[metric] for group in per_group) / len(per_group)
        for metric in per_group[0]
    } | {"group_count": float(len(per_group))}


def train(
    csv_path: Path,
    model_path: Path,
    metrics_path: Path,
    seed: int = 0,
    test_fraction: float = 0.2,
) -> dict:
    try:
        from sklearn.ensemble import HistGradientBoostingClassifier
    except ImportError as exc:
        raise RuntimeError(
            "Stage 1 training requires scikit-learn; install it in a dedicated "
            "training environment, not in the production runtime"
        ) from exc

    with csv_path.open(newline="") as handle:
        raw_rows = list(csv.DictReader(handle))
    trainable_statuses = {"generated", "refinement_accepted"}
    labeled_rows = [
        row for row in raw_rows
        if row.get("status") in trainable_statuses
        and row.get("label_status", "labeled") == "labeled"
        and row.get("dockq") not in (None, "")
    ]
    if not labeled_rows:
        raise ValueError("candidate table contains no labeled rows")
    rows = validate_training_rows(labeled_rows)
    train_indices, test_indices = grouped_split(rows, test_fraction=test_fraction, seed=seed)
    features = list(REQUIRED_FEATURES)
    x_train = [[float(rows[i][name]) for name in features] for i in train_indices]
    x_test = [[float(rows[i][name]) for name in features] for i in test_indices]
    y_train = [rows[i]["native_like"] for i in train_indices]
    y_test = [rows[i]["native_like"] for i in test_indices]
    if len(set(y_train)) < 2:
        raise ValueError("training split must contain both native-like classes")
    model = HistGradientBoostingClassifier(random_state=seed)
    model.fit(x_train, y_train)
    probabilities = model.predict_proba(x_test)[:, 1]
    predictions = [int(value >= 0.5) for value in probabilities]
    accuracy = sum(a == b for a, b in zip(y_test, predictions)) / len(y_test)
    test_rows = [rows[i] for i in test_indices]
    learned_quality = _ranking_quality(test_rows, list(probabilities))
    baseline_scores = [biological_baseline_score(row) or 0.0 for row in test_rows]
    baseline_quality = _ranking_quality(test_rows, baseline_scores)
    accepted_statuses = {"generated", "refinement_accepted"}
    metrics = {
        "features": features,
        "seed": seed,
        "train_rows": len(train_indices),
        "test_rows": len(test_indices),
        "test_accuracy": accuracy,
        "test_native_like_prevalence": sum(y_test) / len(y_test),
        "test_native_complexes": sorted({str(row["native_complex_id"]) for row in test_rows}),
        "test_ranking_group_count": int(learned_quality["group_count"]),
        "test_fraction": test_fraction,
        "baseline_test_top1_native_like": baseline_quality["top1_native_like"],
        "learned_test_top1_native_like": learned_quality["top1_native_like"],
        "baseline_test_top1_dockq": baseline_quality["top1_dockq"],
        "learned_test_top1_dockq": learned_quality["top1_dockq"],
        "baseline_test_top10pct_dockq": baseline_quality["top10pct_dockq"],
        "learned_test_top10pct_dockq": learned_quality["top10pct_dockq"],
        "accepted_candidate_coverage": sum(
            row.get("status") in accepted_statuses for row in raw_rows
        ) / len(raw_rows) if raw_rows else 0.0,
        "supervised_label_coverage": len(labeled_rows) / len(raw_rows) if raw_rows else 0.0,
    }
    model_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    with model_path.open("wb") as handle:
        pickle.dump(model, handle)
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    args = parser.parse_args()
    print(json.dumps(train(args.csv, args.model, args.metrics, args.seed, args.test_fraction), indent=2))


if __name__ == "__main__":
    main()
