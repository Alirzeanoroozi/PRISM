#!/usr/bin/env python3
"""Attach canonical DockQ/iRMSD labels to audited candidate rows."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.ranking_data import native_like_label


def _number(value: str):
    try:
        return float(value)
    except (TypeError, ValueError):
        return ""


def _read_table(path: Path) -> list[dict[str, str]]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _model_path(row: dict[str, str]) -> str:
    for field in ("source_model_path", "model_pdb", "staged_model_path"):
        if row.get(field):
            return str(Path(row[field]).resolve())
    return ""


def attach_labels(candidate_path: Path, scores_path: Path, output_path: Path) -> int:
    candidates = _read_table(candidate_path)
    scores = _read_table(scores_path)
    canonical_scores = "dockq_cross_mean" in (scores[0] if scores else {})
    by_key: dict[tuple[str, str], dict[str, str]] = {}
    duplicates: set[tuple[str, str]] = set()
    for score in scores:
        if score.get("score_status") not in (None, "", "scored", "scored_cross_only"):
            continue
        if score.get("source_gate_status") == "audit_only":
            continue
        model = _model_path(score)
        model_hash = score.get("source_model_sha256", "")
        dataset_row_id = score.get("dataset_row_id", "")
        if canonical_scores and (not dataset_row_id or not model_hash):
            continue
        identity = model_hash or model
        if not identity:
            continue
        key = (dataset_row_id, identity)
        if key in by_key:
            duplicates.add(key)
        by_key[key] = score
    labeled = 0
    output = []
    for row in candidates:
        score_rows = []
        for field in ("model_complex", "model_left", "model_right"):
            model = row.get(field, "")
            if model:
                model_path = Path(model)
                dataset_row_id = row.get("dataset_row_id", "")
                if canonical_scores:
                    model_hash = row.get("source_model_sha256") or row.get("model_sha256") or ""
                    if not model_hash and model_path.is_file():
                        model_hash = _sha256(model_path)
                    identity = model_hash
                else:
                    identity = str(model_path.resolve())
                score = by_key.get((dataset_row_id, identity))
                if score is not None:
                    score_rows.append(score)
        result = dict(row)
        result["dockq"] = ""
        result["irmsd"] = ""
        result["native_like"] = ""
        result["label_status"] = "unlabeled"
        result["label_metric"] = ""
        result["label_model_sha256"] = ""
        result["label_score_scope"] = ""
        if canonical_scores and not row.get("dataset_row_id"):
            result["label_status"] = "missing_dataset_row_id"
        elif len(score_rows) == 1:
            score = score_rows[0]
            key = (
                score.get("dataset_row_id", ""),
                score.get("source_model_sha256", "") or _model_path(score),
            )
            if key in duplicates:
                result["label_status"] = "ambiguous"
                output.append(result)
                continue
            dockq_field = "dockq_cross_mean" if canonical_scores else "dockq"
            irmsd_field = "irmsd_grouped_min" if canonical_scores else "irmsd"
            dockq = _number(score.get(dockq_field, ""))
            irmsd = _number(score.get(irmsd_field, ""))
            if isinstance(dockq, float):
                result["dockq"] = f"{dockq:.8f}"
                result["native_like"] = str(native_like_label(dockq))
                result["label_status"] = "labeled"
                result["label_metric"] = dockq_field
                result["label_model_sha256"] = score.get("source_model_sha256", "")
                result["label_score_scope"] = score.get("score_scope", "")
                labeled += 1
            if isinstance(irmsd, float):
                result["irmsd"] = f"{irmsd:.8f}"
        elif len(score_rows) > 1:
            result["label_status"] = "ambiguous"
        output.append(result)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output[0]) if output else ["label_status"]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    return labeled


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidates", type=Path)
    parser.add_argument("scores", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(f"attached {attach_labels(args.candidates, args.scores, args.output)} native labels")


if __name__ == "__main__":
    main()
