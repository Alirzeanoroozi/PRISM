#!/usr/bin/env python3
"""Select a small, deterministic high-TM-score candidate set from BM55 results."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def find_native(attempt_root: Path, split: str, case_id: str) -> Path | None:
    del split  # case_id already carries the split prefix in the BM55 manifest.
    matches = sorted((attempt_root / "processed" / "pdbs").glob(f"native_{case_id}.pdb"))
    return matches[0] if len(matches) == 1 else None


def load_pair_metadata(
    attempt_root: Path,
    row: dict[str, object],
) -> tuple[dict[str, object], Path] | None:
    manifest_path = attempt_root / "dataset_manifest.json"
    if not manifest_path.is_file():
        return None
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    for pair in payload.get("pairs", []):
        if not isinstance(pair, dict):
            continue
        if pair.get("case_id") != row.get("case_id"):
            continue
        native_value = pair.get("native_pdb")
        if not isinstance(native_value, str):
            return None
        native = Path(native_value)
        if not native.is_absolute():
            native = attempt_root / native
        if not native.is_file():
            return None
        return pair, native.resolve()
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-csv", required=True, type=Path)
    parser.add_argument("--case-results-csv", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--selection-manifest", required=True, type=Path)
    parser.add_argument("--per-pipeline", type=int, default=20)
    parser.add_argument(
        "--skip-input-hashes",
        action="store_true",
        help="Do not hash every selected PDB during manifest construction; workers hash inputs before execution.",
    )
    parser.add_argument(
        "--omit-selected-rows",
        action="store_true",
        help="Keep the CSV authoritative but omit the full row list from the JSON manifest.",
    )
    parser.add_argument(
        "--skip-path-validation",
        action="store_true",
        help="Defer transformed/native file existence checks to the resumable worker.",
    )
    args = parser.parse_args()
    if args.per_pipeline < 1:
        raise SystemExit("--per-pipeline must be positive")

    cases: dict[tuple[str, int], dict[str, str]] = {}
    with args.case_results_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("pipeline") not in {"tmalign", "multiprot"}:
                continue
            cases[(row["pipeline"], int(row["case_index"]))] = row

    valid: list[dict[str, object]] = []
    rejected: list[dict[str, object]] = []
    with args.candidate_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pipeline = row.get("pipeline", "")
            if pipeline not in {"tmalign", "multiprot"} or row.get("status") != "generated":
                continue
            left_tm = number(row.get("tm_score_left", ""))
            right_tm = number(row.get("tm_score_right", ""))
            case = cases.get((pipeline, int(row["case_index"])))
            reason = ""
            if left_tm is None or right_tm is None:
                reason = "missing_tm_score"
            elif case is None:
                reason = "missing_case_result"
            if reason:
                rejected.append({"pipeline": pipeline, "case_index": row.get("case_index", ""), "template": row.get("template", ""), "reason": reason})
                continue
            minimum_tm = min(left_tm, right_tm)
            mean_tm = (left_tm + right_tm) / 2.0
            valid.append({
                **row,
                "attempt_root": str(Path(case["attempt_root"]).resolve()),
                "confidence_tm_min": minimum_tm,
                "confidence_tm_mean": mean_tm,
                "match_count_total": int(row.get("match_count_left", 0)) + int(row.get("match_count_right", 0)),
            })

    valid.sort(key=lambda item: (
        item["pipeline"], -float(item["confidence_tm_min"]),
        -float(item["confidence_tm_mean"]), -int(item["match_count_total"]),
        int(item["case_index"]), str(item["template"]), str(item["orientation"]),
    ))
    selected: list[dict[str, object]] = []
    pair_cache: dict[tuple[Path, str], tuple[dict[str, object], Path] | None] = {}
    for pipeline in ("tmalign", "multiprot"):
        rank = 0
        for row in (item for item in valid if item["pipeline"] == pipeline):
            attempt_root = Path(str(row["attempt_root"]))
            stem = f"{row['template']}_{row['query_left']}_{row['query_right']}_{row['orientation']}"
            left = attempt_root / "processed" / "transformation" / f"{stem}_L.pdb"
            right = attempt_root / "processed" / "transformation" / f"{stem}_R.pdb"
            pair_key = (attempt_root, str(row.get("case_id", "")))
            if pair_key not in pair_cache:
                pair_cache[pair_key] = load_pair_metadata(attempt_root, row)
            pair_data = pair_cache[pair_key]
            native = pair_data[1] if pair_data is not None else find_native(
                attempt_root, str(row.get("split", "")), str(row.get("case_id", ""))
            )
            if not args.skip_path_validation and (not left.is_file() or not right.is_file()):
                rejected.append({"pipeline": pipeline, "case_index": row.get("case_index", ""), "template": row.get("template", ""), "reason": "missing_transformed_pair"})
                continue
            if native is None or (not args.skip_path_validation and not native.is_file()):
                rejected.append({"pipeline": pipeline, "case_index": row.get("case_index", ""), "template": row.get("template", ""), "reason": "missing_native_pdb"})
                continue
            if pair_data is None:
                rejected.append({"pipeline": pipeline, "case_index": row.get("case_index", ""), "template": row.get("template", ""), "reason": "missing_dataset_pair"})
                continue
            pair, native = pair_data
            rank += 1
            selected.append({
                "selection_rank": rank,
                **row,
                "left": str(left.resolve()),
                "right": str(right.resolve()),
                "native_pdb": str(native.resolve()),
                "native_manifest_path": str((attempt_root / "dataset_manifest.json").resolve()),
                "receptor": str(pair.get("receptor", "")),
                "ligand": str(pair.get("ligand", "")),
                "native_receptor_chains": str(pair.get("native_receptor_chains", "")),
                "native_ligand_chains": str(pair.get("native_ligand_chains", "")),
                "benchmark_irmsd_A": str(pair.get("benchmark_irmsd_A", "")),
            })
            if rank >= args.per_pipeline:
                break
    if not args.skip_input_hashes:
        for row in selected:
            row["left_sha256"] = sha256(Path(str(row["left"])))
            row["right_sha256"] = sha256(Path(str(row["right"])))
            row["native_sha256"] = sha256(Path(str(row["native_pdb"])))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in selected for key in row})
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(selected)
    manifest = {
        "status": "ready",
        "candidate_csv": str(args.candidate_csv.resolve()),
        "candidate_csv_sha256": sha256(args.candidate_csv),
        "case_results_csv": str(args.case_results_csv.resolve()),
        "case_results_csv_sha256": sha256(args.case_results_csv),
        "selection_rule": "status=generated; existing transformed L/R and native PDB; descending min(side TM-score), mean TM-score, total match count",
        "per_pipeline_limit": args.per_pipeline,
        "selected_count": len(selected),
        "selected_by_pipeline": {pipeline: sum(row["pipeline"] == pipeline for row in selected) for pipeline in ("tmalign", "multiprot")},
        "rejected_candidate_rows": len(rejected),
        "output_csv": str(args.output_csv.resolve()),
        "selected": [] if args.omit_selected_rows else selected,
        "status_contract": "selected rows are not scientific success claims; they are a high-confidence bounded refinement probe",
    }
    args.selection_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: manifest[key] for key in ("status", "selected_count", "selected_by_pipeline", "rejected_candidate_rows", "output_csv")}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
