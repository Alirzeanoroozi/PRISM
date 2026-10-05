from __future__ import annotations

import csv
import math
from pathlib import Path

from benchmark.scripts.aggregate_final_matched_comparison import aggregate


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_final_join_requires_exact_pair_for_delta(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.tsv"
    transformed = tmp_path / "transformed.tsv"
    refined = tmp_path / "refined.tsv"
    _write(candidate, [{
        "method": "tmalign", "status": "generated", "case_id": "c1",
        "template": "1abcAB", "orientation": "o1", "query_left": "qA",
        "query_right": "qB", "chain_left": "A", "chain_right": "B",
        "tm_score_left": "0.7", "tm_score_right": "0.8",
    }])
    _write(transformed, [{
        "case_id": "c1", "template": "1abcAB", "orientation": "o1",
        "query_left": "qA", "query_right": "qB", "chain_left": "A",
        "chain_right": "B", "score_status": "scored", "dockq_global": "0.4",
    }])
    _write(refined, [{
        "case_id": "c1", "template": "1abcAB", "orientation": "o1",
        "query_left": "qA", "query_right": "qB", "chain_left": "A",
        "chain_right": "B", "external_rosetta_status": "scored",
        "external_rosetta_global_dockq": "0.6",
    }, {
        "case_id": "c1", "template": "other", "orientation": "o1",
        "query_left": "qA", "query_right": "qB", "chain_left": "A",
        "chain_right": "B", "external_rosetta_status": "scored",
        "external_rosetta_global_dockq": "0.99",
    }])
    out = tmp_path / "out"
    result = aggregate(
        output_dir=out,
        candidate_specs=[("tmalign", candidate)],
        transformed_specs=[("tmalign", transformed)],
        refined_specs=[("tmalign", refined)],
        case_list=None,
    )
    assert result["paired_rows"] == 1
    with (out / "paired_refinement.tsv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert math.isclose(float(row["delta_dockq"]), 0.2)
    with (out / "case_results.tsv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["paired_candidate_count"] == "1"


def test_final_join_does_not_short_match_declared_chain_mismatch(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.tsv"
    refined = tmp_path / "refined.tsv"
    common = {
        "case_id": "c1", "template": "1abcAB", "orientation": "o1",
        "query_left": "qA", "query_right": "qB",
        "external_rosetta_status": "scored", "external_rosetta_global_dockq": "0.6",
    }
    _write(candidate, [{**common, "status": "generated", "chain_left": "A", "chain_right": "B"}])
    _write(refined, [{**common, "chain_left": "X", "chain_right": "Y"}])
    out = tmp_path / "out"
    result = aggregate(
        output_dir=out,
        candidate_specs=[("tmalign", candidate)],
        transformed_specs=[],
        refined_specs=[("tmalign", refined)],
        case_list=None,
    )
    assert result["paired_rows"] == 0


def test_final_ranking_reports_no_rank_topk_prodigy_and_oracle(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.tsv"
    transformed = tmp_path / "transformed.tsv"
    candidates = []
    scores = []
    for index, (template, tm, affinity, dockq, cross) in enumerate((
        ("a", "0.8", "-7.0", "0.4", "0.3"),
        ("b", "0.7", "-8.0", "0.6", "0.5"),
        ("c", "0.6", "-9.0", "0.7", "0.6"),
    )):
        common = {
            "case_id": "c1", "template": template, "orientation": "o1",
            "query_left": "qA", "query_right": "qB", "chain_left": "A",
            "chain_right": "B", "status": "generated", "tm_score_left": tm,
            "tm_score_right": tm, "prodigy_affinity_kcal_mol": affinity,
        }
        candidates.append(common)
        scores.append({
            "case_id": "c1", "template": template, "orientation": "o1",
            "query_left": "qA", "query_right": "qB", "chain_left": "A",
            "chain_right": "B", "score_status": "scored", "dockq_global": dockq,
            "dockq_cross_best": cross, "cross_irmsd_mean": str(index + 1),
            "cross_lrmsd_mean": str(index + 2),
        })
    _write(candidate, candidates)
    _write(transformed, scores)
    out = tmp_path / "out"
    result = aggregate(
        output_dir=out,
        candidate_specs=[("tmalign", candidate)],
        transformed_specs=[("tmalign", transformed)],
        refined_specs=[],
        case_list=None,
    )

    assert result["ranking_rows"] == 7
    with (out / "ranking_topk.tsv").open(newline="", encoding="utf-8") as handle:
        rows = {row["strategy"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert rows["all"]["selected_count"] == "3"
    assert rows["all"]["candidates_avoided"] == "0"
    assert math.isclose(float(rows["top3"]["best_selected_dockq"]), 0.7)
    assert math.isclose(float(rows["prodigy"]["selected_dockq"]), 0.7)
    assert math.isclose(float(rows["oracle"]["selected_dockq"]), 0.7)

    with (out / "case_results.tsv").open(newline="", encoding="utf-8") as handle:
        case = next(csv.DictReader(handle, delimiter="\t"))
    assert case["transformed_count"] == "3"
    assert case["transformed_scoreable_count"] == "3"
    assert math.isclose(float(case["best_transformed_cross_dockq"]), 0.6)
    with (out / "method_summary.tsv").open(newline="", encoding="utf-8") as handle:
        summary = next(row for row in csv.DictReader(handle, delimiter="\t") if row["split"] == "all")
    assert summary["candidate_total"] == "3"
    assert summary["top3_cases"] == "1"
    assert math.isclose(float(summary["top3_good_model_retention"]), 1.0)
