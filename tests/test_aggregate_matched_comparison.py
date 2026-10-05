from pathlib import Path
import csv

from benchmark.scripts.aggregate_matched_comparison import aggregate


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_aggregate_emits_case_quality_and_ranking_tables(tmp_path: Path):
    candidates = write(
        tmp_path / "tmalign.csv",
        "pipeline,case_id,template,orientation,status,query_left,query_right,tm_score_left,tm_score_right\n"
        "tmalign,c1,t1,o1,generated,q1,q2,0.8,0.7\n"
        "tmalign,c1,t2,o1,generated,q1,q2,0.6,0.6\n"
        "tmalign,c2,t3,o1,clash_rejected,q1,q2,0.9,0.9\n",
    )
    scores = write(
        tmp_path / "tmalign_scores.tsv",
        "pipeline\tcase_id\ttemplate\torientation\tscore_status\tdockq_global\tdockq_cross_best\tirmsd_grouped_min\n"
        "tmalign\tc1\tt1\to1\tscored\t0.4\t0.5\t2.0\n"
        "tmalign\tc1\tt2\to1\tscored\t0.2\t0.3\t4.0\n",
    )
    case_list = write(tmp_path / "cases.csv", "case_id,split\nc1,rigid\nc2,medium\n")

    summary = aggregate(
        output_dir=tmp_path / "out",
        candidate_specs=[("TMalign", candidates)],
        score_specs=[("TMalign", scores)],
        case_list=case_list,
    )

    assert summary["status"] == "validated_compacted"
    with (tmp_path / "out" / "case_results.tsv").open(newline="", encoding="utf-8") as handle:
        case_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {(row["method"], row["case_id"]) for row in case_rows} == {
        ("TMalign", "c1"),
        ("TMalign", "c2"),
    }
    assert next(row for row in case_rows if row["case_id"] == "c1")["best_dockq_global"] == "0.4"
    with (tmp_path / "out" / "ranking_topk.tsv").open(newline="", encoding="utf-8") as handle:
        ranking_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert "oracle_best_dockq_global" in ranking_rows[0]
    assert "ranking_regret" in ranking_rows[0]


def test_missing_score_is_not_coerced_to_zero(tmp_path: Path):
    candidates = write(
        tmp_path / "gtalign.csv",
        "pipeline,case_id,template,orientation,status\n"
        "gtalign,c1,t1,o1,generated\n",
    )
    summary = aggregate(
        output_dir=tmp_path / "out",
        candidate_specs=[("GTalign", candidates)],
        score_specs=[],
        case_list=None,
    )
    assert summary["status"] == "validated_compacted"
    with (tmp_path / "out" / "case_results.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["best_dockq_global"] == ""
