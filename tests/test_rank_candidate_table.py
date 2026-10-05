import csv

from benchmark.scripts.rank_candidate_table import rank_csv


def test_rank_csv_keeps_failed_rows(tmp_path):
    source = tmp_path / "candidates.csv"
    source.write_text(
        "query_left,query_right,template,orientation,status,tm_score_left,tm_score_right,match_count_left,match_count_right\n"
        "a,b,good,o1,generated,0.8,0.8,40,40\n"
        "a,b,bad,o1,alignment_failed,0.9,0.9,50,50\n"
    )
    output = tmp_path / "ranked.csv"
    assert rank_csv(source, output) == 1
    rows = list(csv.DictReader(output.open()))
    assert rows[0]["baseline_rank"] == "1"
    assert rows[1]["baseline_rank"] == ""
    assert rows[1]["status"] == "alignment_failed"


def test_rank_csv_ranks_each_dataset_row_independently(tmp_path):
    source = tmp_path / "candidates.csv"
    source.write_text(
        "dataset_row_id,query_left,query_right,template,orientation,status,tm_score_left,tm_score_right,match_count_left,match_count_right\n"
        "rigid:000001,a,b,weak,o1,generated,0.4,0.4,20,20\n"
        "rigid:000001,a,b,strong,o1,generated,0.9,0.9,50,50\n"
        "medium:000001,c,d,only,o1,generated,0.5,0.5,25,25\n"
    )
    output = tmp_path / "ranked.csv"

    assert rank_csv(source, output) == 3
    rows = list(csv.DictReader(output.open()))

    assert [row["baseline_rank"] for row in rows] == ["2", "1", "1"]
    assert rows[0]["ranking_group"] == "dataset_row_id=rigid:000001"
    assert rows[2]["ranking_group"] == "dataset_row_id=medium:000001"
