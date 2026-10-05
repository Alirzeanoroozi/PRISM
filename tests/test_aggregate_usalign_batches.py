import json
from pathlib import Path

from benchmark.scripts.aggregate_usalign_batches import aggregate


def make_batch(root: Path, number: int, candidate: str, dockq: str) -> None:
    batch = root / "current" / f"batch_{number:04d}" / "status"
    batch.mkdir(parents=True)
    (batch / "transformed_dockq_status.json").write_text(
        json.dumps({"status": "validated_compacted", "candidate_rows": 1})
    )
    (batch / "candidate_generated.csv").write_text(
        "pipeline,case_id,template,status\nusalign,c1," + candidate + ",generated\n"
    )
    (batch / "transformed_dockq.tsv").write_text(
        "pipeline\tcase_id\ttemplate\tscore_status\tdockq_global\nusalign\tc1\t"
        + candidate
        + "\tscored\t"
        + dockq
        + "\n"
    )
    (batch / "transformed_dockq_interfaces.tsv").write_text(
        "pipeline\tcase_id\ttemplate\tDockQ\nusalign\tc1\t" + candidate + "\t" + dockq + "\n"
    )


def test_aggregate_validates_all_batches_and_compacts_tables(tmp_path: Path):
    make_batch(tmp_path, 1, "t1", "0.4")
    make_batch(tmp_path, 2, "t2", "0.2")
    result = aggregate(tmp_path, tmp_path / "aggregate", 2)
    assert result["status"] == "validated_compacted"
    assert result["candidate_rows"] == 2
    assert result["score_rows"] == 2
    assert (tmp_path / "aggregate" / "aggregation_status.json").is_file()
