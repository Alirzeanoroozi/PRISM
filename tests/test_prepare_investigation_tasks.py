import pytest

from benchmark.scripts.prepare_investigation_tasks import build_task_rows


def test_builds_one_explicit_task_per_row_id():
    ids = [f"rigid:{index:06d}" for index in range(1, 11)]
    rows = build_task_rows(
        ids,
        "python validate.py --dataset-row-id {dataset_row_id} --index {array_index}",
        input_paths=["source_manifest.tsv"],
        output_paths=["manifest/source_manifest.tsv", "manifest/structure_validation.tsv"],
    )
    assert [row["array_index"] for row in rows] == [str(index) for index in range(1, 11)]
    assert rows[0]["task_id"] == "rigid:000001:task-0001"
    assert "rigid:000001" in rows[0]["command"]
    assert rows[-1]["array_index"] == "10"


@pytest.mark.parametrize("ids", [[], ["rigid:000001"] * 10])
def test_rejects_non_unique_ten_row_selection(ids):
    with pytest.raises(ValueError):
        build_task_rows(ids, "true")
