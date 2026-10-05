import csv

import pytest

from benchmark.scripts.merge_bijective_score_retry import merge
from benchmark.scripts.prepare_bijective_score_retry import prepare


def write_rows(path, rows, delimiter):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def test_prepare_and_merge_replace_only_failed_rows(tmp_path):
    base_models = tmp_path / "base-models.tsv"
    base_interfaces = tmp_path / "base-interfaces.tsv"
    stage = tmp_path / "stage.csv"
    retry_manifest = tmp_path / "retry.csv"
    retry_models = tmp_path / "retry-models.tsv"
    retry_interfaces = tmp_path / "retry-interfaces.tsv"
    output_models = tmp_path / "merged-models.tsv"
    output_interfaces = tmp_path / "merged-interfaces.tsv"
    write_rows(base_models, [
        {"dataset_row_id": "r:1", "source_model_sha256": "a", "staged_model_path": "one", "score_status": "scored"},
        {"dataset_row_id": "r:2", "source_model_sha256": "b", "staged_model_path": "two", "score_status": "score_failed"},
    ], "\t")
    write_rows(stage, [
        {"dataset_row_id": "r:1", "source_model_sha256": "a", "staged_model_path": "one"},
        {"dataset_row_id": "r:2", "source_model_sha256": "b", "staged_model_path": "two"},
    ], ",")
    write_rows(base_interfaces, [
        {"staged_model_path": "one", "record_type": "global"},
        {"staged_model_path": "two", "record_type": "stale"},
    ], "\t")
    assert [row["dataset_row_id"] for row in prepare(base_models, stage, retry_manifest)] == ["r:2"]
    write_rows(retry_models, [
        {"dataset_row_id": "r:2", "source_model_sha256": "b", "staged_model_path": "two", "score_status": "scored"},
    ], "\t")
    write_rows(retry_interfaces, [{"staged_model_path": "two", "record_type": "interface"}], "\t")

    assert merge(base_models, base_interfaces, retry_models, retry_interfaces, output_models, output_interfaces) == (2, 2, 1)
    with output_models.open(newline="", encoding="utf-8") as handle:
        models = list(csv.DictReader(handle, delimiter="\t"))
    with output_interfaces.open(newline="", encoding="utf-8") as handle:
        interfaces = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["score_status"] for row in models] == ["scored", "scored"]
    assert [row["record_type"] for row in interfaces] == ["global", "interface"]


def test_merge_refuses_to_replace_successful_base_row(tmp_path):
    paths = [tmp_path / name for name in ("bm", "bi", "rm", "ri", "om", "oi")]
    write_rows(paths[0], [{"dataset_row_id": "r:1", "source_model_sha256": "a", "staged_model_path": "one", "score_status": "scored"}], "\t")
    write_rows(paths[1], [{"staged_model_path": "one", "record_type": "global"}], "\t")
    write_rows(paths[2], [{"dataset_row_id": "r:1", "source_model_sha256": "a", "staged_model_path": "one", "score_status": "scored"}], "\t")
    write_rows(paths[3], [{"staged_model_path": "one", "record_type": "global"}], "\t")
    with pytest.raises(ValueError, match="non-failed"):
        merge(*paths)
