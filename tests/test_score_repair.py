import csv

import pytest

from benchmark.scripts.build_score_repair_manifest import build
from benchmark.scripts.merge_score_repair import merge


def _write(path, rows, delimiter):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def test_build_repair_manifest_uses_failed_durable_identities(tmp_path):
    stage = tmp_path / "stage.csv"
    scores = tmp_path / "scores.tsv"
    output = tmp_path / "repair.csv"
    stages = [
        {"dataset_row_id": "r:1", "source_model_sha256": "a", "status": "staged_symlink"},
        {"dataset_row_id": "r:2", "source_model_sha256": "b", "status": "staged_symlink"},
    ]
    _write(stage, stages, ",")
    _write(scores, [
        {**stages[0], "score_status": "scored"},
        {**stages[1], "score_status": "score_failed"},
    ], "\t")

    assert build(stage, scores, output) == 1
    assert next(csv.DictReader(output.open()))["dataset_row_id"] == "r:2"


def test_build_repair_manifest_can_select_auxiliary_failures(tmp_path):
    stage = tmp_path / "stage.csv"
    scores = tmp_path / "scores.tsv"
    output = tmp_path / "repair.csv"
    stages = [{"dataset_row_id": "r:1", "source_model_sha256": "a", "status": "staged_symlink"}]
    _write(stage, stages, ",")
    _write(scores, [{**stages[0], "score_status": "scored", "irmsd_status": "failed_auxiliary"}], "\t")

    from benchmark.scripts.prepare_bijective_score_retry import prepare

    assert len(prepare(
        scores, stage, output, score_status="scored",
        required_field="irmsd_status", required_value="failed_auxiliary",
    )) == 1


def test_merge_replaces_only_failed_rows_and_their_interfaces(tmp_path):
    base_models = tmp_path / "base_models.tsv"
    base_interfaces = tmp_path / "base_interfaces.tsv"
    repair_models = tmp_path / "repair_models.tsv"
    repair_interfaces = tmp_path / "repair_interfaces.tsv"
    base = [
        {"dataset_row_id": "r:1", "source_model_sha256": "a", "score_status": "scored"},
        {"dataset_row_id": "r:2", "source_model_sha256": "b", "score_status": "score_failed"},
    ]
    repair = [{**base[1], "score_status": "scored", "score_scope": "requested_cross_interfaces_only"}]
    _write(base_models, base, "\t")
    _write(base_interfaces, [
        {"source_model_sha256": "a", "interface": "AB"},
        {"source_model_sha256": "b", "interface": "old"},
    ], "\t")
    _write(repair_models, repair, "\t")
    _write(repair_interfaces, [{"source_model_sha256": "b", "interface": "BA"}], "\t")

    manifest = merge(base_models, base_interfaces, repair_models, repair_interfaces, tmp_path / "merged")

    rows = list(csv.DictReader((tmp_path / "merged/scores_models.tsv").open(), delimiter="\t"))
    interfaces = list(csv.DictReader((tmp_path / "merged/scores_interfaces.tsv").open(), delimiter="\t"))
    assert [row["score_status"] for row in rows] == ["scored", "scored"]
    assert [row["interface"] for row in interfaces] == ["AB", "BA"]
    assert manifest["replaced_rows"] == 1


def test_merge_rejects_replacement_of_accepted_base_row(tmp_path):
    base_models = tmp_path / "base_models.tsv"
    base_interfaces = tmp_path / "base_interfaces.tsv"
    repair_models = tmp_path / "repair_models.tsv"
    repair_interfaces = tmp_path / "repair_interfaces.tsv"
    row = {"dataset_row_id": "r:1", "source_model_sha256": "a", "score_status": "scored"}
    _write(base_models, [row], "\t")
    _write(base_interfaces, [{"source_model_sha256": "a", "interface": "AB"}], "\t")
    _write(repair_models, [row], "\t")
    _write(repair_interfaces, [{"source_model_sha256": "a", "interface": "AB"}], "\t")

    with pytest.raises(ValueError, match="disallowed"):
        merge(base_models, base_interfaces, repair_models, repair_interfaces, tmp_path / "merged")


def test_merge_can_replace_scored_row_with_only_auxiliary_failure(tmp_path):
    base_models = tmp_path / "base_models.tsv"
    base_interfaces = tmp_path / "base_interfaces.tsv"
    repair_models = tmp_path / "repair_models.tsv"
    repair_interfaces = tmp_path / "repair_interfaces.tsv"
    row = {
        "dataset_row_id": "r:1", "source_model_sha256": "a",
        "score_status": "scored", "irmsd_status": "failed_auxiliary",
    }
    _write(base_models, [row], "\t")
    _write(base_interfaces, [{"source_model_sha256": "a", "interface": "AB"}], "\t")
    _write(repair_models, [{**row, "irmsd_status": "scored"}], "\t")
    _write(repair_interfaces, [{"source_model_sha256": "a", "interface": "AB"}], "\t")

    manifest = merge(
        base_models, base_interfaces, repair_models, repair_interfaces, tmp_path / "merged",
        allowed_base_score_statuses={"scored"},
        required_base_field="irmsd_status", required_base_value="failed_auxiliary",
    )
    assert manifest["replaced_rows"] == 1
