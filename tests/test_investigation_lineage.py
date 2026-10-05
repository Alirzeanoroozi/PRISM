import hashlib
import json

import pytest

from benchmark.scripts.investigation_lineage import (
    LINEAGE_RECORD_FIELDS,
    AppendOnlyLineage,
    LineageRecord,
    aggregate_pair_summary,
    copy_pdb_immutable,
    create_pose_record,
    validate_ranking_keys,
)


def _record(**overrides):
    values = {
        "record_id": "event-1",
        "stage": "alignment",
        "pair_id": "pair-1",
        "template_id": "template-1",
        "side": "left",
        "alignment_id": "alignment-1",
        "orientation": "forward",
        "filter_stage": "pending",
        "pose_id": "pose-1",
        "refinement_id": "refinement-1",
        "energy_id": "energy-1",
        "scoring_id": "scoring-1",
        "status": "running",
    }
    values.update(overrides)
    return LineageRecord(**values)


def test_lineage_record_has_fixed_schema_and_terminal_failure_reason():
    record = _record(status="alignment_failed", failure_reason="no_match")

    assert tuple(record.to_dict()) == LINEAGE_RECORD_FIELDS
    assert record.is_terminal
    assert record.failure_reason == "no_match"
    with pytest.raises(AttributeError):
        record.status = "scored"


def test_append_is_idempotent_and_terminal_records_cannot_be_extended(tmp_path):
    log = AppendOnlyLineage(tmp_path / "lineage.jsonl")
    pending = _record()
    terminal = _record(
        record_id="event-2",
        status="alignment_failed",
        failure_reason="no_match",
    )

    assert log.append(pending) is True
    assert log.append(pending) is False
    assert log.append(terminal) is True
    assert log.append(terminal) is False
    with pytest.raises(ValueError, match="terminal"):
        log.append(_record(record_id="event-3", status="running"))

    assert len(log.read()) == 2
    assert len(log.path.read_text().splitlines()) == 2


def test_pose_record_contains_coordinate_sha256(tmp_path):
    pose = tmp_path / "pose.pdb"
    pose.write_bytes(b"ATOM\nMODEL 1\n")

    record = create_pose_record(
        pose,
        pair_id="pair-1",
        template_id="template-1",
        pose_id="pose-1",
        status="pose_created",
    )

    assert record.stage == "pose"
    assert record.coordinate_sha256 == hashlib.sha256(pose.read_bytes()).hexdigest()
    assert record.artifact_path == str(pose)


def test_immutable_pdb_copy_accepts_identical_bytes_and_refuses_overwrite(tmp_path):
    source = tmp_path / "source.pdb"
    destination = tmp_path / "destination.pdb"
    source.write_bytes(b"ATOM  identical\n")

    assert copy_pdb_immutable(source, destination) == destination
    assert copy_pdb_immutable(source, destination) == destination

    source.write_bytes(b"ATOM  changed\n")
    with pytest.raises(FileExistsError, match="bytes differ"):
        copy_pdb_immutable(source, destination)
    assert destination.read_bytes() == b"ATOM  identical\n"


def test_pair_summary_uses_unconditional_zero_without_scoreable_model():
    rows = [
        {
            "pair_id": "pair-no-score",
            "model_id": "model-1",
            "status": "alignment_failed",
            "failure_reason": "no_alignment",
        },
        {
            "pair_id": "pair-score",
            "model_id": "model-2",
            "status": "scored",
            "scoreable": True,
            "dockq": 0.4,
            "irmsd": 1.5,
        },
    ]

    summaries = aggregate_pair_summary(rows)

    no_score = next(row for row in summaries if row["pair_id"] == "pair-no-score")
    assert no_score["scoreable_model_count"] == 0
    assert no_score["secondary_model_count"] == 0
    assert no_score["primary_dockq"] == 0.0
    assert no_score["primary_irmsd"] == 0.0
    assert no_score["dockq_mean"] == 0.0
    assert no_score["irmsd_best"] == 0.0


def test_native_derived_fields_are_rejected_as_ranking_inputs():
    with pytest.raises(ValueError, match="native-derived"):
        validate_ranking_keys(["tm_score_left", "dockq"])

    assert validate_ranking_keys(["tm_score_left", "clash_count"]) == (
        "tm_score_left",
        "clash_count",
    )


def test_string_false_scoreability_is_not_counted_and_ranking_is_native_independent():
    summaries = aggregate_pair_summary(
        [
            {"pair_id": "pair-1", "model_id": "bad", "status": "scored", "scoreable": "False", "dockq": 0.9},
            {"pair_id": "pair-1", "model_id": "good", "status": "scored", "scoreable": "True", "dockq": 0.1, "tm_score_left": 0.8},
        ],
        ranking_keys=["tm_score_left"],
    )
    assert summaries[0]["scoreable_model_count"] == 1
    assert summaries[0]["primary_model_id"] == "good"
