import json

from src import transformation


def test_opt_in_audit_serializes_terminal_status(tmp_path, monkeypatch):
    audit_path = tmp_path / "candidates.jsonl"
    monkeypatch.setattr(transformation, "AUDIT_PATH", str(audit_path))
    transformation.write_audit_record(
        "1abcAB", "left", "right", "A", "B", "o1",
        {"match_count": 2, "tm_score": 0.1},
        {"match_count": 3, "tm_score": 0.2},
        "alignment_failed",
    )
    row = json.loads(audit_path.read_text())
    assert row["status"] == "alignment_failed"


def test_explicit_audit_path_overrides_import_time_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("PRISM_CANDIDATE_AUDIT_PATH", str(tmp_path / "stale.jsonl"))
    explicit_path = tmp_path / "current.jsonl"
    transformation.write_audit_record(
        "1abcAB", "left", "right", "A", "B", "o1",
        {"match_count": 2, "tm_score": 0.1},
        {"match_count": 3, "tm_score": 0.2},
        "alignment_failed",
        audit_path=str(explicit_path),
    )
    assert explicit_path.exists()
    assert not (tmp_path / "stale.jsonl").exists()


def test_audit_records_resolved_transformation_thresholds(tmp_path):
    audit_path = tmp_path / "candidates.jsonl"
    transformation.write_audit_record(
        "1abcAB", "left", "right", "A", "B", "o1",
        {"match_count": 2, "tm_score": 0.1},
        {"match_count": 3, "tm_score": 0.2},
        "alignment_failed",
        audit_path=str(audit_path),
        thresholds={"minimum_residue_match_count": 12},
    )
    row = json.loads(audit_path.read_text())
    assert row["metadata"]["transformation_thresholds"]["minimum_residue_match_count"] == 12


def test_audit_records_real_template_coverage(tmp_path, monkeypatch):
    audit_path = tmp_path / "candidates.jsonl"
    monkeypatch.setattr(transformation, "template_size", {"1abcAB_A": 40, "1abcAB_B": 20})
    transformation.write_audit_record(
        "1abcAB", "left", "right", "A", "B", "o1",
        {"match_count": 20, "tm_score": 0.7},
        {"match_count": 15, "tm_score": 0.6},
        "generated",
        audit_path=str(audit_path),
    )

    row = json.loads(audit_path.read_text())
    assert row["match_coverage_left"] == 50.0
    assert row["match_coverage_right"] == 75.0
