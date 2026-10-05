import json

import pytest

from src.candidate_audit import (
    CandidateAudit,
    CandidateRecord,
    alignment_features,
    record_alignment_pair,
)


def test_alignment_features_are_numeric_and_reproducible():
    result = alignment_features(
        {"match_count": 20, "tm_score": 0.71, "match_dict": {"A.1": "B.2"}},
        template_residue_count=40,
    )
    assert result == {
        "match_count": 20,
        "tm_score": pytest.approx(0.71),
        "match_coverage": pytest.approx(50.0),
        "mapping_count": 1,
    }


def test_audit_retains_one_record_with_both_orientations(tmp_path):
    path = tmp_path / "candidates.jsonl"
    audit = CandidateAudit(str(path))
    record_alignment_pair(
        audit,
        query_left="1abcA",
        query_right="2defB",
        template="3ghiAB",
        chain_left="A",
        chain_right="B",
        orientation="o1",
        left_alignment={"match_count": 18, "tm_score": 0.62, "match_dict": {}},
        right_alignment={"match_count": 16, "tm_score": 0.58, "match_dict": {}},
        template_size_left=30,
        template_size_right=32,
    )
    record_alignment_pair(
        audit,
        query_left="1abcA",
        query_right="2defB",
        template="3ghiAB",
        chain_left="B",
        chain_right="A",
        orientation="o2",
        left_alignment={"match_count": 2, "tm_score": 0.12},
        right_alignment={"match_count": 3, "tm_score": 0.10},
    )

    rows = [json.loads(line) for line in path.read_text().splitlines()]
    assert [row["orientation"] for row in rows] == ["o1", "o2"]
    assert rows[1]["status"] == "generated"
    assert rows[1]["match_count_left"] == 2


def test_unknown_status_is_rejected():
    with pytest.raises(ValueError, match="unknown candidate status"):
        CandidateRecord(
            query_left="a",
            query_right="b",
            template="t",
            chain_left="A",
            chain_right="B",
            orientation="o1",
            status="dropped",
        )
