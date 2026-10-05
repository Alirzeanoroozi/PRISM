"""Tests for candidate_selector module."""

import json
import os
import tempfile
from pathlib import Path

import pytest

from src.candidate_selector import (
    select_top_candidates,
    select_top_candidates_from_audit,
    _parse_output_pdb_path,
)


def test_parse_output_pdb_path():
    """Test parsing of output PDB filenames."""
    rec, lig, tpl, orient = _parse_output_pdb_path("processed/transformation/1abc_1xyzAB_2defCD_o1_L.pdb")
    assert tpl == "1abc"
    assert rec == "1xyzAB"
    assert lig == "2defCD"
    assert orient == "o1"

    rec, lig, tpl, orient = _parse_output_pdb_path("processed/transformation/3broAD_3i6eE_3i6eF_o2_R.pdb")
    assert tpl == "3broAD"
    assert rec == "3i6eE"
    assert lig == "3i6eF"
    assert orient == "o2"

    # Invalid path
    rec, lig, tpl, orient = _parse_output_pdb_path("invalid.pdb")
    assert rec is None and lig is None and tpl is None and orient is None


def test_select_top_candidates_from_audit(tmp_path):
    """Test selecting top-K candidates directly from audit file."""
    audit_path = tmp_path / "audit.jsonl"

    # Write audit records with varying scores
    records = [
        {
            "query_left": "1abcA",
            "query_right": "1abcB",
            "template": "1xyzAB",
            "chain_left": "A",
            "chain_right": "B",
            "orientation": "o1",
            "status": "generated",
            "match_count_left": 30,
            "match_count_right": 25,
            "tm_score_left": 0.7,
            "tm_score_right": 0.6,
            "match_coverage_left": 60.0,
            "match_coverage_right": 50.0,
            "clash_count": 0,
        },
        {
            "query_left": "1abcA",
            "query_right": "1abcB",
            "template": "2xyzCD",
            "chain_left": "C",
            "chain_right": "D",
            "orientation": "o1",
            "status": "generated",
            "match_count_left": 20,
            "match_count_right": 15,
            "tm_score_left": 0.5,
            "tm_score_right": 0.4,
            "match_coverage_left": 40.0,
            "match_coverage_right": 30.0,
            "clash_count": 2,
        },
        {
            "query_left": "1abcA",
            "query_right": "1abcB",
            "template": "3xyzEF",
            "chain_left": "E",
            "chain_right": "F",
            "orientation": "o1",
            "status": "alignment_failed",  # Should be excluded
            "match_count_left": 10,
            "match_count_right": 10,
            "tm_score_left": 0.3,
            "tm_score_right": 0.3,
            "match_coverage_left": 20.0,
            "match_coverage_right": 20.0,
            "clash_count": 0,
        },
    ]

    with open(audit_path, "w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")

    # Create dummy PDB files for the first two templates
    for tpl in ["1xyzAB", "2xyzCD"]:
        left = f"processed/transformation/{tpl}_1abcA_1abcB_o1_L.pdb"
        right = f"processed/transformation/{tpl}_1abcA_1abcB_o1_R.pdb"
        Path(left).parent.mkdir(parents=True, exist_ok=True)
        Path(left).write_text("END\n")
        Path(right).write_text("END\n")

    selected = select_top_candidates_from_audit(str(audit_path), top_k=1)
    assert len(selected) == 1
    # Should pick the highest-scoring one (1xyzAB)
    assert "1xyzAB" in selected[0][0]


def test_select_top_candidates_with_passed_pairs(tmp_path):
    """Test selecting top-K from passed_pairs using audit."""
    audit_path = tmp_path / "audit.jsonl"

    records = [
        {
            "query_left": "1abcA",
            "query_right": "1abcB",
            "template": "1xyzAB",
            "chain_left": "A",
            "chain_right": "B",
            "orientation": "o1",
            "status": "generated",
            "match_count_left": 30,
            "match_count_right": 25,
            "tm_score_left": 0.7,
            "tm_score_right": 0.6,
            "match_coverage_left": 60.0,
            "match_coverage_right": 50.0,
            "clash_count": 0,
        },
        {
            "query_left": "1abcA",
            "query_right": "1abcB",
            "template": "2xyzCD",
            "chain_left": "C",
            "chain_right": "D",
            "orientation": "o1",
            "status": "generated",
            "match_count_left": 20,
            "match_count_right": 15,
            "tm_score_left": 0.5,
            "tm_score_right": 0.4,
            "match_coverage_left": 40.0,
            "match_coverage_right": 30.0,
            "clash_count": 2,
        },
    ]

    with open(audit_path, "w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")

    # Create dummy PDB files
    for tpl in ["1xyzAB", "2xyzCD"]:
        left = f"processed/transformation/{tpl}_1abcA_1abcB_o1_L.pdb"
        right = f"processed/transformation/{tpl}_1abcA_1abcB_o1_R.pdb"
        Path(left).parent.mkdir(parents=True, exist_ok=True)
        Path(left).write_text("END\n")
        Path(right).write_text("END\n")

    passed_pairs = [
        (f"processed/transformation/1xyzAB_1abcA_1abcB_o1_L.pdb",
         f"processed/transformation/1xyzAB_1abcA_1abcB_o1_R.pdb"),
        (f"processed/transformation/2xyzCD_1abcA_1abcB_o1_L.pdb",
         f"processed/transformation/2xyzCD_1abcA_1abcB_o1_R.pdb"),
    ]

    selected = select_top_candidates(passed_pairs, str(audit_path), top_k=1)
    assert len(selected) == 1
    assert "1xyzAB" in selected[0][0]


def test_select_top_candidates_missing_audit(tmp_path):
    """Test graceful handling when audit file is missing."""
    passed_pairs = [
        ("processed/transformation/1xyzAB_1abcA_1abcB_o1_L.pdb",
         "processed/transformation/1xyzAB_1abcA_1abcB_o1_R.pdb"),
    ]
    selected = select_top_candidates(passed_pairs, "/nonexistent/audit.jsonl", top_k=1)
    # Should return original pairs when audit missing
    assert selected == passed_pairs


def test_select_top_candidates_preserves_unmatched_pairs(tmp_path):
    """A partial/stale audit must not empty an otherwise valid pipeline run."""
    audit_path = tmp_path / "audit.jsonl"
    audit_path.write_text(json.dumps({
        "query_left": "otherA", "query_right": "otherB", "template": "staleAB",
        "orientation": "o1", "status": "generated", "tm_score_left": 0.8,
        "tm_score_right": 0.8, "match_count_left": 30, "match_count_right": 30,
    }) + "\n")
    passed_pairs = [(
        "processed/transformation/1xyzAB_1abcA_1abcB_o1_L.pdb",
        "processed/transformation/1xyzAB_1abcA_1abcB_o1_R.pdb",
    )]

    assert select_top_candidates(passed_pairs, str(audit_path), top_k=1) == passed_pairs


def test_select_top_candidates_preserves_unparseable_pairs(tmp_path):
    audit_path = tmp_path / "audit.jsonl"
    audit_path.write_text("{}\n")
    passed_pairs = [("processed/transformation/legacy-name.pdb", "right.pdb")]

    assert select_top_candidates(passed_pairs, str(audit_path), top_k=1) == passed_pairs


def test_select_top_candidates_empty_input():
    """Test with empty passed_pairs."""
    selected = select_top_candidates([], "/tmp/audit.jsonl", top_k=5)
    assert selected == []


def test_select_top_candidates_rejects_nonpositive_top_k(tmp_path):
    with pytest.raises(ValueError, match="top_k must be positive"):
        select_top_candidates(
            [("left.pdb", "right.pdb")], str(tmp_path / "audit.jsonl"), top_k=0,
        )


def test_audit_only_selector_groups_pairs_and_applies_min_score(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    audit_path = tmp_path / "audit.jsonl"
    records = [
        {
            "query_left": "leftA", "query_right": "rightA", "template": "1abcAB",
            "orientation": "o1", "status": "generated", "tm_score_left": 0.8,
            "tm_score_right": 0.8, "match_count_left": 20, "match_count_right": 20,
        },
        {
            "query_left": "leftB", "query_right": "rightB", "template": "2defCD",
            "orientation": "o1", "status": "generated", "tm_score_left": 0.4,
            "tm_score_right": 0.4, "match_count_left": 20, "match_count_right": 20,
        },
    ]
    audit_path.write_text("".join(json.dumps(row) + "\n" for row in records))
    for row in records:
        prefix = (
            f"processed/transformation/{row['template']}_{row['query_left']}_"
            f"{row['query_right']}_{row['orientation']}"
        )
        Path(prefix + "_L.pdb").parent.mkdir(parents=True, exist_ok=True)
        Path(prefix + "_L.pdb").write_text("END\n")
        Path(prefix + "_R.pdb").write_text("END\n")

    assert len(select_top_candidates_from_audit(str(audit_path), top_k=1)) == 2
    selected = select_top_candidates_from_audit(str(audit_path), top_k=1, min_score=0.6)
    assert len(selected) == 1
    assert "leftA_rightA" in selected[0][0]


def test_audit_only_selector_rejects_nonpositive_top_k(tmp_path):
    audit_path = tmp_path / "audit.jsonl"
    audit_path.write_text(json.dumps({
        "query_left": "left", "query_right": "right", "template": "1abcAB",
        "orientation": "o1", "status": "generated",
    }) + "\n")
    with pytest.raises(ValueError, match="top_k must be positive"):
        select_top_candidates_from_audit(str(audit_path), top_k=0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
