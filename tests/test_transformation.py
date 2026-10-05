"""Tests for transformation stage helpers matching the prescript API."""

import os
from pathlib import Path

import pytest

from src.transformation import (
    apply_tm_transform,
    pair_has_acceptable_clashes,
    alignment_passes_thresholds,
)


def test_alignment_passes_thresholds():
    assert alignment_passes_thresholds("1abc_A", {"match_count": 30, "tm_score": 0.6, "len_template": 100}) is True
    assert alignment_passes_thresholds("1abc_A", {"match_count": 5, "tm_score": 0.5, "len_template": 100}) is False


def test_apply_tm_transform_identity(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA B   1       4.000   5.000   6.000  1.00  0.00           C  \n"
    )
    out = tmp_path / "out.pdb"
    apply_tm_transform(str(src), str(out), [0.0, 0.0, 0.0], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    content = out.read_text()
    lines = [l for l in content.splitlines() if l.startswith("ATOM")]
    assert len(lines) == 2
    assert "1.000   2.000   3.000" in lines[0]


def test_apply_tm_transform_translation(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n")
    out = tmp_path / "out.pdb"
    apply_tm_transform(str(src), str(out), [10.0, -5.0, 0.5], [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    content = out.read_text()
    assert "10.000  -5.000   0.500" in content


def test_clash_check(tmp_path):
    r = tmp_path / "r.pdb"
    l = tmp_path / "l.pdb"
    r.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n")
    l.write_text("ATOM      1  CA  ALA B   1      10.000   0.000   0.000  1.00  0.00           C  \n")
    assert pair_has_acceptable_clashes(str(r), str(l)) is True
    # 6 CA atoms within 3Å of each other exceed MAX_CLASHING_COUNT=5
    r.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA A   2       2.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      3  CA  ALA A   3       0.000   2.000   0.000  1.00  0.00           C  \n"
    )
    l.write_text(
        "ATOM      1  CA  ALA B   1       0.500   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA B   2       0.500   2.000   0.000  1.00  0.00           C  \n"
    )
    assert pair_has_acceptable_clashes(str(r), str(l)) is False
