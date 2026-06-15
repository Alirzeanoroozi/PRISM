import json
import os

import pytest

from src.alignment import (
    check_alignment_passes_thresholds,
    _build_match_dict,
    _parse_matrix,
    _parse_tmalign_output,
)


def test_check_alignment_thresholds_pass():
    row = {"match_count": 30, "tm_score": 0.6, "len_template": 100}
    assert check_alignment_passes_thresholds(row) is True


def test_check_alignment_thresholds_fail_low_score():
    row = {"match_count": 30, "tm_score": 0.05, "len_template": 100}
    assert check_alignment_passes_thresholds(row) is False


def test_check_alignment_thresholds_fail_low_count():
    row = {"match_count": 5, "tm_score": 0.5, "len_template": 100}
    assert check_alignment_passes_thresholds(row) is False


def test_check_alignment_thresholds_short_template():
    row = {"match_count": 8, "tm_score": 0.5, "len_template": 40}
    assert check_alignment_passes_thresholds(row) is False
    row = {"match_count": 20, "tm_score": 0.5, "len_template": 40}
    assert check_alignment_passes_thresholds(row) is True


def test_parse_matrix(tmp_path):
    mtx = tmp_path / "m.out"
    mtx.write_text(
        "Header\n"
        " 0     1.234   0.99   0.00   0.00\n"
        " 1     5.678   0.00   0.99   0.00\n"
        " 2     9.000   0.00   0.00   0.99\n"
    )
    translation, rotation = _parse_matrix(str(mtx))
    assert translation == pytest.approx([1.234, 5.678, 9.0])
    for actual_row, expected_row in zip(rotation, [[0.99, 0, 0], [0, 0.99, 0], [0, 0, 0.99]]):
        assert actual_row == pytest.approx(expected_row)


def test_parse_tmalign_output(tmp_path):
    out = tmp_path / "x.tm"
    out.write_text(
        "Header\n"
        "Length of Chain_1:   100 residues\n"
        "Length of Chain_2:    80 residues\n"
        "Aligned length= 75, RMSD= 1.50, Seq_ID=n_identical/n_aligned=0.500\n"
        "TM-score= 0.700 (if normalized by length of Chain_1, by Lref)\n"
        "TM-score= 0.800 (if normalized by length of Chain_2, by Lref)\n"
        '(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)\n'
        "AC-DEF\n"
        "::.:::\n"
        "ACGDEF\n"
    )
    parsed = _parse_tmalign_output(str(out))
    assert parsed["match_count"] == 75
    assert parsed["tm_score"] == pytest.approx(0.8)
    assert parsed["len_target"] == 100
    assert parsed["len_template"] == 80
    assert parsed["seq1"] == "AC-DEF"
    assert parsed["seq2"] == "ACGDEF"


def test_build_match_dict_basic():
    target_residues = [("A", 1, "ALA"), ("A", 2, "CYS"), ("A", 3, "ASP"), ("A", 4, "GLU"), ("A", 5, "PHE")]
    template_residues = [("X", 10, "ALA"), ("X", 11, "CYS"), ("X", 12, "GLY"), ("X", 13, "ASP"), ("X", 14, "GLU"), ("X", 15, "PHE")]
    seq1 = "AC-DEF"
    seq2 = "ACGDEF"
    match = _build_match_dict(seq1, seq2, target_residues, template_residues)
    assert match == {
        "X.A.10": "A.A.1",
        "X.C.11": "A.C.2",
        "X.D.13": "A.D.3",
        "X.E.14": "A.E.4",
        "X.F.15": "A.F.5",
    }


def test_build_match_dict_empty():
    assert _build_match_dict("", "", [], []) == {}
