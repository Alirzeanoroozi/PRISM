import json
import os

import pytest

from src.alignment import (
    build_alignment_command,
    parse_tmalign,
    _has_ca_atoms,
    _valid_tmalign_outputs,
    extract_chain_and_res_ids,
)
from src.transformation import alignment_passes_thresholds


def test_alignment_passes_thresholds_pass():
    """alignment_passes_thresholds lives in transformation.py now."""
    row = {"match_count": 30, "tm_score": 0.6, "len_template": 100}
    assert alignment_passes_thresholds("1abc_A", row) is True


def test_alignment_passes_thresholds_fail_low_score():
    row = {"match_count": 30, "tm_score": 0.05, "len_template": 100}
    assert alignment_passes_thresholds("1abc_A", row) is False


def test_alignment_passes_thresholds_fail_low_count():
    row = {"match_count": 5, "tm_score": 0.5, "len_template": 100}
    assert alignment_passes_thresholds("1abc_A", row) is False


def test_alignment_passes_thresholds_short_template():
    row = {"match_count": 8, "tm_score": 0.5, "len_template": 40}
    assert alignment_passes_thresholds("1abc_A", row) is False
    row = {"match_count": 20, "tm_score": 0.5, "len_template": 40}
    assert alignment_passes_thresholds("1abc_A", row) is True


def test_has_ca_atoms(tmp_path):
    p = tmp_path / "x.pdb"
    p.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \nEND\n")
    assert _has_ca_atoms(str(p)) is True
    empty = tmp_path / "empty.pdb"
    empty.write_text("END\n")
    assert _has_ca_atoms(str(empty)) is False


def test_valid_tmalign_outputs(tmp_path):
    mtx = tmp_path / "m.out"
    mtx.write_text("0 0 1 0 0\n1 0 0 1 0\n2 0 0 0 1\n")
    tm = tmp_path / "t.out"
    tm.write_text('(": denotes residue pairs...\nTM-score= 0.500\nAligned length= 20\n')
    assert _valid_tmalign_outputs(str(mtx), str(tm)) is True
    assert _valid_tmalign_outputs(str(mtx), "/no/such/file") is False


def test_parse_tmalign_reads_usalign_structure_score_labels(tmp_path):
    protein = tmp_path / "query.pdb"
    interface = tmp_path / "interface.pdb"
    pdb = "".join(
        f"ATOM  {index:5d}  CA  ALA A{index:4d}    {float(index):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C  \n"
        for index in range(1, 4)
    ) + "END\n"
    protein.write_text(pdb)
    interface.write_text(pdb)

    matrix = tmp_path / "matrix.out"
    matrix.write_text("0 0 1 0 0\n1 0 0 1 0\n2 0 0 0 1\n")
    output = tmp_path / "usalign.out"
    output.write_text(
        "Aligned length= 3, RMSD= 0.00, Seq_ID=n_identical/n_aligned= 1.000\n"
        "TM-score= 0.40000 (normalized by length of Structure_1: L=3, d0=1.00)\n"
        "TM-score= 0.80000 (normalized by length of Structure_2: L=3, d0=1.00)\n"
        "(\":\" denotes residue pairs)\n"
        "AAA\n"
        ":::\n"
        "AAA\n"
    )

    output_dir = tmp_path / "parsed"
    parse_tmalign(
        str(protein),
        str(interface),
        "query",
        "templAB",
        "A",
        str(matrix),
        str(output),
        str(output_dir),
        aligner_name="USalign",
    )

    row = json.loads((output_dir / "query_templAB_A.json").read_text())
    assert row["tm_score"] == pytest.approx(0.8)
    assert row["aligner"] == "USalign"


def test_parse_tmalign_preserves_query_and_reference_normalizations(tmp_path):
    protein = tmp_path / "query.pdb"
    interface = tmp_path / "interface.pdb"
    pdb = "".join(
        f"ATOM  {index:5d}  CA  ALA A{index:4d}    {float(index):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C  \n"
        for index in range(1, 4)
    ) + "END\n"
    protein.write_text(pdb)
    interface.write_text(pdb)

    matrix = tmp_path / "matrix.out"
    matrix.write_text("0 0 1 0 0\n1 0 0 1 0\n2 0 0 0 1\n")
    output = tmp_path / "usalign.out"
    output.write_text(
        "Aligned length= 3, RMSD= 0.00, Seq_ID=n_identical/n_aligned= 1.000\n"
        "TM-score= 0.40000 (normalized by length of Structure_1: L=3, d0=1.00)\n"
        "TM-score= 0.80000 (normalized by length of Structure_2: L=3, d0=1.00)\n"
        "(\":\" denotes residue pairs)\n"
        "AAA\n"
        ":::\n"
        "AAA\n"
    )

    output_dir = tmp_path / "parsed"
    parse_tmalign(
        str(protein),
        str(interface),
        "query",
        "templAB",
        "A",
        str(matrix),
        str(output),
        str(output_dir),
        aligner_name="USalign",
    )

    row = json.loads((output_dir / "query_templAB_A.json").read_text())
    assert row["tm_score_query"] == pytest.approx(0.4)
    assert row["tm_score_ref"] == pytest.approx(0.8)
    assert row["tm_score"] == pytest.approx(0.8)
    assert row["tm_score_contract"] == "reference_normalized_structure_2"


def test_build_alignment_command_adds_usalign_options_without_changing_tmalign(monkeypatch):
    monkeypatch.setenv("PRISM_ALIGNMENT_NAME", "USalign")
    monkeypatch.setenv("PRISM_USALIGN_FAST", "true")
    assert build_alignment_command("USalign", "q.pdb", "r.pdb", "m.out") == [
        "USalign", "q.pdb", "r.pdb", "-fast", "-outfmt", "-1", "-m", "m.out"
    ]

    monkeypatch.setenv("PRISM_ALIGNMENT_NAME", "TMalign")
    monkeypatch.delenv("PRISM_USALIGN_FAST")
    assert build_alignment_command("TMalign", "q.pdb", "r.pdb", "m.out") == [
        "TMalign", "q.pdb", "r.pdb", "-m", "m.out"
    ]


def test_extract_chain_and_res_ids(tmp_path):
    p = tmp_path / "x.pdb"
    p.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  CYS A   2       1.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      3  CA  ASP B   1       5.000   0.000   0.000  1.00  0.00           C  \n"
    )
    residue_ids, chain_id = extract_chain_and_res_ids("test", str(p))
    assert len(residue_ids) == 3
    assert len(chain_id) == 3
    assert chain_id == ["A", "A", "B"]
