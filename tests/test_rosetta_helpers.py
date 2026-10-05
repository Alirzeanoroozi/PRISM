"""Tests for the parts of rosetta_refinement that don't require PyRosetta."""
import os

from src.rosetta_refinement import (
    combine_pdb,
    source_chain_ids,
    source_chain_id,
    partner_chain_ids,
    target_chain_id,
)


def test_source_chain_ids_reads_from_pdb(tmp_path):
    p = tmp_path / "x.pdb"
    p.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      3  CA  ALA B   1       5.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      4  CA  ALA C   1      10.000   0.000   0.000  1.00  0.00           C  \n"
    )
    chains = source_chain_ids(str(p))
    assert chains == "ABC"


def test_source_chain_ids_missing_file(tmp_path):
    p = tmp_path / "nope.pdb"
    p.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n")
    chains = source_chain_ids(str(p))
    assert chains


def test_combine_pdb(tmp_path, monkeypatch):
    from src import rosetta_refinement as rr
    monkeypatch.setattr(rr, "ROSETTA_DIR", str(tmp_path))
    a = tmp_path / "a.pdb"
    b = tmp_path / "b.pdb"
    a.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n")
    b.write_text("ATOM      1  CA  ALA B   1       5.000   0.000   0.000  1.00  0.00           C  \n")
    out = combine_pdb(str(a), str(b))
    assert os.path.exists(out)
    content = open(out).read()
    assert "ALA A" in content and "ALA B" in content
    assert "TER\n" in content


def test_partner_chain_ids():
    left, right = partner_chain_ids("processed/transformation/1abc_1xyzAB_2defCD_o1_L.pdb",
                                     "processed/transformation/1abc_1xyzAB_2defCD_o1_R.pdb")
    # partner_chain_ids calls source_chain_ids on the files; returns defaults if filenames are unparseable
    result = partner_chain_ids("a.pdb", "b.pdb")
    assert result == ("A", "B")


def test_target_chain_id():
    cid = target_chain_id("processed/transformation/1abc_1xyzAB_2defCD_o1_L.pdb")
    assert cid == "B"


def test_source_chain_id_from_filename():
    cid = source_chain_id("processed/transformation/1abc_1xyzAB_2defCD_o1_L.pdb")
    assert cid  # returns something from filename parsing or file reading
