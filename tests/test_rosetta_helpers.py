"""Tests for the parts of rosetta_refinement that don't require PyRosetta."""
import os

from src.rosetta_refinement import (
    combine_pdb,
    _extract_chain_ids,
    _extract_atom_lines_by_partners,
    _resolve_entry,
)


def test_extract_chain_ids(tmp_path):
    p = tmp_path / "x.pdb"
    p.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      3  CA  ALA B   1       5.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      4  CA  ALA C   1      10.000   0.000   0.000  1.00  0.00           C  \n"
    )
    chains = _extract_chain_ids(str(p))
    assert chains == ["A", "B", "C"]


def test_extract_chain_ids_missing_file():
    assert _extract_chain_ids("/no/such/file.pdb") == []


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


def test_extract_atom_lines_by_partners(tmp_path):
    p = tmp_path / "x.pdb"
    p.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA B   1       5.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      3  CA  ALA C   1      10.000   0.000   0.000  1.00  0.00           C  \n"
    )
    left, right = _extract_atom_lines_by_partners(str(p), ["A", "B"], ["C"])
    assert len(left) == 2 and len(right) == 1


def test_resolve_entry_2tuple(tmp_path):
    a = tmp_path / "a.pdb"; b = tmp_path / "b.pdb"
    a.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n")
    b.write_text("ATOM      1  CA  ALA B   1       5.000   0.000   0.000  1.00  0.00           C  \n")
    rec, lig, lc, rc, combined = _resolve_entry((str(a), str(b)))
    assert rec == str(a) and lig == str(b)
    assert lc == ["A"] and rc == ["B"]
    assert combined is None


def test_resolve_entry_4tuple():
    res = _resolve_entry(("3i6eE", "3i6eF", "2ai9AB", "processed/output/file.pdb"))
    rec, lig, lc, rc, combined = res
    assert "3i6eE" in rec and "3i6eF" in lig
    assert combined == "processed/output/file.pdb"
