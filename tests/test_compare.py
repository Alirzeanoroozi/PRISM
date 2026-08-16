"""Tests for output-vs-native comparison helpers."""

from pathlib import Path

import pytest

from sandbox.compare_outputs import parse_output_filename
from src.compare import ca_rmsd, compare_pair, get_trimmed_native_pdb
from tests.conftest import _atom


def _write_two_chain_native(path: Path):
  lines = [
      _atom(1, "CA", "ALA", "E", 1, 0.0, 0.0, 0.0),
      _atom(2, "CA", "ALA", "E", 2, 4.0, 0.0, 0.0),
      _atom(3, "CA", "ALA", "E", 3, 8.0, 0.0, 0.0),
      _atom(4, "CA", "ALA", "F", 1, 0.0, 5.0, 0.0),
      _atom(5, "CA", "ALA", "F", 2, 4.0, 5.0, 0.0),
      "END\n",
  ]
  path.write_text("".join(lines))


def test_parse_output_filename():
    tpl, rec, lig = parse_output_filename(Path("3broAD_3i6eE_3i6eF.pdb"))
    assert tpl == "3broAD"
    assert rec == "3i6eE"
    assert lig == "3i6eF"


def test_ca_rmsd_identical(tmp_path, repo_root_chdir, monkeypatch):
    native = tmp_path / "processed" / "pdbs" / "3i6e.pdb"
    native.parent.mkdir(parents=True)
    _write_two_chain_native(native)
    model = tmp_path / "model.pdb"
    model.write_text(native.read_text())

    rmsd, n = ca_rmsd(str(model), str(native), ["E"], ["E"])
    assert n == 3
    assert rmsd == pytest.approx(0.0, abs=1e-6)


def test_trimmed_native_same_pdb(tmp_path, repo_root_chdir, monkeypatch):
    native = tmp_path / "processed" / "pdbs" / "3i6e.pdb"
    native.parent.mkdir(parents=True)
    _write_two_chain_native(native)
    monkeypatch.chdir(tmp_path)

    path, rec, lig = get_trimmed_native_pdb("3i6eE", "3i6eF")
    assert path.endswith("3i6e_3i6e_EF.pdb")
    assert rec == "E"
    assert lig == "F"
    assert Path(path).exists()


def test_compare_pair_missing_output(repo_root_chdir):
    row = compare_pair("3i6eE", "3i6eF", "3broAD", "processed/output/missing.pdb")
    assert "not found" in row["error_dockq"].lower()
