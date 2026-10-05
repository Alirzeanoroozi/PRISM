"""Sanity tests for the benchmark scoring helpers (no DockQ/iRMSD required)."""
import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def _import(path, module_name):
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def score_module():
    return _import(REPO_ROOT / "benchmark/scripts/score_single_prism_pair.py", "score_single_prism_pair")


def test_infer_two_chains(tmp_path, score_module):
    pdb = tmp_path / "two.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA B   1       5.000   0.000   0.000  1.00  0.00           C  \n"
    )
    assert score_module.infer_chain_order(pdb) == ("A", "B")


def test_infer_multi_chains_split_on_ter(tmp_path, score_module):
    pdb = tmp_path / "multi.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA B   1       1.000   0.000   0.000  1.00  0.00           C  \n"
        "TER\n"
        "ATOM      3  CA  ALA C   1       2.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      4  CA  ALA D   1       3.000   0.000   0.000  1.00  0.00           C  \n"
    )
    receptor, ligand = score_module.infer_chain_order(pdb)
    assert receptor == "AB"
    assert ligand == "CD"


def test_infer_multi_chains_even_split(tmp_path, score_module):
    pdb = tmp_path / "four.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA B   1       1.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      3  CA  ALA C   1       2.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      4  CA  ALA D   1       3.000   0.000   0.000  1.00  0.00           C  \n"
    )
    receptor, ligand = score_module.infer_chain_order(pdb)
    assert receptor == "AB"
    assert ligand == "CD"


def test_infer_too_few_chains(tmp_path, score_module):
    pdb = tmp_path / "single.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
    )
    with pytest.raises(ValueError):
        score_module.infer_chain_order(pdb)


def test_parse_dockq_short(score_module):
    text = "DockQ 0.85 iRMSD 1.2 LRMSD 3.4 fnat 0.7 fnonnat 0.1 F1 0.8 clashes 0"
    result = score_module.parse_dockq_short(text)
    assert result["DockQ"] == 0.85
    assert result["iRMSD"] == 1.2
    assert result["clashes"] == 0


def test_maybe_float(score_module):
    assert score_module.maybe_float("1.5") == 1.5
    assert score_module.maybe_float("not a number") is None
