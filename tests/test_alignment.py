"""Tests for alignment parsing utilities in src/alignment.py."""
import os
import json
import tempfile
import pytest

from src.alignment import extract_chain_and_res_ids


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_MINI_PDB = """\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C  
ATOM      2  CA  GLY A   2       4.000   5.000   6.000  1.00  0.00           C  
ATOM      3  CA  TRP B   1       7.000   8.000   9.000  1.00  0.00           C  
END
"""


class TestExtractChainAndResIds:
    def test_extracts_residue_ids(self):
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
            f.write(_MINI_PDB)
            path = f.name
        try:
            res_ids, chain_ids = extract_chain_and_res_ids("test", path)
            assert "1" in res_ids
            assert "2" in res_ids
        finally:
            os.unlink(path)

    def test_extracts_chain_ids(self):
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
            f.write(_MINI_PDB)
            path = f.name
        try:
            res_ids, chain_ids = extract_chain_and_res_ids("test", path)
            assert "A" in chain_ids
            assert "B" in chain_ids
        finally:
            os.unlink(path)

    def test_returns_equal_length_lists(self):
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
            f.write(_MINI_PDB)
            path = f.name
        try:
            res_ids, chain_ids = extract_chain_and_res_ids("test", path)
            assert len(res_ids) == len(chain_ids)
        finally:
            os.unlink(path)

    def test_counts_ca_atoms_only(self):
        pdb_with_non_ca = """\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C  
ATOM      2  CB  ALA A   1       1.500   2.500   3.500  1.00  0.00           C  
ATOM      3  CA  GLY A   2       4.000   5.000   6.000  1.00  0.00           C  
END
"""
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
            f.write(pdb_with_non_ca)
            path = f.name
        try:
            res_ids, chain_ids = extract_chain_and_res_ids("test", path)
            # Only 2 residues have CA atoms
            assert len(res_ids) == 2
        finally:
            os.unlink(path)

    def test_empty_pdb_returns_empty_lists(self):
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
            f.write("END\n")
            path = f.name
        try:
            res_ids, chain_ids = extract_chain_and_res_ids("test", path)
            assert res_ids == []
            assert chain_ids == []
        finally:
            os.unlink(path)
