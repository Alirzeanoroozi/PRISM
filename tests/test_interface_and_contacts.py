"""Tests for interface generation and contact detection (prescript API)."""

import os

import pytest

from src.contact import get_contacts_from_atom_lines


def test_get_contacts_from_atom_lines(tmp_path):
    lines_a = [
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n",
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C  \n",
    ]
    lines_b = [
        "ATOM      3  CA  ALA B   1       2.000   0.000   0.000  1.00  0.00           C  \n",
        "ATOM      4  N   ALA B   2      20.000   0.000   0.000  1.00  0.00           N  \n",
    ]
    out = tmp_path / "contacts.txt"
    contacts = get_contacts_from_atom_lines(str(tmp_path / "ignored.pdb"), str(out), lines_a, lines_b)
    # get_contacts_from_atom_lines returns list of (res0, res1) tuples
    assert len(contacts) >= 1
    assert out.exists()


@pytest.mark.skip(reason="generate_interface and get_contacts require templates/pdbs/ directory with real PDB files")
def test_generate_interface_two_chains():
    pass


@pytest.mark.skip(reason="get_contacts requires templates/pdbs/ directory with real PDB files")
def test_get_contacts_two_chains():
    pass
