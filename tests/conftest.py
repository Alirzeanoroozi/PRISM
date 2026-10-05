"""Pytest fixtures: synthetic mini PDB structures for offline unit tests."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _atom(serial, atom_name, res_name, chain, res_seq, x, y, z, element="C", occupancy=1.0, bfactor=0.0):
    return (
        f"ATOM  {serial:5d} {atom_name:>4s} {res_name:3s} {chain}{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}           {element:>2s}  \n"
    )


def _ala_residue(serial0, chain, res_seq, base):
    bx, by, bz = base
    return [
        _atom(serial0 + 0, "N", "ALA", chain, res_seq, bx + 0.0, by + 0.0, bz, "N"),
        _atom(serial0 + 1, "CA", "ALA", chain, res_seq, bx + 1.5, by + 0.0, bz, "C"),
        _atom(serial0 + 2, "C", "ALA", chain, res_seq, bx + 2.5, by + 1.1, bz, "C"),
        _atom(serial0 + 3, "O", "ALA", chain, res_seq, bx + 3.7, by + 1.1, bz, "O"),
        _atom(serial0 + 4, "CB", "ALA", chain, res_seq, bx + 1.5, by - 1.5, bz, "C"),
    ]


@pytest.fixture
def tiny_two_chain_pdb(tmp_path):
    """Two interacting chains A and B with mostly CB-CB contacts < 5 Å."""
    serial = 1
    lines = []
    for i, x in enumerate([0.0, 4.0, 8.0, 12.0, 16.0]):
        lines.extend(_ala_residue(serial, "A", i + 1, (x, 0.0, 0.0)))
        serial += 5
    for i, x in enumerate([0.0, 4.0, 8.0, 12.0, 16.0]):
        lines.extend(_ala_residue(serial, "B", i + 1, (x, 3.5, 0.0)))
        serial += 5
    lines.append("END\n")
    pdb_path = tmp_path / "pdbs" / "test.pdb"
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text("".join(lines))
    return tmp_path


@pytest.fixture
def repo_root_chdir(monkeypatch):
    monkeypatch.chdir(REPO_ROOT)
    return REPO_ROOT
