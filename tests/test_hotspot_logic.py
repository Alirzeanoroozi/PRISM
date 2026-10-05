"""End-to-end-ish tests for hotspot generation on a synthetic mini PDB."""
import json
import os

import pytest
from Bio.PDB import PDBParser

from src import hotspot as hs_mod
from src.hotspot import hotspot_creator, center_of_mass


def _write_test_pdb(pdbs_dir, pdb_id="test"):
    """Write the conftest fixture PDB to a directory so hardcoded paths can find it."""
    pdbs_dir.mkdir(parents=True, exist_ok=True)
    # Reconstruct the mini PDB from conftest's _atom/_ala_residue helpers
    from tests.conftest import _ala_residue
    serial = 1
    lines = []
    for i, x in enumerate([0.0, 4.0, 8.0, 12.0, 16.0]):
        lines.extend(_ala_residue(serial, "A", i + 1, (x, 0.0, 0.0)))
        serial += 5
    for i, x in enumerate([0.0, 4.0, 8.0, 12.0, 16.0]):
        lines.extend(_ala_residue(serial, "B", i + 1, (x, 3.5, 0.0)))
        serial += 5
    lines.append("END\n")
    (pdbs_dir / f"{pdb_id}.pdb").write_text("".join(lines))


def test_hotspot_creator_writes_json(monkeypatch, tmp_path):
    hot_dir = tmp_path / "hotspots"
    hot_dir.mkdir()
    monkeypatch.setattr(hs_mod, "HOTSPOT_DIR", str(hot_dir))
    monkeypatch.setattr(hs_mod, "RELATIVE_ASA_THRESHOLD", 100.0)
    monkeypatch.setattr(hs_mod, "CONTACT_POTENTIAL_THRESHOLD", -1.0)
    monkeypatch.setattr(hs_mod, "get_asa_complex", lambda t, d: {"GLY_1_A": 50.0, "GLY_2_B": 60.0})
    monkeypatch.setattr(hs_mod, "get_contact_potentials", lambda p, c1, c2: {"GLY_1_A": 5.0, "GLY_2_B": 3.0})
    hotspot_creator("testAB")
    out_path = hot_dir / "testAB.json"
    assert out_path.exists()
    loaded = json.loads(out_path.read_text())
    assert "A" in loaded and "B" in loaded


def test_center_of_mass(monkeypatch):
    # fetch_all_atoms_coordinates returns {residue_key: {atom_name: (x,y,z), ...}}
    monkeypatch.setattr(hs_mod, "fetch_all_atoms_coordinates",
                        lambda p, c1, c2: {"GLY_1_A": {"CA": (0.0, 0.0, 0.0), "CB": (1.5, -1.5, 0.0)},
                                            "GLY_2_B": {"CA": (3.5, 0.0, 0.0), "CB": (5.0, -1.5, 0.0)}})
    centers = center_of_mass("test", "A", "B")
    assert centers
    for key, coord in centers.items():
        assert len(coord) == 3
