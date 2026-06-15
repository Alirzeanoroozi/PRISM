"""End-to-end-ish tests for hotspot generation on a synthetic mini PDB."""
import json
import os

import pytest

from src import hotspot as hs_mod
from src.hotspot import hotspot_creator, _center_of_mass


def test_hotspot_creator_writes_json(monkeypatch, tiny_two_chain_pdb):
    out_dir = tiny_two_chain_pdb / "hotspots"
    out_dir.mkdir()
    monkeypatch.setattr(hs_mod, "HOTSPOT_DIR", str(out_dir))
    monkeypatch.setattr(hs_mod, "RELATIVE_ASA_THRESHOLD", 100.0)
    monkeypatch.setattr(hs_mod, "CONTACT_POTENTIAL_THRESHOLD", -1.0)
    result = hotspot_creator("testAB", templates_root=str(tiny_two_chain_pdb))
    assert set(result.keys()) == {"A", "B"}
    out_path = out_dir / "testAB.json"
    assert out_path.exists()
    loaded = json.loads(out_path.read_text())
    assert "A" in loaded and "B" in loaded


def test_center_of_mass(tiny_two_chain_pdb):
    centers = _center_of_mass("testAB", str(tiny_two_chain_pdb))
    assert centers
    for key, coord in centers.items():
        assert len(coord) == 3
        assert key.endswith("_A") or key.endswith("_B")
