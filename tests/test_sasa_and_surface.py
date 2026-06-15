import os

import pytest

from src import sasa_utils
from src import surface_extract as se
from src.sasa_utils import get_asa_complex, get_asa_flat
from src.surface_extract import extract_surface


def test_get_asa_complex_chain_filter(tiny_two_chain_pdb):
    asa = get_asa_complex("testAB", str(tiny_two_chain_pdb))
    assert set(asa.keys()) <= {"A", "B"}
    assert "A" in asa or "B" in asa
    for chain_dict in asa.values():
        for k in chain_dict.keys():
            assert isinstance(k, int)
        for v in chain_dict.values():
            assert v >= 0


def test_get_asa_complex_subset_chain(tiny_two_chain_pdb):
    asa = get_asa_complex("testA", str(tiny_two_chain_pdb))
    assert set(asa.keys()) == {"A"}


def test_get_asa_flat_format(tiny_two_chain_pdb):
    flat = get_asa_flat("testAB", str(tiny_two_chain_pdb))
    assert flat
    for key in flat.keys():
        parts = key.split("_")
        assert len(parts) == 3
        assert parts[2] in ("A", "B")


def test_sasa_missing_pdb_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        get_asa_complex("missingA", str(tmp_path))


def test_extract_surface(monkeypatch, tiny_two_chain_pdb):
    out_dir = tiny_two_chain_pdb / "surfaces"
    out_dir.mkdir()
    monkeypatch.setattr(se, "SURFACE_EXTRACTION_DIR", str(out_dir))
    monkeypatch.setattr(se, "RSATHRESHOLD", 5.0)
    ok = extract_surface("testAB", pdb_root=str(tiny_two_chain_pdb))
    assert ok is True
    out_pdb = out_dir / "testAB.asa.pdb"
    assert out_pdb.exists()
    content = out_pdb.read_text()
    assert "CA" in content


def test_extract_surface_subset_chain(monkeypatch, tiny_two_chain_pdb):
    out_dir = tiny_two_chain_pdb / "surfaces"
    out_dir.mkdir()
    monkeypatch.setattr(se, "SURFACE_EXTRACTION_DIR", str(out_dir))
    monkeypatch.setattr(se, "RSATHRESHOLD", 5.0)
    ok = extract_surface("testA", pdb_root=str(tiny_two_chain_pdb))
    assert ok is True
    out_pdb = out_dir / "testA.asa.pdb"
    content = out_pdb.read_text()
    for line in content.splitlines():
        if line.startswith("ATOM"):
            assert line[21] == "A"
