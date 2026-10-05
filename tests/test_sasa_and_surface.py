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


def test_extract_surface_returns_dict(monkeypatch, tmp_path):
    """extract_surface returns a dict when asa_complex is empty."""
    monkeypatch.setattr(se, "RSATHRESHOLD", 999999.0)  # threshold above any possible RSA
    monkeypatch.setattr(se, "SURFACE_EXTRACTION_DIR", str(tmp_path / "surfaces"))
    monkeypatch.setattr(se, "get_asa_complex_target", lambda p, d: {})
    result = extract_surface("testAB")
    assert result == {}


def test_scaffold_threshold_resolves_explicit_and_environment_values(monkeypatch):
    monkeypatch.setenv("PRISM_SCFF_THRESHOLD", "4.25")
    assert se._resolve_scaffold_threshold() == pytest.approx(4.25)
    assert se._resolve_scaffold_threshold(3.5) == pytest.approx(3.5)


@pytest.mark.skip(reason="Requires writeable processed/pdbs/ directory with test PDB")
def test_extract_surface_with_residues(monkeypatch):
    """extract_surface returns CA coord dict when residues above threshold."""
    monkeypatch.setattr(se, "RSATHRESHOLD", 5.0)
    monkeypatch.setattr(se, "get_asa_complex_target", lambda p, d: {"GLY_1_A": 50.0, "GLY_2_B": 60.0})
    monkeypatch.setattr(se, "SURFACE_EXTRACTION_DIR", "/tmp/surfaces")
    result = extract_surface("testAB")
    assert isinstance(result, dict)
