import pytest

from src.naccess_utils import _relative_areas, _surface_backend


def test_surface_backend_defaults_to_naccess(monkeypatch):
    monkeypatch.delenv("PRISM_SURFACE_BACKEND", raising=False)
    assert _surface_backend() == "naccess"


def test_surface_backend_rejects_unknown_backend(monkeypatch):
    monkeypatch.setenv("PRISM_SURFACE_BACKEND", "unknown")
    with pytest.raises(ValueError, match="PRISM_SURFACE_BACKEND"):
        _surface_backend()


def test_freesasa_areas_use_prism_relative_asa_convention():
    areas = {
        "A": {"10": {"residue_name": "ALA", "total": 107.95}},
        "B": {"11": {"residue_name": "GLY", "total": 40.05}},
    }
    result = _relative_areas(areas, ("A", "B"))
    assert result["ALA_10_A"] == pytest.approx(100.0)
    assert result["GLY_11_B"] == pytest.approx(50.0)
