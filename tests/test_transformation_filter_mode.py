import importlib

import pytest


def test_unknown_filter_mode_fails_closed(monkeypatch):
    monkeypatch.setenv("PRISM_FILTER_MODE", "unknown")
    with pytest.raises(ValueError, match="PRISM_FILTER_MODE"):
        import src.transformation as transformation
        importlib.reload(transformation)


def test_geometry_only_mode_is_explicit_and_published_mode_requires_assets(monkeypatch):
    monkeypatch.setenv("PRISM_FILTER_MODE", "geometry_only_experimental")
    import src.transformation as transformation
    importlib.reload(transformation)
    assert transformation.FILTER_MODE == "geometry_only_experimental"
    assert transformation.hotspot_analysis({}) is True
    monkeypatch.setenv("PRISM_FILTER_MODE", "published_protocol")
    importlib.reload(transformation)
    assert transformation.hotspot_analysis({}) is False
    monkeypatch.delenv("PRISM_FILTER_MODE", raising=False)
    importlib.reload(transformation)


def test_published_thresholds_use_protocol_hotspots(monkeypatch):
    monkeypatch.setenv("PRISM_FILTER_MODE", "published_protocol")
    monkeypatch.setenv("PRISM_MINIMUM_RESIDUE_MATCH_COUNT", "1")
    monkeypatch.setenv("PRISM_TM_SCORE_THRESHOLD", "0.1")
    import src.transformation as transformation
    importlib.reload(transformation)
    alignment = {
        "match_count": 1,
        "tm_score": 0.2,
        "match_dict": {"A.K.1": "A.K.7"},
    }
    assert transformation.alignment_passes_thresholds(
        "1abcA", alignment, protocol_hotspots=["A.K.1"]
    )
    monkeypatch.delenv("PRISM_FILTER_MODE", raising=False)
    monkeypatch.delenv("PRISM_MINIMUM_RESIDUE_MATCH_COUNT", raising=False)
    monkeypatch.delenv("PRISM_TM_SCORE_THRESHOLD", raising=False)
    importlib.reload(transformation)


def test_protocol_assets_are_selected_per_chain_and_orientation(monkeypatch):
    monkeypatch.setenv("PRISM_FILTER_MODE", "published_protocol")
    import src.transformation as transformation
    importlib.reload(transformation)
    assets = {
        "hotspots": [["A", "A", "1"], ["B", "B", "2"]],
        "hotspots_by_chain": {
            "A": [["A", "A", "1"]],
            "B": [["B", "B", "2"]],
        },
        "contacts": [["A..1", "B..2"]],
    }
    assert transformation._protocol_hotspots(assets, "A") == [["A", "A", "1"]]
    assert transformation._protocol_hotspots(assets, "B") == [["B", "B", "2"]]
    assert transformation._protocol_contacts(assets, "A", "B") == [["A..1", "B..2"]]
    assert transformation._protocol_contacts(assets, "B", "A") == [["B..2", "A..1"]]
    monkeypatch.delenv("PRISM_FILTER_MODE", raising=False)
    importlib.reload(transformation)
