import importlib


def test_clash_threshold_can_be_overridden(monkeypatch):
    monkeypatch.setenv("PRISM_MAX_CLASHING_COUNT", "100")
    import src.transformation as transformation
    transformation = importlib.reload(transformation)
    assert transformation.MAX_CLASHING_COUNT == 100
    monkeypatch.delenv("PRISM_MAX_CLASHING_COUNT")
    importlib.reload(transformation)
