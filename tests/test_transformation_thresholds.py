import importlib


def test_transformation_thresholds_can_be_overridden(monkeypatch):
    monkeypatch.setenv("PRISM_TM_SCORE_THRESHOLD", "0.2")
    monkeypatch.setenv("PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE", "30")
    import src.transformation as transformation
    transformation = importlib.reload(transformation)
    assert transformation.TM_SCORE_THRESHOLD == 0.2
    assert transformation.MINIMUM_RESIDUE_MATCH_PERCENTAGE == 30.0
    monkeypatch.delenv("PRISM_TM_SCORE_THRESHOLD")
    monkeypatch.delenv("PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE")
    importlib.reload(transformation)

def test_multiprot_uses_native_gate_not_tmalign_proxy():
    import src.transformation as transformation

    alignment = {
        "aligner": "MultiProt",
        "match_count": 15,
        "tm_score": 0.0,
        "rmsd": 11.8,
    }
    assert transformation.alignment_score_passes(alignment) is True
    assert transformation.alignment_passes_thresholds("1abc_A", alignment) is True


def test_tmalign_still_requires_tm_score_threshold():
    import src.transformation as transformation

    alignment = {"aligner": "TMalign", "match_count": 15, "tm_score": 0.0}
    assert transformation.alignment_score_passes(alignment) is False
    assert transformation.alignment_passes_thresholds("1abc_A", alignment) is False


def test_environment_selects_common_alignment_gate_mode(monkeypatch):
    from src.transformation_config import TransformationThresholds

    monkeypatch.setenv("PRISM_ALIGNMENT_GATE_MODE", "common_match_coverage")
    thresholds = TransformationThresholds.from_environment()
    assert thresholds.alignment_gate_mode == "common_match_coverage"
