import pytest
import json
from argparse import Namespace

import prism
from prism import parse_bool


def test_parse_bool_handles_false_string():
    assert parse_bool("false") is False
    assert parse_bool("true") is True


def test_parse_bool_rejects_unknown_value():
    with pytest.raises(Exception):
        parse_bool("maybe")


def test_main_records_terminal_stage_events(tmp_path, monkeypatch):
    stage_path = tmp_path / "status" / "stages.jsonl"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PRISM_STAGE_STATUS_PATH", str(stage_path))
    monkeypatch.setattr(prism, "pdb_downloader", lambda: ([], []))
    monkeypatch.setattr(prism, "run_analysis", lambda: ([], 0, 0))
    monkeypatch.setattr(prism, "template_generator", lambda: [])
    monkeypatch.setattr(prism, "extract_surfaces", lambda targets: None)
    monkeypatch.setattr(prism, "align", lambda targets, templates, **kwargs: None)
    monkeypatch.setattr(
        prism, "transformer", lambda templates, alignment_dir, inputs_csv=None: []
    )
    monkeypatch.setattr(prism, "refiner", lambda pairs: None)

    prism.main(
        Namespace(
            surface_backend="naccess",
            freesasa_python=None,
            generate_templates=True,
            template_limit=None,
            aligner="tmalign",
            refiner="external_rosetta",
            gtalign_path="gtalign",
            gtalign_dev_min_length=3,
            gtalign_pre_score=0.0,
            gtalign_speed=0,
            gtalign_refinement=3,
        )
    )

    events = [json.loads(line) for line in stage_path.read_text().splitlines()]
    terminal = {"completed", "skipped"}
    terminal_stages = {
        event["stage"] for event in events if event["event"] in terminal
    }
    assert {"input", "alignment", "transformation", "refinement"}.issubset(terminal_stages)
    assert any(
        event["stage"] == "refinement"
        and event["event"] == "skipped"
        and "no candidates" in event["detail"]
        for event in events
    )
    assert all(event["timestamp"] for event in events)


def test_ranked_main_allocates_fresh_audit_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("PRISM_CANDIDATE_AUDIT_PATH", raising=False)
    monkeypatch.setattr(prism, "pdb_downloader", lambda: ([], []))
    monkeypatch.setattr(prism, "run_analysis", lambda: ([], 0, 0))
    monkeypatch.setattr(prism, "template_generator", lambda: [])
    monkeypatch.setattr(prism, "extract_surfaces", lambda targets: None)
    monkeypatch.setattr(prism, "align", lambda targets, templates, **kwargs: None)
    captured = {}
    def fake_transformer(templates, alignment_dir, audit_path, inputs_csv=None):
        captured["audit_path"] = audit_path
        return []

    monkeypatch.setattr(prism, "transformer", fake_transformer)
    monkeypatch.setattr(prism, "select_top_candidates", lambda pairs, **kwargs: pairs)
    monkeypatch.setattr(prism, "refiner", lambda pairs: None)

    prism.main(Namespace(
        surface_backend="naccess", freesasa_python=None, generate_templates=True,
        template_limit=None, aligner="tmalign", refiner="external_rosetta",
        gtalign_path="gtalign", gtalign_dev_min_length=3, gtalign_pre_score=0.0,
        gtalign_speed=0, gtalign_refinement=3, rank=True, top_k=5,
        rank_min_score=0.0, candidate_audit_path=None,
    ))

    assert captured["audit_path"].startswith("processed/candidate_audit/")
    assert captured["audit_path"].endswith(".jsonl")


def test_main_forwards_fixed_orientation_to_transformation(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(prism, "pdb_downloader", lambda: ([], []))
    monkeypatch.setattr(prism, "run_analysis", lambda: ([], 0, 0))
    monkeypatch.setattr(prism, "template_generator", lambda: [])
    monkeypatch.setattr(prism, "extract_surfaces", lambda targets: None)
    monkeypatch.setattr(prism, "align", lambda targets, templates, **kwargs: None)
    captured = {}

    def fake_transformer(templates, alignment_dir, orientation, inputs_csv=None):
        captured["orientation"] = orientation
        return []

    monkeypatch.setattr(prism, "transformer", fake_transformer)

    prism.main(Namespace(
        surface_backend="naccess", freesasa_python=None, generate_templates=True,
        template_limit=None, aligner="tmalign", refiner="external_rosetta",
        gtalign_path="gtalign", gtalign_dev_min_length=3, gtalign_pre_score=0.0,
        gtalign_speed=0, gtalign_refinement=3, orientation="o2", refine=False,
    ))

    assert captured["orientation"] == "o2"


def test_no_refine_records_terminal_refinement_event(tmp_path, monkeypatch):
    stage_path = tmp_path / "status" / "stages.jsonl"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PRISM_STAGE_STATUS_PATH", str(stage_path))
    monkeypatch.setattr(prism, "pdb_downloader", lambda: ([], []))
    monkeypatch.setattr(prism, "run_analysis", lambda: ([], 0, 0))
    monkeypatch.setattr(prism, "template_generator", lambda: [])
    monkeypatch.setattr(prism, "extract_surfaces", lambda targets: None)
    monkeypatch.setattr(prism, "align", lambda targets, templates, **kwargs: None)
    monkeypatch.setattr(
        prism, "transformer", lambda templates, alignment_dir, inputs_csv=None: []
    )

    prism.main(Namespace(
        surface_backend="naccess", freesasa_python=None, generate_templates=True,
        template_limit=None, aligner="tmalign", refiner="external_rosetta",
        gtalign_path="gtalign", gtalign_dev_min_length=3, gtalign_pre_score=0.0,
        gtalign_speed=0, gtalign_refinement=3, refine=False,
    ))

    events = [json.loads(line) for line in stage_path.read_text().splitlines()]
    assert any(
        event["stage"] == "refinement"
        and event["event"] == "skipped"
        and "disabled" in event["detail"]
        for event in events
    )
