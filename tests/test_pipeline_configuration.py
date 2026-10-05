from argparse import Namespace

import pytest

import prism
from src.pipeline_inputs import select_templates
from src.template_filtering import evaluate_protocol_candidate
from src.transformation_config import TransformationThresholds


def test_select_templates_supports_explicit_ids_and_limit(tmp_path):
    template_list = tmp_path / "templates.txt"
    template_list.write_text("1abcAB\n2defCD\n3ghiEF\n", encoding="utf-8")

    assert select_templates(
        ["9zzzXY"], template_list_path=template_list, template_limit=2
    ) == ["1abcAB", "2defCD"]
    assert select_templates(
        ["9zzzXY"], explicit_templates=["2defCD", "1abcAB"], template_limit=1
    ) == ["2defCD"]


def test_select_templates_rejects_invalid_ids_and_limits():
    with pytest.raises(ValueError, match="template ID"):
        select_templates([], explicit_templates=["not-a-template"])
    with pytest.raises(ValueError, match="template-limit"):
        select_templates([], explicit_templates=["1abcAB"], template_limit=0)


def test_parser_exposes_template_and_threshold_controls():
    args = prism.build_parser().parse_args([
        "--inputs-csv", "pairs.csv",
        "--templates", "1abcAB", "2defCD",
        "--minimum-residue-match-count", "12",
        "--minimum-residue-match-percentage", "35",
        "--minimum-hotspot-match-number", "2",
        "--diff-percentage", "25",
        "--template-residue-count", "60",
        "--contact-count-threshold", "4",
        "--clashing-distance", "2.5",
        "--max-clashing-count", "8",
        "--scaffold-threshold", "4.0",
        "--tm-score-threshold", "0.4",
        "--multiprot-minimum-residue-match-count", "9",
        "--multiprot-minimum-residue-match-percentage", "28",
        "--alignment-gate-mode", "common_match_coverage",
    ])

    assert args.inputs_csv == "pairs.csv"
    assert args.templates == ["1abcAB", "2defCD"]
    assert prism.threshold_overrides_from_args(args) == {
        "minimum_residue_match_count": 12,
        "minimum_residue_match_percentage": 35.0,
        "minimum_hotspot_match_number": 2,
        "diff_percentage": 25.0,
        "template_residue_count": 60.0,
        "contact_count_threshold": 4,
        "clashing_distance": 2.5,
        "max_clashing_count": 8,
        "scaffold_threshold": 4.0,
        "tm_score_threshold": 0.4,
        "multiprot_minimum_residue_match_count": 9,
        "multiprot_minimum_residue_match_percentage": 28.0,
        "alignment_gate_mode": "common_match_coverage",
    }


def test_thresholds_preserve_current_defaults_and_accept_overrides():
    defaults = TransformationThresholds.from_environment({})
    assert defaults.minimum_residue_match_count == 15
    assert defaults.minimum_residue_match_percentage == 50.0
    assert defaults.scaffold_threshold == 5.0
    assert defaults.tm_score_threshold == 0.5
    assert defaults.alignment_gate_mode == "native"

    configured = defaults.with_overrides({
        "minimum_residue_match_count": 12,
        "max_clashing_count": 8,
    })
    assert configured.minimum_residue_match_count == 12
    assert configured.max_clashing_count == 8
    assert configured.tm_score_threshold == defaults.tm_score_threshold


def test_explicit_thresholds_reach_transformation_gates(monkeypatch):
    import src.transformation as transformation

    monkeypatch.setattr(transformation, "template_size", {"1abc_A": 40})
    alignment = {
        "aligner": "TMalign",
        "match_count": 12,
        "tm_score": 0.4,
        "match_dict": {},
    }
    thresholds = TransformationThresholds(
        minimum_residue_match_count=12,
        minimum_residue_match_percentage=20.0,
        tm_score_threshold=0.4,
    )

    assert transformation.alignment_passes_thresholds(
        "1abc_A", alignment, thresholds=thresholds
    ) is True


def test_multiprot_uses_its_native_count_and_coverage_contract(monkeypatch):
    import src.transformation as transformation

    monkeypatch.setattr(transformation, "template_size", {"1abc_A": 100})
    alignment = {
        "aligner": "MultiProt",
        "match_count": 10,
        "tm_score": 0.0,
        "match_dict": {f"A.A.{index}": f"Q.Q.{index}" for index in range(10)},
    }

    assert transformation.alignment_passes_thresholds("1abc_A", alignment) is False
    assert transformation.alignment_passes_thresholds(
        "1abc_A", alignment,
        thresholds=TransformationThresholds(
            multiprot_minimum_residue_match_percentage=5.0
        ),
    ) is True
    assert transformation.alignment_passes_thresholds(
        "1abc_A", alignment,
        thresholds=TransformationThresholds(multiprot_minimum_residue_match_count=11),
    ) is False


def test_common_match_coverage_contract_is_identical_for_tmalign_and_multiprot(monkeypatch):
    import src.transformation as transformation

    monkeypatch.setattr(transformation, "template_size", {"1abc_A": 40})
    thresholds = TransformationThresholds(alignment_gate_mode="common_match_coverage")
    tmalign = {
        "aligner": "TMalign",
        "match_count": 20,
        "tm_score": 0.0,
        "match_dict": {},
    }
    multiprot = {
        "aligner": "MultiProt",
        "match_count": 20,
        "tm_score": 0.0,
        "rmsd": 99.0,
        "match_dict": {f"A.A.{index}": f"Q.Q.{index}" for index in range(20)},
    }

    assert transformation.alignment_passes_thresholds(
        "1abc_A", tmalign, thresholds=thresholds
    ) is True
    assert transformation.alignment_passes_thresholds(
        "1abc_A", multiprot, thresholds=thresholds
    ) is True

    below_count = dict(tmalign, match_count=14)
    assert transformation.alignment_passes_thresholds(
        "1abc_A", below_count, thresholds=thresholds
    ) is False


def test_common_match_coverage_uses_same_large_interface_boundary(monkeypatch):
    import src.transformation as transformation

    monkeypatch.setattr(transformation, "template_size", {"1abc_A": 100})
    thresholds = TransformationThresholds(alignment_gate_mode="common_match_coverage")
    at_boundary = {
        "aligner": "TMalign",
        "match_count": 30,
        "tm_score": 0.0,
        "match_dict": {},
    }
    assert transformation.alignment_passes_thresholds(
        "1abc_A", at_boundary, thresholds=thresholds
    ) is True


def test_hotspot_threshold_reaches_published_protocol_gate():
    left = {"A.A.1": "Q.A.1", "A.A.2": "Q.A.2"}
    right = {"B.B.1": "R.B.1", "B.B.2": "R.B.2"}
    contacts = [("A.A.1", "B.B.1")]

    result = evaluate_protocol_candidate(
        left, right, [(1, "A"), (2, "A")], [(1, "B")], contacts,
        minimum_contacts=1, minimum_hotspots=2,
    )

    assert result.reason == "hotspot_threshold_failed"


def test_main_forwards_explicit_input_csv_and_thresholds(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(prism, "pdb_downloader", lambda path=None: ([], []))
    monkeypatch.setattr(prism, "run_analysis", lambda: ([], 0, 0))
    monkeypatch.setattr(prism, "template_generator", lambda: [])
    monkeypatch.setattr(prism, "extract_surfaces", lambda targets: None)
    monkeypatch.setattr(prism, "align", lambda targets, templates, **kwargs: None)
    captured = {}

    def fake_transformer(templates, alignment_dir, inputs_csv, thresholds):
        captured["inputs_csv"] = inputs_csv
        captured["thresholds"] = thresholds
        return []

    monkeypatch.setattr(prism, "transformer", fake_transformer)

    args = Namespace(
        surface_backend="naccess", freesasa_python=None, inputs_csv="pairs.csv",
        generate_templates=True, template_limit=None, templates=["1abcAB"],
        template_list=None, aligner="tmalign", refiner="external_rosetta",
        gtalign_path="gtalign", gtalign_dev_min_length=3, gtalign_pre_score=0.0,
        gtalign_speed=0, gtalign_refinement=3, orientation="native", refine=False,
        minimum_residue_match_count=12, minimum_residue_match_percentage=None,
        minimum_hotspot_match_number=None, diff_percentage=None,
        template_residue_count=None, contact_count_threshold=None,
        clashing_distance=None, max_clashing_count=None, tm_score_threshold=None,
        multiprot_minimum_residue_match_count=None,
        multiprot_minimum_residue_match_percentage=None,
    )

    prism.main(args)

    assert captured["inputs_csv"] == "pairs.csv"
    assert captured["thresholds"].minimum_residue_match_count == 12


def test_main_forwards_explicit_scaffold_threshold(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(prism, "pdb_downloader", lambda path=None: ([], []))
    monkeypatch.setattr(prism, "extract_surfaces", lambda targets, **kwargs: captured.update(kwargs))
    monkeypatch.setattr(prism, "align", lambda targets, templates, **kwargs: None)
    monkeypatch.setattr(prism, "transformer", lambda templates, alignment_dir, inputs_csv=None, thresholds=None: [])
    captured = {}

    args = prism.build_parser().parse_args([
        "--inputs-csv", "pairs.csv",
        "--templates", "1abcAB",
        "--scaffold-threshold", "4.25",
        "--no-refine",
    ])
    prism.main(args)

    assert captured["scaffold_threshold"] == 4.25


def test_main_explicit_template_list_does_not_require_default_manifest(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    template_list = tmp_path / "subset.txt"
    template_list.write_text("1abcAB\n", encoding="utf-8")
    monkeypatch.setattr(prism, "pdb_downloader", lambda path=None: ([], []))
    monkeypatch.setattr(prism, "extract_surfaces", lambda targets: None)
    monkeypatch.setattr(prism, "align", lambda targets, templates, **kwargs: None)
    captured = {}

    def fake_transformer(templates, alignment_dir, inputs_csv=None):
        captured["templates"] = templates
        return []

    monkeypatch.setattr(prism, "transformer", fake_transformer)
    prism.main(Namespace(
        surface_backend="naccess", freesasa_python=None, inputs_csv=None,
        generate_templates=False, template_limit=None, templates=None,
        template_list=str(template_list), aligner="tmalign",
        refiner="external_rosetta", gtalign_path="gtalign",
        gtalign_dev_min_length=3, gtalign_pre_score=0.0, gtalign_speed=0,
        gtalign_refinement=3, orientation="native", refine=False,
    ))

    assert captured["templates"] == ["1abcAB"]


def test_main_always_forwards_literal_default_input_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PRISM_INPUTS_CSV", "environment.csv")
    monkeypatch.setattr(prism, "pdb_downloader", lambda path=None: ([], []))
    monkeypatch.setattr(prism, "run_analysis", lambda: ([], 0, 0))
    monkeypatch.setattr(prism, "template_generator", lambda: [])
    monkeypatch.setattr(prism, "extract_surfaces", lambda targets: None)
    monkeypatch.setattr(prism, "align", lambda targets, templates, **kwargs: None)
    captured = {}

    def fake_transformer(templates, alignment_dir, inputs_csv):
        captured["inputs_csv"] = inputs_csv
        return []

    monkeypatch.setattr(prism, "transformer", fake_transformer)
    prism.main(Namespace(
        surface_backend="naccess", freesasa_python=None, inputs_csv="inputs.csv",
        generate_templates=True, template_limit=None, templates=["1abcAB"],
        template_list=None, aligner="tmalign", refiner="external_rosetta",
        gtalign_path="gtalign", gtalign_dev_min_length=3, gtalign_pre_score=0.0,
        gtalign_speed=0, gtalign_refinement=3, orientation="native", refine=False,
    ))

    assert captured["inputs_csv"] == "inputs.csv"
