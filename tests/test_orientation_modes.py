import json

import pytest

import prism
import src.transformation as transformation
from src.transformation import select_orientations


def test_native_orientation_mode_preserves_both_implicit_branches():
    args = prism.build_parser().parse_args([])

    assert args.orientation == "native"
    assert select_orientations(args.orientation) == ("o1", "o2")


def test_no_orientation_option_is_native_even_if_environment_requests_fixed_mode(monkeypatch):
    monkeypatch.setenv("PRISM_ORIENTATION", "o1")

    args = prism.build_parser().parse_args([])

    assert args.orientation == "native"


@pytest.mark.parametrize(
    ("mode", "expected"),
    [("o1", ("o1",)), ("o2", ("o2",))],
)
def test_fixed_orientation_mode_selects_one_branch(mode, expected):
    args = prism.build_parser().parse_args(["--orientation", mode])

    assert args.orientation == mode
    assert select_orientations(args.orientation) == expected


def test_unknown_orientation_mode_fails_closed():
    with pytest.raises(ValueError, match="orientation"):
        select_orientations("invalid")


@pytest.mark.parametrize(
    ("mode", "expected_suffixes"),
    [
        ("native", {"o1", "o2"}),
        ("o1", {"o1"}),
        ("o2", {"o2"}),
    ],
)
def test_transformation_mode_controls_generated_orientation_suffixes(
    monkeypatch, mode, expected_suffixes
):
    alignment = {
        "match_count": 30,
        "tm_score": 0.6,
        "match_dict": {"A.K.1": "A.K.1"},
        "translation": [0.0, 0.0, 0.0],
        "rotation_mat": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    }
    generated_suffixes = []

    monkeypatch.setattr(transformation, "template_size", {
        "1abc_A": 40,
        "1abc_B": 40,
    })
    monkeypatch.setattr(transformation, "passed_pairs", [])
    monkeypatch.setattr(
        transformation,
        "load_alignment",
        lambda *args, **kwargs: dict(alignment),
    )
    monkeypatch.setattr(
        transformation,
        "create_transformed_pair",
        lambda *args: generated_suffixes.append(args[-1]) or "generated",
    )

    transformation.process_pair_for_template(
        "1abcAB", "A", "B", "left", "right", orientation=mode,
    )

    assert set(generated_suffixes) == expected_suffixes


def test_audit_distinguishes_alignment_threshold_rejection(tmp_path, monkeypatch):
    alignment = {
        "aligner": "TMalign",
        "match_count": 1,
        "tm_score": 0.1,
        "match_dict": {},
    }
    monkeypatch.setattr(transformation, "template_size", {
        "1abc_A": 40,
        "1abc_B": 40,
    })
    monkeypatch.setattr(transformation, "passed_pairs", [])
    monkeypatch.setattr(
        transformation, "load_alignment", lambda *args, **kwargs: dict(alignment)
    )
    audit_path = tmp_path / "audit.jsonl"

    transformation.process_pair_for_template(
        "1abcAB", "A", "B", "left", "right",
        orientation="o1", audit_path=str(audit_path),
    )

    row = json.loads(audit_path.read_text())
    assert row["status"] == "alignment_threshold_rejected"
    assert row["error_reason"] == "alignment_threshold_rejected"
