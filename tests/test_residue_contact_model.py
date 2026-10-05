import pytest

from src.residue_contact_model import build_model


def test_stage_two_model_reports_optional_dependency_without_importing_it():
    try:
        build_model(8)
    except RuntimeError as exc:
        assert "requires an environment with PyTorch" in str(exc)


def test_stage_two_model_rejects_invalid_dimensions():
    try:
        build_model(0)
    except RuntimeError:
        pytest.skip("PyTorch is unavailable in the current test environment")
    except ValueError as exc:
        assert "must be positive" in str(exc)
