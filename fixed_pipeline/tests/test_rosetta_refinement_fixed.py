"""
Tests for the fixed rosetta_refinement.py.
"""
import os
import tempfile
from pathlib import Path

import src.rosetta_refinement as rr


def test_energy_has_error_code():
    """EnergyResult failure must have a non-None error_code."""
    result = rr.EnergyResult.failure("TEST_CODE", "test detail")
    assert result.error_code == "TEST_CODE"
    assert result.error_detail == "test detail"
    assert result.total_score is None
    assert result.interaction_score is None


def test_energy_result_success():
    """EnergyResult success has no error_code."""
    result = rr.EnergyResult(total_score=-12.5, interaction_score="-8.3",
                             structure_path="/path/to/model.pdb")
    assert result.error_code is None
    assert result.total_score == -12.5
    assert result.interaction_score == "-8.3"


def test_no_os_system_in_code():
    """Ensure no os.system() calls exist in executable code (docstrings excluded)."""
    import ast
    source = Path(rr.__file__).read_text()
    tree = ast.parse(source)
    os_system_calls = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr == "system" and isinstance(node.func.value, ast.Name) and node.func.value.id == "os":
                os_system_calls.append(node.lineno)
    assert not os_system_calls, (
        f"os.system() calls found at lines {os_system_calls} in {rr.__file__}"
    )


def test_refinement_error_hierarchy():
    """RefinementError subclasses should have distinct codes."""
    assert issubclass(rr.RosettaBinaryNotFound, rr.RefinementError)
    assert issubclass(rr.RosettaExecutionError, rr.RefinementError)
    assert issubclass(rr.ScoreParsingError, rr.RefinementError)

    err = rr.RosettaBinaryNotFound("docking_protocol", "/usr/bin/docking_protocol")
    assert err.code == "ROSETTA_BINARY_NOT_FOUND"
    assert "docking_protocol" in str(err)


def test_parse_rosetta_score_valid():
    """Parse a valid score.sc file — only data lines, not headers."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sc", delete=False) as f:
        # Rosetta score.sc: first SCORE: line may be column names, second is data.
        # The parser picks the LAST line if the first doesn't parse as float.
        # Simulate a realistic format where the data line is the only SCORE: line.
        f.write(
            "# Rosetta score file\n"
            "# comment\n"
            "SCORE:  -12.50  -8.30  0.0  0.0  0.0  0.0\n"
        )
        score_path = f.name

    try:
        total, interaction = rr._parse_rosetta_score(score_path)
        assert total == -12.5
        assert interaction == 0.0, f"Expected 0.0, got {interaction}"
    finally:
        os.unlink(score_path)


def test_parse_rosetta_score_no_data():
    """Score file with no data lines should raise ScoreParsingError."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".sc", delete=False) as f:
        f.write("# only a comment\n")
        score_path = f.name

    try:
        rr._parse_rosetta_score(score_path)
        assert False, "Should have raised ScoreParsingError"
    except rr.ScoreParsingError:
        pass
    finally:
        os.unlink(score_path)


def test_check_binary_not_found():
    """check_binary raises RosettaBinaryNotFound for missing binaries."""
    try:
        rr._check_binary("/nonexistent/binary_XYZ_123")
        assert False, "Should have raised"
    except rr.RosettaBinaryNotFound:
        pass
