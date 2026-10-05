"""
Tests for the fixed transformation.py.
"""
import os
import tempfile

import src.transformation as tf


def test_transform_result_ok():
    """TransformResult.ok() should be truthy with no error."""
    result = tf.TransformResult.ok()
    assert result.success is True
    assert result.error_code is None
    assert bool(result) is True


def test_transform_result_failure():
    """TransformResult.failure() should be falsy with error_code."""
    result = tf.TransformResult.failure("INPUT_MISSING", "file not found")
    assert result.success is False
    assert result.error_code == "INPUT_MISSING"
    assert bool(result) is False


def test_transform_missing_file():
    """apply_tm_transform should return INPUT_MISSING for nonexistent file."""
    result = tf.apply_tm_transform(
        "/nonexistent/file.pdb",
        "/tmp/output.pdb",
        [0.0, 0.0, 0.0],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    )
    assert result.success is False
    assert result.error_code == "INPUT_MISSING"


def test_apply_transform_valid():
    """apply_tm_transform should succeed on a valid ATOM PDB."""
    pdb_content = (
        "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.500   2.500   3.500  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.000   3.000   4.000  1.00  0.00           C\n"
        "END\n"
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        f.write(pdb_content)
        input_path = f.name

    output_path = input_path + ".transformed.pdb"

    try:
        result = tf.apply_tm_transform(
            input_path, output_path,
            [0.0, 0.0, 0.0],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        )
        assert result.success is True, f"Transform failed: {result.error_code}"
        assert os.path.exists(output_path)

        with open(output_path) as f:
            content = f.read()
            assert "1.000" in content
            assert "2.000" in content
            assert "3.000" in content
    finally:
        for p in [input_path, output_path]:
            if os.path.exists(p):
                os.unlink(p)


def test_translation_applied():
    """apply_tm_transform should translate coordinates."""
    pdb_content = (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        f.write(pdb_content)
        input_path = f.name

    output_path = input_path + ".transformed.pdb"

    try:
        result = tf.apply_tm_transform(
            input_path, output_path,
            [10.0, 20.0, 30.0],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        )
        assert result.success is True

        with open(output_path) as f:
            content = f.read()
            assert "10.000" in content
            assert "20.000" in content
            assert "30.000" in content
    finally:
        for p in [input_path, output_path]:
            if os.path.exists(p):
                os.unlink(p)


def test_nan_detected():
    """NaN coordinates should be detected (logged) and handled gracefully."""
    import logging
    pdb_content = (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
        f.write(pdb_content)
        input_path = f.name

    output_path = input_path + ".transformed.pdb"

    try:
        result = tf.apply_tm_transform(
            input_path, output_path,
            [float("nan"), 0.0, 0.0],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        )
        # NaN coords should cause failure with NAN_COORDS code
        assert result.success is False
        assert result.error_code == "NAN_COORDS"
    finally:
        for p in [input_path, output_path]:
            if os.path.exists(p):
                os.unlink(p)
