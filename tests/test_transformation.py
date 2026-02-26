"""Tests for src/transformation.py pipeline functions."""
import os
import json
import tempfile
import pytest

from src.transformation import (
    alignment_passes_thresholds,
    apply_tm_transform,
    pair_has_acceptable_clashes,
    MINIMUM_RESIDUE_MATCH_COUNT,
    MINIMUM_RESIDUE_MATCH_PERCENTAGE,
    TM_SCORE_THRESHOLD,
    CLASHING_DISTANCE,
    MAX_CLASHING_COUNT,
    TEMPLATE_RESIDUE_COUNT,
    DIFF_PERCENTAGE,
)
import src.transformation as transformation_module


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_alignment(match_count=20, tm_score=0.6, match_dict=None):
    return {
        "match_count": match_count,
        "tm_score": tm_score,
        "match_dict": match_dict or {},
    }


def _set_template_size(key, size):
    transformation_module.template_size[key] = size


def _write_mini_pdb(path, ca_coords):
    """Write a tiny CA-only PDB with the given list of (x,y,z) tuples."""
    with open(path, "w") as f:
        for i, (x, y, z) in enumerate(ca_coords, start=1):
            f.write(
                f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"
            )
        f.write("END\n")


# ---------------------------------------------------------------------------
# alignment_passes_thresholds
# ---------------------------------------------------------------------------

class TestAlignmentPassesThresholds:
    def setup_method(self):
        # Reset module-level dict before each test
        transformation_module.template_size.clear()

    def test_passes_with_good_alignment(self):
        key = "tmpl_A"
        _set_template_size(key, 30)
        align = _make_alignment(match_count=20, tm_score=0.7)
        assert alignment_passes_thresholds(key, align) is True

    def test_fails_low_match_count(self):
        key = "tmpl_A"
        _set_template_size(key, 30)
        align = _make_alignment(match_count=5, tm_score=0.7)
        assert alignment_passes_thresholds(key, align) is False

    def test_fails_low_tm_score(self):
        key = "tmpl_A"
        _set_template_size(key, 30)
        align = _make_alignment(match_count=20, tm_score=0.1)
        assert alignment_passes_thresholds(key, align) is False

    def test_fails_low_match_percentage_small_template(self):
        key = "tmpl_A"
        _set_template_size(key, 40)  # <= TEMPLATE_RESIDUE_COUNT (50)
        # match_count=16 / 40 * 100 = 40% < MINIMUM_RESIDUE_MATCH_PERCENTAGE (50%)
        align = _make_alignment(match_count=16, tm_score=0.7)
        assert alignment_passes_thresholds(key, align) is False

    def test_passes_relaxed_percentage_large_template(self):
        key = "tmpl_A"
        _set_template_size(key, 60)  # > TEMPLATE_RESIDUE_COUNT (50)
        # match_count=19 / 60 * 100 = 31.67% > (50 - 20) = 30%
        align = _make_alignment(match_count=19, tm_score=0.7)
        assert alignment_passes_thresholds(key, align) is True

    def test_no_template_size_falls_back_to_count_and_tm(self):
        key = "nonexistent"
        # key not in template_size => protein_size == 0 -> fallback
        align = _make_alignment(match_count=20, tm_score=0.7)
        assert alignment_passes_thresholds(key, align) is True

    def test_no_template_size_fails_fallback(self):
        key = "nonexistent"
        align = _make_alignment(match_count=5, tm_score=0.3)
        assert alignment_passes_thresholds(key, align) is False


# ---------------------------------------------------------------------------
# apply_tm_transform
# ---------------------------------------------------------------------------

class TestApplyTmTransform:
    def test_identity_transform_leaves_coords_unchanged(self):
        translation = [0.0, 0.0, 0.0]
        rotation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            in_pdb = os.path.join(tmpdir, "in.pdb")
            out_pdb = os.path.join(tmpdir, "out.pdb")
            _write_mini_pdb(in_pdb, [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)])
            apply_tm_transform(in_pdb, out_pdb, translation, rotation)
            with open(out_pdb) as f:
                lines = [l for l in f if l.startswith("ATOM")]
            assert len(lines) == 2
            x = float(lines[0][30:38])
            y = float(lines[0][38:46])
            z = float(lines[0][46:54])
            assert x == pytest.approx(1.0, abs=0.01)
            assert y == pytest.approx(2.0, abs=0.01)
            assert z == pytest.approx(3.0, abs=0.01)

    def test_translation_applied(self):
        translation = [10.0, 20.0, 30.0]
        rotation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            in_pdb = os.path.join(tmpdir, "in.pdb")
            out_pdb = os.path.join(tmpdir, "out.pdb")
            _write_mini_pdb(in_pdb, [(0.0, 0.0, 0.0)])
            apply_tm_transform(in_pdb, out_pdb, translation, rotation)
            with open(out_pdb) as f:
                lines = [l for l in f if l.startswith("ATOM")]
            x = float(lines[0][30:38])
            y = float(lines[0][38:46])
            z = float(lines[0][46:54])
            assert x == pytest.approx(10.0, abs=0.01)
            assert y == pytest.approx(20.0, abs=0.01)
            assert z == pytest.approx(30.0, abs=0.01)

    def test_missing_input_does_not_raise(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_pdb = os.path.join(tmpdir, "out.pdb")
            # Should print an error but not raise
            apply_tm_transform("/nonexistent/path.pdb", out_pdb, [0, 0, 0],
                                [[1, 0, 0], [0, 1, 0], [0, 0, 1]])


# ---------------------------------------------------------------------------
# pair_has_acceptable_clashes
# ---------------------------------------------------------------------------

class TestPairHasAcceptableClashes:
    def test_no_clash_distant_structures(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            left = os.path.join(tmpdir, "left.pdb")
            right = os.path.join(tmpdir, "right.pdb")
            # Place structures far apart (1000 Å away)
            _write_mini_pdb(left, [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
            _write_mini_pdb(right, [(1000.0, 0.0, 0.0), (1001.0, 0.0, 0.0)])
            assert pair_has_acceptable_clashes(left, right) is True

    def test_clash_detected_overlapping_structures(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            left = os.path.join(tmpdir, "left.pdb")
            right = os.path.join(tmpdir, "right.pdb")
            # Place many atoms at the exact same location to force many clashes
            coords = [(float(i), 0.0, 0.0) for i in range(MAX_CLASHING_COUNT + 2)]
            _write_mini_pdb(left, coords)
            _write_mini_pdb(right, coords)  # identical coords -> all pairs clash
            assert pair_has_acceptable_clashes(left, right) is False

    def test_missing_file_returns_true(self):
        # read_ca_coordinates returns [] for missing files -> no clashes
        with tempfile.TemporaryDirectory() as tmpdir:
            left = os.path.join(tmpdir, "left.pdb")
            _write_mini_pdb(left, [(0.0, 0.0, 0.0)])
            result = pair_has_acceptable_clashes(left, "/nonexistent/right.pdb")
            assert result is True


# ---------------------------------------------------------------------------
# Constants sanity checks
# ---------------------------------------------------------------------------

class TestTransformationConstants:
    def test_minimum_match_count_positive(self):
        assert MINIMUM_RESIDUE_MATCH_COUNT > 0

    def test_tm_score_threshold_in_range(self):
        assert 0.0 < TM_SCORE_THRESHOLD < 1.0

    def test_clashing_distance_positive(self):
        assert CLASHING_DISTANCE > 0

    def test_max_clashing_count_positive(self):
        assert MAX_CLASHING_COUNT > 0

    def test_match_percentage_in_range(self):
        assert 0 < MINIMUM_RESIDUE_MATCH_PERCENTAGE <= 100

    def test_diff_percentage_non_negative(self):
        assert DIFF_PERCENTAGE >= 0
