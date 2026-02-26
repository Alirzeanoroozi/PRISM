"""Tests for src/utils.py utility functions."""
import math
import pytest

from src.utils import (
    distance_calculator,
    three2one,
    vdw_radii,
    vdw_radii_extended,
    standard_data,
    STANDARD_AA,
    PAIR_POTENTIAL,
    ATOM_DICT,
    DEFAULT_VDW,
)


class TestDistanceCalculator:
    def test_same_point_is_zero(self):
        assert distance_calculator([0, 0, 0], [0, 0, 0]) == 0.0

    def test_unit_x_axis(self):
        assert distance_calculator([0, 0, 0], [1, 0, 0]) == pytest.approx(1.0)

    def test_known_3d_distance(self):
        # sqrt(1^2 + 2^2 + 2^2) = sqrt(9) = 3
        assert distance_calculator([0, 0, 0], [1, 2, 2]) == pytest.approx(3.0)

    def test_negative_coords(self):
        # distance between (-1,-1,-1) and (1,1,1) = sqrt(12)
        assert distance_calculator([-1, -1, -1], [1, 1, 1]) == pytest.approx(math.sqrt(12))

    def test_symmetric(self):
        a, b = [1.5, 2.5, 3.5], [4.0, 5.0, 6.0]
        assert distance_calculator(a, b) == pytest.approx(distance_calculator(b, a))

    def test_string_coords(self):
        # function uses float() internally
        assert distance_calculator(["0", "0", "0"], ["3", "4", "0"]) == pytest.approx(5.0)


class TestThreeToOne:
    def test_standard_amino_acids(self):
        mapping = {
            "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
            "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
            "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
            "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
        }
        for three, one in mapping.items():
            assert three2one(three) == one, f"Failed for {three}"

    def test_unknown_returns_X(self):
        assert three2one("UNK") == "X"
        assert three2one("HOH") == "X"
        assert three2one("XYZ") == "X"

    def test_all_20_covered(self):
        assert len({three2one(aa) for aa in STANDARD_AA if three2one(aa) != "X"}) == 20


class TestVdwRadii:
    def test_returns_dict(self):
        radii = vdw_radii()
        assert isinstance(radii, dict)

    def test_carbon_radius(self):
        assert vdw_radii()["C"] == pytest.approx(1.76)

    def test_nitrogen_radius(self):
        assert vdw_radii()["N"] == pytest.approx(1.65)

    def test_oxygen_radius(self):
        assert vdw_radii()["O"] == pytest.approx(1.40)


class TestVdwRadiiExtended:
    def test_ala_has_CA(self):
        ala = vdw_radii_extended("ALA")
        assert "CA" in ala

    def test_unknown_returns_default(self):
        result = vdw_radii_extended("XXX")
        assert isinstance(result, dict)
        assert result == vdw_radii()

    def test_all_standard_aa_covered(self):
        # All standard amino acids should have entries
        for aa in STANDARD_AA:
            radii = vdw_radii_extended(aa)
            assert isinstance(radii, dict), f"No entry for {aa}"
            assert len(radii) > 0


class TestStandardData:
    def test_known_values(self):
        assert standard_data("ALA") == pytest.approx(107.95)
        assert standard_data("GLY") == pytest.approx(80.1)
        assert standard_data("TRP") == pytest.approx(249.36)

    def test_unknown_returns_minus_one(self):
        assert standard_data("UNK") == -1
        assert standard_data("HOH") == -1

    def test_all_standard_aa_positive(self):
        for aa in STANDARD_AA:
            val = standard_data(aa)
            assert val > 0, f"standard_data({aa}) should be positive, got {val}"


class TestConstants:
    def test_standard_aa_count(self):
        assert len(STANDARD_AA) == 20

    def test_pair_potential_keys_are_sorted(self):
        for key in PAIR_POTENTIAL:
            parts = key.split("-")
            assert len(parts) == 2
            assert parts[0] <= parts[1], f"Key {key} is not in sorted order"

    def test_default_vdw_positive(self):
        assert DEFAULT_VDW > 0

    def test_atom_dict_keys_subset_of_standard_aa(self):
        for aa in ATOM_DICT:
            if aa != "ACE":
                assert aa in STANDARD_AA, f"{aa} not in STANDARD_AA"
