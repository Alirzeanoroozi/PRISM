from pathlib import Path

import numpy as np

from src.alignment_multiprot import (
    _legacy_solution_to_alignment,
    _parse_multiprot_solutions,
)
from src.transformation import _alignment_variants


def test_legacy_solution_parser_preserves_reference_and_multiple_solutions(tmp_path: Path):
    solution_file = tmp_path / "2_sol.res"
    solution_file.write_text(
        "Molecules:\n"
        "\nResults:\n\n"
        "Solution Num : 0\n\n"
        "Mult Corres Score : 4\n"
        "Reference Molecule : 1\n"
        "Molecule : 0\n"
        "Trans : 0.1 0.2 0.3 0.4 0.5 0.6\n"
        "RMSD : 1.25\n\n"
        "Match List (Chain_ID.Res_Type.Res_Num) : 4\n"
        "C.A.10 I.G.20\n"
        "C.G.11 I.A.21\n"
        "C.L.12 I.V.22\n"
        "C.V.13 I.L.23\n"
        "End of Match List\n\n"
        "Solution Num : 1\n\n"
        "Mult Corres Score : 3\n"
        "Reference Molecule : 0\n"
        "Molecule : 1\n"
        "Trans : 0 0 0 0 0 0\n"
        "RMSD : 2.5\n\n"
        "Match List (Chain_ID.Res_Type.Res_Num) : 3\n"
        "C.A.10 I.G.20\n"
        "C.G.11 I.A.21\n"
        "C.L.12 I.V.22\n"
        "End of Match List\n"
    )

    solutions = _parse_multiprot_solutions(str(solution_file), max_solutions=3)

    assert len(solutions) == 2
    assert solutions[0]["match_count"] == 4
    assert solutions[0]["reference_molecule"] == 1
    assert solutions[0]["trans"] == [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    assert solutions[0]["match_dict"] == {
        "I.G.20": "C.A.10",
        "I.A.21": "C.G.11",
        "I.V.22": "C.L.12",
        "I.L.23": "C.V.13",
    }


def test_legacy_transform_matches_identity_translation_contract():
    alignment = _legacy_solution_to_alignment(
        {
            "reference_molecule": 0,
            "trans": [0.0, 0.0, 0.0, 1.0, 2.0, 3.0],
        }
    )
    np.testing.assert_allclose(alignment["rotation_mat"], [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ], atol=1e-12)
    assert alignment["translation"] == [1.0, 2.0, 3.0]

    inverse = _legacy_solution_to_alignment(
        {
            "reference_molecule": 1,
            "trans": [0.0, 0.0, 0.0, 1.0, 2.0, 3.0],
        }
    )
    np.testing.assert_allclose(inverse["translation"], [-1.0, -2.0, -3.0], atol=1e-12)


def test_legacy_transform_preserves_nontrivial_reference_inversion():
    trans = [0.2, -0.4, 0.7, 1.0, 2.0, 3.0]
    forward = _legacy_solution_to_alignment({"reference_molecule": 0, "trans": trans})
    inverse = _legacy_solution_to_alignment({"reference_molecule": 1, "trans": trans})

    np.testing.assert_allclose(
        inverse["rotation_mat"], np.asarray(forward["rotation_mat"]).T, atol=1e-12
    )
    np.testing.assert_allclose(
        inverse["translation"],
        -np.asarray(inverse["rotation_mat"]) @ np.asarray(trans[3:]),
        atol=1e-12,
    )


def test_transformation_promotes_retained_legacy_solutions():
    base = {
        "aligner": "MultiProt",
        "multiprot_mode": "legacy_compatible",
        "tm_score": 0.0,
        "multiprot_solutions": [
            {"match_count": 12, "match_dict": {"A.A.1": "B.A.2"}},
            {"match_count": 9, "match_dict": {"A.G.3": "B.G.4"}},
        ],
    }
    variants = _alignment_variants(base)
    assert [variant["match_count"] for variant in variants] == [12, 9]
    assert all(variant["aligner"] == "MultiProt" for variant in variants)
    assert all(variant["multiprot_mode"] == "legacy_compatible" for variant in variants)
