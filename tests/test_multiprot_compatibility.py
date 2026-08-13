import csv
import json
from pathlib import Path

import numpy as np

from src import alignment_multiprot
from src import transformation


def _write_ca_pdb(path, coordinates):
    rows = []
    for index, (x, y, z) in enumerate(coordinates, 1):
        rows.append(
            f"ATOM  {index:5d}  CA  ALA A{index:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
        )
    path.write_text("".join(rows) + "END\n")


def test_current_transform_maps_query_into_template_frame(tmp_path):
    query = tmp_path / "query.pdb"
    interface = tmp_path / "interface.pdb"
    query_points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    expected_rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    expected_translation = np.array([4.0, -2.0, 3.0])
    interface_points = (expected_rotation @ query_points.T).T + expected_translation
    _write_ca_pdb(query, query_points)
    _write_ca_pdb(interface, interface_points)
    matches = {f"I.A.{i}": f"Q.A.{i}" for i in range(1, 4)}

    translation, rotation, rmsd = alignment_multiprot._compute_transform_from_matches(
        matches, query, interface
    )

    transformed = (np.asarray(rotation) @ query_points.T).T + np.asarray(translation)
    np.testing.assert_allclose(transformed, interface_points, atol=1e-6)
    assert rmsd < 1e-6


def test_legacy_parser_preserves_reference_and_solutions(tmp_path):
    result = tmp_path / "2_sol.res"
    result.write_text(
        "Solution Num : 0\n"
        "Mult Corres Score : 3\n"
        "Reference Molecule : 1\n"
        "Trans : 0.1 0.2 0.3 0.4 0.5 0.6\n"
        "RMSD : 1.25\n"
        "Match List (Chain_ID.Res_Type.Res_Num) : 3\n"
        "C.A.10 I.G.20\nC.G.11 I.A.21\nC.L.12 I.V.22\n"
        "End of Match List\n"
        "Solution Num : 1\n"
        "Mult Corres Score : 3\n"
        "Reference Molecule : 0\n"
        "Trans : 0 0 0 0 0 0\n"
        "RMSD : 2.5\n"
        "Match List (Chain_ID.Res_Type.Res_Num) : 3\n"
        "C.A.10 I.G.20\nC.G.11 I.A.21\nC.L.12 I.V.22\n"
        "End of Match List\n"
    )

    solutions = alignment_multiprot._parse_multiprot_solutions(result, 3)

    assert len(solutions) == 2
    assert solutions[0]["reference_molecule"] == 1
    assert solutions[0]["match_dict"] == {
        "I.G.20": "C.A.10",
        "I.A.21": "C.G.11",
        "I.V.22": "C.L.12",
    }


def test_legacy_transform_inverts_reference_one():
    transform = alignment_multiprot._legacy_solution_to_alignment({
        "reference_molecule": 1,
        "trans": [0.0, 0.0, 0.0, 1.0, 2.0, 3.0],
    })
    np.testing.assert_allclose(transform["rotation_mat"], np.eye(3), atol=1e-12)
    np.testing.assert_allclose(transform["translation"], [-1.0, -2.0, -3.0])


def test_align_multiprot_writes_transformation_summary(monkeypatch, tmp_path):
    output = tmp_path / "alignment"

    def fake_align(task):
        query, template, chain = task[:3]
        return {
            "protein": query,
            "template": template,
            "chain": chain,
            "match_count": 12,
            "tm_score": 0.0,
            "len_target": 20,
            "len_template": 20,
            "translation": json.dumps([0.0, 0.0, 0.0]),
            "rotation_mat": json.dumps(np.eye(3).tolist()),
        }

    monkeypatch.setattr(alignment_multiprot, "_align_one", fake_align)
    alignment_multiprot.align_multiprot(
        ["1abcA"], ["2defBC"], output_dir=output, max_workers=1,
        multiprot_mode="legacy_compatible", multiprot_solutions=3,
    )

    rows = list(csv.DictReader((output / "1abcA.csv").open()))
    assert [(row["template"], row["chain"]) for row in rows] == [
        ("2defBC", "B"), ("2defBC", "C")
    ]


def test_transformation_expands_retained_solutions(monkeypatch, tmp_path):
    payload = {
        "aligner": "MultiProt",
        "tm_score": 0.0,
        "multiprot_solutions": [
            {"match_count": 12, "translation": [1, 0, 0], "match_dict": {"A.A.1": "X.A.1"}},
            {"match_count": 11, "translation": [2, 0, 0], "match_dict": {"A.A.2": "X.A.2"}},
        ],
    }
    (tmp_path / "1abcA_2defBC_B.json").write_text(json.dumps(payload))
    monkeypatch.setattr(transformation, "ALIGNMENT_DIR", str(tmp_path))

    variants = transformation._alignment_variants(
        "1abcA", "2defBC", "B", {"tm_score": 0.0}
    )

    assert [variant["match_count"] for variant in variants] == [12, 11]
    assert [variant["translation"] for variant in variants] == [[1, 0, 0], [2, 0, 0]]
