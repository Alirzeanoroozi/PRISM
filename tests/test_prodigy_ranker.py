"""Focused tests for the opt-in PRODIGY ranking adapter."""

from pathlib import Path
import sys

import pytest

from src import prodigy_ranker


def _pdb(chain: str, x: float) -> str:
    return (
        f"ATOM      1  CA  ALA {chain}   1      {x:8.3f}  0.000   0.000  1.00 20.00           C  \n"
        f"TER       2      ALA {chain}   1\nEND\n"
    )


def test_parse_affinity_reads_quiet_output():
    assert prodigy_ranker.parse_affinity("candidate  -9.373\n") == pytest.approx(-9.373)


def test_combine_pair_pdbs_namespaces_overlapping_chains(tmp_path):
    left = tmp_path / "left.pdb"
    right = tmp_path / "right.pdb"
    combined = tmp_path / "combined.pdb"
    left.write_text(_pdb("A", 0.0))
    right.write_text(_pdb("A", 1.0))

    selection_left, selection_right = prodigy_ranker.combine_pair_pdbs(
        str(left), str(right), str(combined)
    )

    assert (selection_left, selection_right) == ("A", "B")
    chains = {
        line[21]
        for line in combined.read_text().splitlines()
        if line.startswith("ATOM")
    }
    assert chains == {"A", "B"}


def test_score_candidate_places_input_before_selection_for_prodigy_argparse(tmp_path):
    left = tmp_path / "1xyzAB_1abcA_1abcB_o1_L.pdb"
    right = tmp_path / "1xyzAB_1abcA_1abcB_o1_R.pdb"
    left.write_text(_pdb("A", 0.0))
    right.write_text(_pdb("B", 1.0))
    fake = tmp_path / "fake_prodigy.py"
    fake.write_text(
        "import pathlib, sys\n"
        "argv = sys.argv[1:]\n"
        "input_path = argv[-4]\n"
        "assert pathlib.Path(input_path).is_file(), argv\n"
        "assert argv[-3:] == ['--selection', 'A', 'B'], argv\n"
        "print('-7.5')\n"
    )

    score = prodigy_ranker.score_candidate(
        str(left),
        str(right),
        executable=f"{sys.executable} {fake}",
        output_dir=str(tmp_path / "scores"),
    )

    assert score.status == "scored"
    assert score.affinity_kcal_mol == pytest.approx(-7.5)


def test_prodigy_ranking_uses_lowest_affinity_and_preserves_failed_group(tmp_path, monkeypatch):
    pairs = []
    for template, receptor, ligand in (
        ("1xyzAB", "1abcA", "1abcB"),
        ("2xyzCD", "1abcA", "1abcB"),
        ("3xyzEF", "2abcA", "2abcB"),
    ):
        left = tmp_path / f"{template}_{receptor}_{ligand}_o1_L.pdb"
        right = tmp_path / f"{template}_{receptor}_{ligand}_o1_R.pdb"
        left.write_text(_pdb("A", 0.0))
        right.write_text(_pdb("B", 1.0))
        pairs.append((str(left), str(right)))

    def fake_score(left, right, **kwargs):
        value = {"1xyz": -8.0, "2xyz": -9.0}.get(next((key for key in ("1xyz", "2xyz") if key in left), ""))
        return prodigy_ranker.ProdigyScore(
            left, right, "combined.pdb", "scored" if value is not None else "failed",
            value, 0, ("prodigy",), "stdout", "stderr", "hash",
            None if value is not None else "failed",
        )

    monkeypatch.setattr(prodigy_ranker, "score_candidate", fake_score)
    selected = prodigy_ranker.select_top_candidates_with_prodigy(
        pairs, top_k=1, executable="prodigy", output_dir=str(tmp_path / "scores")
    )

    # The failed candidate keeps the complete group visible; no partial
    # external-tool result is promoted to a scientific ranking decision.
    assert pairs[1] in selected
    assert pairs[0] not in selected
    assert pairs[2] in selected
