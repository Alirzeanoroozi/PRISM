from pathlib import Path

import pytest

import prism
from src import fiberdock_refinement
from src import pyrosetta_refinement


def test_cli_preserves_existing_backend_defaults():
    args = prism.build_parser().parse_args([])
    assert args.aligner == "tmalign"
    assert args.refiner == "external_rosetta"
    assert args.refine is False
    assert args.surface_backend == "freesasa"
    assert args.multiprot_mode == "current"
    assert args.multiprot_solutions == 3


def test_pyrosetta_merged_adapter_preserves_tuple(monkeypatch, tmp_path):
    candidate = tmp_path / "candidate.pdb"
    candidate.write_text("END\n")

    class Adapter:
        def __init__(self, **kwargs):
            pass

        def refine(self, input_pdb, output_pdb, partners):
            Path(output_pdb).write_text("END\n")
            assert partners == "A_B"
            return {"status": "success"}

    monkeypatch.setattr(pyrosetta_refinement, "PyRosettaRefinementAdapter", Adapter)
    result = pyrosetta_refinement.refine_merged_candidates(
        [("1abcA", "1abcB", "tpl", str(candidate))], output_root=tmp_path / "out",
    )
    assert result[0][:3] == ("1abcA", "1abcB", "tpl")
    assert result[0][3].endswith("_pyrosetta.pdb")


def test_fiberdock_split_uses_target_chain_groups(tmp_path):
    merged = tmp_path / "merged.pdb"
    merged.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ATOM      2  CA  ALA B   1       4.000   0.000   0.000  1.00 20.00           C  \nEND\n"
    )
    receptor = tmp_path / "r.pdb"
    ligand = tmp_path / "l.pdb"
    fiberdock_refinement._split_merged_candidate(
        merged, ["A"], ["B"], receptor, ligand,
    )
    assert " A   1" in receptor.read_text()
    assert " B   1" not in receptor.read_text()
    assert " B   1" in ligand.read_text()


def test_fiberdock_split_rejects_missing_partner_chain(tmp_path):
    merged = tmp_path / "merged.pdb"
    merged.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n"
    )
    with pytest.raises(ValueError, match="ligand atoms=0"):
        fiberdock_refinement._split_merged_candidate(
            merged, ["A"], ["B"], tmp_path / "r.pdb", tmp_path / "l.pdb",
        )
