import csv
from pathlib import Path

from benchmark.scripts.run_common_refinement_adapter import split_combined_model


def test_split_combined_model_writes_two_chain_inputs(tmp_path: Path):
    source = tmp_path / "model.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       0.0   0.0   0.0  1.00  0.00           C\n"
        "ATOM      2  CA  GLY B   1       1.0   0.0   0.0  1.00  0.00           C\nEND\n"
    )
    row = {"model": str(source), "chain_left": "A", "chain_right": "B", "source_row": "7"}
    left, right = split_combined_model(row, tmp_path / "split")
    assert left.is_file() and right.is_file()
    assert " A " in left.read_text()
    assert " B " in right.read_text()
    assert " B " not in left.read_text()


def test_split_paths_are_stable_for_resume(tmp_path: Path):
    source = tmp_path / "model.pdb"
    source.write_text("ATOM      1  CA  ALA A   1       0.0   0.0   0.0  1.00  0.00           C\nEND\n")
    row = {"model": str(source), "chain_left": "A", "chain_right": "B", "source_row": "7"}
    first = split_combined_model(row, tmp_path / "split")
    second = split_combined_model(row, tmp_path / "split")
    assert first == second
