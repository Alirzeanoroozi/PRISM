from __future__ import annotations

import csv
import json
from pathlib import Path

from benchmark.scripts.prepare_usalign_refinement_manifest import prepare


def _write_pdb(path: Path) -> None:
    path.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        "END\n",
        encoding="ascii",
    )


def test_prepare_selects_refinable_rows_and_records_rejections(tmp_path: Path) -> None:
    left = tmp_path / "left.pdb"
    right = tmp_path / "right.pdb"
    native = tmp_path / "native.pdb"
    for path in (left, right, native):
        _write_pdb(path)
    source = tmp_path / "scores.tsv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "status", "score_status", "transformed_left", "transformed_right",
                "native_pdb", "native_receptor_chains", "native_ligand_chains",
                "case_id", "template", "candidate_index", "raw_dockq_json_sha256",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({
            "status": "generated", "score_status": "scored", "transformed_left": str(left),
            "transformed_right": str(right), "native_pdb": str(native),
            "native_receptor_chains": "A", "native_ligand_chains": "B",
            "case_id": "case_1", "template": "1abc", "candidate_index": "0",
            "raw_dockq_json_sha256": "abc",
        })
        writer.writerow({
            "status": "generated", "score_status": "score_failed", "transformed_left": str(left),
            "transformed_right": str(right), "native_pdb": str(native),
            "native_receptor_chains": "A", "native_ligand_chains": "B",
            "case_id": "case_2", "template": "2abc", "candidate_index": "1",
            "raw_dockq_json_sha256": "def",
        })
    selected = tmp_path / "selected.csv"
    rejected = tmp_path / "rejected.csv"
    result = prepare(source, selected, rejected)
    assert result["selected_count"] == 1
    assert result["rejected_count"] == 1
    with selected.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["pipeline"] == "usalign"
    assert row["left"] == str(left.resolve())
    with rejected.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["reason"] == "score_status:score_failed"

