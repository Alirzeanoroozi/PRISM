import json
from pathlib import Path

from benchmark.scripts.prepare_gtalign_refinement_manifest import prepare


def test_prepare_gtalign_manifest_preserves_model_and_failure_status(tmp_path: Path):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text("ATOM      1  CA  ALA A   1       0.0   0.0   0.0  1.00  0.00           C\nEND\n")
    native.write_text(model.read_text())
    raw = tmp_path / "raw.json"
    raw.write_text(json.dumps({"model": str(model), "native": str(native)}))
    source = tmp_path / "scores.tsv"
    source.write_text(
        "case_id\ttemplate\torientation\tchain_left\tchain_right\tquery_left\tquery_right\tnative_receptor_chains\tnative_ligand_chains\traw_dockq_json\tsource_status\tsplit\n"
        f"c1\tt1\to1\tA\tB\tq1\tq2\tA\tB\t{raw}\tscore_failed\trigid\n"
    )
    output = tmp_path / "selected.csv"
    rejected = tmp_path / "rejected.csv"
    manifest = prepare(source, output, rejected)
    assert manifest["selected_count"] == 1
    row = output.read_text().splitlines()[1]
    assert str(model) in row
    assert "score_failed" in row
    assert manifest["rejected_count"] == 0
