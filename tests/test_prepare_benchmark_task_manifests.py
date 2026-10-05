import csv
from pathlib import Path

from benchmark.scripts.prepare_benchmark_task_manifests import prepare_manifests


HEADER = [
    "Complex",
    "Cat.",
    "PDB ID 1",
    "Protein 1",
    "PDB ID 2",
    "Protein 2",
    "I-RMSD (Å)",
    "ΔASA(Å2)",
    "BM version introduced",
]


def test_prepare_manifests_preserves_dataset_row_identity_and_partial_size(tmp_path):
    root = tmp_path / "repo"
    for dataset in ("T_Rigid.csv", "T_medium.csv", "T_difficult.csv"):
        path = root / "benchmark" / "data" / dataset
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(HEADER)
            writer.writerow(["1ABC_A:B", "AA", "1AAA_A", "r", "2BBB_B", "l", "1", "2", "5.5"])
    records = prepare_manifests(root, tmp_path / "manifests", chunk_size=10)
    assert len(records) == 1
    assert [record["array_size"] for record in records] == ["3"]
    assert [record["first_dataset_row_id"] for record in records] == ["rigid:000001"]
    manifest = Path(records[0]["task_manifest"])
    rows = list(csv.DictReader(manifest.open(newline="", encoding="utf-8")))
    assert rows[0]["task_id"].startswith("rigid:000001")
    assert "--dataset-row-id rigid:000001" in rows[0]["command"]
    assert rows[0]["output_paths"] == (
        "source/source_manifest.tsv;source/structure_validation.tsv;"
        "source/staged/staged_sources.tsv;source/source_gate_summary.json"
    )
