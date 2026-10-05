import csv

from benchmark.scripts.prepare_full_comparison_batches import prepare_batches


def _write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Complex", "PDB ID 1", "PDB ID 2"])
        writer.writeheader()
        writer.writerows(rows)


def test_prepare_batches_preserves_identical_pairs_for_both_consumers(tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    row = {"Complex": "1ABC_A:B", "PDB ID 1": "1ABC_AB", "PDB ID 2": "2DEF_C"}
    for filename in ("T_Rigid.csv", "T_medium.csv", "T_difficult.csv"):
        _write(data_dir / filename, [row])

    batches = prepare_batches(data_dir, tmp_path / "batches", batch_size=2)

    assert len(batches) == 2
    with (batches[0] / "inputs.csv").open(newline="") as handle:
        first = next(csv.DictReader(handle))
    assert (batches[0] / "pair_list").read_text().splitlines()[0] == f"{first['Receptor']} {first['Ligand']}"
    assert sum(1 for _ in csv.DictReader((tmp_path / "batches" / "shared_manifest.csv").open())) == 3
