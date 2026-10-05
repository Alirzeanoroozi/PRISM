import csv
import json
from pathlib import Path

from benchmark.scripts.collect_matched_benchmark import collect


def _write_tsv(path: Path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_collect_fails_closed_and_uses_minimum_irmsd(tmp_path):
    task_root = tmp_path / "task"
    task_root.mkdir()
    (task_root / "exit.json").write_text(json.dumps({"return_code": 0, "scientific_status": "completed", "elapsed_seconds": 4}), encoding="utf-8")
    _write_tsv(
        task_root / "scores_global.tsv",
        ("GlobalDockQ", "iRMSD"),
        [{"GlobalDockQ": "0.2", "iRMSD": "3.0"}, {"GlobalDockQ": "0.8", "iRMSD": "1.0"}],
    )
    _write_tsv(task_root / "scores_interfaces.tsv", ("interface", "DockQ"), [{"interface": "A:B", "DockQ": "0.7"}])
    manifest = tmp_path / "task_manifest.tsv"
    _write_tsv(manifest, ("task_id", "dataset_row_id", "arm", "experiment", "output_root"), [{
        "task_id": "rigid:000001:tmalign:score", "dataset_row_id": "rigid:000001", "arm": "tmalign_external_rosetta", "experiment": "strict_scoring", "output_root": "task",
    }])
    paths = collect(manifest, tmp_path / "collected")
    rows = list(csv.DictReader(paths["pair_summary"].open(), delimiter="\t"))
    assert rows[0]["best_GlobalDockQ_at_20"] == "0.8"
    assert rows[0]["best_iRMSD_at_20"] == "1.0"
    assert rows[0]["interface_rows"] == "1"
    assert list(csv.DictReader(paths["failures"].open(), delimiter="\t")) == []


def test_collect_records_missing_exit_as_failure(tmp_path):
    manifest = tmp_path / "task_manifest.tsv"
    _write_tsv(manifest, ("task_id", "dataset_row_id", "arm", "experiment", "output_root"), [{
        "task_id": "rigid:000001:tmalign:score", "dataset_row_id": "rigid:000001", "arm": "tmalign_external_rosetta", "experiment": "strict_scoring", "output_root": "missing",
    }])
    paths = collect(manifest, tmp_path / "collected")
    failures = list(csv.DictReader(paths["failures"].open(), delimiter="\t"))
    assert failures[0]["reason"] == "missing_exit_record"


def test_collect_rejects_declared_corrupt_complex_output(tmp_path):
    task_root = tmp_path / "task"
    task_root.mkdir()
    (task_root / "exit.json").write_text(json.dumps({
        "return_code": 0, "scientific_status": "completed", "elapsed_seconds": 1,
        "output_pdb": "refined.pdb",
    }), encoding="utf-8")
    (task_root / "refined.pdb").write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "task_manifest.tsv"
    _write_tsv(manifest, ("task_id", "dataset_row_id", "arm", "experiment", "output_root"), [{
        "task_id": "rigid:000001:score", "dataset_row_id": "rigid:000001",
        "arm": "tmalign_external_rosetta", "experiment": "strict_scoring", "output_root": "task",
    }])
    paths = collect(manifest, tmp_path / "collected")
    failures = list(csv.DictReader(paths["failures"].open(), delimiter="\t"))
    assert any(row["reason"] == "invalid_or_missing_output_pdb" for row in failures)
