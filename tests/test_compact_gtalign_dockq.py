import csv
import json

from benchmark.scripts.compact_gtalign_dockq import compact


def test_compact_gtalign_keeps_global_and_requested_cross_scores(tmp_path):
    raw = tmp_path / "raw.json"
    raw.write_text(json.dumps({
        "GlobalDockQ": 0.25,
        "best_dockq": 0.8,
        "best_mapping_str": "AB:AC",
        "best_result": {
            "AC": {"DockQ": 0.25, "fnat": 0.2, "iRMSD": 2.0, "LRMSD": 3.0},
            "AB": {"DockQ": 0.8, "fnat": 0.9, "iRMSD": 1.0, "LRMSD": 2.0},
        },
    }))
    scores = tmp_path / "scores.csv"
    with scores.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["case_id", "template", "orientation", "status", "raw_dockq_json"])
        writer.writeheader()
        writer.writerow({"case_id": "rigid_1abc_000", "template": "1abcAB", "orientation": "1", "status": "scored", "raw_dockq_json": str(raw)})
    manifest = tmp_path / "dataset.json"
    manifest.write_text(json.dumps({"pairs": [{
        "case_id": "rigid_1abc_000", "split": "rigid", "benchmark_complex": "1ABC_A:B",
        "native_receptor_chains": "A", "native_ligand_chains": "C",
    }]}))

    summary = compact(scores, manifest, tmp_path / "out")

    assert summary["status_counts"] == {"scored": 1}
    with (tmp_path / "out/gtalign_transformed_dockq.tsv").open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["dockq_global"] == "0.25"
    assert row["dockq_best_internal_diagnostic"] == "0.8"
    assert row["dockq_cross_best"] == "0.25"


def test_compact_gtalign_does_not_promote_source_score_failure(tmp_path):
    raw = tmp_path / "raw.json"
    raw.write_text(json.dumps({"GlobalDockQ": 0.5, "best_result": {}}))
    scores = tmp_path / "scores.csv"
    with scores.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["case_id", "status", "raw_dockq_json"])
        writer.writeheader()
        writer.writerow({"case_id": "rigid_1abc_000", "status": "score_failed", "raw_dockq_json": str(raw)})
    manifest = tmp_path / "dataset.json"
    manifest.write_text(json.dumps({"pairs": [{"case_id": "rigid_1abc_000"}]}))

    summary = compact(scores, manifest, tmp_path / "out")

    assert summary["status_counts"] == {"score_failed": 1}
