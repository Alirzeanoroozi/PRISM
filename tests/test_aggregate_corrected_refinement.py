import json
from pathlib import Path

from benchmark.scripts.aggregate_corrected_refinement import aggregate


def test_aggregate_uses_global_dockq_and_preserves_cross_components(tmp_path):
    root = tmp_path / "comparison"
    checkpoint_dir = root / "checkpoints"
    checkpoint_dir.mkdir(parents=True)
    raw = tmp_path / "dockq.json"
    raw.write_text(
        json.dumps(
            {
                "GlobalDockQ": 0.3,
                "best_dockq": 0.9,
                "best_mapping_str": "AB:AC",
                "best_result": {
                    "AC": {"DockQ": 0.3, "fnat": 0.2, "iRMSD": 2.0, "LRMSD": 3.0},
                    "AB": {"DockQ": 0.9, "fnat": 0.8, "iRMSD": 1.0, "LRMSD": 2.0},
                },
            }
        )
    )
    checkpoint = {
        "key": "candidate-1",
        "index": 1,
        "status": "completed",
        "candidate": {
            "manifest_index": 1,
            "pipeline": "tmalign",
            "case_id": "rigid_1abc_000",
            "template": "1abcAB",
        },
        "stages": {
            "dockq_fiberdock": {
                "status": "scored",
                "raw_dockq_json": str(raw),
                "native_receptor_chains": "A",
                "native_ligand_chains": "B",
                "dockq": 0.9,
            },
            "dockq_rosetta": {"status": "not_run_no_model", "reason": "no model"},
        },
    }
    (checkpoint_dir / "candidate.json").write_text(json.dumps(checkpoint))

    summary = aggregate(root, tmp_path / "out", 1)

    assert summary["status"] == "complete"
    rows = (tmp_path / "out/refinement_comparison.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    values = dict(zip(header, rows[1].split("\t")))
    assert values["fiberdock_global_dockq"] == "0.3"
    assert values["fiberdock_dockq_best_internal_diagnostic"] == "0.9"
    assert values["fiberdock_cross_best"] == "0.9"
    assert values["external_rosetta_status"] == "not_run_no_model"
    assert Path(values["fiberdock_raw_dockq_json"]).is_file()


def test_aggregate_preserves_nested_failure_reason(tmp_path):
    root = tmp_path / "comparison"
    checkpoint_dir = root / "checkpoints"
    checkpoint_dir.mkdir(parents=True)
    checkpoint = {
        "key": "failed-candidate",
        "index": 7,
        "status": "failed",
        "candidate": {
            "manifest_index": 7,
            "pipeline": "gtalign",
            "case_id": "medium_1wq1_045",
            "template": "1de4AC",
            "orientation": "o1",
        },
        "stages": {
            "input_normalization": {
                "status": "failed",
                "error": "ValueError: native ligand G* requires two chain segments; model exposes one",
            }
        },
    }
    (checkpoint_dir / "failed.json").write_text(json.dumps(checkpoint))

    summary = aggregate(root, tmp_path / "out", 1)

    assert summary["status"] == "complete"
    import csv

    with (tmp_path / "out/refinement_comparison.tsv").open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["candidate_status"] == "failed"
    assert row["failed_stages"] == "input_normalization"
    assert json.loads(row["failure_reasons"])["input_normalization"].startswith("ValueError:")
