import json
from pathlib import Path

import pytest

from benchmark.scripts.collect_observational_score_replay import metrics, validate_tasks


def test_irmsd_best_is_minimum_and_dockq_best_is_maximum():
    result = metrics([
        {"status": "ready", "dockq": "0.2", "irmsd": "5.0", "pair_id": "p1"},
        {"status": "ready", "dockq": "0.8", "irmsd": "1.0", "pair_id": "p2"},
    ])

    assert result["dockq_best"] == 0.8
    assert result["irmsd_best"] == 1.0


def test_validate_tasks_rejects_stale_or_missing_output_hash(tmp_path: Path):
    replay = tmp_path / "replay"
    shard_dir = replay / "shards"
    task_dir = replay / "tasks" / "task-1"
    shard_dir.mkdir(parents=True)
    task_dir.mkdir(parents=True)
    model_manifest = replay / "model_manifest.csv"
    model_manifest.write_text("model_path,pipeline\nmodel.pdb,test\n")
    shard = shard_dir / "shard_01.csv"
    shard.write_text("model_path,pipeline\nmodel.pdb,test\n")
    scored = task_dir / "scored_models.csv"
    scored.write_text("model_path,pipeline,status\nmodel.pdb,test,ready\n")
    (replay / "replay_manifest.json").write_text(json.dumps({
        "shard_count": 1,
        "model_manifest": str(model_manifest),
        "shards": [{"shard": 1, "path": str(shard), "rows": 1, "sha256": "wrong"}],
    }))
    (task_dir / "exit.json").write_text(json.dumps({
        "array_task_id": 1,
        "return_code": 0,
        "scientific_status": "completed",
        "input": str(shard),
        "output": str(scored),
        "output_sha256": "wrong",
    }))

    with pytest.raises(ValueError, match="shard hash changed"):
        validate_tasks(replay)
