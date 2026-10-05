import json
import csv
import hashlib
import pytest

from benchmark.scripts.investigation_artifacts import (
    freeze_pose_artifacts,
    write_dockq_tsv,
    write_pair_summary_tsv,
    write_poses_tsv,
)
from benchmark.scripts.investigation_lineage import LineageRecord


def test_freeze_pose_and_write_pose_schema(tmp_path):
    source = tmp_path / "source.pdb"
    source.write_bytes(b"ATOM\n")
    rows = freeze_pose_artifacts(
        [{"pose_id": "pose-1", "pair_id": "pair-1", "artifact_path": str(source), "status": "pose_created"}],
        tmp_path / "immutable",
    )
    output = write_poses_tsv(tmp_path / "poses.tsv", rows)
    assert (tmp_path / "immutable" / "pose-1.pdb").read_bytes() == b"ATOM\n"
    assert rows[0]["coordinate_sha256"] == hashlib.sha256(b"ATOM\n").hexdigest()
    assert output.read_text().splitlines()[0].startswith("pose_id\tpair_id")
    assert "pose-1" in output.read_text()

    with pytest.raises(ValueError, match="unsafe pose_id"):
        freeze_pose_artifacts(
            [{"pose_id": "../../escape", "artifact_path": str(source)}],
            tmp_path / "immutable-unsafe",
        )


def test_dockq_writer_separates_global_and_interfaces(tmp_path):
    raw = {
        "GlobalDockQ": 0.4,
        "best_result": {
            "A:B": {"DockQ": 0.4, "iRMSD": 1.2, "LRMSD": 2.3, "fnat": 0.5, "F1": 0.6, "clashes": 0},
            "A:C": {"DockQ": 0.2, "iRMSD": 3.2, "LRMSD": 4.3, "fnat": 0.1, "F1": 0.2, "clashes": 1},
        },
    }
    raw_path = tmp_path / "dockq.json"
    raw_path.write_text(json.dumps(raw), encoding="utf-8")
    global_path = tmp_path / "scores_global.tsv"
    interface_path = tmp_path / "scores_interfaces.tsv"
    write_dockq_tsv(raw_path, global_path, interface_path, raw_json_path=raw_path)

    with global_path.open(newline="") as handle:
        global_rows = list(csv.DictReader(handle, delimiter="\t"))
    with interface_path.open(newline="") as handle:
        interface_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(global_rows) == 1
    assert len(interface_rows) == 2
    assert "0.4" in global_path.read_text()
    assert "A:B" in interface_path.read_text()
    assert global_rows[0]["record_type"] == "global"
    assert {row["interface"] for row in interface_rows} == {"A:B", "A:C"}


def test_pair_summary_writer_preserves_unconditional_zero(tmp_path):
    output = write_pair_summary_tsv(
        tmp_path / "pair_summary.tsv",
        [{"pair_id": "pair-1", "status": "alignment_failed", "failure_reason": "no_match"}],
    )
    text = output.read_text()
    assert "primary_dockq" in text.splitlines()[0]
    assert "\t0.0\t0.0\t" in text
