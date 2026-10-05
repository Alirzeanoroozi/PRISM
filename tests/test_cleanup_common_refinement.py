import hashlib
import json
from pathlib import Path

from benchmark.scripts.cleanup_common_refinement import cleanup


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_cleanup_requires_complete_aggregation_and_preserves_audit_files(tmp_path: Path):
    root = tmp_path / "gtalign_common_refinement_19855"
    (root / "results" / "candidate_a").mkdir(parents=True)
    (root / "adapter_inputs" / "task_a").mkdir(parents=True)
    (root / "results" / "candidate_a" / "scratch.pdb").write_text("pdb")
    (root / "adapter_inputs" / "task_a" / "candidate.csv").write_text("csv")
    results = root / "aggregate" / "refinement_comparison.tsv"
    interfaces = root / "aggregate" / "refinement_interfaces.tsv"
    results.parent.mkdir()
    results.write_text("result")
    interfaces.write_text("interface")
    (root / "aggregate" / "aggregation_status.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "cleanup_eligible": True,
                "results_path": str(results),
                "results_sha256": digest(results),
                "interfaces_path": str(interfaces),
                "interfaces_sha256": digest(interfaces),
            }
        )
    )
    (root / "checkpoints").mkdir()
    (root / "checkpoints" / "keep.json").write_text("keep")

    dry = cleanup(root, apply=False)
    assert dry["status"] == "dry_run"
    assert (root / "results" / "candidate_a" / "scratch.pdb").is_file()
    applied = cleanup(root, apply=True)
    assert applied["status"] == "applied"
    assert not (root / "results" / "candidate_a").exists()
    assert not (root / "adapter_inputs" / "task_a").exists()
    assert (root / "checkpoints" / "keep.json").is_file()
