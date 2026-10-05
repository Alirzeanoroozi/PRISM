from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from benchmark.scripts.cleanup_usalign_run import cleanup, planned_paths


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_cleanup_usalign_requires_all_validated_consumers(tmp_path: Path) -> None:
    root = tmp_path / "usalign_production_19855"
    refinement = tmp_path / "usalign_common_refinement_19855"
    final = tmp_path / "final_matched"
    batch = root / "current" / "batch_0001"
    batch.mkdir(parents=True)
    (batch / "raw.json").write_text("raw")
    aggregate = root / "aggregate"
    aggregate.mkdir(parents=True)
    retained = {}
    for name in ("candidate_generated.tsv", "transformed_dockq.tsv", "transformed_dockq_interfaces.tsv"):
        path = aggregate / name
        path.write_text(name)
        stem = name.removesuffix(".tsv")
        key_stem = "transformed_interfaces" if stem == "transformed_dockq_interfaces" else stem
        retained[f"{key_stem}_path"] = str(path.resolve())
        retained[f"{key_stem}_sha256"] = _hash(path)
    (aggregate / "aggregation_status.json").write_text(json.dumps({"status": "validated_compacted", **retained}))
    handoff = root / "usalign_common_refinement_19855"
    handoff.mkdir()
    (handoff / "handoff_status.json").write_text(json.dumps({"status": "validated_manifest_ready_for_resource_review"}))
    refinement_result = refinement / "aggregate" / "refinement_comparison.tsv"
    refinement_interface = refinement / "aggregate" / "refinement_interfaces.tsv"
    refinement_result.parent.mkdir(parents=True)
    refinement_result.write_text("result")
    refinement_interface.write_text("interface")
    (refinement / "aggregate" / "aggregation_status.json").write_text(json.dumps({
        "status": "complete", "cleanup_eligible": True,
        "results_path": str(refinement_result.resolve()), "results_sha256": _hash(refinement_result),
        "interfaces_path": str(refinement_interface.resolve()), "interfaces_sha256": _hash(refinement_interface),
    }))
    (refinement / "results" / "candidate").mkdir(parents=True)
    (refinement / "results" / "candidate" / "scratch").write_text("scratch")
    final.mkdir()
    (final / "aggregate_manifest.json").write_text(json.dumps({"status": "validated_compacted"}))

    dry = cleanup(root, refinement, final, apply=False)
    assert dry["status"] == "dry_run"
    assert (batch / "raw.json").is_file()
    assert json.loads((root / "cleanup_manifest.json").read_text())["status"] == "dry_run"
    applied = cleanup(root, refinement, final, apply=True)
    assert applied["status"] == "applied"
    assert not batch.exists()
    assert not (refinement / "results").exists()
    assert aggregate.exists()
    assert (refinement / "aggregate" / "refinement_comparison.tsv").is_file()
    assert json.loads((root / "cleanup_manifest.json").read_text())["status"] == "applied"


def test_cleanup_rejects_symlinked_batch(tmp_path: Path) -> None:
    root = tmp_path / "usalign_production_19855"
    current = root / "current"
    current.mkdir(parents=True)
    target = tmp_path / "outside"
    target.mkdir()
    (current / "batch_0001").symlink_to(target, target_is_directory=True)

    with pytest.raises(RuntimeError, match="unsafe USalign batch cleanup target"):
        planned_paths(root, tmp_path / "usalign_common_refinement_19855")
