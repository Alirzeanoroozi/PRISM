from pathlib import Path


def test_prepare_usalign_refinement_job_is_fail_closed_and_non_submitting() -> None:
    path = Path(__file__).parents[1] / "benchmark/jobs/prepare_usalign_refinement.sbatch"
    text = path.read_text(encoding="utf-8")
    assert "blocked_upstream_aggregation" in text
    assert "validated_compacted" in text
    assert "sbatch --" not in text
    assert "prepare_usalign_refinement_manifest.py" in text
