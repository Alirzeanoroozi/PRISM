from benchmark.scripts.run_pipeline_completion_canaries import run


def test_completion_canaries_have_expected_statuses(tmp_path):
    results = run(tmp_path / "canaries")
    assert results["cancelled"]["scientific_status"] == "cancelled"
    assert results["no_prediction"]["scientific_status"] == "completed_no_predictions"
    assert results["positive_full_pose"]["scientific_status"] == "completed"
