import json


def _write_stage_events(path, stages):
    with path.open("w", encoding="utf-8") as handle:
        for stage in stages:
            handle.write(json.dumps({"stage": stage, "event": "completed", "return_code": 0}) + "\n")


def _atom(chain, residue):
    return f"ATOM      1  CA  ALA {chain}{residue:4d}    0.000   0.000   0.000  1.00 20.00           C\n"


def test_signal_always_overrides_zero_return_code(tmp_path):
    from benchmark.scripts.pipeline_completion_contract import classify_run

    run = tmp_path / "run"
    status = run / "status"
    status.mkdir(parents=True)
    (status / "pipeline_returned.json").write_text("{}\n", encoding="utf-8")
    _write_stage_events(status / "stages.jsonl", ["input", "alignment", "transformation", "refinement"])
    (run / "processed" / "transformation").mkdir(parents=True)
    (run / "processed" / "transformation" / "pair_L.pdb").write_text("ATOM\n", encoding="utf-8")

    result = classify_run(run, process_return_code=0, termination_signal="SIGTERM")

    assert result["process_status"] == "cancelled"
    assert result["scientific_status"] == "cancelled"
    assert result["reason"] == "termination_signal:SIGTERM"


def test_completed_stages_without_prediction_is_completed_no_predictions(tmp_path):
    from benchmark.scripts.pipeline_completion_contract import classify_run

    run = tmp_path / "run"
    status = run / "status"
    status.mkdir(parents=True)
    (status / "pipeline_returned.json").write_text("{}\n", encoding="utf-8")
    _write_stage_events(status / "stages.jsonl", ["input", "alignment", "transformation", "refinement"])

    result = classify_run(run, process_return_code=0, termination_signal=None)

    assert result["process_status"] == "completed"
    assert result["scientific_status"] == "completed_no_predictions"
    assert result["paired_transformation_count"] == 0
    assert result["refined_model_count"] == 0


def test_skipped_refinement_is_terminal_for_no_prediction_run(tmp_path):
    from benchmark.scripts.pipeline_completion_contract import classify_run

    run = tmp_path / "run"
    status = run / "status"
    status.mkdir(parents=True)
    (status / "pipeline_returned.json").write_text("{}\n", encoding="utf-8")
    _write_stage_events(status / "stages.jsonl", ["input", "alignment", "transformation"])
    with (status / "stages.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps({
            "stage": "refinement",
            "event": "skipped",
            "return_code": 0,
            "detail": "no candidates",
        }) + "\n")

    result = classify_run(run, process_return_code=0, termination_signal=None)

    assert result["scientific_status"] == "completed_no_predictions"
    assert result["reason"] == "no_refined_models"


def test_missing_refinement_terminal_event_is_incomplete_even_with_transformations(tmp_path):
    from benchmark.scripts.pipeline_completion_contract import classify_run

    run = tmp_path / "run"
    status = run / "status"
    status.mkdir(parents=True)
    (status / "pipeline_returned.json").write_text("{}\n", encoding="utf-8")
    _write_stage_events(status / "stages.jsonl", ["input", "alignment", "transformation"])
    transformation = run / "processed" / "transformation"
    transformation.mkdir(parents=True)
    (transformation / "model_L.pdb").write_text("ATOM\n", encoding="utf-8")
    (transformation / "model_R.pdb").write_text("ATOM\n", encoding="utf-8")

    result = classify_run(run, process_return_code=0, termination_signal=None)

    assert result["scientific_status"] == "incomplete"
    assert result["reason"] == "missing_terminal_stages:refinement"
    assert result["paired_transformation_count"] == 1


def test_cli_writes_atomic_exit_status(tmp_path):
    from benchmark.scripts.pipeline_completion_contract import main

    run = tmp_path / "run"
    (run / "status").mkdir(parents=True)
    output = run / "status" / "exit.json"

    assert main([
        "--run-root", str(run),
        "--return-code", "0",
        "--termination-signal", "SIGTERM",
        "--output", str(output),
        "--pipeline", "current",
        "--aligner", "gtalign",
        "--elapsed-seconds", "12",
    ]) == 0

    status = json.loads(output.read_text(encoding="utf-8"))
    assert status["scientific_status"] == "cancelled"
    assert status["pipeline"] == "current"
    assert status["aligner"] == "gtalign"
    assert status["elapsed_seconds"] == 12


def test_invalid_nonempty_refined_model_cannot_claim_completed(tmp_path):
    from benchmark.scripts.pipeline_completion_contract import classify_run

    run = tmp_path / "run"
    status = run / "status"
    status.mkdir(parents=True)
    (status / "pipeline_returned.json").write_text("{}\n", encoding="utf-8")
    _write_stage_events(status / "stages.jsonl", ["input", "alignment", "transformation", "refinement"])
    refined = run / "processed" / "rosetta_refinement"
    refined.mkdir(parents=True)
    (refined / "bad.pdb").write_text(_atom("A", 1), encoding="ascii")

    result = classify_run(run, process_return_code=0, termination_signal=None)

    assert result["refined_model_count"] == 1
    assert result["valid_refined_model_count"] == 0
    assert result["scientific_status"] == "incomplete"
    assert result["reason"] == "no_valid_refined_models"
