import json
import csv


def _make_run(tmp_path, *, scientific_status, scheduler_stderr=""):
    run = tmp_path / "run"
    status_dir = run / "status"
    status_dir.mkdir(parents=True)
    (status_dir / "exit.json").write_text(
        json.dumps({"return_code": 0, "scientific_status": scientific_status}) + "\n",
        encoding="utf-8",
    )
    (run / "slurm-42_1.err").write_text(scheduler_stderr, encoding="utf-8")
    return run


def test_cancelled_scheduler_evidence_overrides_completed_json(tmp_path):
    from benchmark.scripts.build_pipeline_verification_baseline import classify_retained_run

    run = _make_run(
        tmp_path,
        scientific_status="completed",
        scheduler_stderr="slurmstepd: error: *** JOB 42 ON ai03 CANCELLED AT 2026-07-18T16:06:42 ***\n",
    )

    row = classify_retained_run(run)

    assert row["classification"] == "unsupported"
    assert row["reason"] == "scheduler_cancelled_status_contradiction"


def test_completed_json_without_scheduler_contradiction_is_execution_evidence(tmp_path):
    from benchmark.scripts.build_pipeline_verification_baseline import classify_retained_run

    run = _make_run(tmp_path, scientific_status="completed")

    row = classify_retained_run(run)

    assert row["classification"] == "supported"
    assert row["reason"] == "retained_status_completed"


def test_external_scheduler_log_can_override_nested_run_status(tmp_path):
    from benchmark.scripts.build_pipeline_verification_baseline import classify_retained_run

    run = _make_run(tmp_path, scientific_status="completed")
    log = tmp_path / "slurm-42_1.err"
    log.write_text("JOB 42 CANCELLED\n", encoding="utf-8")

    row = classify_retained_run(run, scheduler_logs=[log])

    assert row["classification"] == "unsupported"


def test_build_baseline_writes_claim_and_artifact_manifests(tmp_path):
    from benchmark.scripts.build_pipeline_verification_baseline import build_baseline

    run = _make_run(tmp_path, scientific_status="completed")
    config = tmp_path / "verification.json"
    config.write_text(
        json.dumps(
            {
                "runs": [
                    {
                        "claim_id": "smoke-execution",
                        "claim": "retained smoke executed",
                        "run_root": str(run),
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "baseline"

    build_baseline(config, output)

    with (output / "claims.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["claim_id"] == "smoke-execution"
    assert rows[0]["classification"] == "supported"
    with (output / "artifact_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        artifacts = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["relative_path"] for row in artifacts} == {"slurm-42_1.err", "status/exit.json"}
    assert all(len(row["sha256"]) == 64 for row in artifacts)


def test_build_baseline_hashes_declared_external_scheduler_log(tmp_path):
    from benchmark.scripts.build_pipeline_verification_baseline import build_baseline

    run = _make_run(tmp_path, scientific_status="completed")
    log = tmp_path / "outside-slurm.err"
    log.write_text("JOB 42 CANCELLED\n", encoding="utf-8")
    config = tmp_path / "verification.json"
    config.write_text(
        json.dumps(
            {
                "runs": [
                    {
                        "claim_id": "cancelled",
                        "run_root": str(run),
                        "scheduler_logs": [str(log)],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "baseline"

    build_baseline(config, output)

    with (output / "claims.tsv").open(newline="", encoding="utf-8") as handle:
        claim = next(csv.DictReader(handle, delimiter="\t"))
    assert claim["classification"] == "unsupported"
    with (output / "artifact_manifest.tsv").open(newline="", encoding="utf-8") as handle:
        artifacts = list(csv.DictReader(handle, delimiter="\t"))
    assert any(row["relative_path"] == "external:outside-slurm.err" for row in artifacts)


def test_build_baseline_refuses_to_mix_with_existing_output(tmp_path):
    from benchmark.scripts.build_pipeline_verification_baseline import build_baseline

    run = _make_run(tmp_path, scientific_status="completed")
    config = tmp_path / "verification.json"
    config.write_text(
        json.dumps({"runs": [{"claim_id": "smoke", "run_root": str(run)}]}),
        encoding="utf-8",
    )
    output = tmp_path / "baseline"
    output.mkdir()
    (output / "previous.tsv").write_text("prior evidence\n", encoding="utf-8")

    try:
        build_baseline(config, output)
    except FileExistsError as exc:
        assert "refusing to mix" in str(exc)
    else:
        raise AssertionError("non-empty baseline output was accepted")
