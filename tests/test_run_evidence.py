import json
import os

import pytest

from src.provenance import (
    ArtifactLedgerError,
    ConsumerGateError,
    ArtifactObservation,
    append_artifact_observation,
    build_declared_contract,
    build_execution_attempt,
    closeout_artifact_ledger,
    main,
    observe_artifact,
    validate_artifact_ledger,
    validate_before_consume,
    write_artifact_ledger,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _full_contract(overrides: dict | None = None) -> dict:
    """Return a valid declared-contract dict."""
    kw = {
        "pipeline_version": "1.0.0",
        "stages_enabled": ["alignment", "scoring"],
        "aligner": "tmalign",
        "refiner": "pyrosetta",
        "input_selectors": {"raw": [], "normalized": []},
        "template_inventory": [],
        "parameters": {"threshold": 0.4},
        "resource_request": {"cpus": 1, "memory_gb": 4},
        "source_inventory": {"git_head": "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "git_diff_hash": "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb", "declared_untracked": []},
        "tool_fingerprints": [{"name": "test", "version": "1.0", "sha256": None}],
    }
    if overrides:
        kw.update(overrides)
    return build_declared_contract(**kw)


def _observed(artifact_path, run_root, row_id="row-1", role="query_structure", rel_path="inputs/input.pdb", stg="alignment", produced_by="TMalign") -> ArtifactObservation:
    """Observe an artifact with sensible defaults."""
    return observe_artifact(
        artifact_path,
        run_id="prism-20240601-120000-12345-a1b2c3d4",
        stage=stg,
        dataset_row_id=row_id,
        scientific_role=role,
        run_relative_path=rel_path,
        produced_by=produced_by,
    )


def _expected_inventory(*entries) -> list[tuple[str, str, str]]:
    """Build an expected inventory from *entries* of (row_id, role, path)."""
    return list(entries)


# ===================================================================
# Declared contract & execution attempt
# ===================================================================

def test_contract_hash_is_deterministic_and_attempts_distinct():
    first = _full_contract({"parameters": {"threshold": 0.4, "api_token": "s3cret"}})
    second = _full_contract({"parameters": {"api_token": "other", "threshold": 0.5}})

    assert first["contract_hash"] != second["contract_hash"]

    # Same parameters (including secrets) → same hash (deterministic)
    duplicate = _full_contract({"parameters": {"threshold": 0.4, "api_token": "s3cret"}})
    assert first["contract_hash"] == duplicate["contract_hash"]

    # Execution attempts point to same contract hash but have distinct run IDs
    attempt_a = build_execution_attempt(first, run_id="run-1")
    attempt_b = build_execution_attempt(
        first,
        run_id="run-2",
        parent_run_id="run-1",
        supersedes_reason="retry with more CPU",
    )
    assert attempt_a["declared_contract_hash"] == attempt_b["declared_contract_hash"]
    assert attempt_a["run_identity"]["run_id"] != attempt_b["run_identity"]["run_id"]
    assert attempt_b["run_identity"]["parent_run_id"] == "run-1"
    assert attempt_b["run_identity"]["supersedes_reason"] == "retry with more CPU"


def test_execution_attempt_records_identity_output_and_safe_environment(tmp_path):
    contract = _full_contract()
    attempt = build_execution_attempt(
        contract,
        run_id="attempt-1",
        command=["python", "prism.py", "--aligner", "tmalign"],
        output_root=tmp_path / "output",
        environment={"CONDA_PREFIX": "/env", "API_TOKEN": "must-not-be-serialized"},
    )

    assert attempt["run_identity"]["attempt_id"] == "attempt-1"
    assert attempt["runtime_context"]["output_root"] == str((tmp_path / "output").resolve())
    assert attempt["runtime_context"]["command_argv"][-1] == "tmalign"
    assert "API_TOKEN" not in json.dumps(attempt)
    assert attempt["runtime_context"]["environment"]["selected"] == {
        "CONDA_PREFIX": "/env"
    }


# ===================================================================
# Artifact observation
# ===================================================================

def test_observe_file_computes_sha256(tmp_path):
    run_root = tmp_path / "run"
    p = run_root / "inputs" / "input.pdb"
    p.parent.mkdir(parents=True)
    p.write_bytes(b"ATOM original\n")

    obs = _observed(p, run_root)
    assert obs.path_kind == "file"
    assert obs.status == "ok"
    assert obs.sha256 is not None
    assert len(obs.sha256) == 64


def test_observe_symlink_hashes_target(tmp_path):
    run_root = tmp_path / "run"
    target = run_root / "source.pdb"
    link = run_root / "inputs" / "input.pdb"
    target.parent.mkdir(parents=True)
    link.parent.mkdir(parents=True)
    target.write_bytes(b"ATOM target\n")
    link.symlink_to("../source.pdb")

    obs = _observed(link, run_root)
    assert obs.path_kind == "symlink"
    assert obs.status == "ok"
    assert obs.link_target == "../source.pdb"
    assert obs.sha256  # hash of the resolved target


def test_observe_missing_file(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    obs = _observed(run_root / "missing.pdb", run_root, rel_path="missing.pdb")
    assert obs.status == "missing"
    assert obs.path_kind == "missing"
    assert obs.sha256 is None


# ===================================================================
# Ledger TSV I/O
# ===================================================================

def test_write_and_read_ledger_roundtrip(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    p = run_root / "input.pdb"
    p.write_bytes(b"ATOM\n")

    obs = _observed(p, run_root)
    ledger = tmp_path / "ledger.tsv"
    write_artifact_ledger(ledger, [obs])

    records = validate_artifact_ledger(ledger)
    assert len(records) == 1
    assert records[0].dataset_row_id == obs.dataset_row_id
    assert records[0].scientific_role == obs.scientific_role
    assert records[0].run_relative_path == obs.run_relative_path
    assert records[0].sha256 == obs.sha256


def test_write_then_append(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    p1 = run_root / "a.pdb"
    p2 = run_root / "b.pdb"
    p1.write_bytes(b"ATOM a\n")
    p2.write_bytes(b"ATOM b\n")

    ledger = tmp_path / "ledger.tsv"
    a = _observed(p1, run_root, rel_path="a.pdb")
    write_artifact_ledger(ledger, [a])
    b = _observed(p2, run_root, rel_path="b.pdb")
    append_artifact_observation(ledger, b)

    records = validate_artifact_ledger(ledger)
    assert len(records) == 2


def test_read_empty_ledger_is_empty_list(tmp_path):
    ledger = tmp_path / "nonexistent.tsv"
    records = validate_artifact_ledger(ledger)
    assert records == []


# ===================================================================
# Consumer gate (validate before consume)
# ===================================================================

def test_consumer_gate_accepts_matching_artifact(tmp_path):
    run_root = tmp_path / "run"
    artifact = run_root / "inputs" / "input.pdb"
    artifact.parent.mkdir(parents=True)
    artifact.write_bytes(b"ATOM original\n")
    obs = _observed(artifact, run_root)
    ledger = tmp_path / "ledger.tsv"
    write_artifact_ledger(ledger, [obs])

    expected = _expected_inventory(("row-1", "query_structure", "inputs/input.pdb"))
    result = validate_before_consume(ledger, expected, run_id="test-run")
    assert result["overall_status"] == "pass"


def test_consumer_gate_and_closeout_detect_mutations(tmp_path):
    run_root = tmp_path / "run"
    artifact = run_root / "inputs" / "input.pdb"
    artifact.parent.mkdir(parents=True)
    artifact.write_bytes(b"ATOM original\n")
    obs = _observed(artifact, run_root)
    ledger = tmp_path / "ledger.tsv"
    write_artifact_ledger(ledger, [obs])

    expected = _expected_inventory(("row-1", "query_structure", "inputs/input.pdb"))
    # Consumer gate passes: artifact exists with OK status
    assert validate_before_consume(ledger, expected)["overall_status"] == "pass"

    # Mutate the artifact on disk — ledger still records original hash
    artifact.write_bytes(b"ATOM changed\n")
    # The consumer gate must reject the mutated bytes before scoring.
    mutated = validate_before_consume(ledger, expected, staging_dir=run_root)
    assert mutated["overall_status"] == "fail"
    # Closeout detects the change
    closed = closeout_artifact_ledger(ledger, run_root)
    assert len(closed) == 1
    assert closed[0].changed is True


def test_missing_expected_artifact_produces_warn(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    ledger = tmp_path / "ledger.tsv"
    write_artifact_ledger(ledger, [])

    expected = _expected_inventory(("row-1", "query_structure", "missing.pdb"))
    result = validate_before_consume(ledger, expected, run_id="test-run")
    assert result["overall_status"] == "warn"
    assert result['overall_status'] == 'warn'


# ===================================================================
# Closeout
# ===================================================================

def test_closeout_writes_records_and_detects_changes(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    present = run_root / "present.pdb"
    present.write_bytes(b"ATOM original\n")
    obs = _observed(present, run_root, rel_path="present.pdb")
    ledger = tmp_path / "ledger.tsv"
    write_artifact_ledger(ledger, [obs])

    # Mutate the file before closeout
    present.write_bytes(b"ATOM changed\n")

    closeout = closeout_artifact_ledger(ledger, run_root)
    assert len(closeout) == 1
    assert closeout[0].changed is True


def test_closeout_sees_new_missing_record(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()

    # Write a ledger with an artifact that will go missing
    p = run_root / "willvanish.pdb"
    p.write_bytes(b"ATOM\n")
    obs = _observed(p, run_root, rel_path="willvanish.pdb")
    ledger = tmp_path / "ledger.tsv"
    write_artifact_ledger(ledger, [obs])

    # Delete the file
    p.unlink()

    closeout = closeout_artifact_ledger(ledger, run_root)
    assert len(closeout) == 1
    assert closeout[0].changed is True


# ===================================================================
# Duplicate key rejection
# ===================================================================

def test_duplicate_key_raises_append_error(tmp_path):
    run_root = tmp_path / "run"
    artifact = run_root / "input.pdb"
    run_root.mkdir()
    artifact.write_bytes(b"ATOM\n")
    obs = _observed(artifact, run_root, rel_path="input.pdb")
    ledger = tmp_path / "ledger.tsv"
    write_artifact_ledger(ledger, [obs])

    with pytest.raises(ArtifactLedgerError, match="duplicate"):
        append_artifact_observation(ledger, obs)


# ===================================================================
# CLI wrapper
# ===================================================================

def test_cli_wrapper_runs_closeout(tmp_path):
    run_root = tmp_path / "run"
    artifact = run_root / "input.pdb"
    run_root.mkdir()
    artifact.write_bytes(b"ATOM\n")
    obs = _observed(artifact, run_root, rel_path="input.pdb")
    ledger = write_artifact_ledger(tmp_path / "ledger.tsv", [obs])

    assert main([
        "--ledger-path", str(ledger),
        "--staging-dir", str(run_root),
        "--closeout",
    ]) == 0


def test_cli_validation_returns_nonzero_for_secret_sentinel(monkeypatch, tmp_path):
    ledger = tmp_path / "ledger.tsv"
    write_artifact_ledger(ledger, [])
    monkeypatch.setenv("PRISM_TEST_SECRET_SENTINEL", "secret-value")

    assert main([
        "--ledger-path", str(ledger),
        "--staging-dir", str(tmp_path),
        "--validate",
    ]) == 2
