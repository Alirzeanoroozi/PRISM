import argparse
import json

from src.provenance import observe_artifact, write_artifact_ledger
from src.validation_gate import validate_gate


RUN_ID = "prism-20260912-120000-12345-abcdef12"


def _args(ledger, run_root, output, expected_inventory=None):
    return argparse.Namespace(
        ledger=str(ledger),
        run_root=str(run_root),
        run_id=RUN_ID,
        declared_contract_path=None,
        expected_inventory=str(expected_inventory) if expected_inventory else None,
        out=str(output),
    )


def test_validation_gate_rejects_mutated_artifact(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    artifact = run_root / "input.pdb"
    artifact.write_text("ATOM original\n", encoding="utf-8")
    observation = observe_artifact(
        artifact,
        run_id=RUN_ID,
        stage="input",
        dataset_row_id="row-1",
        scientific_role="query_structure",
        run_relative_path="input.pdb",
        produced_by="test",
    )
    ledger = run_root / "artifact_ledger.tsv"
    write_artifact_ledger(ledger, [observation])
    artifact.write_text("ATOM changed\n", encoding="utf-8")

    output = tmp_path / "validation.json"
    assert validate_gate(_args(ledger, run_root, output)) == 2
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["overall_status"] == "fail"
    assert result["artifact_checks"][0]["status"] == "mismatch"


def test_validation_gate_rejects_duplicate_primary_key(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    artifact = run_root / "input.pdb"
    artifact.write_text("ATOM original\n", encoding="utf-8")
    observation = observe_artifact(
        artifact,
        run_id=RUN_ID,
        stage="input",
        dataset_row_id="row-1",
        scientific_role="query_structure",
        run_relative_path="input.pdb",
        produced_by="test",
    )
    ledger = run_root / "artifact_ledger.tsv"
    write_artifact_ledger(ledger, [observation])
    lines = ledger.read_text(encoding="utf-8").splitlines()
    ledger.write_text("\n".join(lines + [lines[1]]) + "\n", encoding="utf-8")

    output = tmp_path / "validation.json"
    assert validate_gate(_args(ledger, run_root, output)) == 2
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["overall_status"] == "fail"


def test_validation_gate_uses_declared_inventory_to_detect_missing_artifact(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    ledger = run_root / "artifact_ledger.tsv"
    write_artifact_ledger(ledger, [])
    expected = tmp_path / "expected.tsv"
    expected.write_text(
        "dataset_row_id\tscientific_role\trun_relative_path\n"
        "row-1\tquery_structure\tinput.pdb\n",
        encoding="utf-8",
    )

    output = tmp_path / "validation.json"
    assert validate_gate(_args(ledger, run_root, output, expected)) == 1
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["overall_status"] == "warn"
    assert result["artifact_checks"][0]["status"] == "missing"


def test_validation_gate_rejects_invalid_declared_inventory(tmp_path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    ledger = run_root / "artifact_ledger.tsv"
    write_artifact_ledger(ledger, [])
    expected = tmp_path / "expected.tsv"
    expected.write_text("dataset_row_id\tscientific_role\nrow-1\tquery_structure\n", encoding="utf-8")

    output = tmp_path / "validation.json"
    assert validate_gate(_args(ledger, run_root, output, expected)) == 2
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["overall_status"] == "fail"
    assert "run_relative_path" in result["expected_inventory_error"]
