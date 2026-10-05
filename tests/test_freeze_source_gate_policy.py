import json
from pathlib import Path

from benchmark.scripts.freeze_source_gate_policy import freeze_policy


def test_freeze_policy_requires_matching_failure_identity(tmp_path: Path):
    summary = tmp_path / "summary.json"
    failures = tmp_path / "validation_failures.tsv"
    output = tmp_path / "policy.json"
    summary.write_text(
        json.dumps(
            {
                "expected_row_count": 2,
                "validation_failed_dataset_row_count": 1,
                "validation_failed_dataset_row_ids": ["rigid:000002"],
            }
        )
        + "\n",
        encoding="utf-8",
    )
    failures.write_text(
        "dataset_row_id\tsource_role\tnative_complex\traw_receptor_selector\traw_ligand_selector\tarchive_prefix\tparse_status\texpected_chain_ids\tpolymer_chain_ids\tchain_set_status\terror\n"
        "rigid:000002\tpipeline_receptor\t2ABC_A:B\t2ABC_A\t2ABC_B\t2abc\tok\tA\tA\tok\t\n",
        encoding="utf-8",
    )

    policy = freeze_policy(summary, failures, output)

    assert policy["decision"]["status"] == "blocked_source_authority"
    assert policy["decision"]["strict_source_clean_row_count"] == 1
    assert policy["decision"]["audit_only_row_count"] == 1
    assert policy["decision"]["automatic_orientation_swap"] is False
    assert policy["reason"]["confirmed"].startswith("All 2 curated")
    assert json.loads(output.read_text(encoding="utf-8"))["failed_rows"][0]["dataset_row_id"] == "rigid:000002"


def test_freeze_policy_rejects_duplicate_role_records(tmp_path: Path):
    summary = tmp_path / "summary.json"
    failures = tmp_path / "validation_failures.tsv"
    summary.write_text(
        json.dumps(
            {
                "expected_row_count": 2,
                "validation_failed_dataset_row_count": 1,
                "validation_failed_dataset_row_ids": ["rigid:000002"],
            }
        )
        + "\n",
        encoding="utf-8",
    )
    header = "dataset_row_id\tsource_role\tnative_complex\traw_receptor_selector\traw_ligand_selector\tarchive_prefix\tparse_status\texpected_chain_ids\tpolymer_chain_ids\tchain_set_status\terror\n"
    row = "rigid:000002\tpipeline_receptor\t2ABC_A:B\t2ABC_A\t2ABC_B\t2abc\tok\tA\tA\tok\t\n"
    failures.write_text(header + row + row, encoding="utf-8")

    try:
        freeze_policy(summary, failures, tmp_path / "policy.json")
    except ValueError as exc:
        assert "duplicate" in str(exc)
    else:
        raise AssertionError("duplicate validation failure records must be rejected")
