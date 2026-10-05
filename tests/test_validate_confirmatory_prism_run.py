import csv
import json
from pathlib import Path

from benchmark.scripts.validate_confirmatory_prism_run import validate_confirmatory_preflight


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _policy(path: Path, *, authorized: bool) -> None:
    path.write_text(json.dumps({"decision": {
        "status": "authorized" if authorized else "blocked_source_authority",
        "confirmatory_run_authorized": authorized,
        "excluded_dataset_row_ids": [],
    }}), encoding="utf-8")


def test_blocked_source_policy_prevents_confirmatory_ready_status(tmp_path):
    policy = tmp_path / "policy.json"
    _policy(policy, authorized=False)
    templates = tmp_path / "templates.txt"
    templates.write_text("1abcAB\n", encoding="utf-8")
    eligible = tmp_path / "eligible.tsv"
    _write_tsv(eligible, [{
        "dataset_row_id": "rigid:000001", "template_list_sha256": __import__("hashlib").sha256(templates.read_bytes()).hexdigest(),
        "eligible_template_count": "0", "eligible_templates": "",
    }])

    report = validate_confirmatory_preflight(policy_path=policy, eligible_path=eligible, template_list_path=templates, expected_row_count=1)
    assert report["status"] == "blocked"
    assert report["confirmatory_run_authorized"] is False


def test_authorized_preflight_requires_nonempty_unique_template_exposure(tmp_path):
    policy = tmp_path / "policy.json"
    _policy(policy, authorized=True)
    templates = tmp_path / "templates.txt"
    templates.write_text("1abcAB\n", encoding="utf-8")
    template_hash = __import__("hashlib").sha256(templates.read_bytes()).hexdigest()
    eligible = tmp_path / "eligible.tsv"
    _write_tsv(eligible, [{
        "dataset_row_id": "rigid:000001", "template_list_sha256": template_hash,
        "eligible_template_count": "1", "eligible_templates": "1abcAB",
    }])

    report = validate_confirmatory_preflight(policy_path=policy, eligible_path=eligible, template_list_path=templates, expected_row_count=1)
    assert report["status"] == "ready"


def test_authorized_preflight_rejects_duplicate_templates(tmp_path):
    policy = tmp_path / "policy.json"
    _policy(policy, authorized=True)
    templates = tmp_path / "templates.txt"
    templates.write_text("1abcAB\n", encoding="utf-8")
    template_hash = __import__("hashlib").sha256(templates.read_bytes()).hexdigest()
    eligible = tmp_path / "eligible.tsv"
    _write_tsv(eligible, [{
        "dataset_row_id": "rigid:000001", "template_list_sha256": template_hash,
        "eligible_template_count": "2", "eligible_templates": "1abcAB,1abcAB",
    }])

    report = validate_confirmatory_preflight(policy_path=policy, eligible_path=eligible, template_list_path=templates)
    assert report["status"] == "invalid"
    assert "duplicate_eligible_template:rigid:000001" in report["issues"]
