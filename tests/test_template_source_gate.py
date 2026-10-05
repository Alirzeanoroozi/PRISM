import json

from benchmark.scripts.build_template_source_gate import (
    classify_template,
    confirmatory_row_eligibility,
    load_source_gate_policy,
    sequence_similarity,
)


def test_exact_source_pdb_chain_is_excluded():
    row = classify_template(target_pdb="1abc", target_chains="A", template_id="1abcA", target_sequence="AAAA", template_sequence="AAAT")
    assert row["exclusion_reason"] == "exact_self_hit"


def test_sequence_self_hit_and_close_homolog_are_excluded():
    exact = classify_template(target_pdb="1abc", target_chains="A", template_id="2defB", target_sequence="AAAA", template_sequence="AAAA")
    homolog = classify_template(target_pdb="1abc", target_chains="A", template_id="2defB", target_sequence="AAAAAAAAAA", template_sequence="AAAAAAAATT")
    assert exact["exclusion_reason"] == "exact_self_hit"
    assert homolog["exclusion_reason"] == "homology_excluded_primary"
    assert homolog["excluded_identity_gt_50"] is True


def test_low_coverage_fragment_is_excluded_by_the_preregistered_shorter_coverage_rule():
    row = classify_template(target_pdb="1abc", target_chains="A", template_id="2defB", target_sequence="AAAAAAAAAA", template_sequence="AAAA")
    assert row["exclusion_reason"] == "homology_excluded_primary"
    assert row["query_coverage_percent"] == 40.0
    assert sequence_similarity("AAAA", "AAAA")["identity_percent"] == 100.0


def test_source_policy_blocks_audit_rows_and_unauthorized_confirmatory_runs(tmp_path):
    path = tmp_path / "policy.json"
    path.write_text(json.dumps({"decision": {
        "status": "blocked_source_authority",
        "confirmatory_run_authorized": False,
        "excluded_dataset_row_ids": ["rigid:000064"],
    }}), encoding="utf-8")
    decision = load_source_gate_policy(path)
    assert confirmatory_row_eligibility("rigid:000064", decision) == (False, "source_gate_audit_only")
    assert confirmatory_row_eligibility("rigid:000001", decision) == (False, "source_gate_not_authorized")
