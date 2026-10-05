import csv
import json

from benchmark.scripts.replay_compact_usalign_batch import (
    compact_alignment_records,
    compact_candidate_audit,
    parse_alignment_filename,
)


def test_parse_alignment_filename_recovers_query_template_and_chain():
    assert parse_alignment_filename("1fgnHL_3i6eEF_E.json") == (
        "1fgnHL", "3i6eEF", "E"
    )


def test_compact_alignment_records_writes_case_summary_and_retains_generated_rows(tmp_path):
    alignment = tmp_path / "alignment"
    status = tmp_path / "status"
    alignment.mkdir()
    status.mkdir()
    (alignment / "1fgnHL_3i6eEF_E.json").write_text(json.dumps({
        "status": "success", "return_code": 0, "match_count": 15,
        "tm_score": 0.6, "tm_score_query": 0.4, "tm_score_ref": 0.6,
    }))
    (alignment / "1fgnHL_3i6eEF_F.json").write_text(json.dumps({
        "status": "alignment_unavailable", "return_code": 1,
    }))

    result = compact_alignment_records(alignment, status, aligner="USalign")

    assert result["raw_record_count"] == 2
    assert (status / "alignment_case_summary.csv").is_file()
    with (status / "alignment_case_summary.csv").open(newline="") as handle:
        row = next(csv.DictReader(handle))
    assert float(row["tm_score_query_mean"]) == 0.4
    assert float(row["tm_score_ref_mean"]) == 0.6
    assert float(row["match_count_mean"]) == 15.0
    assert row["tm_score_query_valid_count"] == "1"
    assert row["tm_score_ref_valid_count"] == "1"
    assert row["match_count_valid_count"] == "1"
    assert row["status_alignment_unavailable"] == "1"


def test_compact_alignment_records_leaves_means_blank_when_no_alignment_is_valid(tmp_path):
    alignment = tmp_path / "alignment"
    status = tmp_path / "status"
    alignment.mkdir()
    status.mkdir()
    (alignment / "1fgnHL_3i6eEF_E.json").write_text(json.dumps({
        "status": "alignment_unavailable", "return_code": None,
        "tm_score_query": 0.0, "tm_score_ref": 0.0, "match_count": 0,
    }))

    compact_alignment_records(alignment, status, aligner="USalign")

    with (status / "alignment_case_summary.csv").open(newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["tm_score_query_mean"] == ""
    assert row["tm_score_ref_mean"] == ""
    assert row["match_count_mean"] == ""
    assert row["tm_score_query_valid_count"] == "0"


def test_compact_candidate_audit_retains_explicit_dual_score_provenance(tmp_path):
    audit = tmp_path / "audit.jsonl"
    status = tmp_path / "status"
    status.mkdir()
    audit.write_text(json.dumps({
        "query_left": "1fgnHL", "query_right": "1tfhA", "template": "3i6eEF",
        "chain_left": "E", "chain_right": "F", "orientation": "o1",
        "status": "generated", "match_count_left": 20, "match_count_right": 21,
        "tm_score_left": 0.6, "tm_score_right": 0.7,
        "metadata": {
            "tm_score_query_left": 0.4, "tm_score_ref_left": 0.6,
            "tm_score_contract_left": "reference_normalized_structure_2",
            "tm_score_query_right": 0.5, "tm_score_ref_right": 0.7,
            "tm_score_contract_right": "reference_normalized_structure_2",
        },
    }) + "\n")

    compact_candidate_audit(
        audit,
        status,
        aligner="USalign",
        input_rows=[
            {
                "pair_id": "rigid_0001",
                "benchmark_set": "rigid",
                "source_row": "1",
                "complex": "1AHW_AB:C",
                "Receptor": "1fgnHL",
                "Ligand": "1tfhA",
            }
        ],
    )

    with (status / "candidate_generated.csv").open(newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["tm_score_query_left"] == "0.4"
    assert row["tm_score_ref_right"] == "0.7"
    assert row["pair_id"] == "rigid_0001"
    assert row["complex"] == "1AHW_AB:C"
