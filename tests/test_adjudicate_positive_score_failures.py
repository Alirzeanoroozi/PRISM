import csv

from benchmark.scripts.adjudicate_positive_score_failures import adjudicate, failure_class


def test_failure_class_is_deterministic():
    assert failure_class("ValueError: Buffer has wrong number of dimensions (expected 2, got 1)") == "dockq_runtime_buffer_dimensions"
    assert failure_class("no requested cross interfaces") == "missing_cross_interface"


def test_adjudication_preserves_failures_and_disallows_implicit_retries(tmp_path):
    source = tmp_path / "scores.tsv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("score_status", "pair_id", "score_error"), delimiter="\t")
        writer.writeheader()
        writer.writerow({"score_status": "score_failed", "pair_id": "pair-1", "score_error": "Buffer has wrong number of dimensions"})
        writer.writerow({"score_status": "scored", "pair_id": "pair-2", "score_error": ""})
    output = tmp_path / "adjudication.tsv"
    rows = adjudicate([source], output)
    assert len(rows) == 1
    assert rows[0]["failure_class"] == "dockq_runtime_buffer_dimensions"
    assert rows[0]["retry_authorized"] == "false"
