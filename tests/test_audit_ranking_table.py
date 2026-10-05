import csv

from benchmark.scripts.audit_ranking_table import audit


def _write(path, rows, delimiter=","):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def test_audit_ranking_table_passes_complete_canonical_identity(tmp_path):
    candidate = {
        "dataset_row_id": "rigid:000001", "source_model_sha256": "model-hash",
        "observed_source_model_sha256": "model-hash", "label_model_sha256": "model-hash",
        "status": "refinement_accepted", "label_status": "labeled",
        "label_metric": "dockq_cross_mean", "alignment_left_sha256": "left-json",
        "alignment_right_sha256": "right-json", "alignment_left_raw_output_sha256": "left-raw",
        "alignment_right_raw_output_sha256": "right-raw", "alignment_left_aligner": "GTalign",
        "alignment_right_aligner": "GTalign", "alignment_left_status": "success",
        "alignment_right_status": "success", "match_count_left": "10",
        "match_count_right": "11", "tm_score_left": "0.5", "tm_score_right": "0.6",
        "match_coverage_left": "50.0", "match_coverage_right": "55.0",
        "template_coverage_status": "available", "template_interface_sha256": "template-hash",
        "template_size_left": "20", "template_size_right": "20",
        "batch": "batch_0001",
    }
    score = {
        "dataset_row_id": "rigid:000001", "source_model_sha256": "model-hash",
        "score_status": "scored", "source_gate_status": "strict_clean",
    }
    candidates = tmp_path / "candidates.csv"
    scores = tmp_path / "scores.tsv"
    _write(candidates, [candidate])
    _write(scores, [score], delimiter="\t")

    result = audit(candidates, scores)

    assert result["audit_status"] == "passed"
    assert result["candidate_count"] == 1
    assert result["dataset_row_count"] == 1


def test_audit_ranking_table_accepts_cross_only_score_identity(tmp_path):
    candidate = {
        "dataset_row_id": "rigid:000001", "source_model_sha256": "model-hash",
        "observed_source_model_sha256": "model-hash", "label_model_sha256": "model-hash",
        "status": "refinement_accepted", "label_status": "labeled",
        "label_metric": "dockq_cross_mean", "alignment_left_sha256": "left-json",
        "alignment_right_sha256": "right-json", "alignment_left_raw_output_sha256": "left-raw",
        "alignment_right_raw_output_sha256": "right-raw", "alignment_left_aligner": "GTalign",
        "alignment_right_aligner": "GTalign", "alignment_left_status": "success",
        "alignment_right_status": "success", "match_count_left": "10", "match_count_right": "11",
        "tm_score_left": "0.5", "tm_score_right": "0.6", "match_coverage_left": "50.0",
        "match_coverage_right": "55.0", "template_coverage_status": "available",
        "template_interface_sha256": "template-hash", "template_size_left": "20",
        "template_size_right": "20", "batch": "batch_0001",
    }
    score = {
        "dataset_row_id": "rigid:000001", "source_model_sha256": "model-hash",
        "score_status": "scored_cross_only", "source_gate_status": "strict_clean",
    }
    candidates = tmp_path / "candidates.csv"
    scores = tmp_path / "scores.tsv"
    _write(candidates, [candidate])
    _write(scores, [score], delimiter="\t")

    assert audit(candidates, scores)["audit_status"] == "passed"


def test_audit_ranking_table_retains_explicit_non_scoreable_candidate(tmp_path):
    template = {
        "observed_source_model_sha256": "model-hash", "label_model_sha256": "model-hash",
        "status": "refinement_accepted", "label_status": "labeled",
        "label_metric": "dockq_cross_mean", "alignment_left_sha256": "left-json",
        "alignment_right_sha256": "right-json", "alignment_left_raw_output_sha256": "left-raw",
        "alignment_right_raw_output_sha256": "right-raw", "alignment_left_aligner": "GTalign",
        "alignment_right_aligner": "GTalign", "alignment_left_status": "success",
        "alignment_right_status": "success", "match_count_left": "10", "match_count_right": "11",
        "tm_score_left": "0.5", "tm_score_right": "0.6", "match_coverage_left": "50.0",
        "match_coverage_right": "55.0", "template_coverage_status": "available",
        "template_interface_sha256": "template-hash", "template_size_left": "20",
        "template_size_right": "20", "batch": "batch_0001",
    }
    candidate = {
        **template, "dataset_row_id": "rigid:000023", "source_model_sha256": "model-hash",
        "label_status": "unlabeled", "label_metric": "", "label_model_sha256": "",
    }
    labeled = {**template, "dataset_row_id": "rigid:000024", "source_model_sha256": "other-hash",
               "observed_source_model_sha256": "other-hash", "label_model_sha256": "other-hash"}
    score = {
        "dataset_row_id": "rigid:000024", "source_model_sha256": "other-hash",
        "score_status": "scored", "source_gate_status": "strict_clean",
    }
    not_scoreable = {
        "dataset_row_id": "rigid:000023", "source_model_sha256": "model-hash",
        "score_status": "not_scoreable", "source_gate_status": "strict_clean",
    }
    candidates = tmp_path / "candidates.csv"
    scores = tmp_path / "scores.tsv"
    _write(candidates, [candidate, labeled])
    _write(scores, [score, not_scoreable], delimiter="\t")

    result = audit(candidates, scores)

    assert result["audit_status"] == "passed"
    assert result["checks"]["explicit_not_scoreable_retained"] == 1
    assert result["retained_unrankable_candidate_count"] == 1


def test_audit_ranking_table_fails_unsafe_or_unlabeled_rows(tmp_path):
    candidate = {
        "dataset_row_id": "rigid:000001", "source_model_sha256": "model-hash",
        "observed_source_model_sha256": "different", "label_model_sha256": "",
        "status": "alignment_failed", "label_status": "unlabeled", "label_metric": "",
        "alignment_left_sha256": "", "alignment_right_sha256": "",
        "alignment_left_raw_output_sha256": "", "alignment_right_raw_output_sha256": "",
        "alignment_left_aligner": "", "alignment_right_aligner": "",
        "alignment_left_status": "", "alignment_right_status": "",
        "match_count_left": "", "match_count_right": "", "tm_score_left": "", "tm_score_right": "",
        "match_coverage_left": "", "match_coverage_right": "",
        "template_coverage_status": "missing_template_interface", "template_interface_sha256": "",
        "template_size_left": "", "template_size_right": "",
        "batch": "batch_0001",
    }
    score = {
        "dataset_row_id": "rigid:000001", "source_model_sha256": "model-hash",
        "score_status": "scored", "source_gate_status": "strict_clean",
    }
    candidates = tmp_path / "candidates.csv"
    scores = tmp_path / "scores.tsv"
    _write(candidates, [candidate])
    _write(scores, [score], delimiter="\t")

    result = audit(candidates, scores)

    assert result["audit_status"] == "failed"
    assert "model_hash_mismatch" in result["errors"]
    assert "label_not_attached" in result["errors"]
