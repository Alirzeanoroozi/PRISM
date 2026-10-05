import pytest

from benchmark.scripts.evaluate_ranked_candidates import summarize


def test_summarize_reports_top1_and_oracle_per_group():
    rows = [
        {"dataset_row_id": "a", "baseline_rank": "1", "dockq": "0.1", "native_like": "0", "label_status": "labeled"},
        {"dataset_row_id": "a", "baseline_rank": "2", "dockq": "0.8", "native_like": "1", "label_status": "labeled"},
        {"dataset_row_id": "b", "baseline_rank": "1", "dockq": "0.7", "native_like": "1", "label_status": "labeled"},
        {"dataset_row_id": "b", "baseline_rank": "", "dockq": "", "native_like": "", "label_status": "unlabeled"},
    ]

    result = summarize(rows)

    assert result["ranking_group_count"] == 2
    assert result["oracle_native_like_group_count"] == 2
    assert result["top1_native_like_group_count"] == 1
    assert result["top1_native_like_rate"] == 0.5
    assert result["median_top1_dockq"] == pytest.approx(0.4)
    assert result["baseline_score_version"] == "unversioned"


def test_summarize_reports_one_score_version_and_rejects_mixed_versions():
    base = {
        "dataset_row_id": "a", "baseline_rank": "1", "dockq": "0.5",
        "native_like": "1", "label_status": "labeled",
        "baseline_score_version": "biological-baseline/v2-real-coverage",
    }
    assert summarize([base])["baseline_score_version"] == base["baseline_score_version"]

    with pytest.raises(ValueError, match="mixes baseline score versions"):
        summarize([base, {
            **base, "dataset_row_id": "b", "baseline_score_version": "legacy-v1",
        }])
