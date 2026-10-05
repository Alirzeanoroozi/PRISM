import pytest

from src.ranking_metrics import enrichment_factor, spearman_score, top_k_success


def _rows():
    return [
        {"baseline_score": 0.9, "native_like": 1},
        {"baseline_score": 0.8, "native_like": 0},
        {"baseline_score": 0.2, "native_like": 0},
        {"baseline_score": 0.1, "native_like": 1},
    ]


def test_top_k_and_enrichment_are_computed_from_labels():
    assert top_k_success(_rows(), 0.25) == pytest.approx(1.0)
    assert enrichment_factor(_rows(), 0.25) == pytest.approx(2.0)


def test_spearman_is_positive_when_native_candidate_ranks_high():
    assert spearman_score(_rows()) > 0.0
