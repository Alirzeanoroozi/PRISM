import pytest

from src.ranking_data import grouped_split, native_like_label, validate_training_rows


def _row(group, dockq=0.3):
    return {
        "native_complex_id": group,
        "dockq": dockq,
        "tm_score_left": 0.6,
        "tm_score_right": 0.7,
        "match_count_left": 20,
        "match_count_right": 18,
    }


def test_native_like_label_uses_documented_threshold():
    assert native_like_label(0.23) == 1
    assert native_like_label(0.229) == 0
    with pytest.raises(ValueError):
        native_like_label(1.1)


def test_validation_rejects_unlabeled_rows():
    with pytest.raises(ValueError, match="missing dockq"):
        validate_training_rows([_row("complex-a", dockq=None)])


def test_grouped_split_keeps_all_decoys_together():
    rows = [_row("a"), _row("a", 0.1), _row("b"), _row("c")]
    validated = validate_training_rows(rows)
    train, test = grouped_split(validated, test_fraction=0.34, seed=7)
    train_groups = {validated[i]["native_complex_id"] for i in train}
    test_groups = {validated[i]["native_complex_id"] for i in test}
    assert train_groups.isdisjoint(test_groups)
