import pytest

from benchmark.scripts.train_reranker import _ranking_quality, train


def test_ranking_quality_is_averaged_per_native_complex():
    rows = [
        {"native_complex_id": "a", "native_like": 0, "dockq": 0.1},
        {"native_complex_id": "a", "native_like": 1, "dockq": 0.8},
        {"native_complex_id": "b", "native_like": 1, "dockq": 0.7},
    ]

    quality = _ranking_quality(rows, [0.9, 0.1, 0.2])

    assert quality["group_count"] == 2.0
    assert quality["top1_native_like"] == 0.5
    assert quality["top1_dockq"] == pytest.approx(0.4)


def test_training_entrypoint_reports_missing_optional_dependency(tmp_path):
    data = tmp_path / "labeled.csv"
    data.write_text(
        "native_complex_id,status,dockq,tm_score_left,tm_score_right,match_count_left,match_count_right\n"
        "a,generated,0.3,0.7,0.7,20,20\n"
        "b,generated,0.1,0.4,0.4,10,10\n"
    )
    try:
        train(data, tmp_path / "model.pkl", tmp_path / "metrics.json")
    except RuntimeError as exc:
        assert "requires scikit-learn" in str(exc)
    except ValueError as exc:
        pytest.skip(f"scikit-learn is installed but the tiny fixture cannot train: {exc}")


def test_training_entrypoint_writes_model_and_metrics_when_available(tmp_path):
    data = tmp_path / "labeled.csv"
    rows = [
        ("a", "generated", 0.8, 0.8, 0.8, 40, 40),
        ("b", "generated", 0.7, 0.7, 0.7, 35, 35),
        ("c", "generated", 0.1, 0.2, 0.2, 10, 10),
        ("d", "generated", 0.05, 0.1, 0.1, 5, 5),
        ("e", "generated", 0.9, 0.85, 0.85, 45, 45),
        ("f", "generated", 0.02, 0.15, 0.15, 6, 6),
    ]
    data.write_text(
        "native_complex_id,status,dockq,tm_score_left,tm_score_right,match_count_left,match_count_right\n"
        + "\n".join(",".join(map(str, row)) for row in rows)
        + "\n"
    )
    try:
        metrics = train(data, tmp_path / "model.pkl", tmp_path / "metrics.json", seed=3)
    except RuntimeError:
        pytest.skip("scikit-learn is unavailable")
    assert metrics["train_rows"] > 0
    assert (tmp_path / "model.pkl").is_file()
    assert (tmp_path / "metrics.json").is_file()


def test_training_reports_learned_and_baseline_top1_quality(tmp_path):
    data = tmp_path / "labeled.csv"
    rows = [
        ("a", "generated", 0.8, 0.8, 0.8, 40, 40),
        ("b", "generated", 0.1, 0.2, 0.2, 10, 10),
        ("c", "generated", 0.7, 0.7, 0.7, 35, 35),
        ("d", "generated", 0.05, 0.1, 0.1, 5, 5),
    ]
    data.write_text(
        "native_complex_id,status,dockq,tm_score_left,tm_score_right,match_count_left,match_count_right\n"
        + "\n".join(",".join(map(str, row)) for row in rows)
        + "\n"
    )
    try:
        metrics = train(
            data,
            tmp_path / "model.pkl",
            tmp_path / "metrics.json",
            seed=0,
            test_fraction=0.5,
        )
    except RuntimeError:
        pytest.skip("scikit-learn is unavailable")
    assert "baseline_test_top1_dockq" in metrics
    assert "learned_test_top1_dockq" in metrics
    assert metrics["accepted_candidate_coverage"] == 1.0


def test_training_ignores_retained_unlabeled_rows(tmp_path):
    data = tmp_path / "labeled.csv"
    data.write_text(
        "native_complex_id,status,dockq,label_status,tm_score_left,tm_score_right,match_count_left,match_count_right\n"
        "a,generated,0.3,labeled,0.7,0.7,20,20\n"
        "a,generated,,unlabeled,0.0,0.0,0,0\n"
    )
    try:
        train(data, tmp_path / "model.pkl", tmp_path / "metrics.json")
    except RuntimeError:
        pytest.skip("scikit-learn is unavailable")
    except ValueError as exc:
        assert "two native complexes" in str(exc)
