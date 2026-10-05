from src.candidate_ranker import BASELINE_SCORE_VERSION, biological_baseline_score, rank_candidates


def test_baseline_excludes_failed_candidates():
    assert biological_baseline_score({"status": "alignment_failed"}) is None


def test_baseline_prefers_higher_tm_and_coverage():
    low = {
        "status": "generated", "template": "low", "orientation": "o1",
        "tm_score_left": 0.5, "tm_score_right": 0.5,
        "match_count_left": 20, "match_count_right": 20,
    }
    high = {
        "status": "generated", "template": "high", "orientation": "o1",
        "tm_score_left": 0.8, "tm_score_right": 0.8,
        "match_count_left": 40, "match_count_right": 40,
    }
    assert rank_candidates([low, high])[0]["template"] == "high"
    assert rank_candidates([low, high])[0]["baseline_score_version"] == BASELINE_SCORE_VERSION


def test_baseline_penalizes_known_clashes():
    clean = {"status": "generated", "tm_score_left": 0.7, "tm_score_right": 0.7,
             "match_count_left": 30, "match_count_right": 30, "clash_count": 0}
    clash = {**clean, "clash_count": 10}
    assert biological_baseline_score(clean) > biological_baseline_score(clash)


def test_baseline_does_not_invent_coverage_from_fixed_match_count_denominator():
    short_high_tm = {
        "status": "generated", "tm_score_left": 0.8, "tm_score_right": 0.8,
        "match_count_left": 15, "match_count_right": 15,
    }
    long_lower_tm = {
        "status": "generated", "tm_score_left": 0.7, "tm_score_right": 0.7,
        "match_count_left": 50, "match_count_right": 50,
    }
    assert biological_baseline_score(short_high_tm) > biological_baseline_score(long_lower_tm)


def test_baseline_uses_real_coverage_when_both_denominators_are_known():
    row = {
        "status": "generated", "tm_score_left": 0.6, "tm_score_right": 0.6,
        "match_count_left": 20, "match_count_right": 20,
        "match_coverage_left": 50.0, "match_coverage_right": 50.0,
    }
    assert biological_baseline_score(row) == 0.56
