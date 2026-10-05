import json

from src.template_filtering import count_matched_complementary_contacts, evaluate_hotspots, evaluate_protocol_candidate, load_filter_assets


def _matches(prefix):
    return {f"A.A.{index}": f"Q.Q.{index}" for index in range(1, 7)}


def test_hotspot_requires_asset_and_matching_residue():
    matches = _matches("A")
    assert evaluate_hotspots(matches, None).reason == "hotspot_assets_missing"
    assert evaluate_hotspots(matches, [(1, "A")]).passed
    assert not evaluate_hotspots(matches, [(99, "A")]).passed


def test_contacts_count_only_when_both_partner_matches_exist():
    left, right = _matches("A"), _matches("B")
    contacts = [(("A", "1"), ("A", "1"))]
    assert count_matched_complementary_contacts(left, right, contacts) == 0
    assert count_matched_complementary_contacts({"A.A.1": "Q.Q.1"}, {"A.A.1": "Q.Q.1"}, [("A.A.1", "A.A.1")]) == 1


def test_published_candidate_requires_five_contacts_and_both_hotspot_passes():
    left = {f"A.A.{index}": f"Q.Q.{index}" for index in range(1, 7)}
    right = {f"B.B.{index}": f"R.R.{index}" for index in range(1, 7)}
    contacts = [(f"A.A.{index}", f"B.B.{index}") for index in range(1, 6)]
    passed = evaluate_protocol_candidate(left, right, [(1, "A")], [(1, "B")], contacts)
    assert passed.passed and passed.complementary_contacts == 5
    failed = evaluate_protocol_candidate(left, right, [(1, "A")], [(1, "B")], contacts[:4])
    assert failed.reason == "complementary_contact_threshold_failed"


def test_filter_assets_are_loaded_with_hashes_and_missing_assets_fail_closed(tmp_path):
    (tmp_path / "hotspots").mkdir()
    (tmp_path / "contacts").mkdir()
    (tmp_path / "hotspots" / "1abcAB.json").write_text(json.dumps([["A", "ALA", "1"]]), encoding="utf-8")
    (tmp_path / "contacts" / "1abcAB.json").write_text(json.dumps([["A.1", "B.2"]]), encoding="utf-8")
    loaded = load_filter_assets("1abcAB", tmp_path)
    assert len(loaded["hotspots_sha256"]) == 64
    try:
        load_filter_assets("missing", tmp_path)
    except FileNotFoundError:
        pass
    else:
        raise AssertionError("missing protocol assets must fail closed")


def test_modern_assets_are_normalized_by_chain_and_numeric_contacts(tmp_path):
    (tmp_path / "hotspots").mkdir()
    (tmp_path / "contacts").mkdir()
    (tmp_path / "hotspots" / "1abcAB.json").write_text(
        json.dumps({"A": [["1", "ALA"]], "B": [["2", "GLY"]]}),
        encoding="utf-8",
    )
    (tmp_path / "contacts" / "1abcAB.json").write_text(
        json.dumps([[1, 2]]), encoding="utf-8"
    )
    loaded = load_filter_assets("1abcAB", tmp_path)
    assert loaded["asset_format"] == "modern_json"
    assert loaded["hotspots_by_chain"]["A"] == [["A", "A", "1"]]
    assert loaded["hotspots_by_chain"]["B"] == [["B", "G", "2"]]
    assert loaded["contacts"] == [["A..1", "B..2"]]
    left = {"A.A.1": "Q.A.7"}
    right = {"B.G.2": "R.G.8"}
    assert evaluate_protocol_candidate(
        left, right, loaded["hotspots_by_chain"]["A"],
        loaded["hotspots_by_chain"]["B"], loaded["contacts"], minimum_contacts=1
    ).passed
