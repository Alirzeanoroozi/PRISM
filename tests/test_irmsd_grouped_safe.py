import benchmark.scripts.irmsd_grouped_safe as safe_irmsd


def test_three_identical_chains_expand_without_mutating_groups(monkeypatch):
    monkeypatch.setattr(
        safe_irmsd.legacy,
        "parseChainResiduesFromStructure",
        lambda _path, _chain: ([], "AAAA"),
    )
    monkeypatch.setattr(safe_irmsd.legacy, "almostIdentical", lambda left, right: left == right)

    orders = safe_irmsd.symmetric_chain_orders("model.pdb", "ABC")

    assert len(orders) == 6
    assert set(orders) == {"ABC", "ACB", "BAC", "BCA", "CAB", "CBA"}


def test_independent_symmetry_groups_form_cartesian_product(monkeypatch):
    sequences = {"A": "AAAA", "B": "AAAA", "C": "CCCC", "D": "CCCC"}
    monkeypatch.setattr(
        safe_irmsd.legacy,
        "parseChainResiduesFromStructure",
        lambda _path, chain: ([], sequences[chain]),
    )
    monkeypatch.setattr(safe_irmsd.legacy, "almostIdentical", lambda left, right: left == right)

    orders = safe_irmsd.symmetric_chain_orders("model.pdb", "ABCD")

    assert set(orders) == {"ABCD", "ABDC", "BACD", "BADC"}
