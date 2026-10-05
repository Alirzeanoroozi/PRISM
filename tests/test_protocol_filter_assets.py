from benchmark.scripts.build_protocol_filter_assets import parse_legacy_contacts, parse_legacy_hotspots


def test_legacy_asset_parsers_ignore_comments_and_keep_chain_residue_identity(tmp_path):
    hotspot = tmp_path / "hotspot"
    contact = tmp_path / "contact.txt"
    hotspot.write_text("# header\nB.L.46 L\nD.F.104 F\n", encoding="utf-8")
    contact.write_text("# header\nB.46 D.104\n", encoding="utf-8")
    assert parse_legacy_hotspots(hotspot) == [("B", "L", "46"), ("D", "F", "104")]
    assert parse_legacy_contacts(contact) == [("B.46", "D.104")]
