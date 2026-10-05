import csv
import json

from benchmark.scripts.replay_published_filter import replay


def test_replay_records_missing_alignment_and_asset_hashes(tmp_path):
    template_list = tmp_path / "templates.txt"
    template_list.write_text("1abcAB\n", encoding="utf-8")
    inputs = tmp_path / "inputs.csv"
    inputs.write_text("Receptor,Ligand\n1aaaA,1bbbB\n", encoding="utf-8")
    align = tmp_path / "align"
    align.mkdir()
    payload = {"match_dict": {"A.A.1": "A.A.1"}}
    (align / "1aaaA_1abcAB_A.json").write_text(json.dumps(payload))
    assets = tmp_path / "assets"
    (assets / "hotspots").mkdir(parents=True)
    (assets / "contacts").mkdir()
    (assets / "hotspots" / "1abcAB.json").write_text(json.dumps([["A", "ALA", "1"]]))
    (assets / "contacts" / "1abcAB.json").write_text(json.dumps([["A.A.1", "B.G.2"]]))
    output = tmp_path / "replay.tsv"

    rows = replay(inputs, template_list, align, assets, output)

    assert len(rows) == 2
    assert all(row["status"] == "failed" for row in rows)
    assert any(row["reason"] == "alignment_missing" for row in rows)
    with output.open(newline="", encoding="utf-8") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 2
