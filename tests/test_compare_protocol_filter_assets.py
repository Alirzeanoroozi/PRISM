import csv
import json

from benchmark.scripts.compare_protocol_filter_assets import compare_assets


def test_compare_assets_reports_exact_and_missing(tmp_path):
    modern = tmp_path / "modern"
    legacy = tmp_path / "legacy"
    for root in (modern, legacy):
        (root / "hotspots").mkdir(parents=True)
        (root / "contacts").mkdir(parents=True)
    template_list = tmp_path / "templates.txt"
    template_list.write_text("1abcAB\n1defCD\n", encoding="utf-8")
    (modern / "hotspots" / "1abcAB.json").write_text(json.dumps({"A": [[1, "ALA"]], "B": [[2, "GLY"]]}))
    (legacy / "hotspots" / "1abcAB.json").write_text(json.dumps([["A", "ALA", "1"], ["B", "GLY", "2"]]))
    (modern / "contacts" / "1abcAB.json").write_text(json.dumps([[1, 2]]))
    (legacy / "contacts" / "1abcAB.json").write_text(json.dumps([["A.A.1", "B.G.2"]]))
    (legacy / "hotspots" / "1defCD.json").write_text(json.dumps([["C", "ALA", "1"]]))
    (legacy / "contacts" / "1defCD.json").write_text(json.dumps([["C.A.1", "D.G.2"]]))

    output = tmp_path / "parity.tsv"
    rows = compare_assets(template_list, modern, legacy, output)

    assert rows[0]["parity_status"] == "exact"
    assert rows[1]["parity_status"] == "current_assets_missing"
    with output.open(newline="", encoding="utf-8") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 2
