import csv

from benchmark.scripts.build_asset_coverage import build_coverage


def test_build_coverage_reduces_rows_and_preserves_failures(tmp_path):
    source = tmp_path / "assets.tsv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("template_id", "asset_type", "asset_path", "exists", "format_valid", "valid", "fully_resolvable", "missing", "sha256"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({"template_id": "1abcAB", "asset_type": "contact_json", "asset_path": "c", "exists": "True", "format_valid": "True", "valid": "1", "fully_resolvable": "1", "missing": "", "sha256": "a" * 64})
        writer.writerow({"template_id": "1abcAB", "asset_type": "interface_pdb", "asset_path": "p", "exists": "False", "format_valid": "False", "valid": "0", "fully_resolvable": "0", "missing": "interface_pdb", "sha256": ""})
    output = tmp_path / "coverage.tsv"

    rows = build_coverage(source, output)

    assert len(rows) == 1
    assert rows[0]["template_id"] == "1abcAB"
    assert rows[0]["fully_resolvable"] == "0"
    assert rows[0]["missing"] == "interface_pdb"
    assert len(rows[0]["asset_sha256s"]) == 64
