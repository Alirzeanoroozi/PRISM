import csv

from benchmark.scripts.build_reference_crosswalk import FIELDS, ROWS, write_crosswalk


def test_reference_crosswalk_is_schema_first_and_preserves_blocked_claims(tmp_path):
    path = write_crosswalk(tmp_path / "reference_crosswalk.tsv")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert tuple(rows[0]) == FIELDS
    assert len(rows) == len(ROWS)
    assert any(row["claim_id"] == "paper-unbiased-templates" and row["status"] == "blocked" for row in rows)
    assert any(row["status"] == "implementation-specific" for row in rows)
