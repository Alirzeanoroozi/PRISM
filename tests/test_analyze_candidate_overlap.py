import csv

from benchmark.scripts.analyze_candidate_overlap import analyze


def _write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_overlap_is_casewise_and_fail_closed_for_missing_case_keys(tmp_path):
    fields = {
        "status": "generated", "pair_id": "case1", "template": "1abcAB",
        "query_left": "qA", "query_right": "qB", "orientation": "o1",
        "chain_left": "A", "chain_right": "B",
    }
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    _write(left, [fields, {**fields, "template": "2defAB"}])
    _write(right, [fields, {**fields, "template": "3ghiAB"}, {**fields, "pair_id": "", "template": "ignored"}])

    result = analyze([f"tmalign={left}", f"multiprot={right}"], tmp_path / "out")

    assert result["case_count"] == 1
    rows = list(csv.DictReader((tmp_path / "out/candidate_overlap_summary.tsv").open(), delimiter="\t"))
    assert rows[0]["shared_candidates"] == "1"
    assert rows[0]["unique_left_candidates"] == "1"
    assert rows[0]["unique_right_candidates"] == "1"

