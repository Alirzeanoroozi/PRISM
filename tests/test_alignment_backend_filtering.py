from src.alignment_gtalign import parse_gtalign_hits
from src.alignment_gtalign import parse_gtalign_output_text


def test_gtalign_filter_keeps_only_parsed_hits_above_both_thresholds():
    raw = [
        {"ref_path": "pass", "tm_ref": 0.5, "tm_query": 0.4, "aligned_length": 15},
        {"ref_path": "low-score", "tm_ref": 0.39, "tm_query": 0.9, "aligned_length": 20},
        {"ref_path": "short", "tm_ref": 0.9, "tm_query": 0.9, "aligned_length": 14},
        None,
    ]
    accepted = parse_gtalign_hits(raw)
    assert [hit["ref_path"] for hit in accepted] == ["pass"]
    assert accepted[0]["status"] == "accepted"


def test_gtalign_019_numbered_hit_block_is_parsed():
    raw = """ Query (/queries/query.asa.pdb):
 /queries/query.asa.pdb Chn:A
 Searched: /refs

     1 .../refs/template_A_int.pdb Chn:A 1.0000 1.0000  0.00  4  1-4  1-4  4

1.
>/refs/template_A_int.pdb Chn:A
  Length: Refn. = 4, Query = 4

 TM-score (Refn./Query) = 0.80000 / 0.70000, d0 (Refn./Query) = 2.00 / 2.00,  RMSD = 0.00 A
 Identities = 4/4 (100%), Matched = 4/4 (100%)

Query:     1 AAAA 4
             AAAA
Refn.:     1 AAAA 4
struct         tttt

 Rotation [3,3] and translation [3,1] for Query:
    1.0  0.0  0.0  0.0
    0.0  1.0  0.0  0.0
    0.0  0.0  1.0  0.0
"""
    query, hits = parse_gtalign_output_text(raw)
    assert query == "/queries/query.asa.pdb"
    assert len(hits) == 1
    assert hits[0]["ref_path"] == "/refs/template_A_int.pdb"
    assert hits[0]["tm_ref"] == 0.8
    assert hits[0]["tm_query"] == 0.7
    assert hits[0]["aligned_length"] == 4
    assert hits[0]["query_aln"] == "AAAA"
    assert hits[0]["ref_aln"] == "AAAA"
