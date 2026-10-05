from pathlib import Path

from benchmark.scripts.score_transformed_usalign_batch import (
    candidate_pair_paths,
    combine_transformed_pair,
    cleanup_transformed_halves,
)


def _atom(chain: str, x: float) -> str:
    return (
        f"ATOM      1  CA  ALA {chain}   1       {x:6.3f}   0.000   0.000"
        "  1.00 20.00           C\n"
    )


def test_combine_transformed_pair_renames_overlapping_partner_chains(tmp_path):
    left = tmp_path / "1abcCD_2defAB_3ghiAB_o1_L.pdb"
    right = tmp_path / "1abcCD_2defAB_3ghiAB_o1_R.pdb"
    combined = tmp_path / "combined.pdb"
    left.write_text(_atom("A", 1.0) + _atom("B", 2.0))
    right.write_text(_atom("A", 3.0) + _atom("B", 4.0))

    groups = combine_transformed_pair(left, right, combined)

    assert groups == ("AB", "CD")
    assert {line[21] for line in combined.read_text().splitlines() if line.startswith("ATOM")} == {"A", "B", "C", "D"}


def test_candidate_pair_paths_use_shared_transformation_naming(tmp_path):
    row = {
        "template": "2defCD",
        "query_left": "1abcAB",
        "query_right": "3ghiE",
        "orientation": "o2",
    }
    left, right = candidate_pair_paths(tmp_path, row)
    assert left == Path(tmp_path) / "processed/transformation/2defCD_1abcAB_3ghiE_o2_L.pdb"
    assert right == Path(tmp_path) / "processed/transformation/2defCD_1abcAB_3ghiE_o2_R.pdb"


def test_cleanup_retains_unresolved_transformed_candidates(tmp_path):
    resolved_left = tmp_path / "resolved_L.pdb"
    resolved_right = tmp_path / "resolved_R.pdb"
    failed_left = tmp_path / "failed_L.pdb"
    failed_right = tmp_path / "failed_R.pdb"
    for path in (resolved_left, resolved_right, failed_left, failed_right):
        path.write_text("ATOM\n")

    manifest = cleanup_transformed_halves(
        [
            {
                "status": "generated",
                "score_status": "scored",
                "transformed_left": str(resolved_left),
                "transformed_right": str(resolved_right),
            },
            {
                "status": "generated",
                "score_status": "score_failed",
                "transformed_left": str(failed_left),
                "transformed_right": str(failed_right),
            },
        ],
        tmp_path / "cleanup.json",
    )

    assert manifest["deleted_transformed_half_count"] == 2
    assert manifest["retained_unresolved_candidate_count"] == 1
    assert not resolved_left.exists()
    assert not resolved_right.exists()
    assert failed_left.exists()
    assert failed_right.exists()


def test_cleanup_can_retain_scored_transformed_inputs_for_refinement(tmp_path):
    left = tmp_path / "left.pdb"
    right = tmp_path / "right.pdb"
    left.write_text("ATOM\n")
    right.write_text("ATOM\n")

    manifest = cleanup_transformed_halves(
        [
            {
                "status": "generated",
                "score_status": "scored",
                "transformed_left": str(left),
                "transformed_right": str(right),
            }
        ],
        tmp_path / "cleanup.json",
        preserve_for_downstream=True,
    )

    assert manifest["preserved_for_downstream"] is True
    assert manifest["retention_reason"] == "common_refinement_pending"
    assert manifest["deleted_transformed_half_count"] == 0
    assert left.exists()
    assert right.exists()
