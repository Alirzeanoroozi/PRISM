from pathlib import Path
import json

from benchmark.scripts.build_matched_benchmark_manifest import (
    PILOT_ROW_IDS,
    build_manifests,
    select_rows,
)


ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "tmp/agent/20260713-investigation-implementation/source-gate-aggregate-final/source_gate_policy.json"
PRIOR = ROOT / "tmp/agent/20260712-multichain-full-comparison/scored_models.csv"


def test_pilot_selection_is_frozen_and_balanced():
    rows, policy = select_rows(ROOT, POLICY, "pilot", PRIOR)
    assert [row["dataset_row_id"] for row in rows] == list(PILOT_ROW_IDS)
    assert {row["dataset_row_id"] for row in rows}.isdisjoint(
        policy["decision"]["excluded_dataset_row_ids"]
    )
    for difficulty in ("rigid", "medium", "difficult"):
        subset = [row for row in rows if row["difficulty"] == difficulty]
        assert len(subset) == 4
        assert sum(row["chain_context"] == "single_chain" for row in subset) == 2
        assert sum(row["chain_context"] == "multichain" for row in subset) == 2


def test_strict_clean_selection_excludes_exactly_frozen_rows():
    rows, policy = select_rows(ROOT, POLICY, "strict-clean", PRIOR)
    excluded = set(policy["decision"]["excluded_dataset_row_ids"])
    assert len(rows) == 240
    assert not excluded.intersection(row["dataset_row_id"] for row in rows)


def test_build_manifests_freezes_template_and_task_contract(tmp_path):
    paths = build_manifests(ROOT, tmp_path, cohort="pilot", source_policy=POLICY, prior_models=PRIOR)
    assert set(paths) == {"cohort_manifest", "pipeline_inputs", "template_manifest", "arm_manifest", "task_manifest", "analysis_plan"}
    assert sum(1 for _ in paths["cohort_manifest"].open()) == 13
    assert sum(1 for _ in paths["pipeline_inputs"].open()) == 13
    template_list = [line.strip() for line in (ROOT / "new_template/template/final_list.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
    assert sum(1 for _ in paths["template_manifest"].open()) == len(template_list) + 1  # header + templates
    assert sum(1 for _ in paths["arm_manifest"].open()) == 5
    assert sum(1 for _ in paths["task_manifest"].open()) == 145
    plan = __import__("json").loads(paths["analysis_plan"].read_text())
    assert plan["tm_score_threshold"] == 0.4
    assert plan["tm_score_sensitivity_thresholds"] == [0.30, 0.35]


def test_tmalign_parser_preserves_partial_mapping_for_duplicate_residues(tmp_path):
    from src.alignment import parse_tmalign

    protein = tmp_path / "protein.pdb"
    interface = tmp_path / "interface.pdb"
    # One residue record is intentionally absent from the coordinate list while
    # the alignment text remains well formed; this used to raise IndexError.
    atom = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
    protein.write_text(atom)
    interface.write_text(atom.replace(" A   1", " B   1"))
    matrix = tmp_path / "matrix.out"
    matrix.write_text("0 0 1 0 0\n1 0 0 1 0\n2 0 0 0 1\n")
    tm = tmp_path / "out.tm"
    tm.write_text(
        "Aligned length= 1, RMSD= 0.0, Seq_ID=n_identical/n_aligned= 1.0\n"
        "TM-score= 1.0 (if normalized by length of Chain_1)\n"
        '(":" denotes residue pairs of d <  5.0 Angstrom)\n'
        "A\n:\nA\n"
    )
    output = tmp_path / "alignment"
    parse_tmalign(
        str(protein), str(interface), "q", "tpl", "B", str(matrix), str(tm),
        str(output), raw_output_sha256="raw-hash", return_code=0,
    )
    result = json.loads((output / "q_tpl_B.json").read_text())
    assert result["status"] in {"success", "mapping_truncated"}
    assert result["match_count"] == 1
    assert result["raw_output_sha256"] == "raw-hash"
    assert result["return_code"] == 0
