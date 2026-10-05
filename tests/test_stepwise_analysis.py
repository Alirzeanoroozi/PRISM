import json

from src.stepwise_analysis import (
    alignment_inventory,
    asset_provenance_manifest,
    clash_diagnostics,
    gate_ledger,
    ranking_comparison,
    usalign_preflight,
    validate_usalign_contract,
)
from src.transformation_config import TransformationThresholds


def _ca_pdb(chain, x):
    return (
        f"ATOM      1  CA  ALA {chain}   1    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C  \nEND\n"
    )


def test_alignment_and_gate_ledgers_preserve_missing_gate_evidence(tmp_path):
    root = tmp_path / "run"
    (root / "templates/interfaces").mkdir(parents=True)
    (root / "processed/alignment_tmalign/run1").mkdir(parents=True)
    (root / "templates/interfaces/1abcAB_A_int.pdb").write_text(_ca_pdb("A", 0.0))
    (root / "templates/interfaces/1abcAB_B_int.pdb").write_text(_ca_pdb("B", 0.0))
    (root / "inputs.csv").write_text("Receptor,Ligand\nq1,q2\n")
    payload = {
        "aligner": "TMalign",
        "match_count": 1,
        "tm_score": 0.8,
        "match_dict": {"A.A.1": "Q.A.1"},
        "return_code": 0,
        "raw_output_sha256": "raw-hash",
    }
    for query, chain in (("q1", "A"), ("q2", "B")):
        (root / f"processed/alignment_tmalign/run1/{query}_1abcAB_{chain}.json").write_text(
            json.dumps(payload)
        )
    (root / "processed/transformation").mkdir(parents=True)
    (root / "processed/transformation/1abcAB_q1_q2_o1_L.pdb").write_text(_ca_pdb("A", 0.0))
    (root / "processed/transformation/1abcAB_q1_q2_o1_R.pdb").write_text(_ca_pdb("B", 2.0))

    rows = alignment_inventory(
        root, root / "inputs.csv", ["1abcAB"], arm="o1", orientations=("o1",)
    )
    assert rows[0]["partner_available"] is True
    assert rows[0]["raw_output_hash_complete"] is True
    assert rows[0]["coverage_left"] == 100.0

    thresholds = TransformationThresholds(
        minimum_residue_match_count=1,
        minimum_residue_match_percentage=50.0,
        clashing_distance=3.0,
        max_clashing_count=5,
    ).as_dict()
    gates = gate_ledger(
        root, rows, thresholds=thresholds, compute_clash_diagnostics=True
    )
    assert gates[0]["gate_availability"] is True
    assert gates[0]["gate_score_threshold"] is True
    assert gates[0]["gate_transform_materialization"] is True
    assert gates[0]["gate_clash_filter"] is True
    assert gates[0]["gate_hotspots"] is None
    assert gates[0]["cumulative_pass"] == "unknown"


def test_incomplete_alignment_contract_cannot_produce_gate_pass(tmp_path):
    root = tmp_path / "run"
    (root / "templates/interfaces").mkdir(parents=True)
    (root / "templates/interfaces_lists").mkdir(parents=True)
    (root / "processed/alignment_tmalign/run1").mkdir(parents=True)
    (root / "templates/interfaces/1abcAB_A_int.pdb").write_text(_ca_pdb("A", 0.0))
    (root / "templates/interfaces/1abcAB_B_int.pdb").write_text(_ca_pdb("B", 0.0))
    (root / "templates/interfaces_lists/1abcAB.json").write_text('{"A": [1], "B": [1]}')
    (root / "inputs.csv").write_text("Receptor,Ligand\nq1,q2\n")
    payload = {
        "aligner": "TMalign",
        "match_count": 20,
        "tm_score": 0.8,
        "match_dict": {"A.A.1": "Q.A.1"},
        "return_code": 1,
        "raw_output_sha256": "raw-hash",
    }
    for query, chain in (("q1", "A"), ("q2", "B")):
        (root / f"processed/alignment_tmalign/run1/{query}_1abcAB_{chain}.json").write_text(
            json.dumps(payload)
        )
    (root / "processed/transformation").mkdir(parents=True)
    (root / "processed/transformation/1abcAB_q1_q2_o1_L.pdb").write_text(_ca_pdb("A", 0.0))
    (root / "processed/transformation/1abcAB_q1_q2_o1_R.pdb").write_text(_ca_pdb("B", 5.0))
    rows = alignment_inventory(
        root, root / "inputs.csv", ["1abcAB"], arm="o1", orientations=("o1",)
    )
    assert rows[0]["alignment_evidence_complete"] is False
    gates = gate_ledger(
        root, rows, thresholds=TransformationThresholds().as_dict(),
        compute_clash_diagnostics=False,
    )
    assert gates[0]["gate_score_threshold"] is None
    assert gates[0]["gate_match_count"] is None
    assert gates[0]["cumulative_pass"] == "unknown"
    assert "return_code_nonzero" in gates[0]["gate_score_threshold_status"]


def test_common_gate_replay_uses_shared_match_coverage_contract(tmp_path):
    row = {
        "partner_available": True,
        "aligner_left": "TMalign",
        "aligner_right": "MultiProt",
        "left_payload": {"match_count": 15, "tm_score": 0.0},
        "right_payload": {"match_count": 15, "tm_score": 0.0},
        "coverage_left": 50.0,
        "coverage_right": 50.0,
        "template_size_left": 40,
        "template_size_right": 40,
        "alignment_evidence_complete": True,
    }
    thresholds = TransformationThresholds(
        alignment_gate_mode="common_match_coverage",
        minimum_residue_match_count=15,
        minimum_residue_match_percentage=50.0,
    ).as_dict()
    gates = gate_ledger(
        tmp_path, [row], thresholds=thresholds,
        compute_clash_diagnostics=False,
    )
    assert gates[0]["gate_score_threshold"] is True
    assert gates[0]["gate_score_threshold_status"] == "not_applied_common_match_coverage"
    assert gates[0]["gate_match_count"] is True
    assert gates[0]["gate_coverage"] is True


def test_stepwise_helpers_report_clash_grid_and_contract_boundaries(tmp_path):
    left, right = tmp_path / "left.pdb", tmp_path / "right.pdb"
    left.write_text(_ca_pdb("A", 0.0))
    right.write_text(_ca_pdb("B", 2.0))
    details = clash_diagnostics(left, right, distance_grid=(2.5, 3.0), event_grid=(0, 1, 5))
    assert details["clash_count"] == 1
    assert details["minimum_cross_partner_distance"] == 2.0
    assert len(details["grid"]) == 6

    assert usalign_preflight(str(tmp_path / "missing-USalign"))["available"] is False
    contract = validate_usalign_contract({"aligner": "USalign"})
    assert contract["valid_prism_contract"] is False
    ranked = ranking_comparison(
        {"none": ["a", "b"], "baseline": ["b"], "prodigy": ["a"]},
        native_labels={"a": True, "b": False},
        native_label_source="native_labels.tsv",
        native_label_source_sha256="labels-hash",
    )
    assert all(row["affinity_used_as_label"] is False for row in ranked)
    assert all(row["native_label_provenance_complete"] is True for row in ranked)
    assert ranked[0]["top_k_native_like_recovery"] == 1.0
    assert ranked[1]["top_k_native_like_recovery"] is None


def test_usalign_contract_requires_both_explicit_normalizations():
    contract = validate_usalign_contract(
        {
            "match_dict": {},
            "rotation_mat": [[1.0, 0.0, 0.0]] * 3,
            "translation": [0.0, 0.0, 0.0],
            "tm_score": 0.8,
            "aligner": "USalign",
            "return_code": 0,
            "raw_output_sha256": "raw-hash",
        }
    )
    assert contract["valid_prism_contract"] is False
    assert "tm_score_query" in contract["missing_fields"]
    assert "tm_score_ref" in contract["missing_fields"]
    assert "tm_score_contract" in contract["missing_fields"]


def test_asset_manifest_marks_current_only_run(tmp_path):
    root = tmp_path / "run"
    source = tmp_path / "source"
    (source / "templates/interfaces").mkdir(parents=True)
    (source / "templates/interfaces_lists").mkdir(parents=True)
    (source / "templates/interfaces/1abcAB_A_int.pdb").write_text(_ca_pdb("A", 0.0))
    (source / "templates/interfaces/1abcAB_B_int.pdb").write_text(_ca_pdb("B", 0.0))
    (source / "templates/interfaces_lists/1abcAB.json").write_text('{"A": [1], "B": [1]}')
    root.mkdir()
    (root / "templates").symlink_to(source / "templates", target_is_directory=True)
    inputs = tmp_path / "inputs.csv"
    inputs.write_text("Receptor,Ligand\nq1,q2\n")
    manifest = asset_provenance_manifest(
        root, inputs, ["1abcAB"], source_root=source,
        surface_backend="naccess", filter_mode="geometry_only_experimental",
    )
    assert manifest["asset_mix_status"] == "current_only"
    assert manifest["input_csv"][0]["sha256"]
    assert manifest["template_interface_list_assets"][0]["exists"] is True


def test_asset_manifest_detects_hash_matched_mixed_assets_even_when_renamed(tmp_path):
    root = tmp_path / "run"
    source = tmp_path / "source"
    current = source / "templates/interfaces"
    current_lists = source / "templates/interfaces_lists"
    legacy = source / "working_version/Multiprot-new/prism-fiberdock-cli/template/interfaces"
    current.mkdir(parents=True)
    current_lists.mkdir(parents=True)
    legacy.mkdir(parents=True)
    (current / "current.pdb").write_text(_ca_pdb("A", 0.0))
    (legacy / "legacy.pdb").write_text(_ca_pdb("B", 1.0))
    (current_lists / "1abcAB.json").write_text('{"A": [1], "B": [1]}')
    (root / "templates/interfaces").mkdir(parents=True)
    (root / "templates/interfaces_lists").mkdir(parents=True)
    (root / "templates/interfaces/1abcAB_A_int.pdb").write_text((current / "current.pdb").read_text())
    (root / "templates/interfaces/1abcAB_B_int.pdb").write_text((legacy / "legacy.pdb").read_text())
    (root / "templates/interfaces_lists/1abcAB.json").write_text((current_lists / "1abcAB.json").read_text())
    inputs = tmp_path / "inputs.csv"
    inputs.write_text("Receptor,Ligand\nq1,q2\n")

    manifest = asset_provenance_manifest(
        root, inputs, ["1abcAB"], source_root=source,
        surface_backend="naccess", filter_mode="geometry_only_experimental",
    )
    assert manifest["asset_mix_status"] == "mixed"
    assert set(manifest["consumed_asset_origins"]) == {"current", "legacy"}
