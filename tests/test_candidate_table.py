import json

import pytest

from benchmark.scripts.build_candidate_table import build_rows, build_rows_from_stage_manifest


def test_candidate_table_keeps_missing_alignment_and_both_orientations(tmp_path):
    alignments = tmp_path / "alignment"
    transforms = tmp_path / "transformation"
    alignments.mkdir()
    transforms.mkdir()
    for query, chain, count in (("left", "A", 20), ("right", "B", 18)):
        (alignments / f"{query}_3xyzAB_{chain}.json").write_text(
            json.dumps({"match_count": count, "tm_score": 0.6, "match_dict": {}})
        )
    (transforms / "3xyzAB_left_right_o1_L.pdb").write_text("")
    (transforms / "3xyzAB_left_right_o1_R.pdb").write_text("")

    refinements = tmp_path / "refinements"
    refinements.mkdir()
    (refinements / "3xyzAB_left_right_o1_rosetta.pdb").write_text("")
    rows = build_rows(
        alignments,
        transforms,
        [("left", "right")],
        native_complex_id="native-1",
        refinement_dir=refinements,
    )
    assert [row["orientation"] for row in rows] == ["o1", "o2"]
    assert rows[0]["status"] == "generated"
    assert rows[1]["status"] == "alignment_failed"
    assert rows[1]["error_reason"] == "missing_alignment_json"
    assert rows[0]["native_complex_id"] == "native-1"
    assert rows[0]["model_complex"].endswith("rosetta.pdb")


def test_candidate_table_rejects_nonpositive_template_limit(tmp_path):
    with pytest.raises(ValueError, match="limit must be positive"):
        build_rows(tmp_path, tmp_path, [("left", "right")], limit=0)


def test_candidate_table_builds_durable_row_from_stage_and_gtalign_json(tmp_path):
    batch = tmp_path / "batch_0001"
    processed = batch / "processed"
    alignment = processed / "alignment_gtalign" / "run-1"
    refinement = processed / "rosetta_refinement"
    alignment.mkdir(parents=True)
    refinement.mkdir(parents=True)
    interfaces = batch / "templates" / "interfaces_lists"
    interfaces.mkdir(parents=True)
    (interfaces / "3xyzAB.json").write_text(json.dumps({"A": list(range(30)), "B": list(range(40))}))
    model = refinement / "3xyzAB_leftA_rightB_o2_L_3xyzAB_leftA_rightB_o2_R_rosetta_0001_0001.pdb"
    model.write_text("MODEL\n", encoding="ascii")
    (alignment / "leftA_3xyzAB_B.json").write_text(json.dumps({
        "match_count": 20, "tm_score": 0.7, "tm_score_ref": 0.7,
        "tm_score_query": 0.5, "match_dict": {"x": "y"},
        "raw_output_sha256": "raw-left", "aligner": "GTalign", "status": "success",
    }))
    (alignment / "rightB_3xyzAB_A.json").write_text(json.dumps({
        "match_count": 18, "tm_score": 0.6, "tm_score_ref": 0.6,
        "tm_score_query": 0.4, "match_dict": {"a": "b"},
        "raw_output_sha256": "raw-right", "aligner": "GTalign", "status": "success",
    }))
    stage = tmp_path / "stage.csv"
    stage.write_text(
        "dataset_row_id,status,source_gate_status,source_model_path,source_model_sha256,"
        "template_1,template_2,orientation,receptor,ligand,complex,refinement_backend\n"
        f"rigid:000001,staged_symlink,strict_clean,{model},,3xyzA,3xyzB,2,leftA,rightB,1ABC_A:B,external_rosetta\n"
    )

    rows = build_rows_from_stage_manifest(stage)

    assert len(rows) == 1
    assert rows[0]["dataset_row_id"] == "rigid:000001"
    assert rows[0]["native_complex_id"] == "rigid:000001"
    assert rows[0]["orientation"] == "o2"
    assert rows[0]["chain_left"] == "B"
    assert rows[0]["chain_right"] == "A"
    assert rows[0]["match_count_left"] == 20
    assert rows[0]["tm_score_right"] == 0.6
    assert rows[0]["status"] == "refinement_accepted"
    assert rows[0]["alignment_left_raw_output_sha256"] == "raw-left"
    assert rows[0]["alignment_right_aligner"] == "GTalign"
    assert rows[0]["mapping_count_left"] == 1
    assert rows[0]["alignment_left_sha256"]
    assert rows[0]["match_coverage_left"] == 50.0
    assert rows[0]["match_coverage_right"] == 60.0
    assert rows[0]["template_coverage_status"] == "available"
    assert rows[0]["template_interface_sha256"]


def test_stage_candidate_table_rejects_ambiguous_alignment_runs(tmp_path):
    batch = tmp_path / "batch_0001"
    processed = batch / "processed"
    refinement = processed / "rosetta_refinement"
    refinement.mkdir(parents=True)
    model = refinement / "model.pdb"
    model.write_text("MODEL\n", encoding="ascii")
    payload = json.dumps({
        "match_count": 10, "tm_score": 0.5, "raw_output_sha256": "raw",
        "aligner": "GTalign", "status": "success",
    })
    for run in ("run-1", "run-2"):
        alignment = processed / "alignment_gtalign" / run
        alignment.mkdir(parents=True)
        (alignment / "leftA_3xyzAB_A.json").write_text(payload)
        (alignment / "rightB_3xyzAB_B.json").write_text(payload)
    stage = tmp_path / "stage.csv"
    stage.write_text(
        "dataset_row_id,status,source_gate_status,source_model_path,source_model_sha256,"
        "template_1,template_2,orientation,receptor,ligand,complex,refinement_backend\n"
        f"rigid:000001,staged_symlink,strict_clean,{model},,3xyzA,3xyzB,1,leftA,rightB,1ABC_A:B,external_rosetta\n"
    )

    rows = build_rows_from_stage_manifest(stage)

    assert rows[0]["status"] == "alignment_failed"
    assert "ambiguous_alignment_json" in rows[0]["error_reason"]
