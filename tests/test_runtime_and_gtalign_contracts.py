import json
from pathlib import Path

from benchmark.scripts.validate_runtime_manifest import validate
from benchmark.scripts.stage_legacy_tool_environment import classify_refinement_capability
from src.alignment import _valid_tmalign_outputs
from src.alignment_gtalign import write_alignment_json


def test_runtime_manifest_matches_checked_out_core_binaries():
    assert validate(Path(__file__).resolve().parents[1]) == []


def test_gtalign_alignment_record_retains_both_normalized_scores(tmp_path):
    write_alignment_json(
        tmp_path,
        "1abcA",
        "1xyzAB",
        "A",
        18,
        [1, 2, 3],
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        {"A.A.1": "A.A.1"},
        0.72,
        tm_score_ref=0.61,
        tm_score_query=0.72,
        raw_output_sha256="raw-hash",
    )
    record = json.loads((tmp_path / "1abcA_1xyzAB_A.json").read_text())
    assert record["tm_score"] == 0.72
    assert record["tm_score_ref"] == 0.61
    assert record["tm_score_query"] == 0.72
    assert record["raw_output_sha256"] == "raw-hash"
    assert record["aligner"] == "GTalign"
    assert record["status"] == "success"


def test_gtalign_output_writer_does_not_overwrite_existing_run(tmp_path):
    output = tmp_path / "run"
    output.mkdir()
    (output / "stale.json").write_text("{}")
    # The production stage performs this guard before writing.  This test
    # documents the contract without invoking the external GTalign binary.
    from src.alignment_gtalign import align_gtalign

    try:
        align_gtalign([], [], output_dir=output)
    except RuntimeError as exc:
        assert "non-empty GTalign output directory" in str(exc)
    else:
        raise AssertionError("stale GTalign output directory was accepted")


def test_fiberdock_32_bit_helpers_remain_an_explicit_capability_blocker():
    native = {
        name: {
            "exists": True,
            "file": {"stdout": "ELF 32-bit LSB executable"},
            "ldd": {"dependency_status": "available"},
        }
        for name in (
            "fiberdock/nma",
            "fiberdock/reduce.2.21.030604",
            "fiberdock/reduce.3.23.130521",
            "fiberdock/addHydrogens.pl",
        )
    }
    blockers, architecture = classify_refinement_capability(native)
    assert architecture["fiberdock/nma"] == "ELF 32-bit"
    assert "fiberdock/nma:32-bit-helper" in blockers
    assert "fiberdock/reduce.2.21.030604:32-bit-helper" in blockers
    assert "fiberdock/reduce.3.23.130521:32-bit-helper" in blockers
    assert "fiberdock/full_refinement:not_validated_end_to_end" in blockers


def test_fiberdock_missing_loader_remains_a_capability_blocker():
    native = {
        name: {
            "exists": True,
            "file": {"stdout": "ELF 32-bit LSB executable"},
            "ldd": {"dependency_status": "available" if name.endswith("nma") else "missing"},
        }
        for name in (
            "fiberdock/nma",
            "fiberdock/reduce.2.21.030604",
            "fiberdock/reduce.3.23.130521",
            "fiberdock/addHydrogens.pl",
        )
    }
    blockers, _ = classify_refinement_capability(native)
    assert "fiberdock/reduce.2.21.030604:missing-library" in blockers
    assert "fiberdock/full_refinement:not_validated_end_to_end" in blockers


def test_tmalign_output_contract_rejects_empty_or_incomplete_outputs(tmp_path):
    matrix = tmp_path / "matrix.out"
    output = tmp_path / "out.tm"
    matrix.write_text("0 0 1 0 0\n1 0 0 1 0\n")
    output.write_text("Aligned length=0\n")
    assert not _valid_tmalign_outputs(matrix, output)
    matrix.write_text("0 0 1 0 0\n1 0 0 1 0\n2 0 0 0 1\n")
    output.write_text('Aligned length=3\nTM-score= 0.5\n(":" denotes aligned residues)\nAAA\n:::\nAAA\n')
    assert _valid_tmalign_outputs(matrix, output)
