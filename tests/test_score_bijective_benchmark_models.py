import csv
import sys
from pathlib import Path

import pytest

import benchmark.scripts.score_bijective_benchmark_models as scorer
from benchmark.scripts.score_bijective_benchmark_models import validate_mapping_option, validate_score_candidate


def _atom(chain: str) -> str:
    return (
        f"ATOM      1  CA  ALA {chain}   1    "
        "   0.000   0.000   0.000  1.00 20.00           C\n"
    )


def _complete_residue(chain: str) -> str:
    return "".join(
        f"ATOM  {serial:5d} {atom:>4} ALA {chain}   1    "
        "   0.000   0.000   0.000  1.00 20.00           C\n"
        for serial, atom in enumerate(("N", "CA", "C", "O"), 1)
    )


def test_dockq_command_includes_explicit_cpu_count(tmp_path):
    command = scorer.dockq_command(
        Path("python"),
        None,
        tmp_path / "model.pdb",
        tmp_path / "native.pdb",
        "AB:AB",
        tmp_path / "raw.json",
        n_cpu=2,
    )

    assert command[-2:] == ["--n_cpu", "2"]


def test_raw_json_path_includes_model_identity_digest(tmp_path):
    stage = {"dataset_row_id": "rigid:000001", "pair_id": "pair-1"}
    first = tmp_path / "first" / "model.pdb"
    second = tmp_path / "second" / "model.pdb"
    native = tmp_path / "native.pdb"

    first_path = scorer.raw_json_path_for_model(
        tmp_path / "raw", stage, first, native, "AB:AB"
    )
    second_path = scorer.raw_json_path_for_model(
        tmp_path / "raw", stage, second, native, "AB:AB"
    )

    assert first_path != second_path
    assert "rigid" not in first_path.name or first_path.name != second_path.name


def _stage(model: Path, *, status: str = "staged_symlink") -> dict[str, str]:
    return {
        "status": status,
        "staged_model_path": str(model),
        "model_receptor_chains": "A",
        "model_ligand_chains": "B",
        "complex": "1abc_A:B",
    }


def test_score_candidate_rejects_transformation_half_before_scoring(tmp_path):
    half = tmp_path / "target_1abcA_2defB_o1_L_rosetta.pdb"
    half.write_text(_atom("A") + _atom("B"), encoding="ascii")
    native = tmp_path / "1abc.pdb"
    native.write_text(_atom("A") + _atom("B"), encoding="ascii")

    with pytest.raises(ValueError, match="transformation half"):
        validate_score_candidate(_stage(half), native)


def test_score_candidate_rejects_nonstaged_and_invalid_chain_contract(tmp_path):
    model = tmp_path / "model.pdb"
    model.write_text(_atom("A"), encoding="ascii")
    native = tmp_path / "1abc.pdb"
    native.write_text(_atom("A") + _atom("B"), encoding="ascii")

    with pytest.raises(ValueError, match="not stageable"):
        validate_score_candidate(_stage(model, status="transformation_intermediate"), native)
    with pytest.raises(ValueError, match="model chain contract"):
        validate_score_candidate(_stage(model), native)


def test_score_candidate_accepts_complete_staged_model(tmp_path):
    model = tmp_path / "model.pdb"
    native = tmp_path / "1abc.pdb"
    model.write_text(_atom("A") + _atom("B"), encoding="ascii")
    native.write_text(_atom("A") + _atom("B"), encoding="ascii")

    assert validate_score_candidate(_stage(model), native) == model


def test_scored_model_writes_hashed_global_and_interface_records(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native_root = tmp_path / "native"
    native_root.mkdir()
    native = native_root / "1abc.pdb"
    model.write_text(_complete_residue("A") + _complete_residue("B"), encoding="ascii")
    native.write_text(_complete_residue("A") + _complete_residue("B"), encoding="ascii")
    manifest = tmp_path / "stage.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(_stage(model) | {"pair_id": "pair-1", "benchmark_set": "rigid"}))
        writer.writeheader()
        writer.writerow(_stage(model) | {"pair_id": "pair-1", "benchmark_set": "rigid"})

    def fake_dockq(score_python, executable, model_path, native_path, mapping, raw_json, timeout, no_align=False):
        payload = {
            "GlobalDockQ": 0.5,
            "best_result": {"AB": {"DockQ": 0.5, "iRMSD": 1.0, "LRMSD": 2.0, "fnat": 0.4, "F1": 0.5, "clashes": 0}},
        }
        raw_json.write_text(__import__("json").dumps(payload), encoding="utf-8")
        return payload

    monkeypatch.setattr(scorer, "run_dockq", fake_dockq)
    monkeypatch.setattr(scorer, "run_irmsd", lambda *args: 1.25)
    model_output = tmp_path / "scores_models.tsv"
    interface_output = tmp_path / "scores_interfaces.tsv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "score", "--stage-manifest", str(manifest), "--native-root", str(native_root),
            "--score-python", str(tmp_path / "python"), "--irmsd-script", str(tmp_path / "irmsd.py"),
            "--output", str(model_output), "--interfaces-output", str(interface_output),
            "--raw-json-dir", str(tmp_path / "raw"),
        ],
    )

    assert scorer.main() == 0
    with model_output.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["score_status"] == "scored"
    assert len(row["raw_dockq_json_sha256"]) == 64
    assert len(row["source_model_sha256"]) == 64
    with interface_output.open(newline="", encoding="utf-8") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    assert [record["record_type"] for record in records] == ["global", "interface"]
    assert records[1]["requested_cross_interface"] == "True"


def test_failed_model_retains_native_hash_for_adjudication(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native_root = tmp_path / "native"
    native_root.mkdir()
    native = native_root / "1abc.pdb"
    model.write_text(_complete_residue("A") + _complete_residue("B"), encoding="utf-8")
    native.write_text(_complete_residue("A") + _complete_residue("B"), encoding="utf-8")
    manifest = tmp_path / "stage.csv"
    stage = _stage(model) | {"pair_id": "pair-failure", "benchmark_set": "rigid"}
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(stage))
        writer.writeheader()
        writer.writerow(stage)

    def failing_dockq(*args, **kwargs):
        raise RuntimeError("Buffer has wrong number of dimensions (expected 2, got 1)")

    monkeypatch.setattr(scorer, "run_dockq", failing_dockq)
    monkeypatch.setattr(sys, "argv", [
        "score", "--stage-manifest", str(manifest), "--native-root", str(native_root),
        "--score-python", str(tmp_path / "python"), "--irmsd-script", str(tmp_path / "irmsd.py"),
        "--output", str(tmp_path / "scores.tsv"), "--raw-json-dir", str(tmp_path / "raw"),
    ])

    assert scorer.main() == 0
    with (tmp_path / "scores.tsv").open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["score_status"] == "score_failed"
    assert row["native_pdb_sha256"] == scorer.sha256_file(native)
    assert row["source_model_sha256"] == scorer.sha256_file(model)


def test_empty_interface_crash_recovers_requested_cross_score_without_global(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native_root = tmp_path / "native"
    native_root.mkdir()
    native = native_root / "1abc.pdb"
    model.write_text(_complete_residue("A") + _complete_residue("B"), encoding="ascii")
    native.write_text(_complete_residue("A") + _complete_residue("B"), encoding="ascii")
    stage = _stage(model) | {
        "pair_id": "pair-cross-only",
        "benchmark_set": "rigid",
        "dataset_row_id": "rigid:000001",
        "source_gate_status": "strict_clean",
        "source_model_path": str(model),
    }
    manifest = tmp_path / "stage.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(stage))
        writer.writeheader()
        writer.writerow(stage)

    def fake_dockq(_python, _executable, _model, _native, _mapping, raw_json, _timeout, _no_align=False):
        if ".pair-" not in raw_json.name:
            raise RuntimeError("Buffer has wrong number of dimensions (expected 2, got 1)")
        payload = {
            "GlobalDockQ": 0.4,
            "best_result": {
                "AB": {
                    "DockQ": 0.4, "iRMSD": 1.0, "LRMSD": 2.0,
                    "fnat": 0.3, "F1": 0.4, "clashes": 0,
                }
            },
        }
        raw_json.write_text(__import__("json").dumps(payload), encoding="utf-8")
        return payload

    monkeypatch.setattr(scorer, "run_dockq", fake_dockq)
    monkeypatch.setattr(scorer, "run_irmsd", lambda *args: 1.25)
    output = tmp_path / "scores.tsv"
    interfaces = tmp_path / "interfaces.tsv"
    monkeypatch.setattr(sys, "argv", [
        "score", "--stage-manifest", str(manifest), "--native-root", str(native_root),
        "--score-python", str(tmp_path / "python"), "--irmsd-script", str(tmp_path / "irmsd.py"),
        "--output", str(output), "--interfaces-output", str(interfaces),
        "--raw-json-dir", str(tmp_path / "raw"),
    ])

    assert scorer.main() == 0
    row = next(csv.DictReader(output.open(), delimiter="\t"))
    assert row["score_status"] == "scored_cross_only"
    assert row["score_scope"] == "requested_cross_interfaces_only"
    assert row["dockq_global"] == ""
    assert row["dockq_cross_mean"] == "0.4"
    raw_payload = __import__("json").loads(Path(row["raw_dockq_json"]).read_text())
    assert raw_payload["GlobalDockQ"] is None
    assert raw_payload["evaluation_mode"] == "pairwise_cross_fallback"
    records = list(csv.DictReader(interfaces.open(), delimiter="\t"))
    assert [record["record_type"] for record in records] == ["interface"]
    assert records[0]["requested_cross_interface"] == "True"


def test_no_align_requires_strict_pdb_mapping(tmp_path):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text(_complete_residue("A") + _complete_residue("B"), encoding="ascii")
    native.write_text(_complete_residue("A") + _complete_residue("B"), encoding="ascii")
    assert validate_mapping_option(model, native, "AB:AB", True) == "validated_no_align"
    with pytest.raises(ValueError, match="unsafe --no-align"):
        validate_mapping_option(model, native, "AB:AC", True)


def test_pairwise_cross_fallback_preserves_dockq_when_auxiliary_irmsd_fails(
    tmp_path, monkeypatch
):
    model = tmp_path / "model.pdb"
    native_root = tmp_path / "native"
    native_root.mkdir()
    native = native_root / "1abc.pdb"
    pdb = _complete_residue("A") + _complete_residue("B") + _complete_residue("C")
    model.write_text(pdb, encoding="ascii")
    native.write_text(pdb, encoding="ascii")
    stage = {
        **_stage(model),
        "model_receptor_chains": "AB",
        "model_ligand_chains": "C",
        "complex": "1abc_AB:C",
        "pair_id": "pair-fallback",
        "benchmark_set": "difficult",
    }
    manifest = tmp_path / "stage.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(stage))
        writer.writeheader()
        writer.writerow(stage)

    def fake_dockq(score_python, executable, model_path, native_path, mapping, raw_json, timeout, no_align=False):
        if ".pair-" not in raw_json.name:
            raise RuntimeError("Buffer has wrong number of dimensions (expected 2, got 1)")
        interface = raw_json.stem.rsplit("pair-", 1)[1]
        if interface == "BC":
            raise RuntimeError("Could not find interfaces in the native model")
        payload = {
            "GlobalDockQ": 0.4,
            "best_result": {interface: {
                "DockQ": 0.4, "iRMSD": 2.0, "LRMSD": 3.0,
                "fnat": 0.3, "F1": 0.4, "clashes": 0,
            }},
        }
        raw_json.write_text(__import__("json").dumps(payload), encoding="utf-8")
        return payload

    monkeypatch.setattr(scorer, "run_dockq", fake_dockq)

    def failing_irmsd(*args):
        raise RuntimeError("legacy multichain failure")

    monkeypatch.setattr(scorer, "run_irmsd", failing_irmsd)
    models = tmp_path / "scores.tsv"
    interfaces = tmp_path / "interfaces.tsv"
    monkeypatch.setattr(sys, "argv", [
        "score", "--stage-manifest", str(manifest), "--native-root", str(native_root),
        "--score-python", str(tmp_path / "python"), "--irmsd-script", str(tmp_path / "irmsd.py"),
        "--output", str(models), "--interfaces-output", str(interfaces),
        "--raw-json-dir", str(tmp_path / "raw"),
    ])

    assert scorer.main() == 0
    with models.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["score_status"] == "scored_cross_only"
    assert row["score_scope"] == "requested_cross_interfaces_only"
    assert row["dockq_global"] == ""
    assert row["dockq_cross_mean"] == "0.4"
    assert row["dockq_cross_requested_count"] == "2"
    assert row["dockq_cross_scoreable_count"] == "1"
    assert row["dockq_cross_unscoreable_count"] == "1"
    assert row["irmsd_status"] == "failed_auxiliary"
    assert "legacy multichain failure" in row["irmsd_error"]
    components = __import__("json").loads(row["dockq_component_runs"])
    assert [component["status"] for component in components] == ["scored", "no_native_interface"]
