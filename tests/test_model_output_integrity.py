import json
import sys
from pathlib import Path

import benchmark.scripts.score_single_prism_pair as scorer
from benchmark.scripts.score_single_prism_pair import score_one
from benchmark.scripts.standardized_evaluator import validate_raw_pdb_chain_contract


def _atom(serial, chain, residue):
    return (
        f"ATOM  {serial:5d}  CA  ALA {chain}{residue:4d}    "
        "0.000   0.000   0.000  1.00 20.00           C\n"
    )


def _parseable_atom(serial, chain, residue):
    return (
        f"ATOM  {serial:5d}  CA  ALA {chain}{residue:4d}    "
        f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}           C  \n"
    )


def test_raw_chain_contract_accepts_two_distinct_chains(tmp_path):
    model = tmp_path / "valid.pdb"
    model.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    result = validate_raw_pdb_chain_contract(model, "A", "B")
    assert result.valid


def test_raw_chain_contract_rejects_legacy_chain_collision_and_reset(tmp_path):
    model = tmp_path / "legacy.pdb"
    model.write_text(
        _atom(1, "B", 1)
        + _atom(2, "B", 2)
        + _atom(3, "B", 1)
        + "END\n"
    )
    result = validate_raw_pdb_chain_contract(model, "B", "B")
    assert not result.valid
    assert any("overlap" in error for error in result.errors)
    assert any("residue numbering resets" in error for error in result.errors)


def test_score_one_fails_closed_before_external_scorers(tmp_path):
    model = tmp_path / "legacy.pdb"
    model.write_text(_atom(1, "B", 1) + _atom(2, "B", 1) + "END\n")
    result = score_one(
        model,
        Path("native.pdb"),
        Path("python"),
        timeout_sec=1,
        dockq_no_align=False,
        model_receptor="B",
        model_ligand="B",
        native_receptor="A",
        native_ligand="B",
    )
    assert result["irmsd"] is None
    assert result["dockq"] is None
    assert "output contract rejected" in result["error_dockq"]


def test_score_one_rejects_invalid_irmsd_after_successful_subprocess(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    native.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")

    def fake_run_cmd(cmd, timeout_sec, env=None):
        if str(scorer.REPO_ROOT / "benchmark/scripts/irmsd.py") in cmd:
            return {"ok": True, "stdout": "not-a-number", "stderr": "", "returncode": 0}
        json_path = Path(cmd[cmd.index("--json") + 1])
        json_path.write_text(
            '{"GlobalDockQ": 0.5, "best_result": [{"interface": "A:B", '
            '"DockQ": 0.5, "iRMSD": 1.0, "LRMSD": 2.0, "fnat": 0.4, '
            '"F1": 0.5, "clashes": 0}]}\n',
            encoding="utf-8",
        )
        return {"ok": True, "stdout": "", "stderr": "", "returncode": 0}

    monkeypatch.setattr(scorer, "run_cmd", fake_run_cmd)
    result = score_one(model, native, Path("python"), timeout_sec=1, dockq_no_align=False)

    assert result["irmsd"] is None
    assert "iRMSD" in result["error_irmsd"]
    assert result["dockq"] == 0.5


def test_score_one_does_not_promote_first_interface_to_model_fields(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    native.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")

    def fake_run_cmd(cmd, timeout_sec, env=None):
        if str(scorer.REPO_ROOT / "benchmark/scripts/irmsd.py") in cmd:
            return {"ok": True, "stdout": "1.25", "stderr": "", "returncode": 0}
        json_path = Path(cmd[cmd.index("--json") + 1])
        json_path.write_text(
            '{"GlobalDockQ": 0.5, "best_result": ['
            '{"interface": "A:B", "DockQ": 0.5, "iRMSD": 1.0, "LRMSD": 2.0, "fnat": 0.4, "F1": 0.5, "clashes": 0}, '
            '{"interface": "C:D", "DockQ": 0.3, "iRMSD": 3.0, "LRMSD": 4.0, "fnat": 0.2, "F1": 0.3, "clashes": 1}]}'
            '\n',
            encoding="utf-8",
        )
        return {"ok": True, "stdout": "", "stderr": "", "returncode": 0}

    monkeypatch.setattr(scorer, "run_cmd", fake_run_cmd)
    result = score_one(model, native, Path("python"), timeout_sec=1, dockq_no_align=False)

    assert result["irmsd"] == 1.25
    assert result["dockq"] == 0.5
    assert result["dockq_interface_count"] == 2
    assert result["dockq_irmsd"] is None
    assert result["dockq_lrmsd"] is None
    assert result["dockq_fnat"] is None
    assert result["dockq_f1"] is None
    assert result["dockq_clashes"] is None


def test_score_one_does_not_reuse_stale_dockq_json(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    native.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    json_dir = tmp_path / "dockq"
    json_dir.mkdir()
    stale = json_dir / "model.dockq.json"
    stale.write_text('{"GlobalDockQ": 0.1, "best_result": []}\n')

    def fake_run_cmd(cmd, timeout_sec, env=None):
        if str(scorer.REPO_ROOT / "benchmark/scripts/irmsd.py") in cmd:
            return {"ok": True, "stdout": "1.25", "stderr": "", "returncode": 0}
        json_path = Path(cmd[cmd.index("--json") + 1])
        json_path.write_text(
            '{"GlobalDockQ": 0.7, "best_result": [{"interface": "A:B", '
            '"DockQ": 0.7, "iRMSD": 1.0, "LRMSD": 2.0, "fnat": 0.4, '
            '"F1": 0.5, "clashes": 0}]}\n'
        )
        return {"ok": True, "stdout": "", "stderr": "", "returncode": 0}

    monkeypatch.setattr(scorer, "run_cmd", fake_run_cmd)
    result = score_one(
        model,
        native,
        Path("python"),
        timeout_sec=1,
        dockq_no_align=False,
        dockq_json_dir=json_dir,
    )

    assert result["dockq"] == 0.7
    assert '"GlobalDockQ": 0.1' in stale.read_text()
    assert Path(result["dockq_raw_json_path"]) != stale
    assert Path(result["dockq_raw_json_path"]).is_file()


def test_score_one_rejects_unsafe_no_align_mapping_before_dockq(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text(
        _parseable_atom(1, "A", 1) + _parseable_atom(2, "B", 1) + "END\n"
    )
    native.write_text(
        _parseable_atom(1, "A", 1) + _parseable_atom(2, "C", 2) + "END\n"
    )
    commands = []

    def fake_run_cmd(cmd, timeout_sec, env=None):
        commands.append(cmd)
        if "-m" in cmd and "DockQ" in cmd:
            json_path = Path(cmd[cmd.index("--json") + 1])
            json_path.write_text(
                '{"GlobalDockQ": 0.5, "best_result": []}\n',
                encoding="utf-8",
            )
        return {"ok": True, "stdout": "1.0", "stderr": "", "returncode": 0}

    monkeypatch.setattr(scorer, "run_cmd", fake_run_cmd)
    result = score_one(
        model,
        native,
        Path("python"),
        timeout_sec=1,
        dockq_no_align=True,
        model_receptor="A",
        model_ligand="B",
        native_receptor="A",
        native_ligand="C",
    )

    assert "unsafe no-align" in result["error_dockq"]
    assert not any("DockQ" in command for command in commands)


def test_score_one_does_not_use_multichain_best_sum_without_global_score(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    native.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")

    def fake_run_cmd(cmd, timeout_sec, env=None):
        if str(scorer.REPO_ROOT / "benchmark/scripts/irmsd.py") in cmd:
            return {"ok": True, "stdout": "1.0", "stderr": "", "returncode": 0}
        json_path = Path(cmd[cmd.index("--json") + 1])
        json_path.write_text(
            '{"best_dockq": 1.4, "best_result": ['
            '{"interface": "A:B", "DockQ": 0.8}, '
            '{"interface": "A:C", "DockQ": 0.6}]}\n',
            encoding="utf-8",
        )
        return {"ok": True, "stdout": "", "stderr": "", "returncode": 0}

    monkeypatch.setattr(scorer, "run_cmd", fake_run_cmd)
    result = score_one(model, native, Path("python"), timeout_sec=1, dockq_no_align=False)

    assert result["dockq"] is None
    assert result["dockq_interface_count"] == 2


def test_score_one_rejects_short_output_when_dockq_json_is_missing(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    native.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")

    def fake_run_cmd(cmd, timeout_sec, env=None):
        if str(scorer.REPO_ROOT / "benchmark/scripts/irmsd.py") in cmd:
            return {"ok": True, "stdout": "1.0", "stderr": "", "returncode": 0}
        return {
            "ok": True,
            "stdout": "DockQ 0.9",
            "stderr": "",
            "returncode": 0,
        }

    monkeypatch.setattr(scorer, "run_cmd", fake_run_cmd)
    result = score_one(model, native, Path("python"), timeout_sec=1, dockq_no_align=False)

    assert result["dockq"] is None
    assert "JSON" in result["error_dockq"]


def test_score_json_path_distinguishes_same_stem_from_different_sources(tmp_path):
    first = tmp_path / "first" / "model.pdb"
    second = tmp_path / "second" / "model.pdb"
    native = tmp_path / "native.pdb"
    assert scorer.dockq_json_path(tmp_path / "json", first, native, "AB:AB") != scorer.dockq_json_path(
        tmp_path / "json", second, native, "AB:AB"
    )


def test_score_one_marks_valid_json_without_global_as_unscored(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    native.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")

    def fake_run_cmd(cmd, timeout_sec, env=None):
        if str(scorer.REPO_ROOT / "benchmark/scripts/irmsd.py") in cmd:
            return {"ok": True, "stdout": "1.0", "stderr": "", "returncode": 0}
        json_path = Path(cmd[cmd.index("--json") + 1])
        json_path.write_text(
            '{"best_result": [{"DockQ": 0.5, "iRMSD": 1.0}, '
            '{"DockQ": 0.4, "iRMSD": 2.0}]}\n',
            encoding="utf-8",
        )
        return {"ok": True, "stdout": "", "stderr": "", "returncode": 0}

    monkeypatch.setattr(scorer, "run_cmd", fake_run_cmd)
    result = score_one(model, native, Path("python"), timeout_sec=1, dockq_no_align=False)

    assert result["dockq"] is None
    assert result["dockq_json_status"] == "valid_unscored"
    assert "GlobalDockQ" in result["error_dockq"]


def test_score_one_records_global_score_and_dockq_provenance(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    json_dir = tmp_path / "dockq-json"
    model.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    native.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")

    def fake_run_cmd(cmd, timeout_sec, env=None):
        if str(scorer.REPO_ROOT / "benchmark/scripts/irmsd.py") in cmd:
            return {"ok": True, "stdout": "1.0", "stderr": "", "returncode": 0}
        json_path = Path(cmd[cmd.index("--json") + 1])
        json_path.write_text(
            '{"GlobalDockQ": 0.5, "best_dockq": 0.5, "best_result": ['
            '{"interface": "A:B", "DockQ": 0.5, "iRMSD": 1.0}]}'
            "\n",
            encoding="utf-8",
        )
        return {"ok": True, "stdout": "", "stderr": "", "returncode": 0}

    monkeypatch.setattr(scorer, "run_cmd", fake_run_cmd)
    result = score_one(
        model,
        native,
        Path("python"),
        timeout_sec=1,
        dockq_no_align=False,
        dockq_json_dir=json_dir,
    )

    assert result["dockq"] == 0.5
    assert result["dockq_global"] == 0.5
    assert Path(result["dockq_raw_json_path"]).is_file()
    assert len(result["dockq_raw_json_sha256"]) == 64
    dockq_argv = json.loads(result["dockq_argv"])
    assert "--short" in dockq_argv
    assert dockq_argv[dockq_argv.index("--n_cpu") + 1] == "1"


def test_score_cli_passes_dockq_json_directory(tmp_path, monkeypatch, capsys):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    json_dir = tmp_path / "dockq-json"
    model.write_text(_atom(1, "A", 1) + _atom(2, "B", 1) + "END\n")
    native.write_text(model.read_text())
    captured = {}

    def fake_score_one(**kwargs):
        captured.update(kwargs)
        return {"model_pdb": str(model), "native_pdb": str(native), "dockq": 0.5}

    monkeypatch.setattr(scorer, "score_one", fake_score_one)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "score_single_prism_pair.py",
            str(model),
            str(native),
            "--dockq-json-dir",
            str(json_dir),
            "--json",
        ],
    )

    assert scorer.main() == 0
    assert Path(captured["dockq_json_dir"]) == json_dir
    assert json.loads(capsys.readouterr().out)["dockq"] == 0.5
