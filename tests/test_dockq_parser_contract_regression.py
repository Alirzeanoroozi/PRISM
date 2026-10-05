import json
import hashlib
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark.scripts import dockq as benchmark_dockq_module
from src.eval import dockq as dockq_module


def _run_with_payload(payload, captured):
    def fake_run(command, **kwargs):
        captured["command"] = command
        json_path = Path(command[command.index("--json") + 1])
        json_path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    return fake_run


def test_calculate_dockq_reports_global_average_not_interface_sum(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text("MODEL\nENDMDL\n", encoding="ascii")
    native.write_text("MODEL\nENDMDL\n", encoding="ascii")
    payload = {
        "best_dockq": 1.4,
        "GlobalDockQ": 0.7,
        "best_result": [
            {"interface": "A:B", "DockQ": 0.8, "iRMSD": 1.0},
            {"interface": "A:C", "DockQ": 0.6, "iRMSD": 2.0},
        ],
    }
    monkeypatch.setattr(
        dockq_module.subprocess,
        "run",
        _run_with_payload(payload, {}),
    )
    monkeypatch.setenv("DOCKQ_BIN", "DockQ")

    result = dockq_module.calculate_dockq(model, native, work_dir=tmp_path / "dockq")

    assert result["dockq"] == pytest.approx(0.7)
    assert 0.0 <= result["dockq"] <= 1.0


def test_calculate_dockq_does_not_promote_first_interface_components(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text("MODEL\nENDMDL\n", encoding="ascii")
    native.write_text("MODEL\nENDMDL\n", encoding="ascii")
    payload = {
        "best_dockq": 1.4,
        "GlobalDockQ": 0.7,
        "best_result": [
            {
                "interface": "A:B",
                "DockQ": 0.8,
                "iRMSD": 1.0,
                "LRMSD": 2.0,
                "fnat": 0.4,
                "F1": 0.5,
                "clashes": 0,
            },
            {
                "interface": "A:C",
                "DockQ": 0.6,
                "iRMSD": 2.0,
                "LRMSD": 4.0,
                "fnat": 0.2,
                "F1": 0.3,
                "clashes": 1,
            },
        ],
    }
    monkeypatch.setattr(
        dockq_module.subprocess,
        "run",
        _run_with_payload(payload, {}),
    )
    monkeypatch.setenv("DOCKQ_BIN", "DockQ")

    result = dockq_module.calculate_dockq(model, native, work_dir=tmp_path / "dockq")

    assert result["fnat"] is None
    assert result["irmsd"] is None
    assert result["lrmsd"] is None
    assert result["f1"] is None
    assert result["clashes"] is None


@pytest.mark.parametrize("module", [dockq_module, benchmark_dockq_module])
def test_missing_global_is_unscored_even_for_one_interface(module):
    result = module._parse_dockq_json(
        {"best_result": {"AB": {"DockQ": 0.4, "iRMSD": 1.0}}}
    )

    assert result["dockq"] is None
    assert result["dockq_global"] is None
    assert result["dockq_json_status"] == "valid_unscored"


def test_calculate_dockq_passes_explicit_cpu_count(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text("MODEL\nENDMDL\n", encoding="ascii")
    native.write_text("MODEL\nENDMDL\n", encoding="ascii")
    captured = {}
    monkeypatch.setattr(
        dockq_module.subprocess,
        "run",
        _run_with_payload({"GlobalDockQ": 0.5, "best_result": []}, captured),
    )
    monkeypatch.setenv("DOCKQ_BIN", "DockQ")

    dockq_module.calculate_dockq(
        model,
        native,
        work_dir=tmp_path / "dockq",
        n_cpu=1,
    )

    assert captured["command"][captured["command"].index("--n_cpu") + 1] == "1"


@pytest.mark.parametrize("module", [dockq_module, benchmark_dockq_module])
def test_calculate_dockq_defaults_to_one_cpu(module, tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text("MODEL\nENDMDL\n", encoding="ascii")
    native.write_text("MODEL\nENDMDL\n", encoding="ascii")
    captured = {}
    monkeypatch.setattr(
        module.subprocess,
        "run",
        _run_with_payload({"GlobalDockQ": 0.5, "best_result": []}, captured),
    )
    monkeypatch.setattr(module, "resolve_dockq_executable", lambda: ["DockQ"], raising=False)

    result = module.calculate_dockq(model, native, work_dir=tmp_path / "dockq")

    assert captured["command"][-2:] == ["--n_cpu", "1"]
    assert result["dockq_n_cpu"] == 1


@pytest.mark.parametrize("module", [dockq_module, benchmark_dockq_module])
def test_no_align_rejects_residue_identity_mismatch(module, tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA B   1       5.000   0.000   0.000  1.00  0.00           C  \n",
        encoding="ascii",
    )
    native.write_text(
        "ATOM      1  CA  VAL A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
        "ATOM      2  CA  ALA B   1       5.000   0.000   0.000  1.00  0.00           C  \n",
        encoding="ascii",
    )
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: pytest.fail("DockQ should not run"))

    with pytest.raises(ValueError, match="residue identity"):
        module.calculate_dockq(
            model,
            native,
            mapping="AB:AB",
            work_dir=tmp_path / "dockq",
            no_align=True,
            n_cpu=1,
        )


def test_calculate_dockq_returns_raw_json_and_command_provenance(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text("MODEL\nENDMDL\n", encoding="ascii")
    native.write_text("MODEL\nENDMDL\n", encoding="ascii")
    captured = {}
    monkeypatch.setattr(
        dockq_module.subprocess,
        "run",
        _run_with_payload({"GlobalDockQ": 0.5, "best_result": []}, captured),
    )
    monkeypatch.setenv("DOCKQ_BIN", "DockQ")

    result = dockq_module.calculate_dockq(
        model,
        native,
        mapping="AB:CD",
        work_dir=tmp_path / "dockq",
        n_cpu=1,
    )

    assert Path(result["raw_dockq_json"]).is_file()
    assert result["dockq_argv"] == captured["command"]


def test_benchmark_wrapper_keeps_each_raw_json_and_hash_distinct(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text("MODEL\nENDMDL\n", encoding="ascii")
    native.write_text("MODEL\nENDMDL\n", encoding="ascii")
    calls = []

    def fake_run(command, **kwargs):
        json_path = Path(command[command.index("--json") + 1])
        value = 0.4 + 0.2 * len(calls)
        json_path.write_text(
            json.dumps({"GlobalDockQ": value, "best_result": []}) + "\n",
            encoding="utf-8",
        )
        calls.append(command)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(benchmark_dockq_module.subprocess, "run", fake_run)
    first = benchmark_dockq_module.calculate_dockq(
        model, native, work_dir=tmp_path / "dockq", n_cpu=1
    )
    second = benchmark_dockq_module.calculate_dockq(
        model, native, work_dir=tmp_path / "dockq", n_cpu=1
    )

    first_path = Path(first["raw_dockq_json"])
    second_path = Path(second["raw_dockq_json"])
    assert first_path != second_path
    assert json.loads(first_path.read_text())["GlobalDockQ"] == 0.4
    assert json.loads(second_path.read_text())["GlobalDockQ"] == pytest.approx(0.6)
    assert first["raw_dockq_json_sha256"] == hashlib.sha256(first_path.read_bytes()).hexdigest()
    assert second["raw_dockq_json_sha256"] == hashlib.sha256(second_path.read_bytes()).hexdigest()


@pytest.mark.parametrize("parser_module", [dockq_module, benchmark_dockq_module])
def test_duplicate_dockq_parsers_use_global_average(parser_module):
    result = parser_module._parse_dockq_json(
        {
            "best_dockq": 1.4,
            "GlobalDockQ": 0.7,
            "best_result": [
                {"interface": "A:B", "DockQ": 0.8, "iRMSD": 1.0},
                {"interface": "A:C", "DockQ": 0.6, "iRMSD": 2.0},
            ],
        }
    )

    assert result["dockq"] == pytest.approx(0.7)
