from pathlib import Path
from types import SimpleNamespace

from src.eval import dockq as dockq_module


def test_dockq_python_override_is_used_as_module_entrypoint(tmp_path, monkeypatch):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    model.write_text("MODEL\nENDMDL\n", encoding="ascii")
    native.write_text("MODEL\nENDMDL\n", encoding="ascii")

    override = tmp_path / "dockq-env" / "bin" / "python"
    override.parent.mkdir(parents=True)
    override.write_text("placeholder\n", encoding="ascii")
    monkeypatch.setenv("DOCKQ_PYTHON", str(override))
    monkeypatch.delenv("DOCKQ_BIN", raising=False)

    captured = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        json_path = command[command.index("--json") + 1]
        # The production code reads the exact path supplied to DockQ.
        Path(json_path).write_text(
            '{"GlobalDockQ": 0.5, "best_result": []}\n',
            encoding="utf-8",
        )
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(dockq_module.subprocess, "run", fake_run)

    result = dockq_module.calculate_dockq(
        model,
        native,
        mapping="AB:CD",
        work_dir=tmp_path / "dockq-work",
    )

    assert captured["command"][:3] == [str(override), "-m", "DockQ"]
    mapping_index = captured["command"].index("--mapping")
    assert captured["command"][mapping_index:mapping_index + 2] == ["--mapping", "AB:CD"]
    assert captured["command"][-2:] == ["--n_cpu", "1"]
    assert result["dockq"] == 0.5
