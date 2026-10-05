import os

import pytest

from src import template_generate as tg


def test_process_template_propagates_errors(monkeypatch):
    monkeypatch.setattr(tg, "generate_interface", lambda t: (_ for _ in ()).throw(RuntimeError("boom")))
    name, err = tg.process_template("xyzAB")
    assert name is None
    assert "xyzAB" in err and "boom" in err


def test_process_template_success(monkeypatch):
    monkeypatch.setattr(tg, "generate_interface", lambda t: {"A": [], "B": []})
    monkeypatch.setattr(tg, "hotspot_creator", lambda t: {"A": [], "B": []})
    monkeypatch.setattr(tg, "get_contacts", lambda t: {})
    name, err = tg.process_template("xyzAB")
    assert name == "xyzAB"
    assert err is None


def test_template_generator_missing_input(monkeypatch, tmp_path):
    monkeypatch.setattr("builtins.open", lambda *a, **k: (_ for _ in ()).throw(FileNotFoundError(f"no such file: {a[0]}")))
    with pytest.raises(FileNotFoundError):
        tg.template_generator()


def test_template_generator_runs(monkeypatch, tmp_path):
    inp_path = tmp_path / "checked_templates.txt"
    out_path = tmp_path / "calculated_templates.txt"
    inp_path.parent.mkdir(parents=True, exist_ok=True)
    inp_path.write_text("ok1AB\nfailAB\n")

    orig_open = open
    monkeypatch.setattr("builtins.open", lambda *a, **k: orig_open(str(a[0]).replace("templates/checked_templates.txt", str(inp_path)).replace("templates/calculated_templates.txt", str(out_path)), *a[1:], **k) if a[0].startswith("templates/") else orig_open(*a, **k))

    def fake_process(t):
        if t.startswith("fail"):
            return None, f"{t}: explosion"
        return t, None
    monkeypatch.setattr(tg, "process_template", fake_process)

    result = tg.template_generator()
    assert result == ["ok1AB"]
    assert out_path.read_text().strip() == "ok1AB"
