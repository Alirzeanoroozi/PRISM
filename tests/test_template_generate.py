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
    with pytest.raises(FileNotFoundError):
        tg.template_generator(input_list=str(tmp_path / "nope.txt"),
                              output_list=str(tmp_path / "out.txt"))


def test_template_generator_runs(monkeypatch, tmp_path):
    inp = tmp_path / "in.txt"
    out = tmp_path / "out.txt"
    inp.write_text("ok1AB\nfailAB\n")

    def fake_process(t):
        if t.startswith("fail"):
            return None, f"{t}: explosion"
        return t, None
    monkeypatch.setattr(tg, "process_template", fake_process)

    result = tg.template_generator(input_list=str(inp), output_list=str(out), max_workers=1)
    assert result == ["ok1AB"]
    assert out.read_text().strip() == "ok1AB"
