import hashlib
import importlib
import json
from pathlib import Path

import src.pyrosetta_refinement as pyrosetta_refinement


def test_adapter_module_does_not_import_pyrosetta_at_module_import(monkeypatch):
    calls = []

    def fail_if_imported(name):
        calls.append(name)
        raise AssertionError("PyRosetta was imported while loading the adapter")

    monkeypatch.setattr(pyrosetta_refinement.importlib, "import_module", fail_if_imported)
    importlib.reload(pyrosetta_refinement)

    assert calls == []


def test_probe_returns_fail_closed_unavailable_status_and_import_error(monkeypatch):
    def missing_pyrosetta(name):
        raise ModuleNotFoundError("No module named 'pyrosetta'")

    monkeypatch.setattr(pyrosetta_refinement.importlib, "import_module", missing_pyrosetta)

    report = pyrosetta_refinement.probe_environment()

    assert report["package"] == "pyrosetta"
    assert report["status"] == "unavailable"
    assert report["available"] is False
    assert "No module named 'pyrosetta'" in report["import_error"]
    assert report["version"] is None
    assert "command_metadata" in report
    assert "environment_metadata" in report


def test_unavailable_adapter_does_not_write_structure_or_fallback_to_cli(tmp_path, monkeypatch):
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    input_path.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n")

    def missing_pyrosetta(name):
        raise ImportError("PyRosetta license/runtime unavailable")

    monkeypatch.setattr(pyrosetta_refinement.importlib, "import_module", missing_pyrosetta)
    monkeypatch.setattr(
        pyrosetta_refinement.os,
        "system",
        lambda command: (_ for _ in ()).throw(AssertionError(f"CLI fallback: {command}")),
    )
    adapter = pyrosetta_refinement.PyRosettaRefinementAdapter()

    result = adapter.refine(input_path, output_path)

    assert result["backend"] == "pyrosetta"
    assert result["status"] == "unavailable"
    assert result["available"] is False
    assert result["fallback"] is None
    assert not output_path.exists()
    assert result["input_hashes"][str(input_path)] == hashlib.sha256(input_path.read_bytes()).hexdigest()
    assert result["output_hashes"] == {}
    assert Path(result["metadata_path"]).is_file()
    metadata = json.loads(Path(result["metadata_path"]).read_text())
    assert metadata["import_error"] == result["import_error"]


def test_adapter_rejects_preexisting_output_without_hashing_it(tmp_path, monkeypatch):
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    input_path.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n")
    output_path.write_text("stale output\n")

    monkeypatch.setattr(
        pyrosetta_refinement.importlib,
        "import_module",
        lambda name: (_ for _ in ()).throw(ImportError("PyRosetta unavailable")),
    )
    result = pyrosetta_refinement.PyRosettaRefinementAdapter().refine(input_path, output_path)

    assert result["status"] == "failed"
    assert "already exists" in result["error"]
    assert result["output_hashes"] == {}
    assert output_path.read_text() == "stale output\n"


def test_successful_adapter_records_input_output_hashes_and_runtime_metadata(tmp_path, monkeypatch):
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    input_text = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n"
    input_path.write_text(input_text)

    class FakePose:
        def dump_pdb(self, path):
            Path(path).write_text(input_text.replace("0.000", "1.000", 1))

    class FakeScoreFunction:
        def __call__(self, pose):
            return -3.25

    class FakeDockingProtocol:
        def set_docking_local_refine(self):
            self.local_refine = True

        def set_partners(self, partners):
            self.partners = partners

        def set_scorefxn(self, scorefxn):
            self.scorefxn = scorefxn

        def apply(self, pose):
            return None

    class FakePyRosetta:
        __version__ = "fake-1.0"
        rosetta = type(
            "Rosetta",
            (),
            {"protocols": type("Protocols", (), {"docking": type("Docking", (), {"DockingProtocol": FakeDockingProtocol})})},
        )

        @staticmethod
        def init(options):
            assert "mute" in options

        @staticmethod
        def pose_from_pdb(path):
            assert Path(path) == input_path
            return FakePose()

        @staticmethod
        def get_fa_scorefxn():
            return FakeScoreFunction()

    monkeypatch.setattr(pyrosetta_refinement.importlib, "import_module", lambda name: FakePyRosetta)

    result = pyrosetta_refinement.PyRosettaRefinementAdapter().refine(
        input_path,
        output_path,
        partners="A_B",
    )

    assert result["status"] == "success"
    assert result["available"] is True
    assert result["version"] == "fake-1.0"
    assert result["input_hashes"][str(input_path)] == hashlib.sha256(input_path.read_bytes()).hexdigest()
    assert result["output_hashes"][str(output_path)] == hashlib.sha256(output_path.read_bytes()).hexdigest()
    assert result["total_score"] == -3.25
    assert result["command_metadata"]["executable"]
    assert result["environment_metadata"]["cwd"]
    assert Path(result["metadata_path"]).is_file()
