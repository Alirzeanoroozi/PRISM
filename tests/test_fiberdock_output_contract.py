from pathlib import Path
from types import SimpleNamespace

import pytest

from src import fiberdock_refinement


def test_declared_output_prefix_is_parsed_and_collected(monkeypatch, tmp_path):
    work = tmp_path / "work"
    structures = tmp_path / "structures"
    energies = tmp_path / "energies"
    work.mkdir()
    structures.mkdir()
    energies.mkdir()
    params = work / "fd_params.txt"
    params.write_text("parameters\n")
    monkeypatch.setattr(fiberdock_refinement, "FIBERDOCK_DIR", str(tmp_path))
    monkeypatch.setattr(fiberdock_refinement, "STRUCTURES_DIR", str(structures))
    monkeypatch.setattr(fiberdock_refinement, "ENERGIES_DIR", str(energies))
    def successful_run(*args, **kwargs):
        (work / "fiberdock_energies.ref").write_text("1 | 0.51 | -15.12 | 11.67 |\n")
        (work / "fiberdock_energies_1.ref.pdb").write_text("END\n")
        return SimpleNamespace(stdout="ok", stderr="", returncode=0)

    monkeypatch.setattr(fiberdock_refinement.subprocess, "run", successful_run)

    energy = fiberdock_refinement._run_fiberdock(
        str(params), str(tmp_path), str(work), "rec", "lig", pair_name="pair1"
    )

    assert energy == "0.51"
    assert (energies / "pair1.ref").is_file()
    assert (structures / "pair1_fiberdock.ref.pdb").is_file()
    assert not (work / "fd_params.ref").exists()


def test_failed_run_rejects_and_removes_stale_outputs(monkeypatch, tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    params = work / "fd_params.txt"
    params.write_text("parameters\n")
    stale_energy = work / "fiberdock_energies.ref"
    stale_structure = work / "fiberdock_energies_1.ref.pdb"
    stale_energy.write_text("1 | -99.0 |\n")
    stale_structure.write_text("STALE\n")
    monkeypatch.setattr(fiberdock_refinement, "FIBERDOCK_DIR", str(tmp_path))
    monkeypatch.setattr(
        fiberdock_refinement.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout="", stderr="failed", returncode=7),
    )

    with pytest.raises(RuntimeError, match="exit code 7"):
        fiberdock_refinement._run_fiberdock(
            str(params), str(tmp_path), str(work), "rec", "lig", pair_name="pair1"
        )

    assert not stale_energy.exists()
    assert not stale_structure.exists()
