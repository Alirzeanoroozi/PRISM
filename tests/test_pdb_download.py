import argparse
import os
import pandas as pd
import pytest

from src import pdb_download
from src.pdb_download import split_target_id, pdb_downloader


def test_split_target_id_basic():
    assert split_target_id("3i6eE") == ("3i6e", ["E"])
    assert split_target_id("3i6eEF") == ("3i6e", ["E", "F"])
    assert split_target_id("1fgnHLA") == ("1fgn", ["H", "L", "A"])
    assert split_target_id("3i6e") == ("3i6e", [])


def test_split_target_id_strips_whitespace():
    assert split_target_id("  3i6eE  ") == ("3i6e", ["E"])


def test_split_target_id_too_short():
    with pytest.raises(ValueError):
        split_target_id("abc")


def test_pdb_downloader_multi_chain(monkeypatch, tmp_path):
    csv_path = tmp_path / "inputs.csv"
    pd.DataFrame({"Receptor": ["1abcAB", "2xyzCD"], "Ligand": ["1abcEF", "2xyzGH"]}).to_csv(csv_path, index=False)

    monkeypatch.setattr(pdb_download, "TARGET_DIR", str(tmp_path / "pdbs"))
    os.makedirs(tmp_path / "pdbs", exist_ok=True)

    calls = []
    def fake_download(pdb_id, pdb_dir):
        calls.append(pdb_id)
        open(f"{pdb_dir}/{pdb_id}.pdb", "w").write("")
        return True
    monkeypatch.setattr(pdb_download, "download_pdb_file", fake_download)

    args = argparse.Namespace(inputs_csv=str(csv_path))
    receptors, ligands = pdb_downloader(args)
    assert receptors == ["1abcAB", "2xyzCD"]
    assert ligands == ["1abcEF", "2xyzGH"]
    assert set(calls) == {"1abc", "2xyz"}


def test_pdb_downloader_skip_invalid(monkeypatch, tmp_path):
    csv_path = tmp_path / "inputs.csv"
    pd.DataFrame({"Receptor": ["bad", "1abcA"], "Ligand": ["1abcB", "1abcC"]}).to_csv(csv_path, index=False)

    monkeypatch.setattr(pdb_download, "TARGET_DIR", str(tmp_path / "pdbs"))
    os.makedirs(tmp_path / "pdbs", exist_ok=True)
    monkeypatch.setattr(pdb_download, "download_pdb_file", lambda *a, **k: True)

    args = argparse.Namespace(inputs_csv=str(csv_path))
    receptors, ligands = pdb_downloader(args)
    assert receptors == ["1abcA"]
    assert ligands == ["1abcC"]
