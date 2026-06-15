import os

import pytest

from src import interface as iface_mod
from src import contact as contact_mod
from src.interface import generate_interface
from src.contact import get_contacts, get_contacts_from_atom_lines


def _redirect_outputs(monkeypatch, tmp_path):
    monkeypatch.setattr(iface_mod, "INTERFACE_DIR", str(tmp_path / "interfaces"))
    monkeypatch.setattr(iface_mod, "INTERFACE_LIST_DIR", str(tmp_path / "interfaces_lists"))
    monkeypatch.setattr(contact_mod, "CONTACT_DIR", str(tmp_path / "contacts"))
    for p in (tmp_path / "interfaces", tmp_path / "interfaces_lists", tmp_path / "contacts"):
        os.makedirs(p, exist_ok=True)


def test_generate_interface_two_chains(monkeypatch, tiny_two_chain_pdb):
    pdb_root = tiny_two_chain_pdb
    out_dir = pdb_root / "outputs"
    out_dir.mkdir()
    monkeypatch.setattr(iface_mod, "INTERFACE_DIR", str(out_dir / "interfaces"))
    monkeypatch.setattr(iface_mod, "INTERFACE_LIST_DIR", str(out_dir / "interfaces_lists"))
    os.makedirs(out_dir / "interfaces", exist_ok=True)
    os.makedirs(out_dir / "interfaces_lists", exist_ok=True)

    result = generate_interface("testAB", templates_root=str(pdb_root))
    assert set(result.keys()) == {"A", "B"}
    assert len(result["A"]) >= 1
    assert len(result["B"]) >= 1
    assert os.path.exists(out_dir / "interfaces" / "testAB_A_int.pdb")
    assert os.path.exists(out_dir / "interfaces" / "testAB_B_int.pdb")
    assert os.path.exists(out_dir / "interfaces_lists" / "testAB.json")


def test_get_contacts_two_chains(monkeypatch, tiny_two_chain_pdb):
    pdb_root = tiny_two_chain_pdb
    monkeypatch.setattr(contact_mod, "CONTACT_DIR", str(pdb_root / "contacts"))
    os.makedirs(pdb_root / "contacts", exist_ok=True)

    pair_contacts = get_contacts("testAB", templates_root=str(pdb_root))
    assert ("A", "B") in pair_contacts
    assert len(pair_contacts[("A", "B")]) >= 1


def test_get_contacts_from_atom_lines(tmp_path):
    lines_a = [
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n",
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C  \n",
    ]
    lines_b = [
        "ATOM      3  CA  ALA B   1       2.000   0.000   0.000  1.00  0.00           C  \n",
        "ATOM      4  N   ALA B   2      20.000   0.000   0.000  1.00  0.00           N  \n",
    ]
    out = tmp_path / "contacts.txt"
    contacts = get_contacts_from_atom_lines(str(tmp_path / "ignored.pdb"), str(out), lines_a, lines_b)
    assert ("A", 1, "B", 1) in contacts
    assert all(pair[3] != 2 for pair in contacts)
    assert out.read_text().strip().split("\t") == ["A", "1", "B", "1"]
