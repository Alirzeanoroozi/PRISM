import json

from benchmark.scripts.prepare_distinct_chain_fiberdock_probe import prepare


def _atom(chain):
    return f"ATOM      1  CA  ALA {chain}   1      0.000   0.000   0.000  1.00 20.00           C\nEND\n"


def test_prepare_distinct_chain_probe_isolated_and_hashed(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "smoke_manifest.json").write_text("{}\n")
    (source / "input").mkdir()
    (source / "input" / "pair_list").write_text("pdb1 pdb2\n")
    (source / "input" / "template_list").write_text("old\n")
    receptor = tmp_path / "receptor.pdb"
    ligand = tmp_path / "ligand.pdb"
    receptor.write_text(_atom("B"))
    ligand.write_text(_atom("B"))
    destination = tmp_path / "probe"

    manifest = prepare(source, destination, receptor, ligand)
    assert manifest["status"] == "staged_exploratory_distinct_chain_probe"
    assert (destination / "input" / "pair_list").read_text() == "1rghA 1a19B\n"
    assert (destination / "pdb" / "1rgh.pdb").read_text().splitlines()[0][21] == "A"
    assert (destination / "pdb" / "1a19.pdb").read_text().splitlines()[0][21] == "B"
    assert (destination / "template_default").read_text() == "1b27AD\n"
    if (destination / "run_files" / "compat_runtime.py").is_file():
        assert "for _ in range(3)" in (destination / "run_files" / "compat_runtime.py").read_text()
    assert json.loads((destination / "distinct_chain_probe_manifest.json").read_text())["receptor"]["emitted_chain"] == "A"
