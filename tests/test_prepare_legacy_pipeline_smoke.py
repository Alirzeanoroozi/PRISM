from pathlib import Path

from benchmark.scripts.prepare_legacy_pipeline_smoke import prepare


def test_prepare_preserves_logical_selectors_and_template_chain_assets(tmp_path):
    compat = tmp_path / "compat"
    (compat / "run_files").mkdir(parents=True)
    (compat / "run_files" / "compat_runtime.py").write_text("# test\n")
    (compat / "compatibility_manifest.json").write_text("{}\n")
    tools = tmp_path / "tools"
    (tools / "external_tools" / "multiprot").mkdir(parents=True)
    (tools / "external_tools" / "multiprot" / "params.txt").write_text("params\n")
    templates = tmp_path / "templates"
    (templates / "interfaces").mkdir(parents=True)
    (templates / "contact").mkdir()
    (templates / "hotspot").mkdir()
    for name in ("1b27AD_A.int", "1b27AD_D.int"):
        (templates / "interfaces" / name).write_text(name + "\n")
    (templates / "contact" / "1b27AD.txt").write_text("contact\n")
    (templates / "hotspot" / "hotspot1b27AD").write_text("hotspot\n")
    receptor = tmp_path / "1rghB.pdb"
    ligand = tmp_path / "1a19B.pdb"
    receptor.write_text("ATOM receptor\n")
    ligand.write_text("ATOM ligand\n")

    manifest = prepare(
        tmp_path / "workspace",
        compat,
        tools,
        receptor,
        ligand,
        templates,
        template_id="1b27AD",
        logical_receptor_selector="1RGH_B",
        logical_ligand_selector="1A19_B",
    )

    assert manifest["pair_identifiers"] == {
        "receptor": "pdb1",
        "ligand": "pdb2",
        "logical_receptor": "1RGH_B",
        "logical_ligand": "1A19_B",
    }
    assert manifest["template_chains"] == ["A", "D"]
    assert (tmp_path / "workspace/template/interfaces/1b27AD_A.int").is_file()
    assert (tmp_path / "workspace/template/interfaces/1b27AD_D.int").is_file()
