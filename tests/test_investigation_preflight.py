import json

from benchmark.scripts.run_investigation_preflight import main


def test_preflight_supports_named_current_and_historical_arms(tmp_path):
    current_root = tmp_path / "current"
    legacy_root = tmp_path / "legacy"
    (current_root / "contacts").mkdir(parents=True)
    (current_root / "interfaces_lists").mkdir()
    (current_root / "interfaces").mkdir()
    (current_root / "contacts" / "newAB.json").write_text("{}\n", encoding="utf-8")
    (current_root / "interfaces_lists" / "newAB.json").write_text("{}\n", encoding="utf-8")
    (current_root / "interfaces" / "newAB_A_int.pdb").write_text("ATOM\n", encoding="utf-8")
    (legacy_root / "contact").mkdir(parents=True)
    (legacy_root / "contact" / "oldCD.txt").write_text("1 2\n", encoding="utf-8")
    current_manifest = tmp_path / "current.txt"
    legacy_manifest = tmp_path / "legacy.txt"
    current_manifest.write_text("newAB\n", encoding="utf-8")
    legacy_manifest.write_text("oldCD\n", encoding="utf-8")
    output = tmp_path / "out"

    assert main([
        "--output-dir", str(output),
        "--arm", f"current={current_manifest}={current_root}",
        "--arm", f"historical={legacy_manifest}={legacy_root}",
        "--require-all-templates",
    ]) == 0
    summary = json.loads((output / "preflight_summary.json").read_text())
    assert summary["template_preflight_arms"]["current"]["fully_resolvable"] == 1
    assert summary["template_preflight_arms"]["historical"]["fully_resolvable"] == 1
    assert (output / "template_assets_current.tsv").exists()
    assert (output / "template_assets_historical.tsv").exists()
