import hashlib
import json
import os
from pathlib import Path

from benchmark.scripts.multiprot_compat_adapter import create_snapshot


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_adapter_keeps_pristine_tree_untouched_and_suppresses_side_effect_modules(tmp_path):
    source = Path("working_version/multiprot").resolve()
    output = tmp_path / "compat"
    source_controller_hash = _sha(source / "run_files/mainController.py")
    value = create_snapshot(source, output)

    assert _sha(source / "run_files/mainController.py") == source_controller_hash
    assert value["disabled_behaviors"] == ["database", "mail", "html", "network_download", "destructive_cleanup"]
    effective = "\n".join(path.read_text(encoding="utf-8", errors="replace") for path in (output / "run_files").glob("*.py"))
    assert "MySQLdb" not in effective
    assert "smtplib" not in effective
    assert "FTP(" not in effective
    assert 'os.system("rm' not in effective
    manifest = json.loads((output / "compatibility_manifest.json").read_text(encoding="utf-8"))
    assert manifest["algorithm_modules_boundary_only"]
    assert _sha(output / "run_files/flexibleRefinement.py") == _sha(source / "run_files/flexibleRefinement.py")


def test_compat_runtime_reads_staged_lists_without_network_or_deletion(tmp_path, monkeypatch):
    source = Path("working_version/multiprot").resolve()
    output = tmp_path / "compat"
    create_snapshot(source, output)
    job = tmp_path / "job"
    (job / "lists").mkdir(parents=True)
    (job / "lists/pair_list").write_text("1ABC_A 2DEF_B\n", encoding="utf-8")
    (job / "lists/template_list").write_text("1kcaCH\n", encoding="utf-8")
    events = tmp_path / "events.jsonl"
    monkeypatch.setenv("PRISM_COMPAT_EVENT_LOG", str(events))

    import importlib.util

    spec = importlib.util.spec_from_file_location("compat_runtime_test", output / "run_files/compat_runtime.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    left, right, templates = module.PDBdownload(str(job)).PDBdownloader()
    sentinel = job / "2_sol.res"
    sentinel.write_text("preserve\n", encoding="utf-8")
    module.compat_cleanup("2_sol.res")
    assert left == ["1abcA"]
    assert right == ["2defB"]
    assert templates == ["1kcaCH"]
    assert sentinel.read_text(encoding="utf-8") == "preserve\n"
    lines = events.read_text(encoding="utf-8").splitlines()
    assert any("network_download_suppressed" in line for line in lines)
    assert any("cleanup_suppressed" in line for line in lines)

    (job / "template_default").write_text("1kcaCH\n", encoding="utf-8")
    assert module.TemplateChecker(str(job), ["1kcaCH"]).checker() == [1, ["1kcaCH"]]

    (job / "template_default").unlink()
    (tmp_path / "template_default").write_text("1kcaCH\n", encoding="utf-8")
    assert module.TemplateChecker(str(job), ["1kcaCH"]).checker() == [1, ["1kcaCH"]]

    nested_job = tmp_path / "jobs" / "smoke"
    nested_job.mkdir(parents=True)
    assert module.TemplateChecker(str(nested_job), ["1kcaCH"]).checker() == [1, ["1kcaCH"]]


def test_snapshot_stages_complete_executable_multiprot_runtime(tmp_path):
    source = Path("working_version/multiprot").resolve()
    output = tmp_path / "compat"
    create_snapshot(source, output)
    binary = output / "external_tools/multiprot/multiprot.Linux"
    assert binary.is_file()
    assert binary.stat().st_mode & 0o111
    assert (output / "external_tools/multiprot/params.txt").is_file()
