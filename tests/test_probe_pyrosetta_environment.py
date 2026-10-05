import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmark/scripts/probe_pyrosetta_environment.py"


def test_probe_script_emits_json_status_without_changing_environment_declarations(tmp_path):
    report_path = tmp_path / "pyrosetta_probe.json"
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(report_path)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0
    report = json.loads(report_path.read_text())
    assert report["package"] == "pyrosetta"
    assert report["status"] in {"available", "unavailable"}
    assert isinstance(report["available"], bool)
    assert "command_metadata" in report
    assert "environment_metadata" in report
    assert "pyrosetta" not in (ROOT / "environment.yaml").read_text().lower()
