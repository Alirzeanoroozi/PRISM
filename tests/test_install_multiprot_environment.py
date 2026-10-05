import json
from pathlib import Path

import pytest

from benchmark.scripts.install_multiprot_environment import stage


def test_stage_records_executable_binary_and_dependency_gate(tmp_path):
    python2 = Path("/home/rshadi25/.conda/envs/tmalignRosetta/bin/python2.7")
    if not python2.is_file():
        pytest.skip("the repository host does not provide the validated Python 2.7 runtime")
    output = tmp_path / "multiprot-env"
    manifest = stage(Path("working_version/multiprot/external_tools/multiprot"), output, python2)
    assert manifest["status"].startswith("ready_for_standalone_multiprot")
    assert manifest["dependencies"]["pymysql"]["available"]
    assert manifest["dependencies"]["mysqldb_compat"]["available"]
    assert (output / "multiprot/multiprot.Linux").stat().st_mode & 0o111
    assert (output / "multiprot/params.txt").is_file()
    assert len(manifest["generated_records"]) == 4
    assert json.loads((output / "environment_manifest.json").read_text())["binary_sha256"] == manifest["binary_sha256"]
