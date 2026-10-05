import json
from pathlib import Path

import pytest

from benchmark.scripts.stage_legacy_tool_environment import stage


def test_stage_records_all_tool_payloads_and_explicit_naccess_profile(tmp_path):
    python2 = Path("/home/rshadi25/.conda/envs/tmalignRosetta/bin/python2.7")
    if not python2.is_file():
        pytest.skip("Python 2.7 staging interpreter is unavailable")
    output = tmp_path / "legacy-env"
    manifest = stage(
        Path("working_version/multiprot/external_tools"),
        output,
        python2,
        naccess_root=Path("external_tools/naccess"),
        extra_python_site=Path("tmp/agent/20260713-investigation-implementation/python27-site-numpy"),
        fiberdock_reduce_helper=Path("working_version/multiprot/external_tools/fiberdock/reduce.3.23.130521"),
    )
    assert manifest["naccess_profile"] == "explicit_compatibility_substitution"
    assert manifest["dependencies"]["numpy"]["available"]
    assert manifest["dependencies"]["mysqldb_compat"]["available"]
    for relative in (
        "external_tools/multiprot/multiprot.Linux",
        "external_tools/naccess/accall",
        "external_tools/pops/bin/pops",
        "external_tools/fiberdock/FiberDock",
        "external_tools/fiberdock/nma",
    ):
        assert (output / relative).is_file()
    assert (output / "external_tools/naccess/accall").stat().st_mode & 0o111
    loaded = json.loads((output / "environment_manifest.json").read_text())
    assert loaded["status"] in {"ready_compatibility_naccess_toolchain", "ready_python_but_missing_native_dependency"}
    assert not loaded["capabilities"]["fiberdock_energy_only"]
    assert loaded["capabilities"]["fiberdock_energy_only_loadable"]
    assert not loaded["capabilities"]["fiberdock_full_refinement"]
    assert loaded["capabilities"]["fiberdock_refinement_architecture"]["fiberdock/nma"] == "ELF 32-bit"
    assert "fiberdock/full_refinement:not_validated_end_to_end" in loaded["capabilities"]["fiberdock_full_refinement_blockers"]
    substitution = loaded["fiberdock_reduce_helper_substitution"]
    assert substitution["reason"] == "explicit exploratory FiberDock reduce.3 substitution; not historical-equivalent"
    assert substitution["effective_sha256"] == "6d066f88bff740627d7c1d2fb0200d326978fe70a0c041a2528c8853e682b1ce"
    assert substitution["original_sha256"] != substitution["effective_sha256"]
    executable_record = next(record for record in loaded["source_records"] if record["effective_path"].endswith("external_tools/naccess/accall"))
    assert executable_record["effective_mode"] == "0o755"


def test_stage_records_optional_native_library_root(tmp_path):
    python2 = Path("/home/rshadi25/.conda/envs/tmalignRosetta/bin/python2.7")
    library_root = Path("/opt/ohpc/pub/compiler/gcc/6.5.0/lib64")
    if not python2.is_file():
        pytest.skip("Python 2.7 staging interpreter is unavailable")
    if not (library_root / "libgfortran.so.3.0.0").is_file():
        pytest.skip("historical GCC 6 Fortran runtime is unavailable")
    output = tmp_path / "legacy-env-with-native-libs"
    manifest = stage(
        Path("working_version/multiprot/external_tools"),
        output,
        python2,
        extra_python_site=Path("tmp/agent/20260713-investigation-implementation/python27-site-numpy"),
        native_library_root=library_root,
    )
    assert manifest["native_library_root"] == str(library_root.resolve())
    assert manifest["native_library_files"]["libgfortran.so.3"] == str((library_root / "libgfortran.so.3.0.0").resolve())
    assert str(library_root) in (output / "activate.sh").read_text()


def test_historical_naccess_uses_staged_fortran_library(tmp_path):
    python2 = Path("/home/rshadi25/.conda/envs/tmalignRosetta/bin/python2.7")
    library_root = Path("/opt/ohpc/pub/compiler/gcc/6.5.0/lib64")
    if not python2.is_file():
        pytest.skip("Python 2.7 staging interpreter is unavailable")
    if not (library_root / "libgfortran.so.3.0.0").is_file():
        pytest.skip("historical GCC 6 Fortran runtime is unavailable")
    output = tmp_path / "historical-env-with-native-libs"
    manifest = stage(
        Path("working_version/multiprot/external_tools"),
        output,
        python2,
        extra_python_site=Path("tmp/agent/20260713-investigation-implementation/python27-site-numpy"),
        native_library_root=library_root,
    )
    assert manifest["naccess_profile"] == "historical_working_version"
    assert manifest["native_tools"]["naccess/accall"]["ldd"]["dependency_status"] != "missing"
    assert manifest["capabilities"]["naccess_surface_extraction"]
