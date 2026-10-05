#!/usr/bin/env python3
"""Stage the legacy PRISM external-tool chain in a derived directory.

The checked-out working-version payload is copied without changing it.  A
separate ``--naccess-root`` may be selected when the historical NACCESS binary
cannot load on the host; that substitution is recorded in the manifest and is
never implicit.  Executable bits are repaired only in the derived tree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import shutil
import stat
import subprocess
import zipfile
from pathlib import Path


TOOL_DIRS = ("multiprot", "naccess", "pops", "fiberdock")
SKIP_GENERATED = {"log_multiprot.txt"}
EXECUTABLE_RELATIVE = (
    "multiprot/multiprot.Linux",
    "multiprot/utils/pdb_trans_all_atoms.Linux",
    "multiprot/utils/pdb_trans_frag.Linux",
    "naccess/accall",
    "naccess/naccess",
    "pops/bin/pops",
    "fiberdock/FiberDock",
    "fiberdock/FiberDock.32",
    "fiberdock/nma",
    "fiberdock/reduce.2.21.030604",
    "fiberdock/reduce.3.23.130521",
    "fiberdock/addHydrogens.pl",
    "fiberdock/buildFiberDockParams.pl",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def command_probe(command: list[str], env: dict[str, str] | None = None) -> dict[str, object]:
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=20,
            env=env,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return {"command": command, "available": False, "returncode": None, "stdout": "", "stderr": str(exc)}
    return {
        "command": command,
        "available": result.returncode == 0,
        "returncode": result.returncode,
        "stdout": result.stdout.strip(),
        "stderr": result.stderr.strip(),
    }


def python_probe(python: Path, module: str | None, env: dict[str, str]) -> dict[str, object]:
    code = "import sys; print(sys.version.replace('\\n', ' ')); print(sys.executable)"
    if module:
        code = f"import {module}; print(getattr({module}, '__version__', 'present'))"
    value = command_probe([str(python), "-c", code], env=env)
    value["module"] = module or "python"
    return value


def ldd_probe(path: Path, env: dict[str, str] | None = None) -> dict[str, object]:
    value = command_probe(["ldd", str(path)], env=env)
    missing = [line.strip() for line in value["stdout"].splitlines() if "not found" in line]
    value["missing_libraries"] = missing
    value["dependency_status"] = "missing" if missing else ("available" if value["returncode"] == 0 else "not_dynamic_or_failed")
    return value


def file_probe(path: Path) -> dict[str, object]:
    """Record executable format so 32-bit helper limits are explicit."""
    return command_probe(["file", str(path)])


def classify_refinement_capability(
    native: dict[str, dict[str, object]],
    *,
    end_to_end_validated: bool = False,
) -> tuple[list[str], dict[str, str]]:
    """Classify FiberDock helpers with explicit legacy-runtime blockers.

    Staging metadata establishes presence and loader resolution only. A
    retained positive full-refinement probe is required before capability is
    reported as available. Static 32-bit helpers are recorded as a blocker for
    the historical capability boundary even when one helper happens to run on
    the current host; this prevents a caller from promoting a mixed-width
    toolchain to the primary arm.
    """
    helpers = (
        "fiberdock/nma",
        "fiberdock/reduce.2.21.030604",
        "fiberdock/reduce.3.23.130521",
        "fiberdock/addHydrogens.pl",
    )
    blockers: list[str] = []
    architecture: dict[str, str] = {}
    for relative in helpers:
        record = native.get(relative, {})
        if not record.get("exists"):
            blockers.append(f"{relative}:missing")
            continue
        file_text = str(record.get("file", {}).get("stdout", ""))
        if "ELF 32-bit" in file_text:
            architecture[relative] = "ELF 32-bit"
            if relative in {
                "fiberdock/nma",
                "fiberdock/reduce.2.21.030604",
                "fiberdock/reduce.3.23.130521",
            }:
                blockers.append(f"{relative}:32-bit-helper")
        if record.get("ldd", {}).get("dependency_status") == "missing":
            blockers.append(f"{relative}:missing-library")
    if not end_to_end_validated:
        blockers.append("fiberdock/full_refinement:not_validated_end_to_end")
    return blockers, architecture


def archive_record(path: Path, member_names: tuple[str, ...]) -> dict[str, object]:
    if not path.is_file():
        return {"path": str(path), "exists": False, "sha256": "", "members": {}}
    members: dict[str, object] = {}
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        for member in member_names:
            if member not in names:
                members[member] = None
                continue
            with archive.open(member) as handle:
                digest = hashlib.sha256()
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
                members[member] = digest.hexdigest()
    return {"path": str(path), "exists": True, "sha256": sha256_file(path), "members": members}


def copy_tree(source: Path, target: Path, records: list[dict[str, object]], label: str) -> None:
    if not source.is_dir():
        raise ValueError(f"missing tool directory: {source}")
    for path in sorted(source.rglob("*")):
        relative = path.relative_to(source)
        if path.is_dir():
            (target / relative).mkdir(parents=True, exist_ok=True)
            continue
        if not path.is_file() or path.name in SKIP_GENERATED:
            continue
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        records.append(
            {
                "role": label,
                "source_path": str(path),
                "source_sha256": sha256_file(path),
                "effective_path": str(destination),
                "effective_sha256": sha256_file(destination),
                "source_mode": oct(stat.S_IMODE(path.stat().st_mode)),
                "effective_mode": oct(stat.S_IMODE(destination.stat().st_mode)),
            }
        )


def stage(
    source_tools_root: Path,
    output_root: Path,
    python_executable: Path,
    naccess_root: Path | None = None,
    extra_python_site: Path | None = None,
    native_library_root: Path | None = None,
    fiberdock_reduce_helper: Path | None = None,
) -> dict[str, object]:
    source_tools = source_tools_root.resolve()
    output = output_root.resolve()
    python = python_executable.resolve()
    selected_naccess = (naccess_root or (source_tools / "naccess")).resolve()
    native_libs = native_library_root.resolve() if native_library_root else None
    reduce_helper = fiberdock_reduce_helper.resolve() if fiberdock_reduce_helper else None
    if not source_tools.is_dir():
        raise ValueError(f"tool source root does not exist: {source_tools}")
    if not python.is_file():
        raise ValueError(f"Python executable does not exist: {python}")
    if native_libs and not native_libs.is_dir():
        raise ValueError(f"native library root does not exist: {native_libs}")
    if reduce_helper and not reduce_helper.is_file():
        raise ValueError(f"FiberDock reduce helper does not exist: {reduce_helper}")
    if output == source_tools or source_tools in output.parents:
        raise ValueError("output must not be the tool source root or inside it")
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite non-empty output: {output}")
    output.mkdir(parents=True, exist_ok=True)

    source_records: list[dict[str, object]] = []
    for name in TOOL_DIRS:
        source = selected_naccess if name == "naccess" else source_tools / name
        copy_tree(source, output / "external_tools" / name, source_records, f"{name}_source")

    compatibility_patches: list[dict[str, object]] = []
    if selected_naccess == (source_tools / "naccess").resolve():
        wrapper = output / "external_tools" / "naccess" / "naccess"
        original = wrapper.read_text(encoding="utf-8", errors="replace")
        relocated = original.replace("set EXE_PATH = /cosbi/web/apps/prism/external_tools/naccess", "set EXE_PATH = $0:h")
        if relocated != original:
            wrapper.write_text(relocated, encoding="utf-8")
            compatibility_patches.append(
                {
                    "path": str(wrapper),
                    "source_sha256": next(record["source_sha256"] for record in source_records if record["effective_path"] == str(wrapper)),
                    "effective_sha256": sha256_file(wrapper),
                    "reason": "relocate historical NACCESS wrapper data files into the derived environment",
                }
            )
            for record in source_records:
                if record["effective_path"] == str(wrapper):
                    record["effective_sha256"] = sha256_file(wrapper)

    reduce_substitution: dict[str, object] | None = None
    if reduce_helper:
        destination = output / "external_tools/fiberdock/reduce.2.21.030604"
        original_sha256 = sha256_file(destination)
        shutil.copy2(reduce_helper, destination)
        reduce_substitution = {
            "source_path": str(reduce_helper),
            "source_sha256": sha256_file(reduce_helper),
            "replaced_path": str(destination),
            "original_sha256": original_sha256,
            "effective_sha256": sha256_file(destination),
            "reason": "explicit exploratory FiberDock reduce.3 substitution; not historical-equivalent",
        }
        compatibility_patches.append(reduce_substitution)

    for relative in EXECUTABLE_RELATIVE:
        path = output / "external_tools" / relative
        if not path.is_file():
            continue
        mode = stat.S_IMODE(path.stat().st_mode) | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        path.chmod(mode)
    for record in source_records:
        effective = Path(record["effective_path"])
        if effective.is_file():
            record["effective_mode"] = oct(stat.S_IMODE(effective.stat().st_mode))

    runtime = output / "runtime"
    bin_dir = output / "bin"
    runtime.mkdir(parents=True, exist_ok=True)
    bin_dir.mkdir(parents=True, exist_ok=True)
    if extra_python_site:
        site_source = extra_python_site.resolve()
        if not site_source.is_dir():
            raise ValueError(f"extra Python site does not exist: {site_source}")
        copy_tree(site_source, output / "python-site", source_records, "extra_python_site")
    else:
        site_source = None

    write_text(
        runtime / "sitecustomize.py",
        """# Generated Python 2 compatibility boundary for the staged legacy toolchain.
try:
    import pymysql
    pymysql.install_as_MySQLdb()
except ImportError:
    pass
""",
    )
    write_text(bin_dir / "python2", f"#!/bin/sh\nset -eu\nexec {shlex.quote(str(python))} \"$@\"\n")
    (bin_dir / "python2").chmod(0o755)
    write_text(
        bin_dir / "multiprot",
        "#!/bin/sh\nset -eu\nROOT=$(CDPATH= cd -- \"$(dirname -- \"$0\")/..\" && pwd)\nexec \"$ROOT/external_tools/multiprot/multiprot.Linux\" \"$@\"\n",
    )
    (bin_dir / "multiprot").chmod(0o755)
    activation_native_path = ""
    if native_libs:
        activation_native_path = f'export LD_LIBRARY_PATH="{native_libs}:${{LD_LIBRARY_PATH:+$LD_LIBRARY_PATH}}"\n'
    write_text(
        output / "activate.sh",
        f"""#!/bin/sh
ENV_FILE="${{BASH_SOURCE:-$0}}"
ENV_ROOT=$(CDPATH= cd -- "$(dirname -- "$ENV_FILE")" && pwd)
export PRISM_LEGACY_ENV="$ENV_ROOT"
export PRISM_MULTIPROT_ENV="$ENV_ROOT"
export PRISM_MULTIPROT_ROOT="$ENV_ROOT/external_tools/multiprot"
export PRISM_MULTIPROT_PYTHON="$ENV_ROOT/bin/python2"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$ENV_ROOT/runtime:$ENV_ROOT/python-site${{PYTHONPATH:+:$PYTHONPATH}}"
export PATH="$ENV_ROOT/bin:$PATH"
{activation_native_path}""",
    )
    (output / "activate.sh").chmod(0o755)

    runtime_env = os.environ.copy()
    native_library_files: dict[str, str] = {}
    if native_libs:
        for name in ("libgfortran.so.3", "libgfortran.so.5"):
            path = native_libs / name
            if path.exists():
                native_library_files[name] = str(path.resolve())
        runtime_env["LD_LIBRARY_PATH"] = str(native_libs) + (os.pathsep + runtime_env["LD_LIBRARY_PATH"] if runtime_env.get("LD_LIBRARY_PATH") else "")
    runtime_env.update(
        {
            "PYTHONNOUSERSITE": "1",
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONPATH": str(runtime) + os.pathsep + str(output / "python-site") + (os.pathsep + runtime_env["PYTHONPATH"] if runtime_env.get("PYTHONPATH") else ""),
        }
    )
    dependencies = {
        "python2.7": python_probe(python, None, runtime_env),
        "pymysql": python_probe(python, "pymysql", runtime_env),
        "numpy": python_probe(python, "numpy", runtime_env),
        "mysqldb_compat": python_probe(python, "MySQLdb", runtime_env),
    }

    native = {}
    for relative in (
        "multiprot/multiprot.Linux",
        "naccess/accall",
        "pops/bin/pops",
        "fiberdock/FiberDock",
        "fiberdock/nma",
        "fiberdock/reduce.2.21.030604",
        "fiberdock/reduce.3.23.130521",
        "fiberdock/addHydrogens.pl",
    ):
        path = output / "external_tools" / relative
        native[relative] = {
            "path": str(path),
            "exists": path.is_file(),
            "sha256": sha256_file(path) if path.is_file() else "",
            "mode": oct(stat.S_IMODE(path.stat().st_mode)) if path.is_file() else "",
            "ldd": ldd_probe(path, env=runtime_env) if path.is_file() else {"dependency_status": "missing"},
            "file": file_probe(path) if path.is_file() else {"available": False, "stdout": ""},
        }

    archive_root = source_tools
    archives = {
        "multiprot": archive_record(archive_root / "multiprot.zip", ("multiprot/multiprot.Linux",)),
        "naccess": archive_record(archive_root / "naccess.zip", ("naccess/accall", "naccess/naccess")),
        "pops": archive_record(archive_root / "pops.zip", ("pops/bin/pops",)),
        "fiberdock": archive_record(archive_root / "fiberdock.zip", ("fiberdock/FiberDock", "fiberdock/nma")),
    }
    generated_paths = [runtime / "sitecustomize.py", bin_dir / "python2", bin_dir / "multiprot", output / "activate.sh"]
    generated = [{"path": str(path), "sha256": sha256_file(path), "mode": oct(stat.S_IMODE(path.stat().st_mode))} for path in generated_paths]
    historical_naccess = selected_naccess == (source_tools / "naccess").resolve()
    missing_native = [key for key, value in native.items() if value["ldd"].get("dependency_status") == "missing"]
    refinement_blockers, refinement_architecture = classify_refinement_capability(native)
    required_python_ok = all(dependencies[name]["available"] for name in ("python2.7", "pymysql", "numpy", "mysqldb_compat"))
    if not required_python_ok:
        status = "not_ready_python_dependency_gate"
    elif missing_native:
        status = "ready_python_but_missing_native_dependency"
    elif historical_naccess:
        status = "ready_historical_legacy_toolchain"
    else:
        status = "ready_compatibility_naccess_toolchain"
    manifest = {
        "schema_version": "legacy-prism-tool-environment/v1",
        "status": status,
        "source_tools_root": str(source_tools),
        "output_root": str(output),
        "python_executable": str(python),
        "selected_naccess_source": str(selected_naccess),
        "naccess_profile": "historical_working_version" if historical_naccess else "explicit_compatibility_substitution",
        "native_library_root": str(native_libs) if native_libs else "",
        "native_library_files": native_library_files,
        "source_records": source_records,
        "compatibility_patches": compatibility_patches,
        "fiberdock_reduce_helper_substitution": reduce_substitution,
        "generated_records": generated,
        "native_tools": native,
        "source_archives": archives,
        "dependencies": dependencies,
        "compatibility": {
            "network_access": False,
            "database_access": False,
            "mysql_import": "PyMySQL install_as_MySQLdb shim",
            "extra_python_site": str(site_source) if site_source else "",
        },
        "reference_requirements": {
            "python": "Python 2.5.x tested by PRISM reference",
            "numpy": "required by PRISM reference protocol",
            "multiprot": "1.6",
            "naccess": "2.1-compatible executable/wrapper",
            "pops": "POPS surface extractor",
            "fiberdock": "FiberDock 1.0 tool payload",
        },
        "capabilities": {
            "naccess_surface_extraction": native["naccess/accall"]["ldd"].get("dependency_status") != "missing",
            "pops_surface_extraction": native["pops/bin/pops"]["ldd"].get("dependency_status") != "missing",
            "multiprot_alignment": native["multiprot/multiprot.Linux"]["ldd"].get("dependency_status") != "missing",
            # Loadability is not execution evidence. The direct energy-only
            # probe is retained separately and must not be inferred here.
            "fiberdock_energy_only": False,
            "fiberdock_energy_only_loadable": native["fiberdock/FiberDock"]["ldd"].get("dependency_status") != "missing",
            "fiberdock_full_refinement": not refinement_blockers,
            "fiberdock_full_refinement_blockers": refinement_blockers,
            "fiberdock_refinement_architecture": refinement_architecture,
        },
    }
    (output / "environment_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_text(
        output / "README.md",
        f"""# Staged legacy PRISM tool environment\n\nThis derived environment was copied from `{source_tools}`.\nThe selected NACCESS source is `{selected_naccess}` and is recorded in\n`environment_manifest.json`; no source tree was modified.\n\nActivate it with:\n\n```sh\nsource {output}/activate.sh\npython2 -c 'import numpy, MySQLdb; print(numpy.__version__)'\n```\n\nRun every external tool in a unique scratch directory.  The historical\nNACCESS binary may require a host library that is not present; the explicit\ncompatibility profile can be staged separately and must not be conflated with\nthe historical profile in scientific comparisons.\n""",
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-tools-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--python-executable", type=Path, required=True)
    parser.add_argument("--naccess-root", type=Path)
    parser.add_argument("--extra-python-site", type=Path)
    parser.add_argument("--native-library-root", type=Path)
    parser.add_argument(
        "--fiberdock-reduce-helper",
        type=Path,
        help="Explicit exploratory helper to stage at the historical reduce.2 path; never treated as historical-equivalent.",
    )
    args = parser.parse_args(argv)
    try:
        value = stage(
            args.source_tools_root,
            args.output_root,
            args.python_executable,
            naccess_root=args.naccess_root,
            extra_python_site=args.extra_python_site,
            native_library_root=args.native_library_root,
            fiberdock_reduce_helper=args.fiberdock_reduce_helper,
        )
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_root": value["output_root"], "status": value["status"], "naccess_profile": value["naccess_profile"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
