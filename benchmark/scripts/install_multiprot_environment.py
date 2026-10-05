#!/usr/bin/env python3
"""Stage a reproducible local MultiProt 1.6 runtime without changing inputs.

The historical MultiProt distribution is copied from the supplied working
version into a new derived directory. The source tree is never edited. The
runtime records source and effective hashes, uses an existing Python 2.7
interpreter, and exposes PyMySQL as the legacy ``MySQLdb`` module when present.
NumPy is checked and reported explicitly because the PRISM reference lists it
as a requirement; it is not silently downloaded or substituted.
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_version(executable: str) -> str:
    try:
        result = subprocess.run(
            [executable, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return f"unavailable:{type(exc).__name__}:{exc}"
    output = result.stdout.strip().splitlines()
    return output[0] if output else f"exit:{result.returncode}"


def python_probe(python_executable: Path, module: str | None = None, env: dict[str, str] | None = None) -> dict[str, object]:
    if module:
        code = f"import {module}; print(getattr({module}, '__version__', 'present'))"
    else:
        code = "import sys; print(sys.version.replace('\\n', ' ')); print(sys.executable)"
    result = subprocess.run(
        [str(python_executable), "-c", code],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
        timeout=20,
        env=env,
    )
    return {
        "available": result.returncode == 0,
        "returncode": result.returncode,
        "stdout": result.stdout.strip(),
        "stderr": result.stderr.strip(),
    }


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def archive_entry_sha256(archive: Path, member: str) -> str | None:
    if not archive.is_file():
        return None
    with zipfile.ZipFile(archive) as handle:
        try:
            with handle.open(member) as payload:
                digest = hashlib.sha256()
                for chunk in iter(lambda: payload.read(1024 * 1024), b""):
                    digest.update(chunk)
                return digest.hexdigest()
        except KeyError:
            return None


def stage(source_root: Path, output_root: Path, python_executable: Path) -> dict[str, object]:
    source = source_root.resolve()
    output = output_root.resolve()
    python = python_executable.resolve()
    if not source.is_dir():
        raise ValueError(f"MultiProt source directory does not exist: {source}")
    if not (source / "multiprot.Linux").is_file():
        raise ValueError(f"MultiProt binary is missing: {source / 'multiprot.Linux'}")
    if not (source / "params.txt").is_file():
        raise ValueError(f"MultiProt parameter file is missing: {source / 'params.txt'}")
    if not python.is_file():
        raise ValueError(f"Python executable does not exist: {python}")
    version = python_probe(python)
    if not version["available"] or not version["stdout"].startswith("2.7"):
        raise ValueError(f"expected Python 2.7, got: {version}")
    if output == source or source in output.parents:
        raise ValueError("output must not be the source directory or inside it")
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite non-empty output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    archive = source.parent / "multiprot.zip"

    source_records: list[dict[str, object]] = []
    excluded: list[str] = []
    for path in sorted(source.rglob("*")):
        relative = path.relative_to(source)
        if path.is_dir():
            (output / "multiprot" / relative).mkdir(parents=True, exist_ok=True)
            continue
        if not path.is_file():
            continue
        if path.name == "log_multiprot.txt":
            excluded.append(str(relative))
            continue
        target = output / "multiprot" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        mode = stat.S_IMODE(path.stat().st_mode)
        if path.name == "multiprot.Linux" or path.suffix == ".Linux":
            mode |= stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        target.chmod(mode)
        source_records.append(
            {
                "relative_path": str(relative),
                "source_path": str(path),
                "source_sha256": sha256_file(path),
                "effective_path": str(target),
                "effective_sha256": sha256_file(target),
                "effective_mode": oct(stat.S_IMODE(target.stat().st_mode)),
            }
        )

    runtime = output / "runtime"
    bin_dir = output / "bin"
    runtime.mkdir(parents=True, exist_ok=True)
    bin_dir.mkdir(parents=True, exist_ok=True)
    write_text(
        runtime / "sitecustomize.py",
        """# Generated compatibility shim for the historical Python 2 runtime.
try:
    import pymysql
    pymysql.install_as_MySQLdb()
except ImportError:
    pass
""",
    )
    write_text(
        bin_dir / "multiprot",
        """#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
exec "$ROOT/multiprot/multiprot.Linux" "$@"
""",
    )
    (bin_dir / "multiprot").chmod(0o755)
    write_text(bin_dir / "python2", f"#!/bin/sh\nset -eu\nexec {shlex.quote(str(python))} \"$@\"\n")
    (bin_dir / "python2").chmod(0o755)
    write_text(
        output / "activate.sh",
        """#!/bin/sh
ENV_FILE="${BASH_SOURCE:-$0}"
ENV_ROOT=$(CDPATH= cd -- "$(dirname -- "$ENV_FILE")" && pwd)
export PRISM_MULTIPROT_ENV="$ENV_ROOT"
export PRISM_MULTIPROT_ROOT="$ENV_ROOT/multiprot"
export PRISM_MULTIPROT_PYTHON="$ENV_ROOT/bin/python2"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$ENV_ROOT/runtime${PYTHONPATH:+:$PYTHONPATH}"
export PATH="$ENV_ROOT/bin:$PATH"
""",
    )
    (output / "activate.sh").chmod(0o755)

    runtime_env = os.environ.copy()
    runtime_env["PYTHONNOUSERSITE"] = "1"
    runtime_env["PYTHONDONTWRITEBYTECODE"] = "1"
    runtime_env["PYTHONPATH"] = str(runtime) + (os.pathsep + runtime_env["PYTHONPATH"] if runtime_env.get("PYTHONPATH") else "")
    dependencies = {
        "python2.7": {
            "available": True,
            "path": str(python),
            "sha256": sha256_file(python),
            "version": version["stdout"],
        },
        "pymysql": python_probe(python, "pymysql"),
        "numpy": python_probe(python, "numpy"),
        "mysqldb_compat": python_probe(python, "MySQLdb", env=runtime_env),
    }
    status = "ready_for_standalone_multiprot"
    if not dependencies["pymysql"]["available"] or not dependencies["mysqldb_compat"]["available"]:
        status = "not_ready_missing_mysqldb_compatibility"
    elif not dependencies["numpy"]["available"]:
        status = "ready_for_standalone_multiprot_missing_reference_numpy"
    generated_paths = (runtime / "sitecustomize.py", bin_dir / "multiprot", bin_dir / "python2", output / "activate.sh")
    generated_records = [
        {
            "path": str(path),
            "sha256": sha256_file(path),
            "mode": oct(stat.S_IMODE(path.stat().st_mode)),
        }
        for path in generated_paths
    ]
    selected_binary_hash = sha256_file(output / "multiprot/multiprot.Linux")
    archive_binary_hash = archive_entry_sha256(archive, "multiprot/multiprot.Linux") or ""
    manifest = {
        "schema_version": "multiprot-environment/v1",
        "status": status,
        "source_root": str(source),
        "output_root": str(output),
        "selected_binary_source": "checked_out_working_version_directory",
        "source_archive_comparison": {
            "archive_path": str(archive) if archive.is_file() else "",
            "archive_sha256": sha256_file(archive) if archive.is_file() else "",
            "archive_binary_member": "multiprot/multiprot.Linux",
            "archive_binary_sha256": archive_binary_hash,
            "archive_binary_matches_selected": bool(archive_binary_hash) and archive_binary_hash == selected_binary_hash,
        },
        "source_records": source_records,
        "generated_records": generated_records,
        "excluded_generated_files": excluded,
        "binary_sha256": selected_binary_hash,
        "binary_mode": oct(stat.S_IMODE((output / "multiprot/multiprot.Linux").stat().st_mode)),
        "dependencies": dependencies,
        "compatibility": {
            "python_major_minor": "2.7",
            "mysql_import": "PyMySQL install_as_MySQLdb shim",
            "network_access": False,
            "database_access": False,
        },
        "reference_requirements": {
            "python": "2.5.x tested by the PRISM paper/protocol",
            "numpy": "required by the PRISM protocol",
            "multiprot": "1.6",
        },
    }
    (output / "environment_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_text(
        output / "README.md",
        f"""# Staged MultiProt environment

This derived runtime was staged from `{source}` without modifying the source tree.
The Python interpreter is `{python}`.

```sh
source {output}/activate.sh
python2 -c 'import sys; print(sys.version)'
cp "$PRISM_MULTIPROT_ROOT/params.txt" .
mult iprot /path/to/receptor.pdb /path/to/ligand.pdb
```

MultiProt writes `log_multiprot.txt`, `2_sol.res`, and related files in the current
working directory. Run each call in a unique scratch directory. The staged runtime
does not download structures or write to MySQL.

The reference protocol also lists NumPy. Its availability is recorded in
`environment_manifest.json`; this staging step does not download packages.
""".replace("mult iprot", "multiprot"),
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--python-executable", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        manifest = stage(args.source_root, args.output_root, args.python_executable)
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_root": manifest["output_root"], "status": manifest["status"], "binary_sha256": manifest["binary_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
