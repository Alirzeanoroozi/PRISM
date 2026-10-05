#!/usr/bin/env python3
"""Small, deterministic provenance and template-asset preflight helpers.

The module intentionally uses only the Python standard library.  It is useful
both from tests and from a shell command such as::

    python benchmark/scripts/investigation_provenance.py \
        --repo-root . --output provenance.json --file prism.py \
        --executable python3 --package numpy

No timestamps or machine-specific random values are added to the manifest, so
the same inputs produce byte-identical JSON.  Environment values and command
arguments are redacted when their names look secret-bearing.
"""

from __future__ import annotations

import argparse
import configparser
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from collections.abc import Iterable, Mapping
from typing import Any

try:
    from src.provenance.run_evidence import canonical_json as _canonical_json
    from src.provenance.run_evidence import sha256_file
except ModuleNotFoundError as exc:
    if exc.name != "src":
        raise
    # Preserve direct ``python benchmark/scripts/...py`` execution.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from src.provenance.run_evidence import canonical_json as _canonical_json
    from src.provenance.run_evidence import sha256_file


REDACTED = "[REDACTED]"
SECRET_NAME_PARTS = (
    "password",
    "passwd",
    "secret",
    "token",
    "credential",
    "private_key",
    "privatekey",
    "api_key",
    "apikey",
    "access_key",
    "accesskey",
)

DEFAULT_ENVIRONMENT_KEYS = (
    "PATH",
    "PYTHONPATH",
    "PYTHONHASHSEED",
    "CONDA_DEFAULT_ENV",
    "CONDA_PREFIX",
    "PRISM_ALIGNER",
    "PRISM_INPUTS_CSV",
    "PRISM_SCFF_THRESHOLD",
)
SEED_ENVIRONMENT_KEYS = (
    "PYTHONHASHSEED",
    "PRISM_SEED",
    "RANDOM_SEED",
    "NUMPY_SEED",
    "TORCH_SEED",
)
SLURM_RESOURCE_KEYS = (
    "SLURM_ACCOUNT",
    "SLURM_CPUS_PER_TASK",
    "SLURM_GPUS",
    "SLURM_GPUS_PER_NODE",
    "SLURM_JOB_ID",
    "SLURM_JOB_NAME",
    "SLURM_MEM_PER_CPU",
    "SLURM_MEM_PER_NODE",
    "SLURM_NTASKS",
    "SLURM_NTASKS_PER_NODE",
    "SLURM_NODELIST",
    "SLURM_NNODES",
    "SLURM_PARTITION",
    "SLURM_QOS",
    "SLURM_TIME_LIMIT",
)

TSV_FIELDS = (
    "template_id",
    "manifest_indices",
    "listed_count",
    "duplicate_count",
    "unique",
    "valid",
    "fully_resolvable",
    "missing",
    "asset_type",
    "asset_path",
    "exists",
    "format_valid",
    "validation_error",
    "sha256",
    "size",
)


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def hash_files(
    paths: Iterable[str | os.PathLike[str]],
    root: str | os.PathLike[str] | None = None,
) -> list[dict[str, Any]]:
    """Hash files and return stable records, including explicit missing rows.

    Paths are sorted by their normalized spelling.  When *root* is supplied,
    record paths are relative to it when possible; the files are still opened
    using their resolved filesystem paths.
    """

    root_path = Path(root).resolve() if root is not None else None
    normalized = sorted({Path(path).expanduser() for path in paths}, key=lambda p: str(p))
    records = []
    for path in normalized:
        resolved = path.resolve()
        try:
            display = str(resolved.relative_to(root_path)) if root_path else str(resolved)
        except ValueError:
            display = str(resolved)
        record: dict[str, Any] = {"path": display, "exists": resolved.is_file(), "sha256": "", "size": 0}
        if resolved.is_file():
            record["sha256"] = sha256_file(resolved)
            record["size"] = resolved.stat().st_size
        records.append(record)
    return records


def _is_secret_name(name: str) -> bool:
    normalized = name.lower().replace("-", "_")
    return any(part in normalized for part in SECRET_NAME_PARTS)


def _redact_value(name: str, value: Any) -> Any:
    if _is_secret_name(name):
        return REDACTED
    if isinstance(value, Mapping):
        return {str(key): _redact_value(str(key), item) for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))}
    if isinstance(value, (list, tuple)):
        return [_redact_value(name, item) for item in value]
    return value


def capture_environment(
    names: Iterable[str] | None = None,
    environ: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Capture selected environment variables, sorted and secret-redacted."""

    source = os.environ if environ is None else environ
    selected = DEFAULT_ENVIRONMENT_KEYS if names is None else tuple(names)
    return {
        str(name): str(_redact_value(str(name), source[str(name)]))
        for name in sorted(set(str(name) for name in selected))
        if str(name) in source
    }


def capture_seeds(
    seeds: Mapping[str, Any] | Iterable[str] | None = None,
    environ: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Capture explicit seed values or seed-named environment variables."""

    if seeds is None:
        source = os.environ if environ is None else environ
        return {name: source[name] for name in SEED_ENVIRONMENT_KEYS if name in source}
    if isinstance(seeds, Mapping):
        return {str(key): _redact_value(str(key), value) for key, value in sorted(seeds.items(), key=lambda pair: str(pair[0]))}
    source = os.environ if environ is None else environ
    names = sorted(set(str(name) for name in seeds))
    return {name: _redact_value(name, source[name]) for name in names if name in source}


def capture_slurm_resources(environ: Mapping[str, str] | None = None) -> dict[str, str]:
    """Capture the present Slurm job/resource variables without other env data."""

    source = os.environ if environ is None else environ
    return {name: str(source[name]) for name in SLURM_RESOURCE_KEYS if name in source}


def package_versions(packages: Iterable[str] | None) -> dict[str, str | None]:
    """Return installed distribution versions, using ``None`` when absent."""

    if packages is None:
        return {}
    versions: dict[str, str | None] = {}
    for package in sorted(set(str(package) for package in packages)):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def resolve_executable(value: str | os.PathLike[str]) -> str | None:
    """Resolve a command name with ``shutil.which`` or an explicit ``Path``."""

    candidate = Path(value).expanduser()
    found = shutil.which(str(value))
    if found:
        return str(Path(found).resolve())
    if candidate.is_file():
        return str(candidate.resolve())
    return None


def resolve_executables(values: Iterable[str | os.PathLike[str]] | None) -> dict[str, str | None]:
    """Resolve command names/paths into a stable mapping."""

    if values is None:
        return {}
    names = sorted(set(str(value) for value in values))
    return {name: resolve_executable(name) for name in names}


def _run_git(repo_root: Path, *args: str) -> bytes:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=str(repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
            env={**os.environ, "GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"},
        )
    except OSError:
        return b""
    return result.stdout if result.returncode == 0 else b""


def capture_git_provenance(repo_root: str | os.PathLike[str]) -> dict[str, str]:
    """Capture git HEAD, status text, and hashes of status and ``diff HEAD``."""

    root = Path(repo_root).resolve()
    head = _run_git(root, "rev-parse", "HEAD").decode("utf-8", "replace").strip()
    status_bytes = _run_git(root, "status", "--porcelain=v1", "--untracked-files=all")
    diff_bytes = _run_git(root, "diff", "--no-ext-diff", "--binary", "HEAD", "--")
    return {
        "head": head,
        "status": status_bytes.decode("utf-8", "replace"),
        "status_sha256": _sha256_bytes(status_bytes),
        "diff_sha256": _sha256_bytes(diff_bytes),
        "status_hash": _sha256_bytes(status_bytes),
        "diff_hash": _sha256_bytes(diff_bytes),
    }


def _redact_argv(argv: Iterable[Any]) -> list[str]:
    result: list[str] = []
    redact_next = False
    for raw in argv:
        argument = str(raw)
        if redact_next:
            result.append(REDACTED)
            redact_next = False
            continue
        key = argument.split("=", 1)[0].lstrip("-").replace("-", "_")
        if _is_secret_name(key):
            if "=" in argument:
                result.append(argument.split("=", 1)[0] + "=" + REDACTED)
            else:
                result.append(argument)
                redact_next = True
        else:
            result.append(argument)
    return result


def _load_config_file(path: Path) -> Mapping[str, Any]:
    if path.suffix.lower() == ".json":
        loaded = json.loads(path.read_text(encoding="utf-8"))
        return loaded if isinstance(loaded, Mapping) else {"value": loaded}
    parser = configparser.ConfigParser()
    parser.read(path, encoding="utf-8")
    if parser.sections():
        return {f"{section}.{key}": value for section in parser.sections() for key, value in parser.items(section)}
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith(("#", ";")) and "=" in stripped:
            key, value = stripped.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def effective_config(config: Mapping[str, Any] | str | os.PathLike[str] | None) -> dict[str, Any]:
    """Normalize effective configuration key/value pairs deterministically."""

    if config is None:
        return {}
    source: Mapping[str, Any]
    if isinstance(config, (str, os.PathLike)):
        source = _load_config_file(Path(config))
    else:
        source = config
    return {str(key): _redact_value(str(key), source[key]) for key in sorted(source, key=str)}


def build_provenance_manifest(
    repo_root: str | os.PathLike[str] = ".",
    *,
    files: Iterable[str | os.PathLike[str]] | None = None,
    executables: Iterable[str | os.PathLike[str]] | None = None,
    environment: Iterable[str] | None = None,
    packages: Iterable[str] | None = None,
    command: Iterable[Any] | None = None,
    seeds: Mapping[str, Any] | Iterable[str] | None = None,
    config: Mapping[str, Any] | str | os.PathLike[str] | None = None,
    environ: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Build a JSON-serializable provenance manifest for an investigation."""

    root = Path(repo_root).resolve()
    file_paths = [root / Path(path) if not Path(path).is_absolute() else Path(path) for path in (files or ())]
    return {
        "schema_version": 1,
        "repo_root": str(root),
        "git": capture_git_provenance(root),
        "files": hash_files(file_paths, root=root),
        "executables": resolve_executables(executables),
        "environment": capture_environment(environment, environ=environ),
        "packages": package_versions(packages),
        "command": _redact_argv(sys.argv if command is None else command),
        "seeds": capture_seeds(seeds, environ=environ),
        "slurm": capture_slurm_resources(environ=environ),
        "config": effective_config(config),
    }


def write_provenance_manifest(path: str | os.PathLike[str], manifest: Mapping[str, Any]) -> Path:
    """Write a stable, human-readable JSON manifest and return its path."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_canonical_json(manifest) + "\n", encoding="utf-8")
    return output


def _template_id_from_value(value: Any) -> str:
    text = str(value).strip()
    if not text:
        return ""
    text = Path(text).name
    text = re.sub(r"_[A-Za-z0-9]+_int\.pdb$", "", text)
    for suffix in ("_int.pdb", ".txt", ".json", ".pdb"):
        if text.endswith(suffix):
            text = text[: -len(suffix)]
            break
    return text.strip()


def read_template_ids(
    manifest_or_list: str | os.PathLike[str] | Iterable[Any],
    template_column: str | None = None,
) -> list[str]:
    """Read template IDs from CSV/TSV, JSON, a plain list, or row objects."""

    if not isinstance(manifest_or_list, (str, os.PathLike)):
        values = manifest_or_list
    else:
        path = Path(manifest_or_list)
        suffix = path.suffix.lower()
        if suffix in {".csv", ".tsv"}:
            with path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t" if suffix == ".tsv" else ",")
                candidates = (template_column, "template_id", "template", "Template", "interface")
                values = [next((row[name] for name in candidates if name and name in row), "") for row in reader]
        elif suffix == ".json":
            loaded = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(loaded, Mapping):
                values = next((loaded[name] for name in ("templates", "template_ids", "rows", "data") if name in loaded), [loaded])
            else:
                values = loaded if isinstance(loaded, list) else [loaded]
        else:
            values = [line.split()[0] for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.lstrip().startswith("#")]

    ids: list[str] = []
    for value in values:
        if isinstance(value, Mapping):
            names = (template_column, "template_id", "template", "Template", "interface")
            value = next((value[name] for name in names if name and name in value), "")
        elif isinstance(value, (list, tuple)):
            value = value[0] if value else ""
        template_id = _template_id_from_value(value)
        if template_id:
            ids.append(template_id)
    return ids


def _asset_candidates(root: Path, template_id: str) -> dict[str, list[Path]]:
    contacts = (root / "contacts", root / "contact", root)
    interfaces = (root / "interfaces", root / "interface")
    interface_lists = (root / "interfaces_lists", root / "interface_lists", root / "interfaces")
    contact_txt = [directory / f"{template_id}.txt" for directory in contacts]
    contact_json = [directory / f"{template_id}.json" for directory in contacts]
    interface_json = [directory / f"{template_id}.json" for directory in interface_lists]
    interface_pdb = [directory / f"{template_id}.pdb" for directory in interfaces]
    for directory in interfaces:
        interface_pdb.extend(sorted(directory.glob(f"{template_id}_*_int.pdb")))
    return {
        "legacy_contact_txt": contact_txt,
        "modern_contact_json": contact_json,
        "modern_interface_json": interface_json,
        "modern_interface_pdb": interface_pdb,
    }


def _index_template_assets(root: Path) -> dict[str, dict[str, list[Path]]]:
    """Index the supported asset conventions once per preflight run.

    The historical manifest has 21,072 entries and the legacy asset tree has
    tens of thousands of files.  Repeating ``glob`` for every template turns
    preflight into an accidental quadratic scan, so directory traversal is
    deliberately centralized here.
    """

    contact_dirs = tuple(directory for directory in (root / "contacts", root / "contact", root) if directory.is_dir())
    interface_dirs = tuple(directory for directory in (root / "interfaces", root / "interface") if directory.is_dir())
    interface_list_dirs = tuple(directory for directory in (root / "interfaces_lists", root / "interface_lists", root / "interfaces") if directory.is_dir())
    index: dict[str, dict[str, list[Path]]] = {
        "legacy_contact_txt": {},
        "modern_contact_json": {},
        "modern_interface_json": {},
        "modern_interface_pdb": {},
    }

    def add(kind: str, template_id: str, path: Path) -> None:
        index[kind].setdefault(template_id, []).append(path)

    for directory in contact_dirs:
        for path in directory.iterdir():
            if not path.is_file():
                continue
            if path.suffix.lower() == ".txt":
                add("legacy_contact_txt", path.stem, path)
            elif path.suffix.lower() == ".json":
                add("modern_contact_json", path.stem, path)

    for directory in interface_list_dirs:
        for path in directory.iterdir():
            if path.is_file() and path.suffix.lower() == ".json":
                add("modern_interface_json", path.stem, path)

    for directory in interface_dirs:
        for path in directory.iterdir():
            if not path.is_file() or path.suffix.lower() != ".pdb":
                continue
            name = path.name
            if name.endswith("_int.pdb"):
                template_id = name[: -len("_int.pdb")].rsplit("_", 1)[0]
            else:
                template_id = path.stem
            add("modern_interface_pdb", template_id, path)

    for values in index.values():
        for template_id in values:
            values[template_id].sort(key=str)
    return index


def _asset_record(root: Path, asset_type: str, path: Path) -> dict[str, Any]:
    exists = path.is_file()
    try:
        display = str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        display = str(path.resolve())
    format_valid = False
    validation_error = "missing"
    if exists:
        try:
            if asset_type.endswith("txt"):
                format_valid = bool(path.read_text(encoding="utf-8", errors="replace").strip())
                validation_error = "" if format_valid else "empty text asset"
            elif asset_type.endswith("json"):
                loaded = json.loads(path.read_text(encoding="utf-8"))
                format_valid = isinstance(loaded, (Mapping, list))
                validation_error = "" if format_valid else "JSON root must be object or list"
            elif asset_type.endswith("pdb") or asset_type.startswith("interface_pdb"):
                lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
                format_valid = any(line.startswith(("ATOM", "HETATM")) for line in lines)
                validation_error = "" if format_valid else "no ATOM/HETATM records"
            else:
                format_valid = True
                validation_error = ""
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            validation_error = str(exc)
    return {
        "asset_type": asset_type,
        "asset_path": display,
        "exists": exists,
        "format_valid": format_valid,
        "validation_error": validation_error,
        "sha256": sha256_file(path) if exists else "",
        "size": path.stat().st_size if exists else 0,
    }


def _first_existing(paths: Iterable[Path]) -> Path | None:
    return next((path for path in paths if path.is_file()), None)


def _valid_template_id(template_id: str) -> bool:
    return bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", template_id))


def preflight_template_assets(
    manifest_or_list: str | os.PathLike[str] | Iterable[Any],
    asset_root: str | os.PathLike[str],
    template_column: str | None = None,
) -> dict[str, Any]:
    """Check listed templates against legacy and modern asset conventions.

    ``listed`` counts input rows, ``unique`` counts distinct IDs, ``valid``
    counts recognized/syntactically safe template profiles, and
    ``fully_resolvable`` counts templates for which every required asset in
    the selected profile exists.  ``missing`` counts templates that are not
    fully resolvable.  The returned ``rows`` contain one stable row per asset;
    missing required assets have an empty hash.
    """

    listed_ids = read_template_ids(manifest_or_list, template_column=template_column)
    counts: dict[str, int] = {}
    indices: dict[str, list[int]] = {}
    for index, template_id in enumerate(listed_ids, 1):
        counts[template_id] = counts.get(template_id, 0) + 1
        indices.setdefault(template_id, []).append(index)
    root = Path(asset_root).resolve()
    indexed_assets = _index_template_assets(root)
    templates: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []

    for template_id in sorted(counts):
        legacy_contact = _first_existing(indexed_assets["legacy_contact_txt"].get(template_id, ()))
        modern_contact = _first_existing(indexed_assets["modern_contact_json"].get(template_id, ()))
        modern_interface = _first_existing(indexed_assets["modern_interface_json"].get(template_id, ()))
        modern_pdbs = indexed_assets["modern_interface_pdb"].get(template_id, [])

        if modern_contact or modern_interface or modern_pdbs:
            profile = "modern"
            requirements = {
                "contact_json": modern_contact,
                "interface_json": modern_interface,
            }
            if modern_pdbs:
                # Validate every chain-specific interface PDB.  Checking only
                # the first path can silently admit a template whose other
                # chain file is empty or malformed; GTalign then skips it.
                requirements["interface_pdb"] = modern_pdbs[0]
                for extra_index, path in enumerate(modern_pdbs[1:], 2):
                    requirements[f"interface_pdb_{extra_index}"] = path
            else:
                requirements["interface_pdb"] = None
        elif legacy_contact:
            profile = "legacy"
            requirements = {"contact_txt": legacy_contact}
        else:
            profile = "unknown"
            requirements = {"contact_txt_or_json": None, "interface_json": None, "interface_pdb": None}

        valid = _valid_template_id(template_id) and profile != "unknown"
        missing = sorted(asset_type for asset_type, path in requirements.items() if path is None)
        asset_records: list[dict[str, Any]] = []
        for asset_type, path in requirements.items():
            if path is not None:
                asset_records.append(_asset_record(root, asset_type, path))
            else:
                asset_records.append({"asset_type": asset_type, "asset_path": "", "exists": False, "format_valid": False, "validation_error": "missing", "sha256": "", "size": 0})
        missing.extend(
            f"{asset['asset_type']}_invalid"
            for asset in asset_records
            if asset["exists"] and not asset["format_valid"]
        )
        missing = sorted(set(missing))
        valid = valid and all(asset["format_valid"] for asset in asset_records if asset["exists"])
        fully_resolvable = valid and not missing

        template = {
            "template_id": template_id,
            "manifest_indices": ",".join(str(index) for index in indices[template_id]),
            "listed_count": counts[template_id],
            "duplicate_count": counts[template_id] - 1,
            "unique": len(indices[template_id]) == 1,
            "valid": valid,
            "fully_resolvable": fully_resolvable,
            "missing": missing,
            "profile": profile,
            "assets": asset_records,
        }
        templates.append(template)
        for asset in asset_records:
            rows.append({
                "template_id": template_id,
                "manifest_indices": ",".join(str(index) for index in indices[template_id]),
                "listed_count": counts[template_id],
                "duplicate_count": counts[template_id] - 1,
                "unique": len(indices[template_id]) == 1,
                "valid": int(valid),
                "fully_resolvable": int(fully_resolvable),
                "missing": ";".join(missing),
                **asset,
            })

    rows.sort(key=lambda row: (row["template_id"], row["asset_type"], row["asset_path"]))
    valid_count = sum(int(template["valid"]) for template in templates)
    resolvable_count = sum(int(template["fully_resolvable"]) for template in templates)
    return {
        "listed": len(listed_ids),
        "unique": len(templates),
        "valid": valid_count,
        "fully_resolvable": resolvable_count,
        "missing": len(templates) - resolvable_count,
        "missing_assets": sum(len(template["missing"]) for template in templates),
        "templates": templates,
        "rows": rows,
    }


def write_template_preflight_tsv(
    output_path: str | os.PathLike[str],
    report: Mapping[str, Any],
) -> Path:
    """Write stable per-asset preflight rows as TSV and return the path."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = report.get("rows", [])
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in TSV_FIELDS})
    return output


# Descriptive aliases make the helper easy to discover from existing scripts.
build_template_asset_preflight = preflight_template_assets
template_asset_preflight = preflight_template_assets


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--output", required=True, help="Output provenance JSON path")
    parser.add_argument("--file", action="append", default=[], dest="files")
    parser.add_argument("--executable", action="append", default=[], dest="executables")
    parser.add_argument("--package", action="append", default=[], dest="packages")
    parser.add_argument("--env", action="append", default=[], dest="environment")
    parser.add_argument("--seed", action="append", default=[], dest="seeds", help="NAME=VALUE seed pair")
    parser.add_argument("--config", help="JSON/INI/key=value effective config file")
    parser.add_argument("--template-manifest")
    parser.add_argument("--asset-root")
    parser.add_argument("--template-tsv")
    args = parser.parse_args(list(argv) if argv is not None else None)

    seed_values = dict(item.split("=", 1) for item in args.seeds if "=" in item)
    manifest = build_provenance_manifest(
        args.repo_root,
        files=args.files,
        executables=args.executables,
        environment=args.environment or None,
        packages=args.packages,
        seeds=seed_values or None,
        config=args.config,
    )
    write_provenance_manifest(args.output, manifest)
    if args.template_manifest and args.asset_root:
        report = preflight_template_assets(args.template_manifest, args.asset_root)
        write_template_preflight_tsv(args.template_tsv or str(Path(args.output).with_suffix(".templates.tsv")), report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
