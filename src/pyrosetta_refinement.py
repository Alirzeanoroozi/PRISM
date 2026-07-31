"""Opt-in PyRosetta refinement adapter and runtime probe.

PyRosetta is intentionally not a project dependency.  This module does not
import it at module load time; callers must explicitly invoke the probe or
the adapter.  An unavailable or unusable PyRosetta runtime is reported as an
explicit result and is never replaced with the external Rosetta backend.
"""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import json
import os
import platform
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable


PACKAGE_NAME = "pyrosetta"
BACKEND_NAME = "pyrosetta"
SAFE_ENVIRONMENT_KEYS = (
    "CONDA_DEFAULT_ENV",
    "CONDA_PREFIX",
    "CUDA_VISIBLE_DEVICES",
    "LD_LIBRARY_PATH",
    "PATH",
    "PYTHONPATH",
    "SLURM_JOB_ID",
    "SLURM_JOB_NAME",
)


def _command_metadata() -> dict[str, Any]:
    return {
        "argv": list(sys.argv),
        "executable": sys.executable,
        "python_version": sys.version,
    }


def _environment_metadata() -> dict[str, Any]:
    return {
        "cwd": os.getcwd(),
        "platform": platform.platform(),
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "environment": {
            key: os.environ[key]
            for key in SAFE_ENVIRONMENT_KEYS
            if key in os.environ
        },
    }


def _installed_version() -> str | None:
    for distribution in ("PyRosetta", PACKAGE_NAME):
        try:
            return importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            continue
    return None


def _module_version(module: Any) -> str | None:
    value = getattr(module, "__version__", None)
    if value is not None:
        return str(value)
    installed = _installed_version()
    if installed is not None:
        return installed
    version_function = getattr(module, "version", None)
    if callable(version_function):
        try:
            return str(version_function())
        except Exception:
            return None
    return _installed_version()


def _base_report() -> dict[str, Any]:
    return {
        "backend": BACKEND_NAME,
        "package": PACKAGE_NAME,
        "version": None,
        "status": "unavailable",
        "available": False,
        "import_error": None,
        "module_file": None,
        "input_hashes": {},
        "output_hashes": {},
        "command_metadata": _command_metadata(),
        "environment_metadata": _environment_metadata(),
    }


def _import_report() -> tuple[Any | None, dict[str, Any]]:
    report = _base_report()
    try:
        module = importlib.import_module(PACKAGE_NAME)
    except Exception as exc:
        report["import_error"] = f"{type(exc).__name__}: {exc}"
        return None, report

    report.update(
        {
            "version": _module_version(module),
            "status": "available",
            "available": True,
            "module_file": str(getattr(module, "__file__", "")) or None,
        }
    )
    return module, report


def probe_environment() -> dict[str, Any]:
    """Probe the optional PyRosetta runtime without raising on import failure."""

    _, report = _import_report()
    return report


probe_pyrosetta_environment = probe_environment


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_hashes(paths: list[Path]) -> dict[str, str]:
    return {str(path): _sha256_file(path) for path in paths if path.is_file()}


def _metadata_path(output_path: Path) -> Path:
    return Path(f"{output_path}.pyrosetta.json")


def _write_metadata(path: Path, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")


class PyRosettaRefinementAdapter:
    """Explicit PyRosetta-only refinement adapter.

    ``input_pdb`` must already contain the intended receptor/ligand assembly;
    ``partners`` is passed to the PyRosetta docking protocol.  The adapter
    never invokes or imports ``src.rosetta_refinement``.
    """

    def __init__(
        self,
        *,
        init_options: str = "-mute all",
        protocol_factory: Callable[..., Any] | None = None,
    ) -> None:
        self.init_options = init_options
        self.protocol_factory = protocol_factory
        self._initialized_module_ids: set[int] = set()

    def probe(self) -> dict[str, Any]:
        return probe_environment()

    def refine(
        self,
        input_pdb: str | os.PathLike[str],
        output_pdb: str | os.PathLike[str],
        *,
        partners: str = "A_B",
    ) -> dict[str, Any]:
        input_path = Path(input_pdb).resolve()
        output_path = Path(output_pdb).resolve()
        metadata_path = _metadata_path(output_path)
        module, runtime = _import_report()
        record: dict[str, Any] = {
            **runtime,
            "input_path": str(input_path),
            "output_path": str(output_path),
            "metadata_path": str(metadata_path),
            "input_hashes": {},
            "output_hashes": {},
            "parameters": {
                "init_options": self.init_options,
                "partners": partners,
            },
            "fallback": None,
            "total_score": None,
        }

        if input_path.is_file():
            record["input_hashes"] = _file_hashes([input_path])
        else:
            record.update(
                {
                    "status": "failed",
                    "available": bool(module),
                    "error": f"input PDB does not exist: {input_path}",
                }
            )
            _write_metadata(metadata_path, record)
            return record

        if output_path.exists():
            record.update({
                "status": "failed",
                "available": bool(module),
                "error": f"refinement output already exists: {output_path}",
            })
            _write_metadata(metadata_path, record)
            return record

        if module is None:
            # This is deliberately terminal.  No external Rosetta command is
            # attempted when the optional backend is unavailable.
            _write_metadata(metadata_path, record)
            return record

        try:
            module_id = id(module)
            if module_id not in self._initialized_module_ids:
                module.init(self.init_options)
                self._initialized_module_ids.add(module_id)

            pose = module.pose_from_pdb(str(input_path))
            scorefxn = module.get_fa_scorefxn()
            protocol = self._make_protocol(module, pose, scorefxn, partners)
            protocol.apply(pose)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            temporary = tempfile.NamedTemporaryFile(
                mode="w", suffix=".pdb", prefix=f".{output_path.name}.",
                dir=output_path.parent, delete=False,
            )
            temporary_path = Path(temporary.name)
            temporary.close()
            pose.dump_pdb(str(temporary_path))
            if not temporary_path.is_file() or temporary_path.stat().st_size == 0:
                raise RuntimeError(f"PyRosetta did not create a non-empty temporary PDB: {temporary_path}")
            os.replace(temporary_path, output_path)
            record["total_score"] = float(scorefxn(pose))
            record["status"] = "success"
            record["output_hashes"] = _file_hashes([output_path])
        except Exception as exc:
            if "temporary_path" in locals():
                temporary_path.unlink(missing_ok=True)
            record.update(
                {
                    "status": "failed",
                    "error": f"{type(exc).__name__}: {exc}",
                    "output_hashes": {},
                }
            )

        _write_metadata(metadata_path, record)
        return record

    def _make_protocol(
        self,
        module: Any,
        pose: Any,
        scorefxn: Any,
        partners: str,
    ) -> Any:
        if self.protocol_factory is not None:
            return self.protocol_factory(module, pose, scorefxn, partners)

        protocol = module.rosetta.protocols.docking.DockingProtocol()
        try:
            protocol.set_docking_local_refine(True)
        except TypeError:
            # Older/test-double APIs expose this as a no-argument setter.
            protocol.set_docking_local_refine()
        protocol.set_partners(partners)
        if hasattr(protocol, "set_scorefxn"):
            protocol.set_scorefxn(scorefxn)
        elif hasattr(protocol, "set_highres_scorefxn"):
            protocol.set_highres_scorefxn(scorefxn)
        else:
            raise AttributeError("DockingProtocol exposes no compatible score-function setter")
        return protocol


def refine(
    input_pdb: str | os.PathLike[str],
    output_pdb: str | os.PathLike[str],
    *,
    partners: str = "A_B",
    init_options: str = "-mute all",
) -> dict[str, Any]:
    """Refine one explicit PDB through PyRosetta only."""

    return PyRosettaRefinementAdapter(init_options=init_options).refine(
        input_pdb,
        output_pdb,
        partners=partners,
    )


def refine_pairs(
    passed_pairs: list[tuple[str | os.PathLike[str], str | os.PathLike[str]]],
    *,
    output_root: str | os.PathLike[str] = "processed/pyrosetta_refinement",
    init_options: str | None = None,
) -> list[dict[str, Any]]:
    """Refine transformed receptor/ligand pairs through PyRosetta only.

    The input pair is combined using the same chain-preserving helper used by
    the external Rosetta backend.  Each result is written to an isolated
    structure directory and retains the adapter's JSON provenance record.
    """
    from .rosetta_refinement import combine_pdb, partner_chain_ids

    root = Path(output_root).resolve()
    structures = root / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    adapter = PyRosettaRefinementAdapter(init_options=init_options or "-mute all -constant_seed -jran 12345")
    records: list[dict[str, Any]] = []
    for left, right in passed_pairs:
        left_path, right_path = Path(left), Path(right)
        combined = combine_pdb(str(left_path), str(right_path))
        if not combined:
            records.append({"status": "failed", "error": "combine_pdb_failed", "left": str(left_path), "right": str(right_path)})
            continue
        left_chains, right_chains = partner_chain_ids(str(left_path), str(right_path))
        stem = Path(combined).stem
        output = structures / f"{stem}.pdb"
        record = adapter.refine(combined, output, partners=f"{left_chains}_{right_chains}")
        record.update({"left": str(left_path), "right": str(right_path), "combined_input": str(combined), "partners": f"{left_chains}_{right_chains}"})
        records.append(record)
    (root / "refinement_results.json").write_text(json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return records


def refine_merged_candidates(
    candidates: list[tuple[str, str, str, str]],
    *,
    output_root: str | os.PathLike[str] = "processed/pyrosetta_refinement",
    init_options: str = "-mute all -constant_seed -jran 12345",
) -> list[tuple[str, str, str, str]]:
    """Refine PRISM merged candidates without changing their tuple contract."""
    from .pdb_download import split_target_id

    root = Path(output_root).resolve()
    structures = root / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    adapter = PyRosettaRefinementAdapter(init_options=init_options)
    accepted = []
    for receptor, ligand, template, input_pdb in candidates:
        _, receptor_chains = split_target_id(receptor)
        _, ligand_chains = split_target_id(ligand)
        output = structures / f"{Path(input_pdb).stem}_pyrosetta.pdb"
        record = adapter.refine(
            input_pdb,
            output,
            partners=f"{''.join(receptor_chains)}_{''.join(ligand_chains)}",
        )
        if record.get("status") == "success":
            accepted.append((receptor, ligand, template, str(output)))
    return accepted
