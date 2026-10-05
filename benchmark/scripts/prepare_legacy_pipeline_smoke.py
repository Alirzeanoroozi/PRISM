#!/usr/bin/env python3
"""Prepare an isolated one-pair/one-template legacy PRISM smoke workspace."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def copy_file(source: Path, target: Path, records: list[dict[str, str]], role: str) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    records.append({"role": role, "source": str(source), "target": str(target), "sha256": sha256(target)})


def prepare(
    output_root: Path,
    compat_root: Path,
    tool_env: Path,
    receptor: Path,
    ligand: Path,
    template_root: Path,
    template_id: str = "1kcaCH",
    logical_receptor_selector: str | None = None,
    logical_ligand_selector: str | None = None,
    template_chains: tuple[str, str] | None = None,
) -> dict[str, object]:
    output = output_root.resolve()
    compat = compat_root.resolve()
    tools = tool_env.resolve()
    receptor = receptor.resolve()
    ligand = ligand.resolve()
    templates = template_root.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite non-empty workspace: {output}")
    for path in (compat, tools, receptor, ligand, templates):
        if not path.exists():
            raise ValueError(f"missing smoke input: {path}")
    output.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, str]] = []

    for source in sorted(compat.iterdir()):
        if source.name in {"external_tools", "compatibility_files.tsv"}:
            continue
        target = output / source.name
        if source.is_dir():
            shutil.copytree(source, target)
            for path in sorted(target.rglob("*")):
                if path.is_file():
                    records.append({"role": "compatibility_snapshot", "source": str(source / path.relative_to(target)), "target": str(path), "sha256": sha256(path)})
        else:
            copy_file(source, target, records, "compatibility_snapshot")

    if (compat / "compatibility_manifest.json").is_file():
        copy_file(compat / "compatibility_manifest.json", output / "compatibility_manifest.json", records, "compatibility_snapshot_manifest")

    shutil.copytree(tools / "external_tools", output / "external_tools")
    for path in sorted((output / "external_tools").rglob("*")):
        if path.is_file():
            records.append({"role": "staged_external_tool", "source": str(tools / "external_tools" / path.relative_to(output / "external_tools")), "target": str(path), "sha256": sha256(path)})

    if template_chains is None:
        chain_suffix = template_id[4:]
        if len(chain_suffix) != 2:
            raise ValueError(f"template_id must contain exactly two chain IDs: {template_id}")
        template_chains = (chain_suffix[0], chain_suffix[1])
    if len(template_chains) != 2 or any(not chain for chain in template_chains):
        raise ValueError(f"template_chains must contain two non-empty chain IDs: {template_chains}")
    selected_template_files = (
        ("interfaces", f"{template_id}_{template_chains[0]}.int"),
        ("interfaces", f"{template_id}_{template_chains[1]}.int"),
        ("contact", f"{template_id}.txt"),
        ("hotspot", f"hotspot{template_id}"),
    )
    for directory, name in selected_template_files:
        source = templates / directory / name
        if not source.is_file():
            raise ValueError(f"missing template asset: {source}")
        copy_file(source, output / "template" / directory / name, records, "selected_template_asset")

    copy_file(receptor, output / "jobs" / "smoke" / "pdb" / "pdb1.pdb", records, "benchmark_receptor_staged_as_pdb1")
    copy_file(ligand, output / "jobs" / "smoke" / "pdb" / "pdb2.pdb", records, "benchmark_ligand_staged_as_pdb2")
    copy_file(tools / "external_tools" / "multiprot" / "params.txt", output / "run_files" / "params.txt", records, "multiprot_parameters_for_legacy_controller")
    pair_list = output / "input" / "pair_list"
    template_list = output / "input" / "template_list"
    pair_list.parent.mkdir(parents=True, exist_ok=True)
    pair_list.write_text("pdb1 pdb2\n", encoding="utf-8")
    template_list.write_text(template_id + "\n", encoding="utf-8")
    records.append({"role": "smoke_pair_list", "source": "generated:pdb1 pdb2", "target": str(pair_list), "sha256": sha256(pair_list)})
    records.append({"role": "smoke_template_list", "source": f"generated:{template_id}", "target": str(template_list), "sha256": sha256(template_list)})

    job = output / "jobs" / "smoke"
    (job / "template_default").write_text(template_id + "\n", encoding="utf-8")
    (output / "config.inc").write_text("# local compatibility configuration\n'localhost'\n'local'\n'local'\n'local'\n", encoding="utf-8")
    records.append({"role": "local_database_boundary_config", "source": "generated:fake-local-values", "target": str(output / "config.inc"), "sha256": sha256(output / "config.inc")})

    manifest = {
        "schema_version": "legacy-prism-smoke-workspace/v1",
        "workspace_root": str(output),
        "compatibility_root": str(compat),
        "tool_environment_root": str(tools),
        "pair_identifiers": {
            "receptor": "pdb1",
            "ligand": "pdb2",
            "logical_receptor": logical_receptor_selector or receptor.stem,
            "logical_ligand": logical_ligand_selector or ligand.stem,
        },
        "benchmark_inputs": {"receptor": str(receptor), "ligand": str(ligand), "receptor_sha256": sha256(receptor), "ligand_sha256": sha256(ligand)},
        "template_id": template_id,
        "template_chains": list(template_chains),
        "template_assets": [str(output / directory / name) for directory, name in selected_template_files],
        "command": [str(tools / "bin" / "python2"), "prism.py", "input/pair_list", "input/template_list", "smoke"],
        "scope": "plumbing_smoke_only; pair IDs are synthetic aliases for local legacy filter paths",
        "records": records,
    }
    (output / "smoke_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--compat-root", type=Path, required=True)
    parser.add_argument("--tool-env", type=Path, required=True)
    parser.add_argument("--receptor", type=Path, required=True)
    parser.add_argument("--ligand", type=Path, required=True)
    parser.add_argument("--template-root", type=Path, required=True)
    parser.add_argument("--template-id", default="1kcaCH")
    parser.add_argument("--logical-receptor-selector")
    parser.add_argument("--logical-ligand-selector")
    parser.add_argument("--template-chains", help="Two template chain IDs, e.g. AD; defaults to the suffix of template_id")
    args = parser.parse_args(argv)
    try:
        value = prepare(
            args.output_root,
            args.compat_root,
            args.tool_env,
            args.receptor,
            args.ligand,
            args.template_root,
            args.template_id,
            args.logical_receptor_selector,
            args.logical_ligand_selector,
            tuple(args.template_chains) if args.template_chains else None,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps({"workspace_root": value["workspace_root"], "template_id": value["template_id"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
