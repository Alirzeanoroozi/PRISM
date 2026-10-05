#!/usr/bin/env python3
"""Stage an isolated FiberDock probe with explicitly distinct input chains.

The source workspace and benchmark PDBs are never modified. The receptor
chain is rewritten only in the derived probe workspace so that the diagnostic
tests whether FiberDock preserves chain identity at its input/output
boundary. This is not a benchmark repair and is never provenance-equivalent
to the original B/B selector pair.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from benchmark.scripts.multiprot_compat_adapter import RUNTIME as COMPAT_RUNTIME


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rewrite_chain(source, destination, chain):
    source = Path(source)
    destination = Path(destination)
    if len(chain) != 1:
        raise ValueError("PDB chain IDs must be exactly one character")
    lines = source.read_text(encoding="ascii", errors="replace").splitlines(True)
    output = []
    observed = set()
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 22:
            observed.add(line[21].strip() or "_")
            output.append(line[:21] + chain + line[22:])
        else:
            output.append(line)
    destination.write_text("".join(output))
    return {
        "source": str(source.resolve()),
        "source_sha256": sha256(source),
        "destination": str(destination.resolve()),
        "emitted_chain": chain,
        "source_chains": sorted(observed),
        "destination_sha256": sha256(destination),
    }


def prepare(source_workspace, destination, receptor_pdb, ligand_pdb, template="1b27AD"):
    source_workspace = Path(source_workspace).resolve()
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(f"refusing to reuse existing probe workspace: {destination}")
    shutil.copytree(source_workspace, destination, symlinks=True)
    for relative in ("jobs", "fiberdock_output"):
        path = destination / relative
        if path.exists():
            shutil.rmtree(path)
    pdb_dir = destination / "pdb"
    pdb_dir.mkdir(parents=True)
    # Keep this diagnostic aligned with the repository's current compatibility
    # boundary. The controller's job-local path is valid, but the template
    # manifest is intentionally stored at the workspace root; the adapter
    # resolves that parent-level manifest. Record the effective runtime hash
    # below so the probe cannot be mistaken for the older snapshot.
    compat_runtime = destination / "run_files" / "compat_runtime.py"
    if compat_runtime.is_file():
        compat_runtime.write_text(COMPAT_RUNTIME, encoding="utf-8")
    receptor_record = rewrite_chain(receptor_pdb, pdb_dir / "1rgh.pdb", "A")
    ligand_record = rewrite_chain(ligand_pdb, pdb_dir / "1a19.pdb", "B")
    (destination / "input" / "pair_list").write_text("1rghA 1a19B\n")
    (destination / "input" / "template_list").write_text(template + "\n")
    # The compatibility TemplateChecker reads the workspace-level default
    # list, while the controller reads input/template_list. Keep both lists
    # explicit and identical for this one-template diagnostic.
    (destination / "template_default").write_text(template + "\n")
    (destination / "template_list").write_text(template + "\n")

    manifest = {
        "probe": "distinct_input_chains",
        "source_workspace": str(source_workspace),
        "source_workspace_manifest_sha256": sha256(source_workspace / "smoke_manifest.json"),
        "logical_receptor": "1RGH_A",
        "logical_ligand": "1A19_B",
        "template": template,
        "receptor": receptor_record,
        "ligand": ligand_record,
        "input_pair_sha256": sha256(destination / "input" / "pair_list"),
        "input_template_sha256": sha256(destination / "input" / "template_list"),
        "workspace_template_default_sha256": sha256(destination / "template_default"),
        "workspace_template_list_sha256": sha256(destination / "template_list"),
        "compat_runtime_sha256": sha256(compat_runtime) if compat_runtime.is_file() else "",
        "status": "staged_exploratory_distinct_chain_probe",
    }
    (destination / "distinct_chain_probe_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_workspace")
    parser.add_argument("destination")
    parser.add_argument("receptor_pdb")
    parser.add_argument("ligand_pdb")
    parser.add_argument("--template", default="1b27AD")
    args = parser.parse_args()
    print(json.dumps(prepare(args.source_workspace, args.destination, args.receptor_pdb, args.ligand_pdb, args.template), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
