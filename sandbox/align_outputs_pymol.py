#!/usr/bin/env python3
"""Align chain E across PRISM output PDBs with PyMOL and write to a new folder.

Superimposes each structure's chain E onto a reference chain E (default: chain E
from ``processed/pdbs/3i6e.pdb``). The whole model (including chain F) is moved
with the alignment and saved.

Run headless::

    conda activate new_bg   # or any env with PyMOL
    pymol -cq align_outputs_pymol.py

Or with options::

    pymol -cq align_outputs_pymol.py -- --chain E --jobs 1

Requires PyMOL (``pymol`` on PATH or ``pip install pymol``).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def _launch_pymol():
    try:
        import pymol
        from pymol import cmd
    except ImportError as exc:
        raise SystemExit(
            "PyMOL is not available. Run this script with:\n"
            "  pymol -cq align_outputs_pymol.py\n"
            f"Import error: {exc}"
        ) from exc
    pymol.finish_launching(["pymol", "-qc"])
    return cmd


def align_outputs(
    input_dir: Path,
    output_dir: Path,
    reference_pdb: Path,
    chain: str = "E",
    object_name: str = "mobile",
    ref_name: str = "ref",
) -> tuple[int, int]:
    from pymol import cmd

    if not reference_pdb.is_file():
        raise FileNotFoundError(f"Reference PDB not found: {reference_pdb}")
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    pdbs = sorted(input_dir.glob("*.pdb"))
    if not pdbs:
        raise SystemExit(f"No .pdb files in {input_dir}")

    cmd.reinitialize()
    cmd.set("retain_order", 1)
    cmd.load(str(reference_pdb), ref_name)
    ref_chains = cmd.get_chains(ref_name)
    if chain not in ref_chains:
        raise ValueError(
            f"Chain {chain!r} not found in reference {reference_pdb} (chains: {ref_chains})"
        )

    ok = 0
    failed = 0
    ref_sel = f"{ref_name} and chain {chain}"
    total = len(pdbs)

    for i, pdb_path in enumerate(pdbs, start=1):
        out_path = output_dir / pdb_path.name
        try:
            cmd.delete(object_name)
        except Exception:
            pass
        try:
            cmd.load(str(pdb_path), object_name)
            chains = cmd.get_chains(object_name)
            if chain not in chains:
                print(f"  skip {pdb_path.name}: chain {chain} not in {chains}", flush=True)
                failed += 1
                cmd.delete(object_name)
                continue
            mobile_sel = f"{object_name} and chain {chain}"
            result = cmd.align(mobile_sel, ref_sel)
            rmsd = result[0] if result else None
            cmd.save(str(out_path), object_name)
            cmd.delete(object_name)
            ok += 1
            if i == 1 or i % 50 == 0 or i == total:
                rmsd_s = f"{rmsd:.3f}" if rmsd is not None else "NA"
                print(f"  aligned {i}/{total}: {pdb_path.name} (RMSD={rmsd_s})", flush=True)
        except Exception as exc:
            print(f"  FAIL {pdb_path.name}: {exc}", flush=True)
            failed += 1
            try:
                cmd.delete(object_name)
            except Exception:
                pass

    cmd.delete(ref_name)
    return ok, failed


def main(argv=None):
    repo = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Align chain E in all output PDBs to a reference using PyMOL"
    )
    parser.add_argument(
        "--input-dir",
        default=str(repo / "processed" / "output"),
        help="Folder of model PDBs (default: processed/output)",
    )
    parser.add_argument(
        "--output-dir",
        default=str(repo / "processed" / "aligned_output"),
        help="Folder for aligned PDBs (default: processed/aligned_output)",
    )
    parser.add_argument(
        "--reference",
        default=str(repo / "processed" / "pdbs" / "3i6e.pdb"),
        help="Reference PDB for chain E alignment",
    )
    parser.add_argument(
        "--chain",
        default="E",
        help="Chain id to align (default: E)",
    )
    args = parser.parse_args(argv)

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    reference_pdb = Path(args.reference).resolve()

    print(f"input:     {input_dir}")
    print(f"output:    {output_dir}")
    print(f"reference: {reference_pdb}")
    print(f"chain:     {args.chain}", flush=True)

    _launch_pymol()
    ok, failed = align_outputs(
        input_dir=input_dir,
        output_dir=output_dir,
        reference_pdb=reference_pdb,
        chain=args.chain.strip(),
    )
    print(f"done: {ok} aligned, {failed} failed/skipped -> {output_dir}", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    # pymol -cq script.py -- --extra args
    extra = sys.argv[1:]
    if extra and extra[0] == "--":
        extra = extra[1:]
    sys.exit(main(extra))
