#!/usr/bin/env python3
"""Compare PRISM output PDB(s) to downloaded native inputs and score with DockQ."""

import argparse
import json
import os
import sys

from pathlib import Path

from src.compare import SUMMARY_CSV, compare_and_summarize, compare_pair


def parse_output_filename(pdb_path: Path):
    """Parse ``{template}_{receptor}_{ligand}.pdb`` written by transformation."""
    stem = pdb_path.stem
    parts = stem.rsplit("_", 2)
    if len(parts) != 3:
        raise ValueError(f"cannot parse output filename: {pdb_path.name}")
    return parts[0], parts[1], parts[2]


def pairs_from_output_dir(output_dir: str):
    """Build compare tuples from every ``*.pdb`` in *output_dir*."""
    root = Path(output_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"output directory not found: {output_dir}")
    pairs = []
    skipped = []
    for pdb_path in sorted(root.glob("*.pdb")):
        try:
            template, receptor, ligand = parse_output_filename(pdb_path)
        except ValueError as exc:
            skipped.append((pdb_path.name, str(exc)))
            continue
        pairs.append((receptor, ligand, template, str(pdb_path)))
    return pairs, skipped


def main():
    parser = argparse.ArgumentParser(
        description="Compare PRISM model PDB(s) to native processed/pdbs and compute DockQ"
    )
    parser.add_argument(
        "--model",
        help="Single output PDB path (use with --receptor, --ligand, --template)",
    )
    parser.add_argument("--receptor", help="Receptor target id, e.g. 3i6eE")
    parser.add_argument("--ligand", help="Ligand target id, e.g. 3i6eF")
    parser.add_argument("--template", help="Template id, e.g. 3broAD")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Compare every *.pdb in this folder (filenames: template_receptor_ligand.pdb)",
    )
    parser.add_argument(
        "--pairs-csv",
        default=None,
        help="CSV with columns receptor,ligand,template,output_pdb (batch mode)",
    )
    parser.add_argument(
        "--summary-csv",
        default=SUMMARY_CSV,
        help=f"Output summary CSV (default: {SUMMARY_CSV})",
    )
    parser.add_argument(
        "--dockq-no-align",
        action="store_true",
        help="Pass --no_align to DockQ (default: align, recommended for PRISM outputs)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel workers for batch compare (default: 1)",
    )
    parser.add_argument("--json", action="store_true", help="Print JSON to stdout")
    args = parser.parse_args()

    if args.model:
        if not (args.receptor and args.ligand and args.template):
            parser.error("--model requires --receptor, --ligand, and --template")
        row = compare_pair(
            args.receptor,
            args.ligand,
            args.template,
            args.model,
            dockq_no_align=args.dockq_no_align,
        )
        rows = [row]
        os.makedirs(os.path.dirname(args.summary_csv) or ".", exist_ok=True)
        import csv

        fieldnames = list(row.keys())
        with open(args.summary_csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    elif args.output_dir:
        pairs, skipped = pairs_from_output_dir(args.output_dir)
        if skipped:
            print(f"skipped {len(skipped)} file(s) with unparsable names", file=sys.stderr)
            for name, err in skipped[:5]:
                print(f"  {name}: {err}", file=sys.stderr)
        if not pairs:
            raise SystemExit(f"no .pdb files found in {args.output_dir}")
        print(f"comparing {len(pairs)} structure(s) from {args.output_dir}")
        summary_path, rows = compare_and_summarize(
            pairs,
            summary_csv=args.summary_csv,
            dockq_no_align=args.dockq_no_align,
            n_jobs=args.jobs,
        )
        args.summary_csv = summary_path
    elif args.pairs_csv:
        import csv as csvmod

        pairs = []
        with open(args.pairs_csv, newline="") as fh:
            for row in csvmod.DictReader(fh):
                pairs.append(
                    (
                        row["receptor"].strip(),
                        row["ligand"].strip(),
                        row["template"].strip(),
                        row["output_pdb"].strip(),
                    )
                )
        summary_path, rows = compare_and_summarize(
            pairs,
            summary_csv=args.summary_csv,
            dockq_no_align=args.dockq_no_align,
            n_jobs=args.jobs,
        )
        args.summary_csv = summary_path
    else:
        parser.error("Provide --model, --output-dir, or --pairs-csv")

    if args.json:
        print(json.dumps(rows, indent=2))
    else:
        for row in rows:
            print(f"{row['template']}_{row['receptor']}_{row['ligand']}")
            print(f"  output: {row['output_pdb']}")
            print(f"  native: {row['native_pdb']}")
            print(
                f"  receptor CA RMSD: {row['receptor_ca_rmsd']} "
                f"({row['receptor_ca_matched']} atoms)"
            )
            print(f"  iRMSD (backbone): {row.get('irmsd_backbone')}")
            print(f"  DockQ: {row['dockq']} ({row['dockq_capri']})")
            if row.get("error_irmsd_backbone"):
                print(f"  iRMSD error: {row['error_irmsd_backbone']}")
            if row["error_dockq"]:
                print(f"  DockQ error: {row['error_dockq']}")
        print(f"summary -> {args.summary_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
