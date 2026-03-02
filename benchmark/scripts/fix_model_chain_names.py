#!/usr/bin/env python3
"""
Fix PRISM model PDB chain IDs when both partners collapse to one chain label.

Workflow:
1) Scan PRISM model PDB files under rosetta output root.
2) Detect models with a single observed chain ID but >=2 chain segments (TER-separated).
3) Copy those models to a fixed folder and rewrite chain IDs by partner segment:
   - segment 0 -> tpl1 chain (from model filename; fallback A)
   - segment 1 -> tpl2 chain (from model filename; fallback B/C)
   - if tpl1==tpl2, force tpl2 to another chain ID.
4) Save CSV mapping for downstream analysis.
"""

import argparse
import csv
import re
from pathlib import Path


ROSETTA_PDB_RE = re.compile(
    r"^(?P<target>[A-Za-z0-9]+)_(?P<tpl1>[A-Za-z0-9]+)_(?P<idx1>\d+)_(?P<tpl2>[A-Za-z0-9]+)_(?P<idx2>\d+)\.(?P<rest>rosetta.*)\.pdb$"
)


def parse_tpl(token):
    token = (token or "").strip()
    if "_" in token:
        pdb_part, chains = token.split("_", 1)
    else:
        pdb_part, chains = token[:4], token[4:]
    pdb4 = pdb_part[:4].lower()
    chains = "".join(ch for ch in chains if ch.isalnum())
    return pdb4, chains


def choose_chain2(chain1):
    for c in "BCDEFGHIJKLMNOPQRSTUVWXYZ0123456789":
        if c != chain1:
            return c
    return "Z"


def scan_chain_info(pdb_path):
    chains = []
    chain_set = set()
    seg_count = 1
    has_atom = False
    with pdb_path.open() as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
                has_atom = True
                c = line[21]
                if c not in chain_set:
                    chain_set.add(c)
                    chains.append(c)
            elif line.startswith("TER"):
                seg_count += 1
    if not has_atom:
        seg_count = 0
    return chains, seg_count


def rewrite_chains(src, dst, chain_seg0, chain_seg1):
    dst.parent.mkdir(parents=True, exist_ok=True)
    seg = 0
    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith(("ATOM  ", "HETATM", "ANISOU")) and len(line) > 21:
                out_chain = chain_seg0 if seg == 0 else chain_seg1
                fout.write(line[:21] + out_chain + line[22:])
            elif line.startswith("TER"):
                out_chain = chain_seg0 if seg == 0 else chain_seg1
                if len(line) > 21:
                    fout.write(line[:21] + out_chain + line[22:])
                else:
                    fout.write(line)
                seg += 1
            else:
                fout.write(line)


def main():
    parser = argparse.ArgumentParser(description="Rewrite problematic PRISM model chain IDs and save chain-map CSV")
    parser.add_argument(
        "--rosetta-root",
        default="benchmark/prism_processed/rosetta_output_1",
        help="Root folder with PRISM model PDB files",
    )
    parser.add_argument(
        "--fixed-root",
        default="benchmark/prism_processed/rosetta_output_1_chainfixed",
        help="Output folder for fixed model PDB files",
    )
    parser.add_argument(
        "--map-csv",
        default="benchmark/prism_processed_results/chain_fix_map.csv",
        help="Output CSV mapping used by analyzer",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process first N models (debug)",
    )
    args = parser.parse_args()

    rosetta_root = Path(args.rosetta_root)
    fixed_root = Path(args.fixed_root)
    map_csv = Path(args.map_csv)

    models = sorted([p for p in rosetta_root.rglob("*.pdb") if not p.name.endswith(".pdb.intRes.txt")])
    if args.limit is not None:
        models = models[: args.limit]

    rows = []
    fixed_count = 0
    skipped = 0
    for p in models:
        m = ROSETTA_PDB_RE.match(p.name)
        if not m:
            skipped += 1
            continue
        _, tpl1_chains = parse_tpl(m.group("tpl1"))
        _, tpl2_chains = parse_tpl(m.group("tpl2"))

        c1 = tpl1_chains[:1] if tpl1_chains else "A"
        c2 = tpl2_chains[:1] if tpl2_chains else choose_chain2(c1)
        if c1 == c2:
            c2 = choose_chain2(c1)

        chains_seen, seg_count = scan_chain_info(p)
        rel = p.relative_to(rosetta_root)
        out_p = fixed_root / rel

        status = "unchanged"
        note = ""
        if seg_count >= 2 and len(chains_seen) == 1:
            rewrite_chains(p, out_p, c1, c2)
            status = "fixed_single_chain_multi_segment"
            note = "rewrote chains using TER segments"
            fixed_count += 1
        else:
            # Keep original path for unchanged rows; no file copy needed.
            note = "no rewrite needed"

        rows.append({
            "original_model_path": str(p),
            "fixed_model_path": str(out_p if status.startswith("fixed_") else p),
            "fixed_tpl1_chain": c1,
            "fixed_tpl2_chain": c2,
            "observed_chains": "".join(chains_seen),
            "segment_count": seg_count,
            "status": status,
            "note": note,
        })

    map_csv.parent.mkdir(parents=True, exist_ok=True)
    with map_csv.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "original_model_path",
                "fixed_model_path",
                "fixed_tpl1_chain",
                "fixed_tpl2_chain",
                "observed_chains",
                "segment_count",
                "status",
                "note",
            ],
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("[INFO] models_scanned =", len(models))
    print("[INFO] models_fixed =", fixed_count)
    print("[INFO] models_skipped_nonpattern =", skipped)
    print("[INFO] map_csv =", map_csv)


if __name__ == "__main__":
    raise SystemExit(main())

