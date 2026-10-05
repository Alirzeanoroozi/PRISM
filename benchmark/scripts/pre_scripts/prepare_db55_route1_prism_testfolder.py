#!/usr/bin/env python3
"""
Prepare PRISM test-folder inputs for Docking Benchmark tables by merging
UNBOUND receptor + UNBOUND ligand into ONE PDB per case.

Requirements implemented:
1) File lookup uses the 4-letter PDB code from the **Complex** column (NOT PDB ID 1/2).
   - Expects unbound files named:
       <BOUND4>_r_u.pdb   and   <BOUND4>_l_u.pdb
     inside --pdb-dir.
2) Chain IDs must be written inside the PDB ATOM/HETATM records (column 22), not just filenames.
3) If an unbound file has blank chain IDs (chain column is ' '):
   - Fill blanks ONLY when the intended chain for that partner is a SINGLE letter
     (taken from the Complex column partner side).
   - If intended chain is multi-letter (e.g., 'LH'), DO NOT guess/split; leave unchanged and warn.
4) Remove internal END/ENDMDL/TER from inputs; output has a single final END.
5) Skip non-data rows such as "Rigid-body (162)".
6) Manifest file includes a chains field in format:
     ABD:DF   (receptor_chains:ligand_chains)
   where chains are the final chain IDs present in the merged file for each partner.

Usage example:
  python prepare_db55_route1_prism_testfolder.py \
    --table benchmark/Table_BM5.5.tsv \
    --pdb-dir benchmark/benchmark5.5/structures \
    --out-dir benchmark/prism_inputs \
    --col-complex "Complex" \
    --sep $'\t'
"""

from __future__ import annotations

import argparse
import csv
import re
import string
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# -----------------------------
# Parsing: Complex field
# -----------------------------

def parse_complex_field(complex_field: str) -> Tuple[str, str, str]:
    """
    Example: '3EOA_LH:I' -> (bound_pdb='3EOA', partner1_chains='LH', partner2_chains='I')
    """
    s = complex_field.strip()
    if ":" not in s:
        raise ValueError(f"Complex field missing ':' separator: {complex_field}")
    left, right = s.split(":", 1)
    right = right.strip()

    if "_" in left:
        bound_pdb, ch1 = left.split("_", 1)
        bound_pdb, ch1 = bound_pdb.strip(), ch1.strip()
    else:
        bound_pdb, ch1 = left.strip(), ""

    ch2 = right.strip()
    return bound_pdb, ch1, ch2


def looks_like_case_row(complex_field: str) -> bool:
    """
    Accept only rows that look like:
      4 letters/digits + '_' + something + ':' + something
    Example: 3EOA_LH:I
    Reject: "Rigid-body (162)", blanks, etc.
    """
    s = complex_field.strip()
    if not s:
        return False
    if ":" not in s or "_" not in s:
        return False
    left = s.split(":", 1)[0]
    bound = left.split("_", 1)[0].strip()
    return len(bound) == 4


# -----------------------------
# Optional: parse model numbers from PDB ID columns if present (NMR)
# (We do NOT use these for file lookup; only for selecting MODEL blocks if needed.)
# -----------------------------

def parse_model_hint(field: str) -> Optional[int]:
    """
    Extract trailing '(N)' if present, else None.
    """
    s = (field or "").strip()
    m = re.search(r"\((\d+)\)\s*$", s)
    if not m:
        return None
    return int(m.group(1))


# -----------------------------
# PDB helpers
# -----------------------------

def get_chain_id(atom_line: str) -> str:
    return atom_line[21] if len(atom_line) > 21 else " "


def set_chain_id(atom_line: str, chain: str) -> str:
    chain = chain[:1]
    if len(atom_line) < 22:
        return atom_line
    return atom_line[:21] + chain + atom_line[22:]


def set_serial(atom_line: str, serial: int) -> str:
    if len(atom_line) < 11:
        return atom_line
    return atom_line[:6] + f"{serial:5d}" + atom_line[11:]


def unique_chains_in_order(atom_lines: List[str]) -> List[str]:
    seen = set()
    order = []
    for ln in atom_lines:
        c = get_chain_id(ln)
        if c not in seen:
            seen.add(c)
            order.append(c)
    return order


def iter_atom_lines(pdb_path: Path, keep_hetatm: bool, model: Optional[int]) -> List[str]:
    """
    Return ATOM (and optionally HETATM) lines.
    - If MODEL blocks exist and model is specified, select that model; else MODEL 1.
    - Strips END/ENDMDL/TER and other non-atom lines from the output.
    """
    lines = pdb_path.read_text(errors="ignore").splitlines()
    has_model = any(ln.startswith("MODEL") for ln in lines)

    if has_model:
        target = model if model is not None else 1
        in_model = False
        current = None
        out: List[str] = []
        for ln in lines:
            if ln.startswith("MODEL"):
                try:
                    current = int(ln.split()[1])
                except Exception:
                    current = None
                in_model = (current == target)
                continue
            if ln.startswith("ENDMDL"):
                in_model = False
                continue
            if not in_model:
                continue
            if ln.startswith("ATOM"):
                out.append(ln)
            elif keep_hetatm and ln.startswith("HETATM"):
                out.append(ln)
        return out

    out: List[str] = []
    for ln in lines:
        if ln.startswith("ATOM"):
            out.append(ln)
        elif keep_hetatm and ln.startswith("HETATM"):
            out.append(ln)
    return out


def write_with_ter_between_chains(atom_lines: List[str]) -> List[str]:
    out: List[str] = []
    prev = None
    for ln in atom_lines:
        c = get_chain_id(ln)
        if prev is not None and c != prev:
            out.append("TER")
        out.append(ln)
        prev = c
    if out and out[-1] != "TER":
        out.append("TER")
    return out


def available_chain_ids(exclude: set[str]) -> str:
    for c in (string.ascii_uppercase + string.digits + string.ascii_lowercase):
        if c not in exclude:
            return c
    return "Z"


def fill_blank_chain_if_single_target(
    atom_lines: List[str],
    target_chain: Optional[str],
    label: str,
    warnings: List[str],
) -> List[str]:
    """
    If chain IDs include blank (' '):
      - Only fill if target_chain is a single letter.
      - Otherwise leave unchanged and warn.
    """
    present = unique_chains_in_order(atom_lines)
    if " " not in present:
        return atom_lines

    if target_chain and len(target_chain) == 1:
        return [set_chain_id(ln, target_chain) if get_chain_id(ln) == " " else ln for ln in atom_lines]

    warnings.append(f"{label}: blank chain present but intended chain is multi/unknown ({target_chain}); left unchanged.")
    return atom_lines


def avoid_chain_collisions(
    rec_lines: List[str],
    lig_lines: List[str],
    warnings: List[str],
    enable: bool = True,
) -> Tuple[List[str], List[str]]:
    """
    If any chain ID appears in both partners, rename ligand chain IDs to unused IDs.
    """
    if not enable:
        return rec_lines, lig_lines

    rec_chains = set(unique_chains_in_order(rec_lines))
    lig_chains = set(unique_chains_in_order(lig_lines))

    overlap = (rec_chains & lig_chains) - {" "}
    if not overlap:
        return rec_lines, lig_lines

    used = set(c for c in rec_chains if c != " ")
    chain_map: Dict[str, str] = {}
    for ch in sorted(overlap):
        newch = available_chain_ids(used)
        used.add(newch)
        chain_map[ch] = newch
        warnings.append(f"chain collision: ligand chain {ch} renamed to {newch}")

    lig_lines2 = [set_chain_id(ln, chain_map.get(get_chain_id(ln), get_chain_id(ln))) for ln in lig_lines]
    return rec_lines, lig_lines2


def chain_string(atom_lines: List[str]) -> str:
    """
    Return chains as concatenated string in appearance order, excluding blank if any remain.
    Example: ['A','B','D'] -> "ABD"
    """
    ch = [c for c in unique_chains_in_order(atom_lines) if c != " "]
    return "".join(ch)


# -----------------------------
# Input file lookup: USE BOUND CODE FROM COMPLEX
# -----------------------------

def find_unbound_files_by_bound_code(pdb_dir: Path, bound_code: str) -> Tuple[Path, Path]:
    """
    Uses bound code from Complex column, expects:
      <BOUND>_r_u.pdb and <BOUND>_l_u.pdb
    """
    rec = pdb_dir / f"{bound_code}_r_u.pdb"
    lig = pdb_dir / f"{bound_code}_l_u.pdb"
    if rec.exists() and lig.exists():
        return rec, lig

    # try case variants
    rec2 = pdb_dir / f"{bound_code.upper()}_r_u.pdb"
    lig2 = pdb_dir / f"{bound_code.upper()}_l_u.pdb"
    if rec2.exists() and lig2.exists():
        return rec2, lig2

    rec3 = pdb_dir / f"{bound_code.lower()}_r_u.pdb"
    lig3 = pdb_dir / f"{bound_code.lower()}_l_u.pdb"
    if rec3.exists() and lig3.exists():
        return rec3, lig3

    raise FileNotFoundError(
        f"Missing unbound files for bound code {bound_code} in {pdb_dir}. "
        f"Expected {bound_code}_r_u.pdb and {bound_code}_l_u.pdb"
    )


# -----------------------------
# Merge one case
# -----------------------------

def merge_one_case(
    bound_pdb: str,
    complex_ch1: str,
    complex_ch2: str,
    rec_file: Path,
    lig_file: Path,
    rec_model: Optional[int],
    lig_model: Optional[int],
    keep_hetatm: bool,
    out_pdb: Path,
    collision_fix: bool,
) -> Tuple[str, str, List[str]]:
    """
    Returns (rec_chain_str, lig_chain_str, warnings)
    """
    warnings: List[str] = []

    rec_atoms = iter_atom_lines(rec_file, keep_hetatm=keep_hetatm, model=rec_model)
    lig_atoms = iter_atom_lines(lig_file, keep_hetatm=keep_hetatm, model=lig_model)

    # Fill blank chains using Complex-side single-letter intended chain ONLY.
    rec_atoms = fill_blank_chain_if_single_target(
        rec_atoms,
        target_chain=complex_ch1 if len(complex_ch1) == 1 else None,
        label="receptor",
        warnings=warnings,
    )
    lig_atoms = fill_blank_chain_if_single_target(
        lig_atoms,
        target_chain=complex_ch2 if len(complex_ch2) == 1 else None,
        label="ligand",
        warnings=warnings,
    )

    # Avoid chain collisions (optional but recommended for PRISM)
    rec_atoms, lig_atoms = avoid_chain_collisions(rec_atoms, lig_atoms, warnings, enable=collision_fix)

    # Prepare blocks with TER inside each partner when chain changes
    rec_block = write_with_ter_between_chains(rec_atoms)
    lig_block = write_with_ter_between_chains(lig_atoms)

    # Compute final chain strings for manifest
    rec_chs = chain_string(rec_atoms)
    lig_chs = chain_string(lig_atoms)

    merged: List[str] = []
    merged.append("REMARK   PRISM MERGED INPUT: unbound receptor then unbound ligand")
    merged.append(f"REMARK   bound_code: {bound_pdb}")
    merged.append(f"REMARK   receptor_source: {rec_file.name}")
    merged.append(f"REMARK   ligand_source:   {lig_file.name}")
    merged.append(f"REMARK   complex_chains_partner1: {complex_ch1}")
    merged.append(f"REMARK   complex_chains_partner2: {complex_ch2}")
    if rec_model is not None:
        merged.append(f"REMARK   receptor_model_selected: {rec_model}")
    if lig_model is not None:
        merged.append(f"REMARK   ligand_model_selected:   {lig_model}")
    if warnings:
        for w in warnings:
            merged.append(f"REMARK   WARNING: {w}")

    serial = 1
    for ln in rec_block:
        if ln == "TER":
            merged.append("TER")
        else:
            merged.append(set_serial(ln, serial))
            serial += 1

    merged.append("TER")  # separator between partners

    for ln in lig_block:
        if ln == "TER":
            merged.append("TER")
        else:
            merged.append(set_serial(ln, serial))
            serial += 1

    merged.append("END")
    out_pdb.write_text("\n".join(merged) + "\n")

    return rec_chs, lig_chs, warnings


# -----------------------------
# Main
# -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Merge benchmark unbound receptor+ligand into single PRISM input PDBs named by Complex bound code (XXXX.pdb)."
    )
    ap.add_argument("--table", required=True, help="Benchmark table TSV/CSV (text).")
    ap.add_argument("--sep", default="\t", help="Delimiter (default: tab).")
    ap.add_argument("--pdb-dir", required=True, help="Directory containing unbound PDBs named XXXX_r_u.pdb and XXXX_l_u.pdb (XXXX from Complex column).")
    ap.add_argument("--out-dir", required=True, help="Output directory for merged PDBs.")
    ap.add_argument("--keep-hetatm", action="store_true", help="Keep HETATM lines (default: ATOM only).")
    ap.add_argument("--no-collision-fix", action="store_true", help="Disable chain-collision renaming (not recommended).")

    ap.add_argument("--col-complex", default="Complex", help="Column name for Complex (e.g., 3EOA_LH:I).")

    # Optional columns: only used to extract (NMR) model numbers like '(6)' if present.
    ap.add_argument("--col-pdb1", default=None, help="Optional column name for PDB ID 1 (for extracting model like (6)). Not used for file lookup.")
    ap.add_argument("--col-pdb2", default=None, help="Optional column name for PDB ID 2 (for extracting model like (6)). Not used for file lookup.")

    args = ap.parse_args()

    table = Path(args.table)
    pdb_dir = Path(args.pdb_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = out_dir / "merge_manifest.tsv"
    with manifest.open("w", newline="") as mf:
        mf.write("Complex\tBound\tInFiles\tChains\tOutFile\tWarnings\n")

        with table.open("r", newline="") as f:
            reader = csv.DictReader(f, delimiter=args.sep)

            for row_i, row in enumerate(reader, start=1):
                complex_field = (row.get(args.col_complex) or "").strip()

                # Skip group header rows / empty rows
                if not looks_like_case_row(complex_field):
                    continue

                bound_pdb, complex_ch1, complex_ch2 = parse_complex_field(complex_field)

                # Find unbound files by BOUND code from Complex column
                rec_file, lig_file = find_unbound_files_by_bound_code(pdb_dir, bound_pdb)

                # Optional: model numbers from PDB ID columns (if you provide them)
                rec_model = parse_model_hint(row.get(args.col_pdb1, "")) if args.col_pdb1 else None
                lig_model = parse_model_hint(row.get(args.col_pdb2, "")) if args.col_pdb2 else None

                # Output file name: ONLY 4-letter bound code
                out_pdb = out_dir / f"{bound_pdb}.pdb"

                rec_chs, lig_chs, warnings = merge_one_case(
                    bound_pdb=bound_pdb,
                    complex_ch1=complex_ch1,
                    complex_ch2=complex_ch2,
                    rec_file=rec_file,
                    lig_file=lig_file,
                    rec_model=rec_model,
                    lig_model=lig_model,
                    keep_hetatm=args.keep_hetatm,
                    out_pdb=out_pdb,
                    collision_fix=(not args.no_collision_fix),
                )

                chains_fmt = f"{rec_chs}:{lig_chs}"
                infiles = f"{rec_file.name},{lig_file.name}"
                warn_str = " | ".join(warnings)

                mf.write(f"{complex_field}\t{bound_pdb}\t{infiles}\t{chains_fmt}\t{out_pdb.name}\t{warn_str}\n")

    print(f"Done. Outputs in: {out_dir}")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
