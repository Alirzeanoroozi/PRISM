#!/usr/bin/env python3
"""Audit BM55 native chain mappings against the selected PRISM manifest.

The benchmark ``Complex`` field is the expected native mapping, while the
selected manifest and assembled native PDBs are the run-specific evidence.
This script never edits either source.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


COMPLEX_RE = re.compile(r"^(?P<pdb>[A-Za-z0-9]{4})_(?P<receptor>[A-Za-z0-9]+):(?P<ligand>[A-Za-z0-9]+)$")
NATIVE_RE = re.compile(r"^native_(?P<split>[A-Za-z0-9]+)_(?P<pdb>[A-Za-z0-9]{4})_(?P<index>\d+)\.pdb$", re.I)


def clean_complex(value: str) -> str:
    return re.sub(r"\s*\*+\s*$", "", value.strip())


def normalize_group(value: object) -> str:
    return re.sub(r"\*+$", "", "".join(str(value or "").split()))


def parse_complex(value: str) -> tuple[str, str, str]:
    cleaned = clean_complex(value)
    match = COMPLEX_RE.fullmatch(cleaned)
    if not match:
        raise ValueError(f"invalid benchmark Complex value: {value!r}")
    return (
        match.group("pdb").lower(),
        match.group("receptor"),
        match.group("ligand"),
    )


def pdb_chain_ids(path: Path) -> list[str]:
    seen: list[str] = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")) or len(line) <= 21:
                continue
            chain = line[21].strip() or "_"
            if chain not in seen:
                seen.append(chain)
    return seen


def read_expected(t_csv_dir: Path) -> tuple[dict[tuple[str, str], list[dict[str, str]]], Counter[str]]:
    expected: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    stats: Counter[str] = Counter()
    for path in sorted(t_csv_dir.glob("T_*.csv")):
        split = path.stem[2:].lower()
        with path.open(newline="", encoding="utf-8") as handle:
            for line_number, row in enumerate(csv.DictReader(handle), start=2):
                stats[split] += 1
                pdb, receptor, ligand = parse_complex(row.get("Complex", ""))
                expected[(split, pdb)].append(
                    {
                        "split": split,
                        "pdb": pdb,
                        "receptor": receptor,
                        "ligand": ligand,
                        "complex": clean_complex(row["Complex"]),
                        "source": str(path),
                        "line": str(line_number),
                    }
                )
    return expected, stats


def case_identity(row: dict[str, str]) -> tuple[str, str]:
    split = (row.get("split") or "").strip().lower()
    native_name = Path(row.get("native_pdb", "")).name
    match = NATIVE_RE.match(native_name)
    if match:
        return match.group("split").lower(), match.group("pdb").lower()
    case_id = (row.get("case_id") or "").strip().lower()
    parts = case_id.split("_")
    pdb = next((part for part in parts if re.fullmatch(r"[a-z0-9]{4}", part)), "")
    return split, pdb


def load_cases(selected_csv: Path) -> dict[str, list[dict[str, str]]]:
    cases: dict[str, list[dict[str, str]]] = defaultdict(list)
    with selected_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cases[row.get("case_id", "")].append(row)
    return cases


def inventory(selected_csv: Path, output: Path) -> None:
    cases = load_cases(selected_csv)
    rows: list[dict[str, str]] = []
    for case_id, case_rows in sorted(cases.items()):
        row = case_rows[0]
        native = Path(row.get("native_pdb", ""))
        record = {
            "case_id": case_id,
            "split": row.get("split", ""),
            "native_pdb": str(native),
            "native_receptor_chains": normalize_group(row.get("native_receptor_chains", "")),
            "native_ligand_chains": normalize_group(row.get("native_ligand_chains", "")),
            "actual_chain_ids": "",
            "actual_chain_count": "",
            "actual_missing_expected": "",
            "actual_extra_chains": "",
            "status": "missing_file",
            "error": "",
        }
        try:
            actual = pdb_chain_ids(native)
            expected = list(
                normalize_group(row.get("native_receptor_chains", ""))
                + normalize_group(row.get("native_ligand_chains", ""))
            )
            missing = sorted(set(expected) - set(actual))
            extra = sorted(set(actual) - set(expected))
            if missing:
                status = "missing_expected_chains"
            elif extra:
                status = "extra_chains"
            else:
                status = "exact"
            record.update(
                actual_chain_ids="".join(actual),
                actual_chain_count=str(len(actual)),
                actual_missing_expected="".join(missing),
                actual_extra_chains="".join(extra),
                status=status,
            )
        except Exception as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
        rows.append(record)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["case_id"])
        writer.writeheader()
        writer.writerows(rows)


def audit(selected_csv: Path, t_csv_dir: Path, actual_tsv: Path | None, output_dir: Path) -> None:
    expected, t_counts = read_expected(t_csv_dir)
    cases = load_cases(selected_csv)
    actual_by_case: dict[str, dict[str, str]] = {}
    if actual_tsv and actual_tsv.is_file():
        with actual_tsv.open(newline="", encoding="utf-8") as handle:
            actual_by_case = {row["case_id"]: row for row in csv.DictReader(handle)}

    rows: list[dict[str, str]] = []
    for case_id, case_rows in sorted(cases.items()):
        first = case_rows[0]
        split, pdb = case_identity(first)
        mappings = {
            (
                normalize_group(r.get("native_receptor_chains", "")),
                normalize_group(r.get("native_ligand_chains", "")),
            )
            for r in case_rows
        }
        candidates = expected.get((split, pdb), [])
        manifest_mapping = next(iter(mappings)) if len(mappings) == 1 else ("", "")
        exact = [r for r in candidates if (r["receptor"], r["ligand"]) == manifest_mapping]
        status = "matched" if len(mappings) == 1 and len(exact) == 1 else "chain_mismatch"
        if not candidates:
            status = "missing_T_row"
        elif len(candidates) > 1 and not exact:
            status = "ambiguous_T_rows"
        actual = actual_by_case.get(case_id, {})
        rows.append(
            {
                "case_id": case_id,
                "split": split,
                "pdb": pdb,
                "candidate_rows": str(len(case_rows)),
                "manifest_receptor_chains": manifest_mapping[0],
                "manifest_ligand_chains": manifest_mapping[1],
                "expected_complexes": ";".join(r["complex"] for r in candidates),
                "expected_receptor_chains": ";".join(r["receptor"] for r in candidates),
                "expected_ligand_chains": ";".join(r["ligand"] for r in candidates),
                "native_pdb": first.get("native_pdb", ""),
                "manifest_consistent": "yes" if len(mappings) == 1 else "no",
                "actual_chain_ids": actual.get("actual_chain_ids", ""),
                "actual_chain_status": actual.get("status", "not_checked"),
                "status": status,
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    table = output_dir / "native_chain_mapping_audit.csv"
    with table.open("w", newline="", encoding="utf-8") as handle:
        fields = list(rows[0]) if rows else ["case_id"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "selected_manifest": str(selected_csv),
        "benchmark_dir": str(t_csv_dir),
        "selected_rows": sum(int(row["candidate_rows"]) for row in rows),
        "unique_cases": len(rows),
        "benchmark_rows": dict(sorted(t_counts.items())),
        "audit_status_counts": dict(sorted(Counter(row["status"] for row in rows).items())),
        "actual_chain_status_counts": dict(sorted(Counter(row["actual_chain_status"] for row in rows).items())),
        "audit_csv": str(table.resolve()),
    }
    (output_dir / "audit_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-csv", required=True, type=Path)
    parser.add_argument("--t-csv-dir", type=Path)
    parser.add_argument("--actual-tsv", type=Path)
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.inventory_only:
        inventory(args.selected_csv.resolve(), args.output.resolve())
    else:
        if args.t_csv_dir is None:
            parser.error("--t-csv-dir is required unless --inventory-only is used")
        audit(args.selected_csv.resolve(), args.t_csv_dir.resolve(), args.actual_tsv.resolve() if args.actual_tsv else None, args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
