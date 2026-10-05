
#!/usr/bin/env python3
"""Replay retained six-pair MultiProt transforms through isolated FiberDock."""

from __future__ import annotations

import concurrent.futures
import csv
import json
import os
import shutil
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def pdb_integrity(path: Path) -> dict:
    from Bio.PDB import PDBParser
    result = {"path": str(path), "exists": path.is_file(), "valid": False}
    if not path.is_file():
        return result
    try:
        structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
        atoms = 0
        residues = 0
        chains = []
        for model in structure:
            for chain in model:
                chains.append(chain.id)
                for residue in chain:
                    residues += residue.id[0] == " "
                    atoms += len(residue.child_list)
        result.update({"valid": bool(atoms and chains), "atoms": atoms,
                       "residues": residues, "chains": "".join(chains)})
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def parse_energy(path: Path | None) -> dict | None:
    if path is None or not path.is_file():
        return None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [item.strip() for item in line.split("|")]
        if len(parts) >= 3 and parts[0].isdigit():
            try:
                energy = float(parts[1])
            except ValueError:
                energy = parts[1]
            return {"solution": int(parts[0]), "global_energy": energy, "raw_line": line}
    return None


def alignment_record(root: Path, query: str, template: str, chain: str) -> dict:
    path = root / "processed" / "alignment" / f"{query}_{template}_{chain}.json"
    if not path.is_file():
        return {"status": "missing", "path": str(path)}
    data = json.loads(path.read_text(encoding="utf-8"))
    return {"status": data.get("status", "missing"),
            "match_count": data.get("match_count", 0),
            "tm_score": data.get("tm_score", 0.0),
            "rmsd": data.get("rmsd"), "path": str(path)}


def discover_candidates(mp_root: Path, tm_root: Path) -> list[dict]:
    inputs = list(csv.DictReader((mp_root / "inputs.csv").open(encoding="utf-8")))
    transformation = mp_root / "processed" / "transformation"
    candidates = []
    for left_path in sorted(transformation.glob("*_L.pdb")):
        stem = left_path.stem[:-2]
        right_path = transformation / f"{stem}_R.pdb"
        if not right_path.is_file():
            continue
        matched = None
        for row in inputs:
            ql, qr = row["Receptor"], row["Ligand"]
            suffix = f"_{ql}_{qr}_"
            if suffix in stem:
                matched = (ql, qr, stem.split(suffix, 1)[0], stem.rsplit("_", 1)[-1])
                break
        if matched is None:
            continue
        ql, qr, template, orientation = matched
        if len(template) != 6 or orientation not in {"o1", "o2"}:
            continue
        cl, cr = (template[4], template[5]) if orientation == "o1" else (template[5], template[4])
        tm_left = tm_root / "processed" / "transformation" / left_path.name
        tm_right = tm_root / "processed" / "transformation" / right_path.name
        candidates.append({
            "label": stem, "query_left": ql, "query_right": qr,
            "template": template, "orientation": orientation,
            "left_chain": cl, "right_chain": cr,
            "mp_left": str(left_path.resolve()), "mp_right": str(right_path.resolve()),
            "tm_left": str(tm_left.resolve()), "tm_right": str(tm_right.resolve()),
            "mp_alignment_left": alignment_record(mp_root, ql, template, cl),
            "mp_alignment_right": alignment_record(mp_root, qr, template, cr),
            "tm_alignment_left": alignment_record(tm_root, ql, template, cl),
            "tm_alignment_right": alignment_record(tm_root, qr, template, cr),
            "tm_transform_exists": tm_left.is_file() and tm_right.is_file(),
        })
    return candidates


def run_candidate(candidate: dict, output_root: Path, fiber_source: Path) -> dict:
    label = candidate["label"]
    root = output_root / label
    input_root, work = root / "input", root / "work"
    tools = root / "external_tools" / "fiberdock"
    input_root.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy2(candidate["mp_left"], input_root / f"{label}_L.pdb")
    shutil.copy2(candidate["mp_right"], input_root / f"{label}_R.pdb")
    shutil.copytree(fiber_source, tools, symlinks=True, dirs_exist_ok=True)
    os.chdir(root)
    os.environ["PRISM_FIBERDOCK_DIR"] = str(tools)
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    import src.fiberdock_refinement as fd
    fd.FIBERDOCK_DIR = str(tools)
    left, right = input_root / f"{label}_L.pdb", input_root / f"{label}_R.pdb"
    start = time.perf_counter()
    result = {
        "label": label, "query_left": candidate["query_left"],
        "query_right": candidate["query_right"], "template": candidate["template"],
        "orientation": candidate["orientation"],
        "mp_alignment_left": candidate["mp_alignment_left"],
        "mp_alignment_right": candidate["mp_alignment_right"],
        "tm_alignment_left": candidate["tm_alignment_left"],
        "tm_alignment_right": candidate["tm_alignment_right"],
        "tm_transform_exists": candidate["tm_transform_exists"],
        "source_left": candidate["mp_left"], "source_right": candidate["mp_right"],
    }
    try:
        r_hb, l_hb = fd._add_hydrogens(str(right), str(work)), fd._add_hydrogens(str(left), str(work))
        r_ca, r_sizes = fd._create_ca_pdb(str(right), str(work))
        l_ca, l_sizes = fd._create_ca_pdb(str(left), str(work))
        result["hydrogens"] = {"right": bool(r_hb), "left": bool(l_hb)}
        result["ca_inputs"] = {"right": str(r_ca), "left": str(l_ca),
                               "right_sizes": r_sizes, "left_sizes": l_sizes}
        if not r_hb or not l_hb:
            raise RuntimeError("hydrogenation failed")
        r_nma, l_nma = fd._run_nma(r_ca, str(work)), fd._run_nma(l_ca, str(work))
        result["nma"] = {"right": bool(r_nma), "left": bool(l_nma)}
        if not r_nma or not l_nma:
            raise RuntimeError("NMA failed")
        if sum(r_sizes.values()) >= sum(l_sizes.values()):
            receptor_base, ligand_base = right.stem, left.stem
        else:
            receptor_base, ligand_base = left.stem, right.stem
            r_hb, l_hb = l_hb, r_hb
        params = fd._build_fiberdock_params(r_hb, l_hb, str(work), receptor_base, ligand_base)
        result["params"] = str(params)
        if not Path(params).is_file():
            raise RuntimeError(f"params file missing: {params}")
        result["current_parser_energy"] = fd._run_fiberdock(
            params, str(tools), str(work), receptor_base, ligand_base
        )
        declared_stem = next(
            (line.split(None, 1)[1].strip()
             for line in Path(params).read_text(encoding="utf-8").splitlines()
             if line.startswith("energiesOutFileName ")),
            None,
        )
        result["declared_energy_stem"] = declared_stem
        declared_ref = Path(declared_stem + ".ref") if declared_stem else None
        result["declared_ref"] = str(declared_ref) if declared_ref else None
        result["declared_solution"] = parse_energy(declared_ref)
        result["refined_pdbs"] = [
            pdb_integrity(path) for path in sorted(work.glob("fiberdock_energies*.ref.pdb"))
        ]
        result["status"], result["return_code"] = "completed", 0
    except Exception as exc:
        result["status"], result["return_code"] = "failed", 1
        result["error"] = f"{type(exc).__name__}: {exc}"
    result["elapsed_seconds"] = time.perf_counter() - start
    return result


def run_task(task: tuple[dict, Path, Path]) -> dict:
    return run_candidate(*task)


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: run_multiprot_fiberdock_broader_replay.py MP_ROOT TM_ROOT OUTPUT_ROOT")
    mp_root, tm_root, output_root = (Path(value).resolve() for value in sys.argv[1:])
    if output_root.exists() and any(output_root.iterdir()):
        raise RuntimeError(f"output root is not empty: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    candidates = discover_candidates(mp_root, tm_root)
    if not candidates:
        raise RuntimeError("no transformed candidates discovered")
    fiber_source = mp_root / "external_tools" / "fiberdock"
    workers = min(len(candidates), max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))))
    tasks = [(candidate, output_root, fiber_source) for candidate in candidates]
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        records = list(pool.map(run_task, tasks))
    records.sort(key=lambda row: row["label"])
    summary = {
        "candidate_count": len(records),
        "completed_count": sum(row["status"] == "completed" for row in records),
        "tm_transform_exists_count": sum(row["tm_transform_exists"] for row in records),
        "declared_output_count": sum(row.get("declared_solution") is not None for row in records),
        "valid_pdb_count": sum(item["valid"] for row in records for item in row.get("refined_pdbs", [])),
        "current_parser_dash_count": sum(row.get("current_parser_energy") == "-" for row in records),
        "corrected_energy_values": [row["declared_solution"]["global_energy"] for row in records if row.get("declared_solution") is not None],
    }
    report = {"experiment": "multiprot_fiberdock_broader_replay",
              "mp_root": str(mp_root), "tm_root": str(tm_root),
              "output_root": str(output_root), "candidate_count": len(candidates),
              "workers": workers, "records": records, "summary": summary}
    (output_root / "multiprot_fiberdock_broader_replay.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if summary["completed_count"] == len(records) else 2


if __name__ == "__main__":
    raise SystemExit(main())
