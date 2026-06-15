"""Compare PRISM output complexes to downloaded native PDBs and score with DockQ."""

from __future__ import annotations

import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Sequence, Tuple

from Bio.PDB import PDBParser
from Bio.PDB.Superimposer import Superimposer

from .eval.dockq import calculate_dockq, dockq_to_capri_class
from .eval.irmsd_backbone import calculate_irmsd_backbone
from .pdb_download import split_target_id

PROCESSED_PDB_DIR = "processed/pdbs"
SUMMARY_CSV = "processed/summary.csv"
TRIMMED_NATIVE_DIR = "processed/compare/native"


def _residue_key(residue) -> Tuple[str, int, str]:
    return (residue.get_parent().id, residue.id[1], residue.id[2].strip() or "")


def _ca_by_key(pdb_path: str, chain_ids: Sequence[str]) -> Dict[Tuple[str, int, str], object]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("cmp", pdb_path)
    model = next(structure.get_models())
    out = {}
    chain_set = set(chain_ids)
    for chain in model:
        if chain.id not in chain_set:
            continue
        for residue in chain:
            if residue.id[0] != " " or "CA" not in residue:
                continue
            out[_residue_key(residue)] = residue["CA"]
    return out


def ca_rmsd(
    model_pdb: str,
    native_pdb: str,
    model_chains: Sequence[str],
    native_chains: Sequence[str],
) -> Tuple[Optional[float], int]:
    """Superimpose model chains onto native and return (RMSD Å, n_matched_CA)."""
    if len(model_chains) != len(native_chains):
        raise ValueError("model_chains and native_chains must have the same length")
    model_map = _ca_by_key(model_pdb, model_chains)
    native_map = _ca_by_key(native_pdb, native_chains)
    ref_atoms = []
    mob_atoms = []
    for mc, nc in zip(model_chains, native_chains):
        if mc == nc:
            keys = sorted(k for k in model_map if k in native_map and k[0] == mc)
            for key in keys:
                ref_atoms.append(native_map[key])
                mob_atoms.append(model_map[key])
        else:
            for key in sorted(k for k in model_map if k[0] == mc):
                native_key = (nc, key[1], key[2])
                if native_key in native_map:
                    ref_atoms.append(native_map[native_key])
                    mob_atoms.append(model_map[key])
    if len(ref_atoms) < 3:
        return None, len(ref_atoms)
    sup = Superimposer()
    sup.set_atoms(ref_atoms, mob_atoms)
    return float(sup.rms), len(ref_atoms)


def _write_trimmed_native(
    rec_path: str,
    rec_chains: Sequence[str],
    lig_path: str,
    lig_chains: Sequence[str],
    output_path: str,
) -> None:
    keep_rec = set(rec_chains)
    keep_lig = set(lig_chains)
    serial = 1
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as out_f:
        for path, keep in ((rec_path, keep_rec), (lig_path, keep_lig)):
            with open(path) as fh:
                for line in fh:
                    if line.startswith(("ATOM", "HETATM")):
                        if len(line) > 21 and line[21] in keep:
                            out_f.write(f"{line[:6]}{serial:5d}{line[11:]}")
                            serial += 1
                    elif line.startswith("TER"):
                        out_f.write(line)
            out_f.write("TER\n")
        out_f.write("END\n")


def get_trimmed_native_pdb(receptor: str, ligand: str) -> Tuple[str, str, str]:
    """Return cached native PDB with only receptor/ligand chains (trimmed)."""
    rec_pdb_id, rec_chains = split_target_id(receptor)
    lig_pdb_id, lig_chains = split_target_id(ligand)
    rec_chain_str = "".join(rec_chains)
    lig_chain_str = "".join(lig_chains)

    rec_path = os.path.join(PROCESSED_PDB_DIR, f"{rec_pdb_id}.pdb")
    lig_path = os.path.join(PROCESSED_PDB_DIR, f"{lig_pdb_id}.pdb")
    if not os.path.exists(rec_path):
        raise FileNotFoundError(f"Native PDB not found: {rec_path}")
    if not os.path.exists(lig_path):
        raise FileNotFoundError(f"Native PDB not found: {lig_path}")

    os.makedirs(TRIMMED_NATIVE_DIR, exist_ok=True)
    trimmed = os.path.join(
        TRIMMED_NATIVE_DIR,
        f"{rec_pdb_id}_{lig_pdb_id}_{rec_chain_str}{lig_chain_str}.pdb",
    )
    if not os.path.exists(trimmed):
        if rec_pdb_id == lig_pdb_id:
            _write_trimmed_native(rec_path, rec_chains, rec_path, lig_chains, trimmed)
        else:
            _write_trimmed_native(rec_path, rec_chains, lig_path, lig_chains, trimmed)
    return trimmed, rec_chain_str, lig_chain_str


def prepare_trimmed_natives(pairs: Sequence[Tuple]) -> None:
    """Build trimmed native PDB cache for all unique receptor/ligand pairs."""
    seen = set()
    for entry in pairs:
        if len(entry) < 2:
            continue
        key = (entry[0], entry[1])
        if key in seen:
            continue
        seen.add(key)
        get_trimmed_native_pdb(key[0], key[1])


def run_dockq(
    model_pdb: str,
    native_pdb: str,
    model_receptor_chains: str,
    model_ligand_chains: str,
    native_receptor_chains: str,
    native_ligand_chains: str,
    no_align: bool = False,
) -> Dict:
    mapping = (
        f"{model_receptor_chains}{model_ligand_chains}:"
        f"{native_receptor_chains}{native_ligand_chains}"
    )
    try:
        result = calculate_dockq(
            model_pdb,
            native_pdb,
            mapping=mapping,
            no_align=no_align,
        )
        result["capri_class"] = dockq_to_capri_class(result.get("dockq"))
        result["error_dockq"] = ""
    except Exception as exc:
        result = {
            "dockq": None,
            "fnat": None,
            "irmsd": None,
            "lrmsd": None,
            "fnonnat": None,
            "f1": None,
            "clashes": None,
            "capri_class": "Unknown",
            "error_dockq": str(exc),
        }
    result["dockq_mapping"] = mapping
    return result


def run_irmsd_backbone(
    model_pdb: str,
    native_pdb: str,
    model_receptor_chains: str,
    model_ligand_chains: str,
    native_receptor_chains: str,
    native_ligand_chains: str,
) -> Tuple[Optional[float], str]:
    try:
        value = calculate_irmsd_backbone(
            model_pdb,
            model_receptor_chains,
            model_ligand_chains,
            native_pdb,
            native_receptor_chains,
            native_ligand_chains,
        )
        return value, ""
    except Exception as exc:
        return None, str(exc)


def compare_pair(
    receptor: str,
    ligand: str,
    template: str,
    output_pdb: str,
    dockq_no_align: bool = False,
) -> Dict:
    """Compare one PRISM output to trimmed native input structures."""
    rec_pdb_id, rec_chains = split_target_id(receptor)
    lig_pdb_id, lig_chains = split_target_id(ligand)
    model_rec = "".join(rec_chains)
    model_lig = "".join(lig_chains)

    row = {
        "receptor": receptor,
        "ligand": ligand,
        "template": template,
        "output_pdb": output_pdb,
        "native_pdb": "",
        "receptor_ca_rmsd": None,
        "receptor_ca_matched": 0,
        "ligand_ca_rmsd": None,
        "ligand_ca_matched": 0,
        "irmsd_backbone": None,
        "dockq": None,
        "dockq_fnat": None,
        "dockq_irmsd": None,
        "dockq_lrmsd": None,
        "dockq_fnonnat": None,
        "dockq_f1": None,
        "dockq_clashes": None,
        "dockq_capri": "",
        "dockq_mapping": "",
        "error_receptor_rmsd": "",
        "error_ligand_rmsd": "",
        "error_irmsd_backbone": "",
        "error_dockq": "",
    }

    if not os.path.exists(output_pdb):
        row["error_dockq"] = f"output PDB not found: {output_pdb}"
        return row

    try:
        native_pdb, native_rec, native_lig = get_trimmed_native_pdb(receptor, ligand)
        row["native_pdb"] = native_pdb
    except Exception as exc:
        row["error_dockq"] = str(exc)
        return row

    try:
        rmsd, n = ca_rmsd(output_pdb, native_pdb, rec_chains, rec_chains)
        row["receptor_ca_rmsd"] = rmsd
        row["receptor_ca_matched"] = n
    except Exception as exc:
        row["error_receptor_rmsd"] = str(exc)

    try:
        rmsd, n = ca_rmsd(output_pdb, native_pdb, lig_chains, lig_chains)
        row["ligand_ca_rmsd"] = rmsd
        row["ligand_ca_matched"] = n
    except Exception as exc:
        row["error_ligand_rmsd"] = str(exc)

    irmsd, irmsd_err = run_irmsd_backbone(
        output_pdb,
        native_pdb,
        model_rec,
        model_lig,
        native_rec,
        native_lig,
    )
    row["irmsd_backbone"] = irmsd
    row["error_irmsd_backbone"] = irmsd_err

    dockq = run_dockq(
        output_pdb,
        native_pdb,
        model_rec,
        model_lig,
        native_rec,
        native_lig,
        no_align=dockq_no_align,
    )
    row["dockq"] = dockq.get("dockq")
    row["dockq_fnat"] = dockq.get("fnat")
    row["dockq_irmsd"] = dockq.get("irmsd")
    row["dockq_lrmsd"] = dockq.get("lrmsd")
    row["dockq_fnonnat"] = dockq.get("fnonnat")
    row["dockq_f1"] = dockq.get("f1")
    row["dockq_clashes"] = dockq.get("clashes")
    row["dockq_capri"] = dockq.get("capri_class", "")
    row["dockq_mapping"] = dockq.get("dockq_mapping", "")
    row["error_dockq"] = dockq.get("error_dockq", "")
    return row


def _compare_task(args: Tuple) -> Dict:
    receptor, ligand, template, output_pdb, dockq_no_align = args
    return compare_pair(
        receptor,
        ligand,
        template,
        output_pdb,
        dockq_no_align=dockq_no_align,
    )


SUMMARY_FIELDNAMES = [
    "receptor",
    "ligand",
    "template",
    "output_pdb",
    "native_pdb",
    "receptor_ca_rmsd",
    "receptor_ca_matched",
    "ligand_ca_rmsd",
    "ligand_ca_matched",
    "irmsd_backbone",
    "dockq",
    "dockq_fnat",
    "dockq_irmsd",
    "dockq_lrmsd",
    "dockq_fnonnat",
    "dockq_f1",
    "dockq_clashes",
    "dockq_capri",
    "dockq_mapping",
    "error_receptor_rmsd",
    "error_ligand_rmsd",
    "error_irmsd_backbone",
    "error_dockq",
]


def compare_and_summarize(
    passed_pairs: Sequence[Tuple],
    summary_csv: str = SUMMARY_CSV,
    dockq_no_align: bool = False,
    n_jobs: int = 1,
) -> Tuple[str, List[Dict]]:
    """Score all accepted complexes; write ``processed/summary.csv``."""
    os.makedirs(os.path.dirname(summary_csv) or ".", exist_ok=True)
    pairs = [e for e in passed_pairs if len(e) >= 4]
    total = len(pairs)
    if total == 0:
        with open(summary_csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDNAMES)
            writer.writeheader()
        return summary_csv, []

    prepare_trimmed_natives(pairs)
    tasks = [
        (entry[0], entry[1], entry[2], entry[3], dockq_no_align) for entry in pairs
    ]

    rows: List[Dict] = [None] * total  # type: ignore
    n_workers = max(1, min(n_jobs, total))

    if n_workers == 1:
        for i, task in enumerate(tasks, start=1):
            rows[i - 1] = _compare_task(task)
            if total > 10 and (i == 1 or i % 50 == 0 or i == total):
                print(f"  compared {i}/{total}: {task[2]}_{task[0]}_{task[1]}", flush=True)
    else:
        print(f"  parallel compare with {n_workers} workers", flush=True)
        done = 0
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_map = {
                executor.submit(_compare_task, task): idx for idx, task in enumerate(tasks)
            }
            for future in as_completed(future_map):
                idx = future_map[future]
                rows[idx] = future.result()
                done += 1
                if done == 1 or done % 50 == 0 or done == total:
                    task = tasks[idx]
                    print(
                        f"  compared {done}/{total}: {task[2]}_{task[0]}_{task[1]}",
                        flush=True,
                    )

    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return summary_csv, rows
