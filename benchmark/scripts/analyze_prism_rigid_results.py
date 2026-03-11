#!/usr/bin/env python3
"""
Analyze PRISM Rosetta outputs against Benchmark 5.0 T_Rigid chain-pair definitions.

This script does not modify existing repository code. It creates analysis CSV files:

1) per_prediction.csv  : one row per PRISM prediction model
2) pair_summary.csv    : aggregated metrics per target chain-pair (across templates)

It maps PRISM filenames like:
    1ahwBC_1fgnH_0_1tfhA_0.rosetta_0001_0001.pdb
to T_Rigid rows like:
    Complex=1AHW_AB:C, PDB ID 1=1FGN_LH, PDB ID 2=1TFH_A

and evaluates the specific pair interface (e.g. 1AHW B:C).

Notes
-----
- Native bound complex PDBs (e.g. 1ahw.pdb) are required for DockQ/iRMSD scoring.
  They are often *not* present in benchmark/data/pdbs/rigid (that directory mostly holds
  unbound constituents). Missing native files are reported in the output.
- The script assumes the user intended to run both benchmark/scripts/dockq.py and
  benchmark/scripts/irmsd.py (the prompt repeated dockq.py twice).
"""

import argparse
import csv
import gzip
import math
import re
import statistics
import subprocess
import shutil
import socket
import sys
import urllib.request
import os
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple


ROSETTA_PDB_RE = re.compile(
    r"^(?P<target>[A-Za-z0-9]+)_(?P<tpl1>[A-Za-z0-9]+)_(?P<idx1>\d+)_(?P<tpl2>[A-Za-z0-9]+)_(?P<idx2>\d+)\.(?P<rest>rosetta.*)\.pdb$"
)


class PdbToken(NamedTuple):
    raw: str
    pdb4: str
    chains: str


class TRigidPairDef(NamedTuple):
    row_index: int
    complex_raw: str
    bound_pdb4: str
    bound_receptor_group: str
    bound_ligand_group: str
    ref_receptor_chain: str
    ref_ligand_chain: str
    pdb1: PdbToken
    pdb2: PdbToken
    pair_target_token: str  # e.g., 1ahwBC


class PredictionRecord(NamedTuple):
    model_path: Path
    joblist: str
    filename: str
    target_token: str
    tpl1: PdbToken
    tpl2: PdbToken
    tpl1_idx: int
    tpl2_idx: int
    rosetta_tag: str


def parse_bound_complex_field(s: str) -> Tuple[str, str, str]:
    s = (s or "").strip()
    if not s or "_" not in s or ":" not in s:
        raise ValueError(f"Unexpected Complex format: {s!r}")
    pdb_part, chains_part = s.split("_", 1)
    rec, lig = chains_part.split(":", 1)
    return pdb_part[:4].lower(), rec.strip(), lig.strip()


def parse_pdb_token(token: str) -> PdbToken:
    token = (token or "").strip()
    if not token:
        return PdbToken(raw="", pdb4="", chains="")
    if "_" in token:
        pdb_part, chains = token.split("_", 1)
    else:
        pdb_part, chains = token[:4], token[4:]
    pdb4 = pdb_part[:4].lower()
    chains = "".join(ch for ch in chains if ch.isalnum())
    return PdbToken(raw=token, pdb4=pdb4, chains=chains)


def _unordered_pair_key(a: str, b: str) -> Tuple[str, str]:
    a = (a or "").lower()
    b = (b or "").lower()
    return (a, b) if a <= b else (b, a)


def load_t_rigid_pairs(csv_path: Path) -> Dict[Tuple[str, str], List[TRigidPairDef]]:
    # Index by unordered (PDB ID 1, PDB ID 2) names only, intentionally ignoring target token.
    pair_index: Dict[Tuple[str, str], List[TRigidPairDef]] = {}
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row_index, row in enumerate(reader, start=2):
            complex_raw = (row.get("Complex") or row.get("complex") or "").strip()
            if not complex_raw:
                continue
            try:
                bound_pdb4, rec_group, lig_group = parse_bound_complex_field(complex_raw)
            except ValueError:
                continue

            pdb1 = parse_pdb_token(row.get("PDB ID 1", ""))
            pdb2 = parse_pdb_token(row.get("PDB ID 2", ""))
            if not pdb1.pdb4 or not pdb2.pdb4:
                continue

            rec = TRigidPairDef(
                row_index=row_index,
                complex_raw=complex_raw,
                bound_pdb4=bound_pdb4,
                bound_receptor_group=rec_group,
                bound_ligand_group=lig_group,
                ref_receptor_chain=rec_group[0] if rec_group else "",
                ref_ligand_chain=lig_group[0] if lig_group else "",
                pdb1=pdb1,
                pdb2=pdb2,
                pair_target_token="",
            )
            key = _unordered_pair_key(pdb1.pdb4, pdb2.pdb4)
            pair_index.setdefault(key, []).append(rec)
    return pair_index


def iter_prism_predictions(rosetta_root: Path) -> Iterable[PredictionRecord]:
    for path in sorted(rosetta_root.rglob("*.pdb")):
        if path.name.endswith(".pdb.intRes.txt"):
            continue
        m = ROSETTA_PDB_RE.match(path.name)
        if not m:
            continue
        yield PredictionRecord(
            model_path=path,
            joblist=path.parent.name,
            filename=path.name,
            target_token=m.group("target").lower(),
            tpl1=parse_pdb_token(m.group("tpl1")),
            tpl2=parse_pdb_token(m.group("tpl2")),
            tpl1_idx=int(m.group("idx1")),
            tpl2_idx=int(m.group("idx2")),
            rosetta_tag=m.group("rest"),
        )


def collect_t_rigid_bound_complex_ids(csv_path: Path) -> List[str]:
    ids = []
    seen = set()
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            complex_raw = (row.get("Complex") or row.get("complex") or "").strip()
            if not complex_raw:
                continue
            pdb4 = complex_raw[:4].lower()
            if len(pdb4) == 4 and pdb4 not in seen:
                seen.add(pdb4)
                ids.append(pdb4)
    return ids


def download_bound_pdb(pdb4: str, out_dir: Path, overwrite: bool = False) -> Tuple[bool, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdb = out_dir / (pdb4.lower() + ".pdb")
    if out_pdb.exists() and out_pdb.stat().st_size > 0 and not overwrite:
        return True, "exists"

    tmp_gz = out_dir / (pdb4.lower() + ".ent.gz")
    # Try PDBj mirror first (same source pattern used elsewhere in repo).
    urls = [
        "https://files.pdbj.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz" % pdb4.lower(),
        "https://files.rcsb.org/download/%s.pdb" % pdb4.upper(),
    ]

    last_err = None
    for url in urls:
        try:
            resp = urllib.request.urlopen(url, timeout=45)
            data = resp.read()
            if url.endswith(".ent.gz"):
                with tmp_gz.open("wb") as fh:
                    fh.write(data)
                with gzip.open(str(tmp_gz), "rb") as f_in, out_pdb.open("wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
                if tmp_gz.exists():
                    tmp_gz.unlink()
            else:
                with out_pdb.open("wb") as fh:
                    fh.write(data)

            if not out_pdb.exists() or out_pdb.stat().st_size == 0:
                raise RuntimeError("downloaded file empty")
            return True, "downloaded"
        except Exception as exc:
            last_err = str(exc)
            try:
                if tmp_gz.exists():
                    tmp_gz.unlink()
            except Exception:
                pass
            try:
                if out_pdb.exists() and out_pdb.stat().st_size == 0:
                    out_pdb.unlink()
            except Exception:
                pass
            continue
    return False, last_err or "download_failed"


def ensure_native_bound_complexes(csv_path: Path, native_dir: Path, verbose: bool = True) -> Dict[str, str]:
    statuses: Dict[str, str] = {}
    for pdb4 in collect_t_rigid_bound_complex_ids(csv_path):
        ok, status = download_bound_pdb(pdb4, native_dir, overwrite=False)
        statuses[pdb4] = status if ok else ("error:" + status)
        if verbose:
            print("[NATIVE] %s -> %s" % (pdb4, statuses[pdb4]))
    return statuses


def choose_pair_definitions(candidates: Sequence[TRigidPairDef]) -> Tuple[List[TRigidPairDef], str]:
    if not candidates:
        return [], "no_t_rigid_match"
    if len(candidates) == 1:
        return [candidates[0]], "unique_pdbid12_pair"
    return list(candidates), "multiple_complex_rows_same_pdbid12_pair"


def _chains_compatible_with_csv(model_or_tpl_chains: str, csv_chains: str) -> bool:
    # No chain listed in CSV means whole PDB / unspecified chain -> accept any.
    if not csv_chains:
        return True
    # Missing chain in filename (e.g. 1fgn_0) handled later via model-chain inference; do not reject here.
    if not model_or_tpl_chains:
        return True
    return set(model_or_tpl_chains).issubset(set(csv_chains))


def determine_orientation(pred: PredictionRecord, pair_def: TRigidPairDef) -> str:
    normal = pred.tpl1.pdb4 == pair_def.pdb1.pdb4 and pred.tpl2.pdb4 == pair_def.pdb2.pdb4
    swapped = pred.tpl1.pdb4 == pair_def.pdb2.pdb4 and pred.tpl2.pdb4 == pair_def.pdb1.pdb4
    if normal and not swapped:
        return "normal"
    if swapped and not normal:
        return "swapped"
    # Fallback if names are missing/unexpected.
    return "normal_assumed"


def read_model_chain_order(model_pdb: Path) -> str:
    seen = set()
    order: List[str] = []
    with model_pdb.open() as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
                ch = line[21]
                if ch not in seen:
                    seen.add(ch)
                    order.append(ch)
    return "".join(order)


def read_pdb_chain_set(pdb_path: Path) -> set:
    chains = set()
    with pdb_path.open() as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
                chains.add(line[21])
    return chains


def infer_model_partner_chains(
    pred: PredictionRecord,
    pair_def: TRigidPairDef,
    orientation: str,
) -> Tuple[Optional[str], Optional[str], str]:
    # Preferred source: template chain suffixes in filename.
    if orientation.startswith("normal"):
        tpl_r, tpl_l = pred.tpl1.chains, pred.tpl2.chains
    else:
        tpl_r, tpl_l = pred.tpl2.chains, pred.tpl1.chains

    if tpl_r and tpl_l:
        return tpl_r, tpl_l, "from_filename_template_chains"

    # Fallback when filenames omit template chains (e.g. 1fgn_0_1tfh_0).
    model_chain_order = read_model_chain_order(pred.model_path)

    # For T_Rigid pair-level targets, native interface is chain-vs-chain, so default 1+1.
    n_r = len(tpl_r) if tpl_r else 1
    n_l = len(tpl_l) if tpl_l else 1
    expected = n_r + n_l
    if len(model_chain_order) < expected:
        return None, None, f"cannot_infer_model_chains_only_{len(model_chain_order)}_chains_found_expected_{expected}"

    # Assume PRISM writes partner1 chains first then partner2 chains.
    first = model_chain_order[:n_r]
    second = model_chain_order[n_r:n_r + n_l]
    if orientation.startswith("normal"):
        return first, second, "inferred_from_model_chain_order"
    return second, first, "inferred_from_model_chain_order_swapped"


def model_chains_match_csv_partner_def(model_r: Optional[str], model_l: Optional[str], pair_def: TRigidPairDef, orientation: str) -> bool:
    if not model_r or not model_l:
        return False
    if orientation.startswith("normal"):
        return (
            _chains_compatible_with_csv(model_r, pair_def.pdb1.chains)
            and _chains_compatible_with_csv(model_l, pair_def.pdb2.chains)
        )
    return (
        _chains_compatible_with_csv(model_r, pair_def.pdb2.chains)
        and _chains_compatible_with_csv(model_l, pair_def.pdb1.chains)
    )


def compact_error_list(err_items: List[str], limit: int = 4) -> str:
    if not err_items:
        return ""
    if len(err_items) <= limit:
        return " | ".join(err_items)
    return " | ".join(err_items[:limit]) + " | ... (%d more)" % (len(err_items) - limit)


def load_model_fix_map(csv_path: Optional[Path]) -> Dict[str, dict]:
    if not csv_path:
        return {}
    if not csv_path.exists():
        return {}
    out = {}
    with csv_path.open(newline="") as fh:
        r = csv.DictReader(fh)
        for row in r:
            k = str(Path(row.get("original_model_path", "")).resolve())
            if not k:
                continue
            out[k] = row
    return out


def run_irmsd(
    irmsd_script: Path,
    python_exec: str,
    model_pdb: Path,
    model_r: str,
    model_l: str,
    native_pdb: Path,
    ref_r: str,
    ref_l: str,
    timeout_sec: int,
) -> Tuple[Optional[float], Optional[str]]:
    py_path = Path(python_exec)
    py_bin_dir = str((Path.cwd() / py_path).parent) if not py_path.is_absolute() else str(py_path.parent)
    cmd = [
        python_exec,
        str(irmsd_script),
        str(model_pdb),
        model_r,
        model_l,
        str(native_pdb),
        ref_r,
        ref_l,
    ]
    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=timeout_sec,
            env=dict(
                os.environ,
                PATH=py_bin_dir + os.pathsep + os.environ.get("PATH", ""),
            ),
        )
    except subprocess.TimeoutExpired:
        return None, f"irmsd_timeout_after_{timeout_sec}s"
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip().replace("\n", " | ")
        return None, f"irmsd_failed_exit_{proc.returncode}:{err}"
    out = (proc.stdout or "").strip().splitlines()
    if not out:
        return None, "irmsd_no_output"
    try:
        return float(out[-1].strip()), None
    except ValueError:
        return None, f"irmsd_parse_error:{out[-1].strip()}"


def run_dockq_script(
    dockq_script: Path,
    python_exec: str,
    model_pdb: Path,
    native_pdb: Path,
    mapping: str,
    timeout_sec: int,
) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    py_path = Path(python_exec)
    py_bin_dir = str((Path.cwd() / py_path).parent) if not py_path.is_absolute() else str(py_path.parent)
    cmd = [
        python_exec,
        str(dockq_script),
        str(model_pdb),
        str(native_pdb),
        "--mapping",
        mapping,
    ]
    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=timeout_sec,
            env=dict(
                os.environ,
                PATH=py_bin_dir + os.pathsep + os.environ.get("PATH", ""),
            ),
        )
    except subprocess.TimeoutExpired:
        return None, None, f"dockq_timeout_after_{timeout_sec}s"
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip().replace("\n", " | ")
        return None, None, "dockq_failed_exit_%s:%s" % (proc.returncode, err)
    dockq_val = None
    capri_val = None
    for line in (proc.stdout or "").splitlines():
        line = line.strip()
        if line.startswith("DockQ:"):
            try:
                dockq_val = float(line.split(":", 1)[1].strip())
            except Exception:
                pass
        elif line.startswith("CAPRI class:"):
            capri_val = line.split(":", 1)[1].strip()
    if dockq_val is None:
        return None, capri_val, "dockq_parse_error"
    return dockq_val, capri_val, None


def safe_float(x) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def summarize(values: List[float], best_mode: str) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if not values:
        return None, None, None
    mean_v = float(sum(values)) / float(len(values))
    var_v = statistics.pvariance(values) if len(values) > 1 else 0.0
    best_v = max(values) if best_mode == "max" else min(values)
    return mean_v, var_v, best_v


def write_csv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description="Analyze PRISM rosetta outputs against T_Rigid pair interfaces")
    parser.add_argument(
        "--rosetta-root",
        default="benchmark/prism_processed/rosetta_output_1",
        help="Root directory containing joblist_* PRISM output PDBs",
    )
    parser.add_argument(
        "--t-rigid-csv",
        default="benchmark/data/T_Rigid.csv",
        help="Path to benchmark T_Rigid CSV",
    )
    parser.add_argument(
        "--native-dir",
        default="benchmark/prism_processed_results/native_bound_complexes_t_rigid",
        help="Directory containing downloaded native bound complex PDBs (e.g. 1ahw.pdb).",
    )
    parser.add_argument(
        "--download-native-missing",
        action="store_true",
        help="Download all T_Rigid bound complex PDBs from public PDB mirrors into --native-dir before scoring.",
    )
    parser.add_argument(
        "--score-python",
        default=sys.executable,
        help="Python interpreter used to run benchmark/scripts/dockq.py and irmsd scripts",
    )
    parser.add_argument(
        "--out-dir",
        default="benchmark/prism_processed_results/prism_rigid_analysis",
        help="Output directory for analysis CSV files",
    )
    parser.add_argument(
        "--metrics",
        default="dockq,irmsd",
        help="Comma-separated metrics: dockq, irmsd, irmsd_backbone, none",
    )
    parser.add_argument(
        "--no-irmsd-reciprocal",
        action="store_true",
        help="Skip reciprocal iRMSD run (default computes both directions and stores min)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only first N predictions (for quick testing)",
    )
    parser.add_argument(
        "--verbose-every",
        type=int,
        default=200,
        help="Progress print frequency",
    )
    parser.add_argument(
        "--score-timeout-sec",
        type=int,
        default=120,
        help="Timeout (seconds) for each DockQ/iRMSD subprocess call",
    )
    parser.add_argument(
        "--model-fix-map-csv",
        default=None,
        help="Optional CSV from fix_model_chain_names.py with fixed model paths and corrected tpl1/tpl2 chains",
    )
    args = parser.parse_args()

    rosetta_root = Path(args.rosetta_root)
    t_rigid_csv = Path(args.t_rigid_csv)
    native_dir = Path(args.native_dir)
    out_dir = Path(args.out_dir)
    scripts_dir = Path(__file__).resolve().parent
    model_fix_map = load_model_fix_map(Path(args.model_fix_map_csv) if args.model_fix_map_csv else None)

    metrics = {m.strip() for m in args.metrics.split(",") if m.strip()}
    if "none" in metrics:
        metrics = set()
    valid_metrics = {"dockq", "irmsd", "irmsd_backbone"}
    unknown = metrics - valid_metrics
    if unknown:
        print(f"[ERROR] Unknown metrics: {sorted(unknown)}", file=sys.stderr)
        return 2

    pair_index = load_t_rigid_pairs(t_rigid_csv)
    predictions = list(iter_prism_predictions(rosetta_root))
    if args.limit is not None:
        predictions = predictions[: args.limit]

    print(f"[INFO] Loaded {sum(len(v) for v in pair_index.values())} T_Rigid pair definitions across {len(pair_index)} PDB-ID pairs")
    print(f"[INFO] Found {len(predictions)} PRISM prediction PDB files under {rosetta_root}")
    if args.download_native_missing:
        print("[INFO] Downloading/ensuring native bound complex PDBs in %s" % native_dir)
        ensure_native_bound_complexes(t_rigid_csv, native_dir, verbose=True)

    irmsd_script = scripts_dir / "irmsd.py"
    irmsd_bb_script = scripts_dir / "irmsd_backbone.py"
    dockq_script = scripts_dir / "dockq.py"

    detailed_rows: List[dict] = []
    pdb_chain_cache: Dict[str, set] = {}

    for i, pred in enumerate(predictions, start=1):
        pair_key = _unordered_pair_key(pred.tpl1.pdb4, pred.tpl2.pdb4)
        pair_candidates = pair_index.get(pair_key, [])
        pair_defs, pair_match_status = choose_pair_definitions(pair_candidates)

        base_row = {
            "model_path": str(pred.model_path),
            "model_path_used": str(pred.model_path),
            "joblist": pred.joblist,
            "filename": pred.filename,
            "target_token": pred.target_token,
            "tpl1_raw": pred.tpl1.raw,
            "tpl1_pdb4": pred.tpl1.pdb4,
            "tpl1_chains": pred.tpl1.chains,
            "tpl1_idx": pred.tpl1_idx,
            "tpl2_raw": pred.tpl2.raw,
            "tpl2_pdb4": pred.tpl2.pdb4,
            "tpl2_chains": pred.tpl2.chains,
            "tpl2_idx": pred.tpl2_idx,
            "rosetta_tag": pred.rosetta_tag,
            "pair_match_status": pair_match_status,
            "pdb_pair_key": "|".join(pair_key),
            "n_complex_rows_for_pdb_pair": len(pair_defs),
            "complex_row_rank": None,
            "t_rigid_row": None,
            "complex_raw": None,
            "bound_pdb4": None,
            "ref_receptor_chain": None,
            "ref_ligand_chain": None,
            "pdb1_raw": None,
            "pdb2_raw": None,
            "orientation": None,
            "model_receptor_chains": None,
            "model_ligand_chains": None,
            "model_chain_inference": None,
            "native_pairs_evaluated": None,
            "dockq_pair_scores": None,
            "irmsd_pair_scores_forward": None,
            "irmsd_pair_scores_reverse": None,
            "irmsd_backbone_pair_scores_forward": None,
            "irmsd_backbone_pair_scores_reverse": None,
            "native_pdb": None,
            "native_exists": None,
            "dockq": None,
            "dockq_capri": None,
            "dockq_error": None,
            "irmsd_forward": None,
            "irmsd_reverse": None,
            "irmsd_min": None,
            "irmsd_error": None,
            "irmsd_backbone_forward": None,
            "irmsd_backbone_reverse": None,
            "irmsd_backbone_min": None,
            "irmsd_backbone_error": None,
        }

        if not pair_defs:
            detailed_rows.append(base_row)
            continue

        for row_rank, pair_def in enumerate(pair_defs, start=1):
            row = dict(base_row)
            row["complex_row_rank"] = row_rank

            # Map model partner chains to CSV PDB ID 1 / PDB ID 2.
            # If model name order is swapped (e.g. ..._1tfhA_..._1fgnH_...), set swapped.
            orientation = determine_orientation(pred, pair_def)
            fixed = model_fix_map.get(str(pred.model_path.resolve()))
            model_pdb_for_scoring = pred.model_path
            if fixed:
                fpath = fixed.get("fixed_model_path", "").strip()
                if fpath:
                    model_pdb_for_scoring = Path(fpath)
                f_tpl1 = (fixed.get("fixed_tpl1_chain") or "").strip()
                f_tpl2 = (fixed.get("fixed_tpl2_chain") or "").strip()
                if f_tpl1 and f_tpl2:
                    if orientation.startswith("normal"):
                        model_r, model_l = f_tpl1, f_tpl2
                    else:
                        model_r, model_l = f_tpl2, f_tpl1
                    chain_inference = "from_model_fix_map_csv"
                else:
                    model_r, model_l, chain_inference = infer_model_partner_chains(pred, pair_def, orientation)
            else:
                model_r, model_l, chain_inference = infer_model_partner_chains(pred, pair_def, orientation)
            native_pdb = native_dir / f"{pair_def.bound_pdb4}.pdb"
            native_pairs = [(r, l) for r in pair_def.bound_receptor_group for l in pair_def.bound_ligand_group]

            row.update(
                {
                    "t_rigid_row": pair_def.row_index,
                    "complex_raw": pair_def.complex_raw,
                    "bound_pdb4": pair_def.bound_pdb4,
                    "ref_receptor_chain": pair_def.ref_receptor_chain,
                    "ref_ligand_chain": pair_def.ref_ligand_chain,
                    "pdb1_raw": pair_def.pdb1.raw,
                    "pdb2_raw": pair_def.pdb2.raw,
                    "orientation": orientation,
                    "model_receptor_chains": model_r,
                    "model_ligand_chains": model_l,
                    "model_chain_inference": chain_inference,
                    "model_path_used": str(model_pdb_for_scoring),
                    "native_pairs_evaluated": ",".join("%s%s" % (r, l) for r, l in native_pairs),
                    "native_pdb": str(native_pdb),
                    "native_exists": native_pdb.exists(),
                }
            )

            if not model_r or not model_l:
                row["dockq_error"] = row["dockq_error"] or "missing_model_chain_assignment"
                row["irmsd_error"] = row["irmsd_error"] or "missing_model_chain_assignment"
                row["irmsd_backbone_error"] = row["irmsd_backbone_error"] or "missing_model_chain_assignment"
                detailed_rows.append(row)
                continue

            if not native_pdb.exists():
                msg = "native_bound_complex_not_found"
                row["dockq_error"] = row["dockq_error"] or msg
                row["irmsd_error"] = row["irmsd_error"] or msg
                row["irmsd_backbone_error"] = row["irmsd_backbone_error"] or msg
                detailed_rows.append(row)
                continue

            model_chain_key = str(model_pdb_for_scoring.resolve())
            native_chain_key = str(native_pdb.resolve())
            model_chain_set = pdb_chain_cache.get(model_chain_key)
            if model_chain_set is None:
                model_chain_set = read_pdb_chain_set(model_pdb_for_scoring)
                pdb_chain_cache[model_chain_key] = model_chain_set
            native_chain_set = pdb_chain_cache.get(native_chain_key)
            if native_chain_set is None:
                native_chain_set = read_pdb_chain_set(native_pdb)
                pdb_chain_cache[native_chain_key] = native_chain_set

            missing_model = [c for c in set((model_r or "") + (model_l or "")) if c not in model_chain_set]
            if missing_model:
                msg = "model_chain_missing_in_pdb:%s" % ",".join(sorted(missing_model))
                row["dockq_error"] = row["dockq_error"] or msg
                row["irmsd_error"] = row["irmsd_error"] or msg
                row["irmsd_backbone_error"] = row["irmsd_backbone_error"] or msg
                detailed_rows.append(row)
                continue

            valid_native_pairs = []
            missing_native_pairs = []
            for ref_r, ref_l in native_pairs:
                if ref_r in native_chain_set and ref_l in native_chain_set:
                    valid_native_pairs.append((ref_r, ref_l))
                else:
                    missing_native_pairs.append((ref_r, ref_l))
            if not valid_native_pairs:
                miss = ",".join("%s%s" % (r, l) for r, l in missing_native_pairs)
                msg = "native_chain_missing_in_pdb:%s" % miss
                row["dockq_error"] = row["dockq_error"] or msg
                row["irmsd_error"] = row["irmsd_error"] or msg
                row["irmsd_backbone_error"] = row["irmsd_backbone_error"] or msg
                detailed_rows.append(row)
                continue

            if "dockq" in metrics:
                dockq_errors = []
                dockq_pair_scores = []
                best_dockq = None
                best_capri = None
                for ref_r, ref_l in valid_native_pairs:
                    mapping = "%s%s:%s%s" % (model_r, model_l, ref_r, ref_l)
                    dval, dcapri, derr = run_dockq_script(
                        dockq_script,
                        args.score_python,
                    model_pdb_for_scoring,
                    native_pdb,
                    mapping,
                    args.score_timeout_sec,
                )
                    if dval is not None:
                        dockq_pair_scores.append("%s%s=%.6f" % (ref_r, ref_l, dval))
                        if best_dockq is None or dval > best_dockq:
                            best_dockq = dval
                            best_capri = dcapri
                    elif derr:
                        dockq_errors.append("%s%s:%s" % (ref_r, ref_l, derr))
                row["dockq_pair_scores"] = ";".join(dockq_pair_scores) if dockq_pair_scores else None
                row["dockq"] = best_dockq
                row["dockq_capri"] = best_capri
                if dockq_errors and best_dockq is None:
                    row["dockq_error"] = compact_error_list(dockq_errors)

            if "irmsd" in metrics:
                irmsd_errors = []
                pair_fwd = []
                pair_rev = []
                best_fwd = None
                best_rev = None
                for ref_r, ref_l in valid_native_pairs:
                    fwd, err = run_irmsd(
                        irmsd_script,
                        args.score_python,
                        model_pdb_for_scoring,
                        model_r,
                        model_l,
                        native_pdb,
                        ref_r,
                        ref_l,
                        args.score_timeout_sec,
                    )
                    if fwd is not None:
                        pair_fwd.append("%s%s=%.3f" % (ref_r, ref_l, fwd))
                        if best_fwd is None or fwd < best_fwd:
                            best_fwd = fwd
                    elif err:
                        irmsd_errors.append("%s%s:fwd:%s" % (ref_r, ref_l, err))

                    if not args.no_irmsd_reciprocal:
                        rev, rev_err = run_irmsd(
                            irmsd_script,
                            args.score_python,
                            model_pdb_for_scoring,
                            model_l,
                            model_r,
                            native_pdb,
                            ref_l,
                            ref_r,
                            args.score_timeout_sec,
                        )
                        if rev is not None:
                            pair_rev.append("%s%s=%.3f" % (ref_r, ref_l, rev))
                            if best_rev is None or rev < best_rev:
                                best_rev = rev
                        elif rev_err:
                            irmsd_errors.append("%s%s:rev:%s" % (ref_r, ref_l, rev_err))

                row["irmsd_pair_scores_forward"] = ";".join(pair_fwd) if pair_fwd else None
                row["irmsd_pair_scores_reverse"] = ";".join(pair_rev) if pair_rev else None
                row["irmsd_forward"] = best_fwd
                row["irmsd_reverse"] = best_rev
                vals = [v for v in [best_fwd, best_rev] if isinstance(v, (int, float))]
                if vals:
                    row["irmsd_min"] = min(vals)
                if irmsd_errors and not vals:
                    row["irmsd_error"] = compact_error_list(irmsd_errors)

            if "irmsd_backbone" in metrics:
                bb_errors = []
                pair_fwd = []
                pair_rev = []
                best_fwd = None
                best_rev = None
                for ref_r, ref_l in valid_native_pairs:
                    fwd, err = run_irmsd(
                        irmsd_bb_script,
                        args.score_python,
                        model_pdb_for_scoring,
                        model_r,
                        model_l,
                        native_pdb,
                        ref_r,
                        ref_l,
                        args.score_timeout_sec,
                    )
                    if fwd is not None:
                        pair_fwd.append("%s%s=%.3f" % (ref_r, ref_l, fwd))
                        if best_fwd is None or fwd < best_fwd:
                            best_fwd = fwd
                    elif err:
                        bb_errors.append("%s%s:fwd:%s" % (ref_r, ref_l, err))

                    if not args.no_irmsd_reciprocal:
                        rev, rev_err = run_irmsd(
                            irmsd_bb_script,
                            args.score_python,
                            model_pdb_for_scoring,
                            model_l,
                            model_r,
                            native_pdb,
                            ref_l,
                            ref_r,
                            args.score_timeout_sec,
                        )
                        if rev is not None:
                            pair_rev.append("%s%s=%.3f" % (ref_r, ref_l, rev))
                            if best_rev is None or rev < best_rev:
                                best_rev = rev
                        elif rev_err:
                            bb_errors.append("%s%s:rev:%s" % (ref_r, ref_l, rev_err))

                row["irmsd_backbone_pair_scores_forward"] = ";".join(pair_fwd) if pair_fwd else None
                row["irmsd_backbone_pair_scores_reverse"] = ";".join(pair_rev) if pair_rev else None
                row["irmsd_backbone_forward"] = best_fwd
                row["irmsd_backbone_reverse"] = best_rev
                vals = [v for v in [best_fwd, best_rev] if isinstance(v, (int, float))]
                if vals:
                    row["irmsd_backbone_min"] = min(vals)
                if bb_errors and not vals:
                    row["irmsd_backbone_error"] = compact_error_list(bb_errors)

            detailed_rows.append(row)

        if args.verbose_every and i % args.verbose_every == 0:
            print(f"[INFO] Processed {i}/{len(predictions)} predictions")

    # Per-prediction output
    detailed_fieldnames = list(detailed_rows[0].keys()) if detailed_rows else []
    detailed_csv = out_dir / "per_prediction.csv"
    write_csv(detailed_csv, detailed_rows, detailed_fieldnames)

    # Pair-level aggregation (across template combinations / joblists / predictions)
    grouped: Dict[Tuple[str, str, str, str, str], List[dict]] = {}
    for row in detailed_rows:
        if not row.get("complex_raw"):
            continue
        key = (
            str(row.get("pdb_pair_key") or ""),
            str(row.get("pdb1_raw") or ""),
            str(row.get("pdb2_raw") or ""),
            str(row.get("complex_raw") or ""),
            str(row.get("bound_pdb4") or ""),
        )
        grouped.setdefault(key, []).append(row)

    summary_rows: List[dict] = []
    for (pdb_pair_key, pdb1_raw, pdb2_raw, complex_raw, bound_pdb4), rows in sorted(grouped.items()):
        scored_dockq = [r["dockq"] for r in rows if isinstance(r.get("dockq"), (int, float))]
        scored_irmsd = [r["irmsd_min"] for r in rows if isinstance(r.get("irmsd_min"), (int, float))]
        scored_irmsd_bb = [r["irmsd_backbone_min"] for r in rows if isinstance(r.get("irmsd_backbone_min"), (int, float))]

        dockq_mean, dockq_var, dockq_best = summarize(scored_dockq, "max")
        irmsd_mean, irmsd_var, irmsd_best = summarize(scored_irmsd, "min")
        irmsd_bb_mean, irmsd_bb_var, irmsd_bb_best = summarize(scored_irmsd_bb, "min")

        best_dockq_row = max((r for r in rows if isinstance(r.get("dockq"), (int, float))), key=lambda r: r["dockq"], default=None)
        best_irmsd_row = min((r for r in rows if isinstance(r.get("irmsd_min"), (int, float))), key=lambda r: r["irmsd_min"], default=None)
        best_irmsd_bb_row = min((r for r in rows if isinstance(r.get("irmsd_backbone_min"), (int, float))), key=lambda r: r["irmsd_backbone_min"], default=None)

        summary_rows.append(
            {
                "pdb_pair_key": pdb_pair_key,
                "pdb1_raw": pdb1_raw,
                "pdb2_raw": pdb2_raw,
                "complex_raw": complex_raw,
                "bound_pdb4": bound_pdb4,
                "n_predictions_total": len(rows),
                "n_predictions_mapped": sum(1 for r in rows if r.get("pair_match_status") not in ("no_t_rigid_match", None)),
                "n_native_missing": sum(1 for r in rows if r.get("native_exists") is False),
                "n_dockq_scored": len(scored_dockq),
                "dockq_mean": dockq_mean,
                "dockq_variance": dockq_var,
                "dockq_best": dockq_best,
                "dockq_best_file": best_dockq_row.get("model_path") if best_dockq_row else None,
                "dockq_best_templates": f"{best_dockq_row.get('tpl1_raw')}|{best_dockq_row.get('tpl2_raw')}" if best_dockq_row else None,
                "n_irmsd_scored": len(scored_irmsd),
                "irmsd_mean": irmsd_mean,
                "irmsd_variance": irmsd_var,
                "irmsd_best_min": irmsd_best,
                "irmsd_best_file": best_irmsd_row.get("model_path") if best_irmsd_row else None,
                "irmsd_best_templates": f"{best_irmsd_row.get('tpl1_raw')}|{best_irmsd_row.get('tpl2_raw')}" if best_irmsd_row else None,
                "n_irmsd_backbone_scored": len(scored_irmsd_bb),
                "irmsd_backbone_mean": irmsd_bb_mean,
                "irmsd_backbone_variance": irmsd_bb_var,
                "irmsd_backbone_best_min": irmsd_bb_best,
                "irmsd_backbone_best_file": best_irmsd_bb_row.get("model_path") if best_irmsd_bb_row else None,
                "irmsd_backbone_best_templates": f"{best_irmsd_bb_row.get('tpl1_raw')}|{best_irmsd_bb_row.get('tpl2_raw')}" if best_irmsd_bb_row else None,
                "example_ref_pair": f"{bound_pdb4}:{rows[0].get('ref_receptor_chain')}:{rows[0].get('ref_ligand_chain')}" if rows else None,
            }
        )

    summary_fieldnames = list(summary_rows[0].keys()) if summary_rows else []
    summary_csv = out_dir / "pair_summary.csv"
    write_csv(summary_csv, summary_rows, summary_fieldnames)

    print(f"[INFO] Wrote {detailed_csv}")
    print(f"[INFO] Wrote {summary_csv}")

    no_match = sum(1 for r in detailed_rows if r.get("pair_match_status") == "no_t_rigid_match")
    native_missing = sum(1 for r in detailed_rows if r.get("native_exists") is False)
    print(f"[INFO] Unmatched PRISM predictions vs T_Rigid PDB-ID pairs: {no_match}")
    print(f"[INFO] Predictions with missing native bound PDBs: {native_missing}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
