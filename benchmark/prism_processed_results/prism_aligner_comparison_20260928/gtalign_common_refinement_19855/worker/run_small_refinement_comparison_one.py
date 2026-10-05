#!/usr/bin/env python3
"""Run one selected transformed pair through FiberDock, Rosetta, and DockQ.

The worker is intentionally one-candidate-per-Slurm-task.  Each stage has an
atomic checkpoint and an event record so an interrupted task can be resumed
without recomputing completed stages or changing the selected input files.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def append_event(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            handle.write(json.dumps(payload, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def chain_ids(path: Path) -> list[str]:
    result: list[str] = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                chain = line[21].strip() or "_"
                if chain not in result:
                    result.append(chain)
    if not result:
        raise ValueError(f"no atom chains in {path}")
    return result


def _line_residue_number(line: str) -> int | None:
    try:
        return int(line[22:26].strip())
    except (TypeError, ValueError):
        return None


def _input_segments(path: Path, expected_count: int) -> list[int]:
    """Return one segment number per atom record in a transformed partner."""
    records: list[tuple[str, int | None, str]] = []
    explicit: list[str] = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")) or len(line) <= 21:
                continue
            chain = line[21].strip()
            if chain and chain not in explicit:
                explicit.append(chain)
            records.append((chain, _line_residue_number(line), line))
    if not records:
        raise ValueError(f"no atom records in transformed input {path}")
    if explicit:
        if len(explicit) != expected_count:
            raise ValueError(f"input {path} has chains {explicit}, expected {expected_count}")
        return [explicit.index(chain) for chain, _, _ in records]
    if expected_count == 1:
        return [0] * len(records)
    segments: list[int] = []
    segment = 0
    previous_residue: int | None = None
    for _, residue, _ in records:
        if residue is not None and previous_residue is not None and residue < previous_residue:
            segment += 1
        segments.append(segment)
        if residue is not None:
            previous_residue = residue
    if segment + 1 != expected_count:
        raise ValueError(
            f"cannot recover {expected_count} chain segments from chainless input {path}; "
            f"observed {segment + 1} residue-number segments"
        )
    return segments


def normalize_partner(path: Path, output: Path, desired_chains: str) -> dict[str, object]:
    """Give a chainless transformed partner deterministic PDB chain IDs."""
    desired_chains = "".join(desired_chains.split())
    if not desired_chains or len(set(desired_chains)) != len(desired_chains):
        raise ValueError(f"invalid desired chain group {desired_chains!r}")
    segments = _input_segments(path, len(desired_chains))
    output.parent.mkdir(parents=True, exist_ok=True)
    atom_index = 0
    previous_segment: int | None = None
    with path.open(encoding="ascii", errors="replace") as source, output.open("w", encoding="ascii") as target:
        for line in source:
            if not line.startswith(("ATOM", "HETATM")) or len(line) <= 21:
                continue
            segment = segments[atom_index]
            if previous_segment is not None and segment != previous_segment:
                target.write("TER\n")
            padded = line.rstrip("\n").ljust(80)
            target.write(padded[:21] + desired_chains[segment] + padded[22:] + "\n")
            previous_segment = segment
            atom_index += 1
        target.write("TER\nEND\n")
    if atom_index == 0 or not output.is_file() or output.stat().st_size == 0:
        raise ValueError(f"normalized input is empty: {output}")
    return {
        "source": str(path.resolve()),
        "source_sha256": sha256(path),
        "normalized": str(output.resolve()),
        "normalized_sha256": sha256(output),
        "chains": desired_chains,
        "atom_records": atom_index,
    }


def input_chain_specs(path: Path) -> list[tuple[str, int]]:
    """Return normalized chain IDs and residue counts in file order."""
    counts: dict[str, int] = {}
    seen: dict[str, set[tuple[str, str, str]]] = {}
    order: list[str] = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM") or len(line) <= 21:
                continue
            chain = line[21].strip() or "_"
            key = (line[17:20].strip(), line[22:26].strip(), line[26].strip())
            if chain not in seen:
                seen[chain] = set()
                order.append(chain)
            if key not in seen[chain]:
                seen[chain].add(key)
                counts[chain] = counts.get(chain, 0) + 1
    if not order:
        raise ValueError(f"no ATOM chains in normalized input {path}")
    return [(chain, counts[chain]) for chain in order]


def prepare_candidate_inputs(row: dict[str, str], candidate_root: Path) -> dict[str, object]:
    """Normalize chainless transformed files for both refinement adapters."""
    native_left = "".join(row["native_receptor_chains"].split())
    native_right = "".join(row["native_ligand_chains"].split())
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    total = len(native_left) + len(native_right)
    if total < 2 or total > len(alphabet):
        raise ValueError(f"unsupported native chain groups: {native_left!r}, {native_right!r}")
    left_chains = alphabet[: len(native_left)]
    right_chains = alphabet[len(native_left) : total]
    input_dir = candidate_root / "inputs"
    basename = f"normalized_srcx{left_chains}_srcy{right_chains}_o1"
    left_out = input_dir / f"{basename}_L.pdb"
    right_out = input_dir / f"{basename}_R.pdb"
    left_info = normalize_partner(Path(row["left"]).resolve(), left_out, left_chains)
    right_info = normalize_partner(Path(row["right"]).resolve(), right_out, right_chains)
    return {
        "prepared_left": str(left_out.resolve()),
        "prepared_right": str(right_out.resolve()),
        "prepared_left_chains": left_chains,
        "prepared_right_chains": right_chains,
        "normalization_left": left_info,
        "normalization_right": right_info,
    }


def model_chain_tokens(path: Path) -> tuple[list[list[str]], int]:
    """Return TER-separated chain occurrences and the number of atom records."""
    blocks: list[list[str]] = []
    current: list[str] = []
    atom_count = 0
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                atom_count += 1
                chain = line[21].strip() or "_"
                if chain not in current:
                    current.append(chain)
            if line.startswith("TER"):
                if current:
                    blocks.append(current)
                    current = []
        if current:
            blocks.append(current)
    return blocks, atom_count


def canonicalize_model(
    model: Path,
    left: Path,
    right: Path,
    native_receptor: str,
    native_ligand: str,
    output: Path,
) -> dict[str, object]:
    """Rewrite a refined model to a unique, explicit receptor/ligand mapping."""
    left_specs = input_chain_specs(left)
    right_specs = input_chain_specs(right)
    left_source = [chain for chain, _ in left_specs]
    right_source = [chain for chain, _ in right_specs]
    native_receptor = "".join(native_receptor.split())
    native_ligand = "".join(native_ligand.split())
    if len(left_source) != len(native_receptor) or len(right_source) != len(native_ligand):
        raise ValueError(
            "transformed/native chain-count mismatch: "
            f"left {left_source}/{native_receptor}, right {right_source}/{native_ligand}"
        )
    blocks, atom_count = model_chain_tokens(model)
    expected_count = len(left_source) + len(right_source)
    tokens = [(block_index, chain) for block_index, block in enumerate(blocks) for chain in block]
    canonical = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    model_all = "".join(canonical[:expected_count])
    native_all = native_receptor + native_ligand
    direct_mapping = len(tokens) == expected_count and not (
        len(tokens) == 1 and tokens[0][1] == "_"
    )
    token_map = {token: canonical[index] for index, token in enumerate(tokens)} if direct_mapping else {}
    residue_bounds: list[tuple[int, str]] = []
    if not direct_mapping:
        residue_offset = 0
        for chain, count in [*left_specs, *right_specs]:
            residue_bounds.append((residue_offset + count, chain))
            residue_offset += count
        observed_residues = 0
        previous_residue: tuple[str, str, str] | None = None
        with model.open(encoding="ascii", errors="replace") as handle:
            for line in handle:
                if not line.startswith("ATOM"):
                    continue
                residue = (line[17:20].strip(), line[22:26].strip(), line[26].strip())
                if residue != previous_residue:
                    observed_residues += 1
                    previous_residue = residue
        expected_residues = sum(count for _, count in [*left_specs, *right_specs])
        if observed_residues < expected_residues:
            raise ValueError(
                f"refined model has {observed_residues} ATOM residues; expected at least {expected_residues}"
            )
    output.parent.mkdir(parents=True, exist_ok=True)
    block_index = 0
    residue_index = -1
    previous_residue = None
    fallback_chain = canonical[expected_count - 1]
    with model.open(encoding="ascii", errors="replace") as source, output.open("w", encoding="ascii") as target:
        for line in source:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                old_chain = line[21].strip() or "_"
                if direct_mapping:
                    replacement = token_map.get((block_index, old_chain))
                    if replacement is None:
                        raise ValueError(f"unmapped refined-model chain occurrence {(block_index, old_chain)}")
                elif line.startswith("ATOM"):
                    residue = (line[17:20].strip(), line[22:26].strip(), line[26].strip())
                    if residue != previous_residue:
                        residue_index += 1
                        previous_residue = residue
                    replacement = canonical[expected_count - 1]
                    for chain_index, (bound, _) in enumerate(residue_bounds):
                        if residue_index < bound:
                            replacement = canonical[chain_index]
                            break
                    fallback_chain = replacement
                else:
                    replacement = fallback_chain
                target.write(line[:21] + replacement + line[22:])
            else:
                target.write(line)
            if line.startswith("TER"):
                block_index += 1
    if atom_count == 0 or not output.is_file() or output.stat().st_size == 0:
        raise ValueError(f"canonical model is empty: {output}")
    return {
        "model_pdb": str(output.resolve()),
        "model_sha256": sha256(output),
        "model_receptor_chains": model_all[: len(left_source)],
        "model_ligand_chains": model_all[len(left_source) :],
        "native_receptor_chains": native_receptor,
        "native_ligand_chains": native_ligand,
        "mapping": f"{model_all}:{native_all}",
        "source_chain_groups": {"left": left_source, "right": right_source},
        "model_chain_blocks": blocks,
        "model_mapping_mode": "direct_chain_blocks" if direct_mapping else "residue_order_fallback",
    }


def stage_record(status: str, started_at: str, started: float, **values: object) -> dict[str, object]:
    return {
        "status": status,
        "started_at": started_at,
        "finished_at": utc_now(),
        "elapsed_seconds": time.perf_counter() - started,
        **values,
    }


def run_fiberdock(row: dict[str, str], candidate_root: Path, repo: Path, fiberdock_dir: Path) -> dict[str, object]:
    work_root = candidate_root / "fiberdock"
    pair_name = "selected_pair"
    pair_work = work_root / "processed" / "fiberdock_refinement" / pair_name
    work_root.mkdir(parents=True, exist_ok=True)
    pair_work.mkdir(parents=True, exist_ok=True)
    left = Path(row["prepared_left"]).resolve()
    right = Path(row["prepared_right"]).resolve()
    previous_cwd = Path.cwd()
    started = time.perf_counter()
    started_at = utc_now()
    try:
        os.chdir(work_root)
        os.environ["PRISM_FIBERDOCK_DIR"] = str(fiberdock_dir)
        sys.path.insert(0, str(repo))
        from src import fiberdock_refinement as fd

        fd.FIBERDOCK_DIR = str(fiberdock_dir)
        # FiberDock itself writes into its current directory.  The bundle is
        # intentionally shared read-only across tasks, so every process must
        # use a distinct output prefix or concurrent tasks can move/parse one
        # another's ``*.ref`` and ``*.ref.pdb`` files.
        fd.FIBERDOCK_OUTPUT_PREFIX = "fd_" + hashlib.sha256(
            str(candidate_root.resolve()).encode("utf-8")
        ).hexdigest()[:20]
        r_hb = fd._add_hydrogens(str(right), str(pair_work))
        l_hb = fd._add_hydrogens(str(left), str(pair_work))
        if not r_hb or not l_hb:
            raise RuntimeError("FiberDock hydrogenation failed")
        r_ca, r_sizes = fd._create_ca_pdb(str(right), str(pair_work))
        l_ca, l_sizes = fd._create_ca_pdb(str(left), str(pair_work))
        r_nma = fd._run_nma(r_ca, str(pair_work))
        l_nma = fd._run_nma(l_ca, str(pair_work))
        if not r_nma or not l_nma:
            raise RuntimeError("FiberDock NMA failed")
        if sum(r_sizes.values()) >= sum(l_sizes.values()):
            receptor_base, ligand_base = right.stem, left.stem
        else:
            receptor_base, ligand_base = left.stem, right.stem
            r_hb, l_hb = l_hb, r_hb
        params = fd._build_fiberdock_params(r_hb, l_hb, str(pair_work), receptor_base, ligand_base)
        if not Path(params).is_file():
            raise RuntimeError(f"FiberDock parameter file missing: {params}")
        energy = fd._run_fiberdock(
            params, str(fiberdock_dir), str(pair_work), receptor_base, ligand_base, pair_name=pair_name
        )
        refined = pair_work / f"{fd.FIBERDOCK_OUTPUT_PREFIX}_1.ref.pdb"
        if not refined.is_file() or refined.stat().st_size == 0:
            raise RuntimeError("FiberDock completed without a refined PDB")
        values = {
            "refined_model": str(refined.resolve()),
            "refined_model_sha256": sha256(refined),
            "energy": energy,
            "output_prefix": fd.FIBERDOCK_OUTPUT_PREFIX,
            "params": str(Path(params).resolve()),
            "receptor_ca_residues": sum(r_sizes.values()),
            "ligand_ca_residues": sum(l_sizes.values()),
        }
        return stage_record("completed", started_at, started, **values)
    except Exception as exc:
        return stage_record("failed", started_at, started, error=f"{type(exc).__name__}: {exc}")
    finally:
        os.chdir(previous_cwd)


def run_rosetta(
    row: dict[str, str],
    candidate_root: Path,
    repo: Path,
    rosetta_prepack: str,
    rosetta_dock: str,
    rosetta_db: str,
) -> dict[str, object]:
    work_root = candidate_root / "rosetta"
    work_root.mkdir(parents=True, exist_ok=True)
    previous_cwd = Path.cwd()
    started = time.perf_counter()
    started_at = utc_now()
    try:
        os.chdir(work_root)
        os.environ.update(
            PRISM_ROSETTA_PREPACK=rosetta_prepack,
            PRISM_ROSETTA_DOCK=rosetta_dock,
            PRISM_ROSETTA_DB=rosetta_db,
        )
        sys.path.insert(0, str(repo))
        from src.rosetta_refinement import calculate_energy

        totalscore, interaction_score, structure = calculate_energy(row["prepared_left"], row["prepared_right"])
        values: dict[str, object] = {
            "totalscore": totalscore,
            "interaction_score": interaction_score,
            "structure_returned": structure,
        }
        if structure == "-":
            return stage_record("no_model", started_at, started, **values)
        model = Path(structure)
        if not model.is_absolute():
            model = (work_root / model).resolve()
        if not model.is_file() or model.stat().st_size == 0:
            return stage_record("no_model", started_at, started, **values, error="Rosetta returned a missing structure")
        values["refined_model"] = str(model)
        values["refined_model_sha256"] = sha256(model)
        return stage_record("completed", started_at, started, **values)
    except Exception as exc:
        return stage_record("failed", started_at, started, error=f"{type(exc).__name__}: {exc}")
    finally:
        os.chdir(previous_cwd)


def run_dockq(
    row: dict[str, str],
    candidate_root: Path,
    repo: Path,
    model_stage: dict[str, object],
    dockq_python: str,
    label: str,
) -> dict[str, object]:
    started = time.perf_counter()
    started_at = utc_now()
    if model_stage.get("status") != "completed" or not model_stage.get("refined_model"):
        return stage_record("not_run_no_model", started_at, started, reason="refiner did not produce a model")
    model = Path(str(model_stage["refined_model"])).resolve()
    native = Path(row["native_pdb"]).resolve()
    canonical = candidate_root / "dockq" / f"{label}_canonical.pdb"
    try:
        canonical_info = canonicalize_model(
            model,
            Path(row["prepared_left"]).resolve(),
            Path(row["prepared_right"]).resolve(),
            row["native_receptor_chains"],
            row["native_ligand_chains"],
            canonical,
        )
        os.environ["DOCKQ_PYTHON"] = dockq_python
        sys.path.insert(0, str(repo))
        from src.eval.dockq import calculate_dockq, dockq_to_capri_class
        from src.eval.irmsd_backbone import calculate_irmsd_backbone

        work_dir = candidate_root / "dockq" / label
        dockq = calculate_dockq(
            str(canonical), str(native), mapping=str(canonical_info["mapping"]), work_dir=str(work_dir)
        )
        raw_candidates = sorted(work_dir.glob("*.json"), key=lambda path: path.stat().st_mtime)
        raw_json = raw_candidates[-1] if raw_candidates else None
        irmsd = calculate_irmsd_backbone(
            str(canonical),
            str(canonical_info["model_receptor_chains"]),
            str(canonical_info["model_ligand_chains"]),
            str(native),
            row["native_receptor_chains"],
            row["native_ligand_chains"],
        )
        return stage_record(
            "scored",
            started_at,
            started,
            **canonical_info,
            native_pdb=str(native),
            native_sha256=sha256(native),
            dockq=dockq.get("dockq"),
            dockq_global=dockq.get("dockq_global"),
            dockq_sum=dockq.get("dockq_sum"),
            dockq_capri=dockq_to_capri_class(dockq.get("dockq")),
            fnat=dockq.get("fnat"),
            dockq_irmsd=dockq.get("irmsd"),
            dockq_lrmsd=dockq.get("lrmsd"),
            irmsd_backbone=irmsd,
            dockq_json=dockq.get("raw_dockq_json"),
            dockq_json_sha256=dockq.get("raw_dockq_json_sha256"),
            dockq_argv=dockq.get("dockq_argv"),
            raw_dockq_json=str(raw_json.resolve()) if raw_json else "",
            raw_dockq_json_sha256=sha256(raw_json) if raw_json else "",
        )
    except Exception as exc:
        return stage_record("score_failed", started_at, started, error=f"{type(exc).__name__}: {exc}")


def load_row(path: Path, index: int) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if index < 0 or index >= len(rows):
        raise IndexError(f"candidate index {index} outside 0..{len(rows) - 1}")
    row = rows[index]
    required = {
        "pipeline", "case_id", "left", "right", "native_pdb",
        "native_receptor_chains", "native_ligand_chains",
    }
    missing = sorted(key for key in required if not row.get(key))
    if missing:
        raise ValueError(f"selected candidate is missing required fields: {missing}")
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-root", required=True, type=Path)
    parser.add_argument("--selected-csv", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int)
    parser.add_argument("--pipeline-repo", required=True, type=Path)
    parser.add_argument("--fiberdock-dir", required=True, type=Path)
    parser.add_argument("--dockq-python", required=True)
    parser.add_argument("--rosetta-prepack", required=True)
    parser.add_argument("--rosetta-dock", required=True)
    parser.add_argument("--rosetta-db", required=True)
    args = parser.parse_args()

    comparison_root = args.comparison_root.resolve()
    selected_csv = args.selected_csv.resolve()
    row = load_row(selected_csv, args.index)
    identity = "\n".join(f"{key}={row.get(key, '')}" for key in sorted(row))
    key = f"v2-chain-normalized-{row['pipeline']}-{row['case_id']}-{hashlib.sha256(identity.encode()).hexdigest()[:16]}"
    candidate_root = comparison_root / "results" / key
    checkpoint = comparison_root / "checkpoints" / f"{key}.json"
    previous = json.loads(checkpoint.read_text(encoding="utf-8")) if checkpoint.is_file() else {}
    if previous.get("status") == "completed":
        print(json.dumps({"status": "resumed", "checkpoint": str(checkpoint), "key": key}, sort_keys=True))
        return 0

    candidate = dict(row)
    stages: dict[str, object] = dict(previous.get("stages", {}))
    preparation = stages.get("input_normalization", {})
    previous_candidate = previous.get("candidate", {})
    if (
        preparation.get("status") == "completed"
        and isinstance(previous_candidate, dict)
        and Path(str(previous_candidate.get("prepared_left", ""))).is_file()
        and Path(str(previous_candidate.get("prepared_right", ""))).is_file()
    ):
        candidate.update(
            prepared_left=str(previous_candidate["prepared_left"]),
            prepared_right=str(previous_candidate["prepared_right"]),
        )
    else:
        prep_started = time.perf_counter()
        prep_started_at = utc_now()
        try:
            candidate.update(prepare_candidate_inputs(candidate, candidate_root))
            stages["input_normalization"] = stage_record(
                "completed", prep_started_at, prep_started,
                prepared_left=candidate["prepared_left"],
                prepared_right=candidate["prepared_right"],
                prepared_left_sha256=sha256(Path(str(candidate["prepared_left"]))),
                prepared_right_sha256=sha256(Path(str(candidate["prepared_right"]))),
                prepared_left_chains=candidate["prepared_left_chains"],
                prepared_right_chains=candidate["prepared_right_chains"],
            )
        except Exception as exc:
            stages["input_normalization"] = stage_record(
                "failed", prep_started_at, prep_started, error=f"{type(exc).__name__}: {exc}"
            )
            failed_record = {
                "status": "failed",
                "key": key,
                "index": args.index,
                "candidate": candidate,
                "stages": stages,
                "finished_at": utc_now(),
            }
            atomic_json(checkpoint, failed_record)
            append_event(comparison_root / "events.jsonl", {"event": "candidate_failed", "key": key, "record": stages["input_normalization"], "at": utc_now()})
            print(json.dumps({"status": "failed", "key": key, "index": args.index, "stage": "input_normalization"}, sort_keys=True))
            return 1

    record: dict[str, object] = {
        "status": "running",
        "key": key,
        "index": args.index,
        "candidate": candidate,
        "input_hashes": {
            "left": sha256(Path(row["left"])),
            "right": sha256(Path(row["right"])),
            "native": sha256(Path(row["native_pdb"])),
            "prepared_left": sha256(Path(str(candidate["prepared_left"]))),
            "prepared_right": sha256(Path(str(candidate["prepared_right"]))),
        },
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
        "started_at": previous.get("started_at", utc_now()),
        "stages": stages,
    }
    atomic_json(checkpoint, record)
    event_path = comparison_root / "events.jsonl"
    append_event(event_path, {"event": "candidate_started", "key": key, "index": args.index, "at": utc_now()})

    for name, function in (
        ("fiberdock", lambda: run_fiberdock(candidate, candidate_root, args.pipeline_repo.resolve(), args.fiberdock_dir.resolve())),
        ("external_rosetta", lambda: run_rosetta(candidate, candidate_root, args.pipeline_repo.resolve(), args.rosetta_prepack, args.rosetta_dock, args.rosetta_db)),
    ):
        if stages.get(name, {}).get("status") == "completed":
            continue
        stages[name] = function()
        atomic_json(checkpoint, record)
        append_event(event_path, {"event": "stage_finished", "key": key, "stage": name, "record": stages[name], "at": utc_now()})

    for name, label in (("dockq_fiberdock", "fiberdock"), ("dockq_rosetta", "rosetta")):
        if stages.get(name, {}).get("status") == "scored":
            continue
        refiner_name = "fiberdock" if name == "dockq_fiberdock" else "external_rosetta"
        stages[name] = run_dockq(
            candidate,
            candidate_root,
            args.pipeline_repo.resolve(),
            stages.get(refiner_name, {}),
            args.dockq_python,
            label,
        )
        atomic_json(checkpoint, record)
        append_event(event_path, {"event": "stage_finished", "key": key, "stage": name, "record": stages[name], "at": utc_now()})

    failed = [name for name, value in stages.items() if value.get("status") in {"failed", "score_failed"}]
    record["status"] = "completed_with_stage_failures" if failed else "completed"
    record["failed_stages"] = failed
    record["finished_at"] = utc_now()
    atomic_json(checkpoint, record)
    append_event(event_path, {"event": "candidate_finished", "key": key, "status": record["status"], "failed_stages": failed, "at": utc_now()})
    print(json.dumps({"status": record["status"], "key": key, "index": args.index, "failed_stages": failed, "checkpoint": str(checkpoint)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
