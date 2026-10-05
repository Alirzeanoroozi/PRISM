#!/usr/bin/env python3
"""Score staged PRISM models with explicit multichain assignments.

DockQ is run from the installed DockQ package with a complete bijective
model:native chain map.  The benchmark result is taken from the requested
native receptor-ligand interfaces, not from DockQ's best internal interface.
iRMSD uses the exact benchmark ``irmsd.py`` command with grouped receptor and
ligand chain strings.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import subprocess
from pathlib import Path

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

try:
    from benchmark.scripts.standardized_evaluator import (
        no_align_is_safe,
        standardize_dockq_json,
        validate_pdb_mapping,
        validate_raw_pdb_chain_contract,
    )
except ModuleNotFoundError:  # Direct `python benchmark/scripts/...` invocation.
    from standardized_evaluator import no_align_is_safe, standardize_dockq_json, validate_pdb_mapping, validate_raw_pdb_chain_contract


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_complex(value: str) -> tuple[str, str, str]:
    pdb, groups = value.strip().split("_", 1)
    receptor, ligand = groups.split(":", 1)
    return pdb[:4].lower(), "".join(receptor), "".join(ligand)


def template_chain_group(value: str) -> str:
    """Extract chains from both ``1ABCAB`` and ``1ABC_AB`` template IDs."""
    token = (value or "").strip()
    if len(token) < 4:
        return ""
    return token[4:].replace("_", "")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_n_cpu(n_cpu: int | None = None) -> int:
    """Resolve and validate the DockQ CPU count used by this scorer."""

    value = n_cpu if n_cpu is not None else os.environ.get("DOCKQ_N_CPU", "1")
    try:
        value = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("DockQ n_cpu must be a positive integer") from exc
    if value < 1:
        raise ValueError("DockQ n_cpu must be a positive integer")
    return value


def raw_json_path_for_model(
    raw_json_dir: Path,
    stage: dict[str, str],
    model: Path,
    native: Path,
    mapping: str,
) -> Path:
    """Create a raw-JSON filename tied to the durable scoring identity."""
    identity = "\0".join(
        (
            stage.get("dataset_row_id", ""),
            stage.get("pair_id", ""),
            str(model.resolve()),
            str(native.resolve()),
            mapping,
        )
    )
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:16]
    return raw_json_dir / f"{model.stem}.{digest}.json"


def dockq_version(score_python: Path) -> str:
    """Return installed DockQ package version without making scoring depend on it."""

    command = [str(score_python), "-c", "from importlib.metadata import version; print(version('DockQ'))"]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"unavailable:{type(exc).__name__}"
    if result.returncode:
        return "unavailable:metadata_lookup_failed"
    return result.stdout.strip() or "unavailable:empty_version"


def is_transformation_half(path: Path) -> bool:
    """Identify an unassembled PRISM transformation half by its filename."""

    return "_L" in path.stem and "_R" not in path.stem


def validate_score_candidate(stage: dict[str, str], native: Path) -> Path:
    """Fail closed before parsing or externally scoring a staged model."""

    if stage.get("status") != "staged_symlink":
        raise ValueError(f"model is not stageable: {stage.get('status', '')}")
    model = Path(stage.get("staged_model_path", ""))
    if not model.is_file():
        raise ValueError(f"staged model is missing: {model}")
    if is_transformation_half(model):
        raise ValueError("transformation half cannot be scored")
    if not native.is_file():
        raise ValueError(f"native PDB is missing: {native}")
    contract = validate_raw_pdb_chain_contract(
        model,
        stage.get("model_receptor_chains", ""),
        stage.get("model_ligand_chains", ""),
    )
    if not contract.valid:
        raise ValueError("model chain contract rejected: " + "; ".join(contract.errors))
    return model


def validate_mapping_option(model: Path, native: Path, mapping: str, no_align: bool) -> str:
    """Require the strict PDB mapping contract before DockQ no-align mode."""

    if not no_align:
        return "aligned_default"
    left, right = mapping.split(":", 1)
    validation = validate_pdb_mapping(model, native, dict(zip(left, right)))
    if not no_align_is_safe(validation):
        raise ValueError("unsafe --no-align mapping: " + "; ".join(validation.errors))
    return "validated_no_align"


def chain_sequence(structure, chain_id: str) -> str:
    chain = structure[0][chain_id]
    return "".join(
        protein_letters_3to1[residue.resname]
        for residue in chain
        if residue.id[0] == " " and all(atom in residue for atom in ("N", "CA", "C", "O"))
    )


def sequence_score(a: str, b: str) -> tuple[int, int, int]:
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    alignment = aligner.align(a, b)[0]
    row_a, row_b = str(alignment[0]), str(alignment[1])
    paired = sum(x != "-" and y != "-" for x, y in zip(row_a, row_b))
    identity = sum(x == y and x != "-" for x, y in zip(row_a, row_b))
    return identity, paired, int(round(alignment.score))


def assign_group(model_structure, native_structure, model_chains: str, native_chains: str) -> tuple[str, str]:
    if len(model_chains) != len(native_chains):
        raise NonBijectiveGroupError(f"non-bijective group: model={model_chains!r} native={native_chains!r}")
    scores = {
        (model_chain, native_chain): sequence_score(
            chain_sequence(model_structure, model_chain),
            chain_sequence(native_structure, native_chain),
        )
        for model_chain in model_chains
        for native_chain in native_chains
    }
    candidates = []
    for permutation in itertools.permutations(model_chains):
        total = tuple(sum(scores[(m, n)][i] for m, n in zip(permutation, native_chains)) for i in range(3))
        candidates.append((total, "".join(permutation)))
    _, best_model_order = max(candidates)
    diagnostics = ";".join(
        f"{m}>{n}:{scores[(m, n)][0]}/{scores[(m, n)][1]}"
        for m, n in zip(best_model_order, native_chains)
    )
    return best_model_order, diagnostics


class NonBijectiveGroupError(ValueError):
    """A declared partner group cannot form a one-to-one chain mapping."""


def run_dockq(score_python: Path, dockq_executable: Path | None, model: Path, native: Path, mapping: str, raw_json: Path, timeout: int, no_align: bool = False, n_cpu: int | None = None) -> dict:
    command = ([str(dockq_executable)] if dockq_executable else [str(score_python), "-m", "DockQ"])
    command.extend([str(model), str(native), "--mapping", mapping, "--json", str(raw_json)])
    command.extend(["--n_cpu", str(resolve_n_cpu(n_cpu))])
    if no_align:
        command.append("--no_align")
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError((result.stderr or result.stdout).strip()[-2000:])
    return json.loads(raw_json.read_text())


def dockq_command(score_python: Path, dockq_executable: Path | None, model: Path, native: Path, mapping: str, raw_json: Path, no_align: bool = False, n_cpu: int | None = None) -> list[str]:
    command = ([str(dockq_executable)] if dockq_executable else [str(score_python), "-m", "DockQ"])
    command.extend([str(model), str(native), "--mapping", mapping, "--json", str(raw_json)])
    command.extend(["--n_cpu", str(resolve_n_cpu(n_cpu))])
    if no_align:
        command.append("--no_align")
    return command


def is_recoverable_complete_mapping_error(exc: Exception) -> bool:
    """Limit pairwise recovery to the observed DockQ empty-array crash."""

    return "Buffer has wrong number of dimensions" in str(exc)


def run_pairwise_cross_dockq(
    score_python: Path,
    dockq_executable: Path | None,
    model: Path,
    native: Path,
    model_r: str,
    model_l: str,
    native_r: str,
    native_l: str,
    raw_json: Path,
    timeout: int,
    no_align: bool = False,
) -> tuple[dict, list[dict[str, object]], list[list[str]]]:
    """Score only requested receptor-ligand interfaces after full DockQ fails.

    DockQ 2.1.3 can crash while evaluating an unrelated internal interface in
    a multichain mapping. Pairwise execution avoids that internal interface.
    Native chain pairs with no interface are explicit skips. The synthetic
    document intentionally has no GlobalDockQ because cross-only scores are
    not a complete-complex metric.
    """

    successful: dict[str, dict] = {}
    pairwise_records: list[dict[str, object]] = []
    commands: list[list[str]] = []
    for receptor_index, native_receptor in enumerate(native_r):
        for ligand_index, native_ligand in enumerate(native_l):
            model_receptor = model_r[receptor_index]
            model_ligand = model_l[ligand_index]
            interface = native_receptor + native_ligand
            pair_mapping = f"{model_receptor}{model_ligand}:{interface}"
            pair_json = raw_json.with_name(f"{raw_json.stem}.pair-{interface}.json")
            try:
                payload = run_dockq(
                    score_python, dockq_executable, model, native,
                    pair_mapping, pair_json, timeout, no_align,
                )
            except RuntimeError as exc:
                if "Could not find interfaces in the native model" not in str(exc):
                    raise
                pairwise_records.append({
                    "interface": interface,
                    "model_interface": model_receptor + model_ligand,
                    "mapping": pair_mapping,
                    "status": "no_native_interface",
                    "error": str(exc),
                })
                continue
            results = payload.get("best_result", {})
            result = results.get(interface)
            if result is None and len(results) == 1:
                result = next(iter(results.values()))
            if result is None:
                raise RuntimeError(
                    f"pairwise DockQ omitted requested interface {interface}: {pair_mapping}"
                )
            successful[interface] = result
            command = dockq_command(
                score_python, dockq_executable, model, native,
                pair_mapping, pair_json, no_align,
            )
            commands.append(command)
            pairwise_records.append({
                "interface": interface,
                "model_interface": model_receptor + model_ligand,
                "mapping": pair_mapping,
                "status": "scored",
                "raw_json": str(pair_json),
                "raw_json_sha256": sha256_file(pair_json),
                "argv": command,
            })
    if not successful:
        raise RuntimeError("pairwise DockQ found no scoreable requested cross interface")
    combined = {
        "schema": "prism-requested-cross-interface-dockq/v1",
        "evaluation_mode": "pairwise_cross_fallback",
        "GlobalDockQ": None,
        "model": str(model),
        "native": str(native),
        "score_scope": "requested_cross_interfaces_only",
        "best_result": successful,
        "component_runs": pairwise_records,
    }
    raw_json.write_text(json.dumps(combined, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return combined, pairwise_records, commands


def run_irmsd(score_python: Path, script: Path, model: Path, model_r: str, model_l: str, native: Path, native_r: str, native_l: str, timeout: int) -> float:
    command = [str(score_python), str(script), str(model), model_r, model_l, str(native), native_r, native_l]
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError((result.stderr or result.stdout).strip()[-2000:])
    return float(result.stdout.strip().splitlines()[-1])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-manifest", type=Path, required=True)
    parser.add_argument("--native-root", type=Path, required=True)
    parser.add_argument("--score-python", type=Path, required=True)
    parser.add_argument("--dockq-executable", type=Path, help="Optional DockQ-compatible executable; production defaults to python -m DockQ")
    parser.add_argument("--irmsd-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--interfaces-output",
        type=Path,
        help="TSV for raw DockQ global/component records; defaults beside --output.",
    )
    parser.add_argument("--raw-json-dir", type=Path, required=True)
    parser.add_argument("--benchmark-set", choices=("rigid", "medium", "difficult"))
    parser.add_argument("--dataset-row-id", action="append", help="Score only this durable row identity; repeatable.")
    parser.add_argument("--limit", type=int, help="Score at most this many manifest rows after filtering.")
    parser.add_argument("--shard-count", type=int, default=1, help="Deterministically split filtered rows across workers.")
    parser.add_argument("--shard-index", type=int, default=0, help="Zero-based shard selected from --shard-count.")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--no-align", action="store_true", help="Use DockQ --no_align only after strict PDB mapping validation.")
    args = parser.parse_args()

    parser_pdb = PDBParser(QUIET=True)
    args.raw_json_dir.mkdir(parents=True, exist_ok=True)
    dockq_cpu = resolve_n_cpu()
    installed_dockq_version = dockq_version(args.score_python)
    output_rows: list[dict[str, str]] = []
    interface_rows: list[dict[str, object]] = []
    stages = read_csv(args.stage_manifest)
    if args.benchmark_set:
        stages = [stage for stage in stages if stage.get("benchmark_set") == args.benchmark_set]
    if args.dataset_row_id:
        requested_ids = set(args.dataset_row_id)
        stages = [stage for stage in stages if stage.get("dataset_row_id") in requested_ids]
    if args.shard_count <= 0 or not 0 <= args.shard_index < args.shard_count:
        parser.error("require --shard-count > 0 and 0 <= --shard-index < --shard-count")
    stages = [stage for index, stage in enumerate(stages) if index % args.shard_count == args.shard_index]
    if args.limit is not None:
        if args.limit <= 0:
            parser.error("--limit must be positive")
        stages = stages[: args.limit]
    for stage in stages:
        row = dict(stage)
        row["dockq_n_cpu"] = dockq_cpu
        if stage.get("status") != "staged_symlink":
            row["score_status"] = "not_scoreable"
            output_rows.append(row)
            continue
        if stage.get("source_gate_status") == "audit_only":
            row["score_status"] = "not_scoreable"
            row["score_error"] = "source_gate_audit_only"
            output_rows.append(row)
            continue
        model = Path(stage["staged_model_path"])
        source_model = Path(stage.get("source_model_path") or stage["staged_model_path"])
        if source_model.is_file():
            row["source_model_sha256"] = stage.get("source_model_sha256") or sha256_file(source_model)
        try:
            native_pdb, native_r, native_l = parse_complex(stage["complex"])
            native = Path(stage["native_pdb_path"]) if stage.get("native_pdb_path") else args.native_root / f"{native_pdb}.pdb"
            # Hash the native input before any later validation or DockQ call.
            # Failure rows must retain the same provenance anchors as scored
            # rows; otherwise a reproducible runtime failure cannot be tied to
            # the exact native structure that was attempted.
            if native.is_file():
                row["native_pdb_sha256"] = sha256_file(native)
            model = validate_score_candidate(stage, native)
            model_structure = parser_pdb.get_structure("model", str(model))
            native_structure = parser_pdb.get_structure("native", str(native))
            model_r_hint = stage.get("model_receptor_chains") or template_chain_group(stage["template_1"])
            model_l_hint = stage.get("model_ligand_chains") or template_chain_group(stage["template_2"])
            if set(model_r_hint).intersection(model_l_hint):
                raise ValueError(f"overlapping model receptor/ligand chains: {model_r_hint}:{model_l_hint}")
            try:
                model_r, receptor_diag = assign_group(model_structure, native_structure, model_r_hint, native_r)
                model_l, ligand_diag = assign_group(model_structure, native_structure, model_l_hint, native_l)
            except NonBijectiveGroupError as exc:
                row["score_status"] = "not_scoreable"
                row["score_error"] = "non_bijective_partner_cardinality"
                row["score_error_detail"] = str(exc)
                output_rows.append(row)
                continue
            mapping = f"{model_r}{model_l}:{native_r}{native_l}"
            mapping_validation_status = validate_mapping_option(model, native, mapping, args.no_align)
            row.update({
                "model_receptor_chains": model_r,
                "model_ligand_chains": model_l,
                "native_receptor_chains": native_r,
                "native_ligand_chains": native_l,
                "dockq_mapping": mapping,
                "mapping_validation_status": mapping_validation_status,
                "receptor_assignment": receptor_diag,
                "ligand_assignment": ligand_diag,
            })
            raw_json = raw_json_path_for_model(args.raw_json_dir, stage, model, native, mapping)
            score_status = "scored"
            score_scope = "complete_bijective_mapping"
            dockq_full_error = ""
            fallback_components: list[dict[str, object]] = []
            fallback_commands: list[list[str]] = []
            try:
                dockq = run_dockq(
                    args.score_python, args.dockq_executable, model, native,
                    mapping, raw_json, args.timeout, args.no_align,
                )
                standardized = standardize_dockq_json(raw_json)
                raw_json_sha256 = standardized.raw_json_sha256
            except RuntimeError as exc:
                if not is_recoverable_complete_mapping_error(exc):
                    raise
                dockq_full_error = str(exc)
                dockq, fallback_components, fallback_commands = run_pairwise_cross_dockq(
                    args.score_python, args.dockq_executable, model, native,
                    model_r, model_l, native_r, native_l,
                    raw_json, args.timeout, args.no_align,
                )
                standardized = None
                raw_json_sha256 = sha256_file(raw_json)
                score_scope = "requested_cross_interfaces_only"
                score_status = "scored_cross_only"
            cross_keys = [f"{r}{l}" for r in native_r for l in native_l]
            cross = [dockq["best_result"][key] for key in cross_keys if key in dockq.get("best_result", {})]
            if not cross:
                raise RuntimeError(
                    f"no requested cross interface in DockQ result: expected={cross_keys} "
                    f"observed={sorted(dockq.get('best_result', {}))}"
                )
            fallback_status_by_interface = {
                str(component.get("interface", "")): str(component.get("status", ""))
                for component in fallback_components
            }
            if score_scope == "requested_cross_interfaces_only" and set(fallback_status_by_interface) != set(cross_keys):
                raise RuntimeError("pairwise fallback did not account for every requested cross interface")
            unscoreable_cross = [
                key for key in cross_keys if fallback_status_by_interface.get(key) == "no_native_interface"
            ]
            grouped = reciprocal = None
            irmsd_status = "scored"
            irmsd_error = ""
            try:
                grouped = run_irmsd(args.score_python, args.irmsd_script, model, model_r, model_l, native, native_r, native_l, args.timeout)
                reciprocal = run_irmsd(args.score_python, args.irmsd_script, model, model_l, model_r, native, native_l, native_r, args.timeout)
            except Exception as exc:
                irmsd_status = "failed_auxiliary"
                irmsd_error = str(exc)
            row.update(
                {
                    "score_status": score_status,
                    "score_scope": score_scope,
                    "dockq_full_error": dockq_full_error,
                    "dockq_global": str(dockq.get("GlobalDockQ", "")) if score_scope == "complete_bijective_mapping" else "",
                    "dockq_global_status": "scored" if score_scope == "complete_bijective_mapping" else "unavailable_cross_only",
                    "dockq_best_internal": str(dockq.get("best_dockq", "")) if score_scope == "complete_bijective_mapping" else "",
                    "dockq_cross_best": str(max(x["DockQ"] for x in cross)),
                    "dockq_cross_mean": str(sum(x["DockQ"] for x in cross) / len(cross)),
                    "dockq_cross_components": json.dumps(cross, sort_keys=True),
                    "dockq_cross_requested_count": str(len(cross_keys)),
                    "dockq_cross_scoreable_count": str(len(cross)),
                    "dockq_cross_unscoreable_count": str(len(unscoreable_cross)),
                    "dockq_cross_unscoreable_interfaces": json.dumps(unscoreable_cross),
                    "irmsd_status": irmsd_status,
                    "irmsd_error": irmsd_error,
                    "irmsd_grouped_forward": "" if grouped is None else str(grouped),
                    "irmsd_grouped_reverse": "" if reciprocal is None else str(reciprocal),
                    "irmsd_grouped_min": "" if grouped is None or reciprocal is None else str(min(grouped, reciprocal)),
                    "raw_dockq_json": str(raw_json),
                    "raw_dockq_json_sha256": raw_json_sha256,
                    "dockq_component_runs": json.dumps(fallback_components, sort_keys=True),
                    "source_model_sha256": sha256_file(model),
                    "native_pdb_sha256": sha256_file(native),
                    "dockq_argv": json.dumps(dockq_command(
                        args.score_python, args.dockq_executable, model, native,
                        mapping, raw_json, args.no_align,
                    )),
                    "dockq_fallback_argv": json.dumps(fallback_commands),
                    "dockq_version": installed_dockq_version,
                }
            )
            interface_context = {
                "pair_id": stage.get("pair_id", ""),
                "benchmark_set": stage.get("benchmark_set", ""),
                "source_model_path": stage.get("source_model_path", ""),
                "staged_model_path": str(model),
                "source_model_sha256": sha256_file(model),
                "native_pdb_path": str(native),
                "native_pdb_sha256": sha256_file(native),
                "dockq_mapping": mapping,
                "mapping_validation_status": mapping_validation_status,
                "raw_dockq_json": str(raw_json),
                "requested_cross_interface": "",
                "dockq_version": installed_dockq_version,
                "score_scope": score_scope,
                "dockq_global_status": "scored" if score_scope == "complete_bijective_mapping" else "unavailable_cross_only",
            }
            if standardized is not None:
                for score_record in standardized:
                    interface_rows.append(
                        {
                            **interface_context,
                            **dict(score_record),
                            "requested_cross_interface": str(score_record.get("interface", "")) in cross_keys,
                        }
                    )
            else:
                for interface in sorted(dockq["best_result"]):
                    interface_rows.append(
                        {
                            **interface_context,
                            "record_type": "interface",
                            "interface": interface,
                            "GlobalDockQ": "",
                            **dockq["best_result"][interface],
                            "grouped_iRMSD": "",
                            "raw_json_sha256": raw_json_sha256,
                            "requested_cross_interface": True,
                        }
                    )
        except Exception as exc:
            row["score_status"] = "score_failed"
            row["score_error"] = str(exc)
        output_rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in output_rows for key in row})
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)
    interfaces_output = args.interfaces_output or args.output.with_name(args.output.stem + "_interfaces.tsv")
    interfaces_output.parent.mkdir(parents=True, exist_ok=True)
    interface_fields = sorted({key for row in interface_rows for key in row})
    with interfaces_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=interface_fields or ["record_type"], delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(interface_rows)
    print(f"wrote {len(output_rows)} model rows to {args.output}; {len(interface_rows)} interface rows to {interfaces_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
