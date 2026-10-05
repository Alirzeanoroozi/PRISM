#!/usr/bin/env python3
"""Resume an isolated PRISM run through Rosetta refinement and DockQ scoring."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


MODEL_RE = re.compile(
    r"^(?P<template>[^_]+)_(?P<receptor>[^_]+)_(?P<ligand>[^_]+)_o(?P<orientation>[12])_(?P<side>[LR])\.pdb$"
)


def now() -> str:
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


def append_event(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def model_fields(left: Path) -> dict[str, str] | None:
    match = MODEL_RE.match(left.name)
    if match is None or match.group("side") != "L":
        return None
    return {key: match.group(key) for key in ("template", "receptor", "ligand", "orientation")}


def score_refined_stage(
    *,
    stage: str,
    entries: list[dict],
    root: Path,
    repo: Path,
    dockq_python: Path,
    downstream: Path,
    script: Path,
) -> tuple[str, float, Path]:
    input_manifest = downstream / f"{stage}_dockq_inputs.json"
    output_root = root / "processed" / f"dockq_irmsd_refined_{stage}"
    summary_path = output_root / "score_summary.json"
    atomic_json(input_manifest, entries)
    if not entries:
        output_root.mkdir(parents=True, exist_ok=True)
        atomic_json(
            summary_path,
            {
                "status": "completed",
                "stage": stage,
                "model_groups": 0,
                "rows": 0,
                "scored": 0,
                "unscored": 0,
                "reason": "no_refined_candidates",
            },
        )
        (output_root / "dockq_irmsd.csv").write_text("stage,status,error\n", encoding="utf-8")
        return "no_candidates", 0.0, summary_path
    started = time.perf_counter()
    environment = dict(os.environ)
    environment.update({"DOCKQ_PYTHON": str(dockq_python), "PYTHONNOUSERSITE": "1"})
    process = subprocess.run(
        [
            str(dockq_python),
            str(script),
            "--run-root",
            str(root),
            "--pipeline-repo",
            str(repo),
            "--input-manifest",
            str(input_manifest),
            "--output-root",
            str(output_root),
        ],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    (downstream / f"{stage}_dockq.console.log").write_text(
        process.stdout + process.stderr, encoding="utf-8"
    )
    status = "completed" if process.returncode == 0 and summary_path.is_file() else "failed"
    return status, time.perf_counter() - started, summary_path


def pair_inputs(run_root: Path) -> list[tuple[Path, Path]]:
    transformation = run_root / "processed" / "transformation"
    pairs = []
    for left in sorted(transformation.glob("*_L.pdb")):
        right = left.with_name(left.name[:-6] + "_R.pdb")
        if right.is_file():
            pairs.append((left, right))
    return pairs


def checkpoint_key(left: Path, right: Path) -> str:
    value = f"{left.name}\n{right.name}".encode("utf-8")
    return hashlib.sha256(value).hexdigest()[:24]


def probe_pyrosetta(python_bin: Path) -> dict:
    started = time.perf_counter()
    if not python_bin.is_file():
        return {"status": "unavailable", "reason": "python_not_found", "python": str(python_bin), "elapsed_seconds": time.perf_counter() - started}
    result = subprocess.run(
        [str(python_bin), "-c", "import pyrosetta; print(getattr(pyrosetta, '__version__', 'available'))"],
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "status": "available" if result.returncode == 0 else "unavailable",
        "python": str(python_bin),
        "return_code": result.returncode,
        "version": result.stdout.strip(),
        "error": result.stderr.strip()[-1000:],
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--pipeline-repo", required=True)
    parser.add_argument("--pipeline-python", required=True)
    parser.add_argument("--dockq-python", required=True)
    parser.add_argument("--pyrosetta-python", required=True)
    args = parser.parse_args()

    root = Path(args.run_root).resolve()
    repo = Path(args.pipeline_repo).resolve()
    python_bin = Path(args.pipeline_python).resolve()
    dockq_python = Path(args.dockq_python).resolve()
    pyrosetta_python = Path(args.pyrosetta_python).resolve()
    downstream = root / "downstream"
    checkpoints = downstream / "checkpoints"
    events = downstream / "events.jsonl"
    timing_path = downstream / "timing.json"
    console_log = downstream / "dockq.console.log"
    started = time.perf_counter()

    summary_path = root / "run_summary.json"
    if not summary_path.is_file():
        atomic_json(downstream / "status.json", {"status": "blocked", "reason": "missing_run_summary", "updated_at": now()})
        return 2
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "completed":
        atomic_json(downstream / "status.json", {"status": "blocked", "reason": "parent_run_not_completed", "parent": summary, "updated_at": now()})
        return 2
    if summary.get("rank") is True or summary.get("rank_method") == "prodigy":
        atomic_json(downstream / "status.json", {"status": "blocked", "reason": "ranking_must_be_disabled", "parent": summary, "updated_at": now()})
        return 2

    atomic_json(
        downstream / "manifest.json",
        {
            "status": "running",
            "started_at": now(),
            "job_id": os.environ.get("SLURM_JOB_ID", ""),
            "parent_run_root": str(root),
            "parent_run_summary_sha256": sha256(summary_path),
            "refiner": "external_rosetta",
            "ranking": False,
            "prodigy": False,
            "dockq_python": str(dockq_python),
            "pyrosetta_python": str(pyrosetta_python),
            "resume_contract": "completed per-candidate checkpoints are skipped; incomplete checkpoints are retried",
        },
    )

    # Import only after the run root is selected: the legacy refiner creates
    # its output directories relative to the current working directory.
    previous_cwd = Path.cwd()
    os.chdir(root)
    try:
        sys.path.insert(0, str(root))
        sys.path.insert(0, str(repo))
        from src.rosetta_refinement import calculate_energy
        from src.rosetta_refinement import partner_chain_ids

        candidates = pair_inputs(root)
        refined_entries = {"external_rosetta": [], "pyrosetta": []}
        pyrosetta_probe = probe_pyrosetta(pyrosetta_python)
        atomic_json(downstream / "pyrosetta_probe.json", pyrosetta_probe)
        append_event(events, {"stage": "pyrosetta_probe", "event": pyrosetta_probe["status"], "at": now(), **pyrosetta_probe})

        pyrosetta_started = time.perf_counter()
        pyrosetta_completed = 0
        pyrosetta_failed = 0
        if pyrosetta_probe["status"] == "available":
            pyrosetta_script = Path(__file__).with_name("pyrosetta_refine_one.py")
            for left, right in candidates:
                key = checkpoint_key(left, right)
                checkpoint = checkpoints / "pyrosetta" / f"{key}.json"
                if checkpoint.is_file() and json.loads(checkpoint.read_text(encoding="utf-8")).get("status") == "completed":
                    pyrosetta_completed += 1
                    append_event(events, {"stage": "pyrosetta", "event": "resumed", "checkpoint": str(checkpoint), "at": now()})
                    continue
                item_started = time.perf_counter()
                output_root = root / "processed" / "pyrosetta_refinement"
                process = subprocess.run(
                    [str(pyrosetta_python), str(pyrosetta_script), "--run-root", str(root), "--pipeline-repo", str(repo), "--left", str(left), "--right", str(right), "--output-root", str(output_root)],
                    cwd=root,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                status = "completed" if process.returncode == 0 else "failed"
                if status == "completed":
                    pyrosetta_completed += 1
                else:
                    pyrosetta_failed += 1
                atomic_json(checkpoint, {"status": status, "left": str(left), "right": str(right), "return_code": process.returncode, "stdout": process.stdout[-2000:], "stderr": process.stderr[-2000:], "elapsed_seconds": time.perf_counter() - item_started, "at": now()})
                append_event(events, {"stage": "pyrosetta", "event": status, "left": str(left), "right": str(right), "elapsed_seconds": time.perf_counter() - item_started, "at": now()})
        pyrosetta_seconds = time.perf_counter() - pyrosetta_started

        refinement_started = time.perf_counter()
        refinement_dir = checkpoints / "refinement"
        completed = 0
        terminal_no_output = 0
        for left, right in candidates:
            key = checkpoint_key(left, right)
            checkpoint = refinement_dir / f"{key}.json"
            if checkpoint.is_file():
                previous = json.loads(checkpoint.read_text(encoding="utf-8"))
                if previous.get("status") in {"completed", "terminal_no_output"}:
                    append_event(events, {"stage": "refinement", "event": "resumed", "checkpoint": str(checkpoint), "at": now()})
                    completed += previous.get("status") == "completed"
                    terminal_no_output += previous.get("status") == "terminal_no_output"
                    continue
            item_started = time.perf_counter()
            status = "terminal_no_output"
            result = {"totalscore": "-", "interaction_score": "-", "structure": "-"}
            error = ""
            try:
                totalscore, interaction_score, structure = calculate_energy(str(left), str(right))
                result = {"totalscore": totalscore, "interaction_score": interaction_score, "structure": structure}
                if structure != "-" and Path(structure).is_file():
                    status = "completed"
                    completed += 1
                else:
                    terminal_no_output += 1
            except Exception as exc:  # checkpoint the failure for resumable retry
                status = "failed"
                error = f"{type(exc).__name__}: {exc}"
            atomic_json(
                checkpoint,
                {
                    "status": status,
                    "left": str(left),
                    "right": str(right),
                    "result": result,
                    "error": error,
                    "started_at": now(),
                    "elapsed_seconds": time.perf_counter() - item_started,
                },
            )
            append_event(events, {"stage": "refinement", "event": status, "left": str(left), "right": str(right), "elapsed_seconds": time.perf_counter() - item_started, "at": now()})
            if status == "failed":
                atomic_json(downstream / "status.json", {"status": "failed", "stage": "refinement", "left": str(left), "right": str(right), "error": error, "updated_at": now()})
                return 1
        refinement_seconds = time.perf_counter() - refinement_started

        # Rebuild the refined candidate list from checkpoints so a rerun after
        # a timeout or interruption scores previously completed candidates.
        for left, right in candidates:
            checkpoint = refinement_dir / f"{checkpoint_key(left, right)}.json"
            if not checkpoint.is_file():
                continue
            record = json.loads(checkpoint.read_text(encoding="utf-8"))
            result = record.get("result", {})
            fields = model_fields(left)
            structure = result.get("structure")
            if record.get("status") == "completed" and fields and structure and structure != "-" and Path(structure).is_file():
                left_chains, right_chains = partner_chain_ids(str(left), str(right))
                refined_entries["external_rosetta"].append(
                    {
                        **fields,
                        "stage": "external_rosetta",
                        "left": str(left),
                        "right": str(right),
                        "model_pdb": str(Path(structure).resolve()),
                        "model_receptor_chains": left_chains,
                        "model_ligand_chains": right_chains,
                    }
                )

        if pyrosetta_probe["status"] == "available":
            for checkpoint in sorted((checkpoints / "pyrosetta").glob("*.json")):
                record = json.loads(checkpoint.read_text(encoding="utf-8"))
                if record.get("status") != "completed":
                    continue
                try:
                    printed = json.loads(record.get("stdout", ""))
                    result = printed[0] if isinstance(printed, list) and printed else {}
                except json.JSONDecodeError:
                    result = {}
                fields = model_fields(Path(record.get("left", "")))
                output = result.get("output_path")
                if fields and output and Path(output).is_file():
                    partners = str(result.get("partners", "A_B"))
                    left_chains, right_chains = partners.split("_", 1)
                    refined_entries["pyrosetta"].append(
                        {
                            **fields,
                            "stage": "pyrosetta",
                            "left": record["left"],
                            "right": record["right"],
                            "model_pdb": str(Path(output).resolve()),
                            "model_receptor_chains": left_chains,
                            "model_ligand_chains": right_chains,
                        }
                    )

        score_checkpoint = checkpoints / "dockq.json"
        score_summary = root / "processed" / "dockq_irmsd" / "score_summary.json"
        score_started = time.perf_counter()
        if score_checkpoint.is_file() and score_summary.is_file():
            score_status = "resumed"
        else:
            scorer = repo.parent / "PRISM" / "notebooks" / "tmp_prism_all_pipelines" / "score_transformed_models.py"
            if not scorer.is_file():
                scorer = root / "score_transformed_models.py"
            environment = dict(os.environ)
            environment.update({"DOCKQ_PYTHON": str(dockq_python), "PYTHONNOUSERSITE": "1"})
            with console_log.open("a", encoding="utf-8") as handle:
                process = subprocess.run(
                    [str(python_bin), str(scorer), "--run-root", str(root), "--pipeline-repo", str(repo)],
                    cwd=root,
                    env=environment,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    text=True,
                )
            score_status = "completed" if process.returncode == 0 and score_summary.is_file() else "failed"
            atomic_json(score_checkpoint, {"status": score_status, "return_code": process.returncode, "scorer": str(scorer), "elapsed_seconds": time.perf_counter() - score_started, "at": now()})
            if score_status == "failed":
                atomic_json(downstream / "status.json", {"status": "failed", "stage": "dockq", "return_code": process.returncode, "updated_at": now()})
                return 1
        dockq_seconds = time.perf_counter() - score_started
        refined_dockq_script = Path(__file__).with_name("score_refined_models.py")
        external_refined_status, external_refined_seconds, external_refined_summary = score_refined_stage(
            stage="external_rosetta",
            entries=refined_entries["external_rosetta"],
            root=root,
            repo=repo,
            dockq_python=dockq_python,
            downstream=downstream,
            script=refined_dockq_script,
        )
        if pyrosetta_probe["status"] == "available":
            pyro_refined_status, pyro_refined_seconds, pyro_refined_summary = score_refined_stage(
                stage="pyrosetta",
                entries=refined_entries["pyrosetta"],
                root=root,
                repo=repo,
                dockq_python=dockq_python,
                downstream=downstream,
                script=refined_dockq_script,
            )
        else:
            pyro_refined_status, pyro_refined_seconds, pyro_refined_summary = "not_run_unavailable", 0.0, None
        timing = {
            "status": "completed",
            "refinement_candidates": len(candidates),
            "pyrosetta_status": pyrosetta_probe["status"],
            "pyrosetta_completed": pyrosetta_completed,
            "pyrosetta_failed": pyrosetta_failed,
            "pyrosetta_seconds": pyrosetta_seconds,
            "refinement_completed": completed,
            "refinement_terminal_no_output": terminal_no_output,
            "refinement_seconds": refinement_seconds,
            "dockq_seconds": dockq_seconds,
            "refined_dockq_external_status": external_refined_status,
            "refined_dockq_external_candidates": len(refined_entries["external_rosetta"]),
            "refined_dockq_external_seconds": external_refined_seconds,
            "refined_dockq_external_summary": str(external_refined_summary),
            "refined_dockq_pyrosetta_status": pyro_refined_status,
            "refined_dockq_pyrosetta_candidates": len(refined_entries["pyrosetta"]),
            "refined_dockq_pyrosetta_seconds": pyro_refined_seconds,
            "refined_dockq_pyrosetta_summary": str(pyro_refined_summary) if pyro_refined_summary else None,
            "total_seconds": time.perf_counter() - started,
            "updated_at": now(),
            "dockq_summary": str(score_summary),
        }
        atomic_json(timing_path, timing)
        atomic_json(downstream / "status.json", timing)
        manifest = json.loads((downstream / "manifest.json").read_text(encoding="utf-8"))
        manifest.update({"status": "completed", "finished_at": now(), "timing": timing})
        atomic_json(downstream / "manifest.json", manifest)
        print(json.dumps(timing, indent=2))
        return 0
    finally:
        os.chdir(previous_cwd)


if __name__ == "__main__":
    raise SystemExit(main())
