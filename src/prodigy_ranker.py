"""Opt-in PRODIGY ranking for PRISM's merged candidate PDB model."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import shutil
import subprocess

from .pdb_download import split_target_id


_AFFINITY_LINE = re.compile(
    r"(?:^|\s)(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$"
)


@dataclass(frozen=True)
class ProdigyScore:
    candidate_pdb: str
    status: str
    affinity_kcal_mol: float | None
    return_code: int | None
    command: tuple[str, ...]
    stdout_path: str
    stderr_path: str
    input_sha256: str | None
    error: str | None = None


def parse_affinity(stdout):
    """Read PRODIGY quiet output or its labelled affinity line."""
    values = []
    for line in stdout.splitlines():
        match = _AFFINITY_LINE.search(line.strip())
        if match:
            values.append(float(match.group(1)))
    if not values:
        raise ValueError("PRODIGY output contains no affinity value")
    affinity = values[-1]
    if not math.isfinite(affinity):
        raise ValueError("PRODIGY affinity is not finite")
    return affinity


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def score_candidate(
    candidate_pdb,
    receptor,
    ligand,
    *,
    executable="prodigy",
    output_dir="processed/ranking/prodigy",
    distance_cutoff=5.5,
    acc_threshold=0.05,
    temperature=25.0,
    timeout=120.0,
):
    """Score one already-merged PRISM candidate and retain command evidence."""
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(candidate_pdb).stem)[:180]
    stdout_path = root / f"{stem}.stdout.txt"
    stderr_path = root / f"{stem}.stderr.txt"
    metadata_path = root / f"{stem}.json"
    _, receptor_chains = split_target_id(receptor)
    _, ligand_chains = split_target_id(ligand)
    command = tuple(
        shlex.split(executable)
        + [
            "-q",
            "--distance-cutoff", str(distance_cutoff),
            "--acc-threshold", str(acc_threshold),
            "--temperature", str(temperature),
            str(candidate_pdb),
            "--selection", "".join(receptor_chains), "".join(ligand_chains),
        ]
    )
    try:
        completed = subprocess.run(
            list(command), check=False, capture_output=True, text=True, timeout=timeout,
        )
        stdout_path.write_text(completed.stdout, encoding="utf-8")
        stderr_path.write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(f"PRODIGY exited with return code {completed.returncode}")
        result = ProdigyScore(
            str(candidate_pdb), "scored", parse_affinity(completed.stdout),
            completed.returncode, command, str(stdout_path), str(stderr_path),
            _sha256(candidate_pdb),
        )
    except Exception as exc:
        stdout_path.touch(exist_ok=True)
        stderr_path.touch(exist_ok=True)
        result = ProdigyScore(
            str(candidate_pdb), "failed", None, None, command,
            str(stdout_path), str(stderr_path),
            _sha256(candidate_pdb) if Path(candidate_pdb).is_file() else None,
            f"{type(exc).__name__}: {exc}",
        )
    metadata_path.write_text(json.dumps(asdict(result), indent=2, sort_keys=True) + "\n")
    return result


def select_top_candidates(
    candidates,
    *,
    top_k=5,
    executable="prodigy",
    output_dir="processed/ranking/prodigy",
    distance_cutoff=5.5,
    acc_threshold=0.05,
    temperature=25.0,
    timeout=120.0,
):
    """Select the most favorable affinities per receptor-ligand pair."""
    if top_k < 1:
        raise ValueError("top_k must be positive")
    executable_path = shlex.split(executable)[0]
    if shutil.which(executable_path) is None and not Path(executable_path).exists():
        raise FileNotFoundError(f"PRODIGY executable not found: {executable!r}")

    groups = {}
    for candidate in candidates:
        groups.setdefault((candidate[0], candidate[1]), []).append(candidate)

    selected = []
    for group in groups.values():
        scores = [
            score_candidate(
                candidate[3], candidate[0], candidate[1], executable=executable,
                output_dir=output_dir, distance_cutoff=distance_cutoff,
                acc_threshold=acc_threshold, temperature=temperature, timeout=timeout,
            )
            for candidate in group
        ]
        if any(score.status != "scored" for score in scores):
            selected.extend(group)
            continue
        ranked = sorted(zip(group, scores), key=lambda item: item[1].affinity_kcal_mol)
        selected.extend(candidate for candidate, _ in ranked[:top_k])
    return selected
