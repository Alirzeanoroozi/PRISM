"""PRODIGY-backed candidate ranking for the opt-in PRISM ranking stage.

PRODIGY scores a two-part complex, while PRISM keeps transformed receptor and
ligand structures in separate PDB files.  This adapter combines each pair in
an isolated, chain-renamed input and retains the exact input, command, and
result under the run's ranking output directory.

The adapter deliberately does not make PRODIGY a default or silently discard
an unscorable candidate.  A group is returned unchanged when any candidate in
that receptor/ligand group cannot be scored.
"""

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


_CHAIN_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
_AFFINITY_LINE = re.compile(r"(?:^|\s)(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$")


@dataclass(frozen=True)
class ProdigyScore:
    """One retained PRODIGY score observation."""

    left_pdb: str
    right_pdb: str
    input_pdb: str
    status: str
    affinity_kcal_mol: float | None
    return_code: int | None
    command: tuple[str, ...]
    stdout_path: str
    stderr_path: str
    input_sha256: str | None
    error: str | None = None


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _chain_ids(path: Path) -> list[str]:
    """Return unique chain identifiers in first-seen order."""

    seen: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                chain = line[21]
                if chain not in seen:
                    seen.append(chain)
    return seen


def _rewrite_chain(line: str, chain: str) -> str:
    if len(line) <= 21:
        return line
    return f"{line[:21]}{chain}{line[22:]}"


def combine_pair_pdbs(left_pdb: str, right_pdb: str, output_pdb: str) -> tuple[str, str]:
    """Combine two transformed PDBs and return PRODIGY selection groups.

    The original files may use overlapping chain identifiers.  Each side is
    therefore rewritten to a unique one-character chain namespace before the
    files are concatenated.  Only coordinate records and TER records are
    copied; the result is a deliberately minimal complex input.
    """

    left = Path(left_pdb)
    right = Path(right_pdb)
    destination = Path(output_pdb)
    if not left.is_file() or not right.is_file():
        raise FileNotFoundError(f"PRODIGY inputs must exist: {left_pdb}, {right_pdb}")

    left_chains = _chain_ids(left)
    right_chains = _chain_ids(right)
    labels = iter(_CHAIN_ALPHABET)
    left_map = {chain: next(labels) for chain in left_chains}
    right_map = {chain: next(labels) for chain in right_chains}
    if not left_chains or not right_chains:
        raise ValueError("PRODIGY requires coordinate records on both pair sides")

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as out:
        out.write("REMARK PRISM PRODIGY ranking input; chain IDs were namespaced by side\n")
        for source, mapping in ((left, left_map), (right, right_map)):
            with source.open("r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith(("ATOM", "HETATM", "TER")):
                        old_chain = line[21] if len(line) > 21 else ""
                        if old_chain in mapping:
                            out.write(_rewrite_chain(line, mapping[old_chain]))
        out.write("END\n")

    left_selection = ",".join(left_map.values())
    right_selection = ",".join(right_map.values())
    return left_selection, right_selection


def parse_affinity(stdout: str) -> float:
    """Parse the numeric affinity emitted by ``prodigy -q``."""

    values: list[float] = []
    for line in stdout.splitlines():
        match = _AFFINITY_LINE.search(line.strip())
        if match:
            values.append(float(match.group(1)))
    if not values:
        raise ValueError(f"PRODIGY quiet output contained no affinity: {stdout!r}")
    affinity = values[-1]
    if not math.isfinite(affinity):
        raise ValueError(f"PRODIGY affinity is not finite: {affinity!r}")
    return affinity


def _safe_stem(left_pdb: str, right_pdb: str) -> str:
    raw = f"{Path(left_pdb).stem}__{Path(right_pdb).stem}"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", raw)[:180]


def score_candidate(
    left_pdb: str,
    right_pdb: str,
    *,
    executable: str = "prodigy",
    output_dir: str = "processed/ranking/prodigy",
    distance_cutoff: float = 5.5,
    acc_threshold: float = 0.05,
    temperature: float = 25.0,
    timeout: float = 120.0,
) -> ProdigyScore:
    """Run PRODIGY once and retain its complete machine-readable evidence."""

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    stem = _safe_stem(left_pdb, right_pdb)
    input_path = root / f"{stem}.pdb"
    stdout_path = root / f"{stem}.stdout.txt"
    stderr_path = root / f"{stem}.stderr.txt"
    metadata_path = root / f"{stem}.json"
    command = tuple(
        shlex.split(executable)
        + [
            "-q",
            "--distance-cutoff",
            str(distance_cutoff),
            "--acc-threshold",
            str(acc_threshold),
            "--temperature",
            str(temperature),
        ]
    )

    result: ProdigyScore
    try:
        left_selection, right_selection = combine_pair_pdbs(left_pdb, right_pdb, str(input_path))
        command = tuple(command + (str(input_path), "--selection", left_selection, right_selection))
        completed = subprocess.run(
            list(command),
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        stdout_path.write_text(completed.stdout, encoding="utf-8")
        stderr_path.write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(f"PRODIGY exited with return code {completed.returncode}")
        affinity = parse_affinity(completed.stdout)
        result = ProdigyScore(
            left_pdb=left_pdb,
            right_pdb=right_pdb,
            input_pdb=str(input_path),
            status="scored",
            affinity_kcal_mol=affinity,
            return_code=completed.returncode,
            command=command,
            stdout_path=str(stdout_path),
            stderr_path=str(stderr_path),
            input_sha256=_sha256_file(input_path),
        )
    except subprocess.TimeoutExpired as exc:
        stdout_path.write_text(str(exc.stdout or ""), encoding="utf-8")
        stderr_path.write_text(str(exc.stderr or ""), encoding="utf-8")
        result = ProdigyScore(
            left_pdb, right_pdb, str(input_path), "timeout", None, None, command,
            str(stdout_path), str(stderr_path), _sha256_file(input_path) if input_path.exists() else None,
            f"PRODIGY timed out after {timeout}s",
        )
    except Exception as exc:
        if not stdout_path.exists():
            stdout_path.write_text("", encoding="utf-8")
        if not stderr_path.exists():
            stderr_path.write_text("", encoding="utf-8")
        result = ProdigyScore(
            left_pdb, right_pdb, str(input_path), "failed", None, None, command,
            str(stdout_path), str(stderr_path), _sha256_file(input_path) if input_path.exists() else None,
            f"{type(exc).__name__}: {exc}",
        )

    metadata_path.write_text(json.dumps(asdict(result), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def _group_key(left_pdb: str, right_pdb: str) -> tuple[str, str]:
    """Group candidates by receptor/ligand parsed from the left output name."""

    from .candidate_selector import _parse_output_pdb_path

    receptor, ligand, _, _ = _parse_output_pdb_path(left_pdb)
    return (receptor or left_pdb, ligand or right_pdb)


def select_top_candidates_with_prodigy(
    passed_pairs: list[tuple[str, str]],
    *,
    top_k: int = 5,
    executable: str = "prodigy",
    output_dir: str = "processed/ranking/prodigy",
    distance_cutoff: float = 5.5,
    acc_threshold: float = 0.05,
    temperature: float = 25.0,
    timeout: float = 120.0,
) -> list[tuple[str, str]]:
    """Select low-affinity (more favorable) candidates per receptor/ligand."""

    if top_k < 1:
        raise ValueError("top_k must be positive")
    if not passed_pairs:
        return []
    if shutil.which(shlex.split(executable)[0]) is None and not Path(shlex.split(executable)[0]).exists():
        raise FileNotFoundError(
            f"PRODIGY executable not found: {executable!r}; install prodigy-prot in a separate environment or set --prodigy-executable"
        )

    groups: dict[tuple[str, str], list[tuple[str, str]]] = {}
    for pair in passed_pairs:
        groups.setdefault(_group_key(*pair), []).append(pair)

    selected: list[tuple[str, str]] = []
    for candidates in groups.values():
        scores = [
            score_candidate(
                left,
                right,
                executable=executable,
                output_dir=output_dir,
                distance_cutoff=distance_cutoff,
                acc_threshold=acc_threshold,
                temperature=temperature,
                timeout=timeout,
            )
            for left, right in candidates
        ]
        if any(score.status != "scored" for score in scores):
            # Preserve all candidates if PRODIGY cannot score the complete
            # group; ranking must not turn a partial external-tool run into a
            # false scientific success.
            selected.extend(candidates)
            continue
        ranked = sorted(
            zip(candidates, scores),
            key=lambda item: (float(item[1].affinity_kcal_mol), item[0][0], item[0][1]),
        )
        selected.extend(pair for pair, _ in ranked[:top_k])
    return selected
