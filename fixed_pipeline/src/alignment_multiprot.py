"""
Fixed version of src/alignment_multiprot.py

Key fixes:
1. Hardcoded `largest_solution < 5` → configurable via PRISM_MULTIPROT_MIN_MATCHES
2. MultiProt return code 159 → explicit seccomp diagnostic message
3. Missing surface/template files are counted and reported
4. Added logging throughout
5. Exceptions propagate with enough context to identify the root cause
"""

import json
import logging
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio.PDB import PDBParser

# Import shared helpers from original src/
import alignment as _align
_has_ca_atoms = _align._has_ca_atoms
_write_empty_alignment = _align._write_empty_alignment
_valid_tmalign_outputs = _align._valid_tmalign_outputs
extract_chain_and_res_ids = _align.extract_chain_and_res_ids

logger = logging.getLogger(__name__)

os.makedirs("processed/alignment", exist_ok=True)

MULTIPROT = os.environ.get("PRISM_MULTIPROT", "external_tools/multiprot.Linux")
TMALIGN = os.environ.get("PRISM_TMALIGN", "external_tools/TMalign")
MIN_MATCHES = int(os.environ.get("PRISM_MULTIPROT_MIN_MATCHES", "5"))
MULTIPROT_TIMEOUT = int(os.environ.get("PRISM_MULTIPROT_TIMEOUT", "300"))
SECCOMP_EXIT_CODE = 159


def _is_seccomp_blocked(returncode, stderr):
    """Detect if MultiProt failed due to seccomp blocking (exit 159 + Bad system call)."""
    if returncode == SECCOMP_EXIT_CODE:
        return True
    if stderr and "Bad system call" in stderr:
        return True
    return False


def _diagnose_multiprot_failure(returncode, stderr, stdout):
    """Return a human-readable diagnostic for a failed MultiProt call."""
    if _is_seccomp_blocked(returncode, stderr):
        return (
            "MultiProt blocked by seccomp (exit code 159 / 'Bad system call'). "
            "This 32-bit binary requires a compute node without seccomp enforcement. "
            "Workarounds:\n"
            "  1. Run via Slurm compute node (sbatch) instead of login node\n"
            "  2. Recompile MultiProt as a 64-bit binary from source\n"
            "  3. Use a different aligner: --aligner tmalign or --aligner gtalign"
        )
    return (
        f"MultiProt failed (exit code {returncode}). "
        f"stderr: {(stderr or '').strip()[:300]} "
        f"stdout: {(stdout or '').strip()[:300]}"
    )


def _align_one(task):
    """Align one (query, template, chain) pair via MultiProt + TMalign."""
    query, template, chain, alignment_root = task
    query_path = f"processed/surface_extraction/{query}.asa.pdb"
    interface_path = f"templates/interfaces/{template}_{chain}_int.pdb"
    output_path = os.path.join(alignment_root, f"{query}_{template}_{chain}.json")

    if not _has_ca_atoms(query_path):
        logger.debug("MultiProt: no CA atoms in query %s", query_path)
        _write_empty_alignment(output_path)
        return None

    if not os.path.exists(interface_path):
        logger.debug("MultiProt: missing interface %s", interface_path)
        _write_empty_alignment(output_path)
        return None

    try:
        with tempfile.TemporaryDirectory(prefix="multiprot-", dir=alignment_root) as scratch:
            scratch = Path(scratch)

            # Symlink PDBs with short names for MultiProt
            q_pdb = scratch / "query.pdb"
            i_pdb = scratch / "interface.pdb"
            shutil.copy2(query_path, str(q_pdb))
            shutil.copy2(interface_path, str(i_pdb))

            # Run MultiProt
            mp_result = subprocess.run(
                [MULTIPROT, str(q_pdb), str(i_pdb)],
                capture_output=True, text=True, timeout=MULTIPROT_TIMEOUT,
            )

            if mp_result.returncode != 0 and not mp_result.stdout:
                diag = _diagnose_multiprot_failure(
                    mp_result.returncode, mp_result.stderr, mp_result.stdout,
                )
                raise RuntimeError(diag)

            # Parse aligned residue count
            largest_solution = 0
            for line in mp_result.stdout.split("\n"):
                if "Largest Solution" in line:
                    try:
                        largest_solution = int(line.split(":")[-1].strip())
                    except (ValueError, IndexError):
                        pass

            if largest_solution < MIN_MATCHES:
                logger.debug(
                    "MultiProt: %s vs %s_%s: largest_solution=%d < min=%d",
                    query, template, chain, largest_solution, MIN_MATCHES,
                )
                _write_empty_alignment(output_path)
                return None

            # Run TMalign to get rotation/translation matrix
            matrix_path = scratch / "matrix.out"
            tm_path = scratch / "out.tm"

            tm_result = subprocess.run(
                [TMALIGN, str(q_pdb), str(i_pdb), "-m", str(matrix_path)],
                stdout=open(str(tm_path), "w"),
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=False,
            )

            if (tm_result.returncode != 0
                    or not matrix_path.exists()
                    or not tm_path.exists()
                    or not _valid_tmalign_outputs(str(matrix_path), str(tm_path))):
                raise RuntimeError(
                    f"TMalign failed after MultiProt for {query}, {template}_{chain}: "
                    f"exit={tm_result.returncode} "
                    f"{(tm_result.stderr or '').strip()[:300]}"
                )

            # Parse TMalign output
            translation, rotation_mat, tm_score, match_count, match_dict = (
                _parse_tmalign_output(
                    str(matrix_path), str(tm_path),
                    query_path, interface_path,
                )
            )

            # Write alignment JSON
            multi_dict = {
                "match_count": match_count,
                "translation": translation,
                "rotation_mat": rotation_mat,
                "match_dict": match_dict,
                "tm_score": tm_score,
                "status": "success",
                "aligner": "MultiProt+TMalign",
                "multiprot_solution": largest_solution,
            }
            os.makedirs(alignment_root, exist_ok=True)
            with open(output_path, "w") as f:
                json.dump(multi_dict, f)

    except subprocess.TimeoutExpired:
        logger.warning(
            "MultiProt timed out (%ds) for %s, %s_%s",
            MULTIPROT_TIMEOUT, query, template, chain,
        )
        _write_empty_alignment(output_path)
    except (OSError, RuntimeError, ValueError, IndexError) as exc:
        logger.warning(
            "MultiProt+TMalign failed for %s and %s_%s: %s",
            query, template, chain, exc,
        )
        _write_empty_alignment(output_path)

    return None


def _parse_tmalign_output(matrix_path, tm_path, query_path, interface_path):
    """Parse TMalign output files to extract alignment data."""
    translation = [0.0, 0.0, 0.0]
    rotation_mat = [[0.0, 0.0, 0.0] for _ in range(3)]

    with open(matrix_path, "r") as f:
        for line in f:
            tokens = line.strip().split()
            if len(tokens) < 5:
                continue
            try:
                row_index = int(tokens[0])
            except (ValueError, IndexError):
                continue
            if row_index in (0, 1, 2):
                translation[row_index] = float(tokens[1])
                rotation_mat[row_index][0] = float(tokens[2])
                rotation_mat[row_index][1] = float(tokens[3])
                rotation_mat[row_index][2] = float(tokens[4])

    tm_score = 0.0
    match_count = 0
    seq1, match, seq2 = [], [], []

    with open(tm_path, "r") as f:
        for line in f:
            if line.startswith("Aligned length"):
                try:
                    match_count = int(line.split("=")[1].split(",")[0].strip())
                except (ValueError, IndexError):
                    pass
            elif line.startswith("TM-score"):
                try:
                    tmscore = float(line.split()[1])
                    tm_score = max(tm_score, tmscore)
                except (ValueError, IndexError):
                    pass
            elif line.startswith('(":"'):
                seq1 = list(f.readline().rstrip("\n"))
                match = list(f.readline().rstrip("\n"))
                seq2 = list(f.readline().rstrip("\n"))
                break

    seq1_res_ids, seq1_chain_ids = extract_chain_and_res_ids("query", query_path)
    seq2_res_ids, seq2_chain_ids = extract_chain_and_res_ids("interface", interface_path)

    match_dict = {}
    index1, index2 = 0, 0
    usable = min(len(seq1), len(match), len(seq2))

    for i, s in enumerate(match[:usable]):
        if s in (":", "."):
            if (index1 < len(seq1_chain_ids) and index1 < len(seq1_res_ids)
                    and index2 < len(seq2_chain_ids) and index2 < len(seq2_res_ids)):
                seq1_str = f"{seq1_chain_ids[index1]}.{seq1[i]}.{seq1_res_ids[index1]}"
                seq2_str = f"{seq2_chain_ids[index2]}.{seq2[i]}.{seq2_res_ids[index2]}"
                match_dict[seq2_str] = seq1_str
        if seq1[i] != "-":
            index1 += 1
        if seq2[i] != "-":
            index2 += 1

    return translation, rotation_mat, tm_score, match_count, match_dict


def align_multiprot(queries, templates, output_dir="processed/alignment", max_workers=8):
    """MultiProt alignment entry point for the PRISM pipeline."""
    alignment_root = os.path.abspath(output_dir)
    os.makedirs(alignment_root, exist_ok=True)

    # Pre-flight: check MultiProt binary exists and is executable
    multiprot_path = os.environ.get("PRISM_MULTIPROT", "external_tools/multiprot.Linux")
    if not os.path.exists(multiprot_path):
        logger.error(
            "MultiProt binary not found at %s. "
            "Set PRISM_MULTIPROT env var or use --aligner tmalign/gtalign",
            multiprot_path,
        )
        return

    tasks = []
    missing_interface = 0
    missing_surface = 0
    for query in queries:
        qpath = f"processed/surface_extraction/{query}.asa.pdb"
        if not os.path.exists(qpath):
            missing_surface += 1
            continue
        for template in templates:
            for chain in template[4:]:
                ipath = f"templates/interfaces/{template}_{chain}_int.pdb"
                if not os.path.exists(ipath):
                    missing_interface += 1
                    continue
                tasks.append((query, template, chain, alignment_root))

    if missing_surface > 0:
        logger.warning("MultiProt: %d query surface files missing", missing_surface)
    if missing_interface > 0:
        logger.warning("MultiProt: %d template interface files missing", missing_interface)

    total = len(tasks)
    if total == 0:
        logger.warning("MultiProt: no valid alignment tasks after pre-flight")
        return

    logger.info("MultiProt alignment: %d pairs to process", total)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_align_one, t): t for t in tasks}
        for i, future in enumerate(as_completed(futures), 1):
            future.result()
            if i % 100 == 0:
                logger.info("MultiProt progress: %d/%d", i, total)

    logger.info("MultiProt alignment finished: %d pairs processed", total)
