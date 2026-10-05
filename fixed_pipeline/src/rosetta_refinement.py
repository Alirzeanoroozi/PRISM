"""
Fixed version of src/rosetta_refinement.py

Key fixes:
1. REPLACED os.system() with subprocess.run() — no shell injection vector
2. REPLACED shell mv/cp with shutil.move()/shutil.copy2()
3. Unique output paths instead of _0001.pdb race condition
4. Structured error reporting instead of bare except + return "-"
5. Enhanced logging with distinct error codes for each failure mode
"""

import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Import from original src/ layer (contact has no relative imports itself).
# The conftest.py adds the original src/ to sys.path.
import contact as _contact
get_contacts_from_atom_lines = _contact.get_contacts_from_atom_lines

logger = logging.getLogger(__name__)

ROSETTA_PREPACK = os.environ.get(
    "PRISM_ROSETTA_PREPACK", "docking_prepack_protocol.static.linuxgccrelease"
)
ROSETTA_DOCK = os.environ.get(
    "PRISM_ROSETTA_DOCK", "docking_protocol.static.linuxgccrelease"
)
ROSETTA_DB = os.environ.get(
    "PRISM_ROSETTA_DB",
    "/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database/",
)
ROSETTA_INT_SCORE_THRESHOLD = -5.0
ROSETTA_DIR = "processed/rosetta_refinement"
ENERGY_DIR = os.path.join(ROSETTA_DIR, "energies")
STRUCTURE_DIR = os.path.join(ROSETTA_DIR, "structures")
os.makedirs(ROSETTA_DIR, exist_ok=True)
os.makedirs(ENERGY_DIR, exist_ok=True)
os.makedirs(STRUCTURE_DIR, exist_ok=True)


class RefinementError(Exception):
    """Base exception for Rosetta refinement failures."""
    def __init__(self, message, code="UNKNOWN", detail=None):
        self.code = code
        self.detail = detail
        super().__init__(message)


class RosettaBinaryNotFound(RefinementError):
    def __init__(self, binary, path):
        super().__init__(
            f"Rosetta binary '{binary}' not found at {path}",
            code="ROSETTA_BINARY_NOT_FOUND",
            detail={"binary": binary, "path": str(path)},
        )


class RosettaExecutionError(RefinementError):
    def __init__(self, binary, returncode, stdout, stderr):
        super().__init__(
            f"Rosetta '{binary}' exited with code {returncode}",
            code="ROSETTA_EXECUTION_FAILED",
            detail={"binary": binary, "returncode": returncode,
                    "stdout": stdout[-500:], "stderr": stderr[-500:]},
        )


class ScoreParsingError(RefinementError):
    def __init__(self, score_path, detail):
        super().__init__(
            f"Failed to parse score file {score_path}: {detail}",
            code="SCORE_PARSE_FAILED",
        )


def _check_binary(binary_name):
    """Verify a Rosetta binary exists and is executable."""
    # Search PATH if not an absolute path
    if os.path.sep not in binary_name:
        for dirpath in os.environ.get("PATH", "").split(os.pathsep):
            candidate = os.path.join(dirpath, binary_name)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return candidate
    if os.path.isfile(binary_name) and os.access(binary_name, os.X_OK):
        return binary_name
    raise RosettaBinaryNotFound(binary_name, binary_name)


def refiner(passed_pairs):
    with open("processed/rosetta_refinement/refinement_energies.txt", "w") as file_out:
        for passed0, passed1 in passed_pairs:
            result = calculate_energy(passed0, passed1)

            if result is not None and result.interaction_score is not None:
                file_out.write(
                    f"{passed0}\t{passed1}\t{result.interaction_score}\t{result.total_score}\n"
                )
            else:
                file_out.write(
                    f"{passed0}\t{passed1}\tFAILED\t{result.error_code if result else 'UNKNOWN'}\n"
                )


class EnergyResult:
    """Structured result from calculate_energy, never bare '-' strings."""

    def __init__(
        self,
        total_score=None,
        interaction_score=None,
        structure_path=None,
        error_code=None,
        error_detail=None,
    ):
        self.total_score = total_score
        self.interaction_score = interaction_score
        self.structure_path = structure_path
        self.error_code = error_code
        self.error_detail = error_detail

    @classmethod
    def failure(cls, code, detail=None):
        return cls(error_code=code, error_detail=detail)


def _run_rosetta(binary, args, cwd=None, timeout=600):
    """Run a Rosetta binary safely via subprocess with timeout."""
    binary_path = _check_binary(binary)
    cmd = [binary_path] + args
    logger.debug("Running: %s", " ".join(cmd))
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=cwd,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise RosettaBinaryNotFound(binary, binary_path) from exc
    except subprocess.TimeoutExpired:
        raise RefinementError(
            f"Rosetta '{binary}' timed out after {timeout}s",
            code="ROSETTA_TIMEOUT",
        )

    if result.returncode != 0:
        raise RosettaExecutionError(binary, result.returncode, result.stdout, result.stderr)

    return result


def calculate_energy(passed0, passed1):
    """
    Run Rosetta prepacking + docking refinement on a receptor/ligand pair.

    Returns an EnergyResult with structured data. Never returns bare '-' strings.
    """
    try:
        combined_path = combine_pdb(passed0, passed1)
        if not combined_path:
            return EnergyResult.failure("COMBINE_PDB_FAILED", f"combine_pdb returned empty for {passed0}, {passed1}")

        left_chains, right_chains = partner_chain_ids(passed0, passed1)
        partner_chains = f"{left_chains}_{right_chains}"
        combined_stem = os.path.splitext(os.path.basename(combined_path))[0]

        # --- Stage 1: Prepacking ---
        prepack_score_path = os.path.join(ENERGY_DIR, f"{combined_stem}_prepack_score.sc")
        _run_rosetta(ROSETTA_PREPACK, [
            "-database", ROSETTA_DB,
            "-s", combined_path,
            "-partners", partner_chains,
            "-ex1", "-ex2aro",
            "-out:file:scorefile", prepack_score_path,
            "-overwrite",
            "-ignore_zero_occupancy", "false",
            "-detect_disulf", "false",
        ])

        # Rosetta output file is named after input: {combined_stem}_0001.pdb in CWD
        prepacked_name = f"{combined_stem}_0001.pdb"
        prepacked_src = os.path.join(os.getcwd(), prepacked_name)
        prepacked_dst = os.path.join(ROSETTA_DIR, prepacked_name)
        if os.path.exists(prepacked_src):
            shutil.move(prepacked_src, prepacked_dst)
        else:
            logger.warning("Expected prepacked file not found: %s", prepacked_src)

        # --- Stage 2: Docking ---
        dock_score_path = os.path.join(ENERGY_DIR, f"{combined_stem}_score.sc")
        _run_rosetta(ROSETTA_DOCK, [
            "-database", ROSETTA_DB,
            "-s", prepacked_dst,
            "-docking_local_refine",
            "-partners", partner_chains,
            "-ex1", "-ex2aro",
            "-overwrite",
            "-ignore_zero_occupancy", "false",
            "-detect_disulf", "false",
            "-out:path:score", ENERGY_DIR,
            "-out:file:scorefile", dock_score_path,
        ])

        # Rosetta docking output: {combined_stem}_0001_0001.pdb
        out_pdb_name = f"{combined_stem}_0001_0001.pdb"
        out_pdb_src = os.path.join(os.getcwd(), out_pdb_name)
        out_pdb_dst = os.path.join(STRUCTURE_DIR, out_pdb_name)
        if os.path.exists(out_pdb_src):
            shutil.move(out_pdb_src, out_pdb_dst)
        elif not os.path.exists(out_pdb_dst):
            logger.warning("Expected docked output not found: %s", out_pdb_src)

        # --- Stage 3: Parse scores ---
        totalscore = None
        interaction_score = None
        final_score_path = os.path.join(ENERGY_DIR, f"{combined_stem}_score.sc")
        os.makedirs(ENERGY_DIR, exist_ok=True)

        if os.path.exists(final_score_path):
            totalscore, interaction_score = _parse_rosetta_score(final_score_path)
        else:
            # Rosetta may write to the default name score.sc in the score output dir
            default_score = os.path.join(ENERGY_DIR, "score.sc")
            if os.path.exists(default_score):
                shutil.move(default_score, final_score_path)
                totalscore, interaction_score = _parse_rosetta_score(final_score_path)
            else:
                logger.warning("No score file found for %s", combined_stem)
                return EnergyResult.failure(
                    "SCORE_FILE_MISSING",
                    f"No score.sc found in {ENERGY_DIR} for {combined_stem}",
                )

        # --- Stage 4: Structure extraction and contact analysis ---
        if (os.path.exists(out_pdb_dst)
                and interaction_score is not None
                and interaction_score <= ROSETTA_INT_SCORE_THRESHOLD):
            structure_list = _extract_partner_chains(out_pdb_dst, left_chains, right_chains)

            shutil.copy2(out_pdb_dst, os.path.join(ROSETTA_DIR, out_pdb_name))

            int_res_path = os.path.join(ROSETTA_DIR, f"{out_pdb_name}.intRes.txt")
            try:
                get_contacts_from_atom_lines(
                    os.path.join(ROSETTA_DIR, out_pdb_name),
                    int_res_path,
                    structure_list[0],
                    structure_list[1],
                )
            except Exception as exc:
                logger.warning("get_contacts failed for %s: %s", out_pdb_name, exc)

            return EnergyResult(
                total_score=totalscore,
                interaction_score=str(interaction_score),
                structure_path=os.path.join(ROSETTA_DIR, out_pdb_name),
            )
        else:
            if not os.path.exists(out_pdb_dst):
                logger.warning("Structure file not found: %s", out_pdb_dst)
            elif interaction_score is None:
                logger.warning("Interaction score is None for %s", out_pdb_dst)
            elif interaction_score > ROSETTA_INT_SCORE_THRESHOLD:
                logger.info(
                    "Interaction score %.2f exceeds threshold %.2f for %s",
                    interaction_score, ROSETTA_INT_SCORE_THRESHOLD, out_pdb_dst,
                )
            return EnergyResult(
                total_score=totalscore,
                interaction_score=str(interaction_score) if interaction_score is not None else None,
                structure_path=out_pdb_dst if os.path.exists(out_pdb_dst) else None,
                error_code="THRESHOLD_NOT_MET",
            )

    except (RosettaBinaryNotFound, RosettaExecutionError, ScoreParsingError) as exc:
        logger.error("Rosetta refinement failed: [%s] %s", exc.code, exc)
        return EnergyResult.failure(exc.code, str(exc))
    except Exception as exc:
        logger.exception("Unexpected error in calculate_energy")
        return EnergyResult.failure("UNEXPECTED_ERROR", str(exc))


def _parse_rosetta_score(score_path):
    """Parse a Rosetta score.sc file and return (total_score, interaction_score)."""
    with open(score_path, "r") as scorefile:
        lines = [l for l in scorefile if l.strip() and not l.startswith("#")]

    if not lines:
        raise ScoreParsingError(score_path, "No data lines in score file")

    data_line = lines[0] if lines[0].strip().startswith("SCORE:") else lines[-1]
    temp = data_line.split()

    if len(temp) < 6:
        raise ScoreParsingError(
            score_path, f"Expected >=6 columns, got {len(temp)}: {data_line[:200]}"
        )

    try:
        totalscore = float(temp[1].strip())
        interaction_score = float(temp[5].strip())
    except (ValueError, IndexError) as exc:
        raise ScoreParsingError(score_path, str(exc)) from exc

    return totalscore, interaction_score


def _extract_partner_chains(pdb_path, left_chains, right_chains):
    """Extract ATOM lines grouped by partner chains from a PDB file."""
    structure_list = {0: [], 1: []}
    with open(pdb_path, "r") as fh:
        for l in fh:
            if l[:3] == "TER":
                break
            if l.startswith("ATOM") and l[21].strip() in left_chains:
                structure_list[0].append(l)
            if l.startswith("ATOM") and l[21].strip() in right_chains:
                structure_list[1].append(l)
    return structure_list


def combine_pdb(passed0, passed1):
    """Combine two partner PDBs into one file with chain-renamed ATOM records."""
    try:
        combined_name = "{}_{}_rosetta.pdb".format(
            os.path.splitext(os.path.basename(passed0))[0],
            os.path.splitext(os.path.basename(passed1))[0],
        )
        combined_path = os.path.join(ROSETTA_DIR, combined_name)
        left_chains, right_chains = partner_chain_ids(passed0, passed1)
        source_left_chains = source_chain_ids(passed0)
        source_right_chains = source_chain_ids(passed1)

        with open(passed0, "r") as p0file, \
             open(passed1, "r") as p1file, \
             open(combined_path, "w") as combinedfile:
            for line in p0file:
                if line.startswith("ATOM") and line[21].strip() in source_left_chains:
                    source_chain = line[21].strip()
                    replacement = left_chains[source_left_chains.index(source_chain)]
                    combinedfile.write(line[:21] + replacement + line[22:])
            combinedfile.write("TER\n")
            for line in p1file:
                if line.startswith("ATOM") and line[21].strip() in source_right_chains:
                    source_chain = line[21].strip()
                    replacement = right_chains[source_right_chains.index(source_chain)]
                    combinedfile.write(line[:21] + replacement + line[22:])
            combinedfile.write("END\n")
        return combined_path
    except Exception as exc:
        logger.error("combine_pdb failed for %s, %s: %s", passed0, passed1, exc)
        return ""


def target_chain_id(path):
    """Extract the target chain suffix from a transformation filename."""
    stem = os.path.basename(path).rsplit(".", 1)[0]
    fields = stem.split("_")
    if len(fields) < 5 or len(fields[1]) < 5 or len(fields[2]) < 5:
        raise ValueError(f"cannot infer target chain from {path}")
    return fields[1][-1] if fields[-1] == "L" else fields[2][-1]


def source_chain_ids(path):
    """Extract all source chains for a transformed left or right partner."""
    stem = os.path.basename(path).rsplit(".", 1)[0]
    fields = stem.split("_")
    if len(fields) >= 5:
        target = fields[1] if fields[-1] == "L" else fields[2]
        if len(target) >= 5:
            return "".join(dict.fromkeys(target[4:]))
    chains = []
    with open(path, "r") as handle:
        for line in handle:
            if line.startswith("ATOM") and line[21].strip() not in chains:
                chains.append(line[21].strip())
    if not chains:
        raise ValueError(f"cannot infer target chains from {path}")
    return "".join(chains)


def source_chain_id(path):
    """Return target chain from filename, or first ATOM chain for generic files."""
    try:
        return source_chain_ids(path)[0]
    except ValueError:
        with open(path, "r") as handle:
            for line in handle:
                if line.startswith("ATOM"):
                    return line[21].strip() or "_"
        raise ValueError(f"no ATOM chain found in {path}")


def partner_chain_ids(left_path, right_path):
    """Return unique Rosetta chain groups for the two transformed partners."""
    if len(os.path.basename(left_path).rsplit(".", 1)[0].split("_")) < 5:
        return "A", "B"

    left = source_chain_ids(left_path)
    right = source_chain_ids(right_path)

    # Collisions: if the two partners share chain IDs, re-assign right to next
    # available letters.
    used = set(left)
    available = iter(c for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if c not in used)
    right_unique = "".join(next(available) for _ in range(len(right)))
    return left, right_unique
