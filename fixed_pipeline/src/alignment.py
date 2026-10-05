"""
Fixed version of src/alignment.py

Key fixes:
1. Mapping truncation status is now exposed as a public constant for
   downstream transformer validation.
2. Alignment JSON includes a human-readable warning when truncated.
3. extract_chain_and_res_ids is re-exported for transformer.py imports.
4. Worker count accepts env var with same convention as other modules.
"""

import json
import logging
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio.PDB import PDBParser

logger = logging.getLogger(__name__)

os.makedirs("processed/alignment", exist_ok=True)

# Public constant: transformer.py can check this to detect truncated mappings.
MAPPING_TRUNCATED_STATUS = "mapping_truncated"
MAPPING_SUCCESS_STATUS = "success"


def iter_bounded_results(tasks, worker, workers, max_pending=None):
    """Yield one result per task while bounding submitted-but-unfinished work."""
    iterator = iter(tasks)
    pending_limit = max_pending or max(1, 2 * int(workers))
    with ThreadPoolExecutor(max_workers=int(workers)) as executor:
        pending = {}
        exhausted = False
        while pending or not exhausted:
            while not exhausted and len(pending) < pending_limit:
                try:
                    task = next(iterator)
                except StopIteration:
                    exhausted = True
                    break
                future = executor.submit(worker, task)
                pending[future] = task
            if not pending:
                break
            future = next(iter(as_completed(tuple(pending))))
            task = pending.pop(future)
            yield task, future.result()


def _align_one(args):
    """Align one (protein, template, chain) pair for parallel dispatch."""
    protein, template, chain, alignment_root, aligner = args
    protein_path = f"processed/surface_extraction/{protein}.asa.pdb"
    interface_path = f"templates/interfaces/{template}_{chain}_int.pdb"
    output_path = os.path.join(alignment_root, f"{protein}_{template}_{chain}.json")
    if not _has_ca_atoms(protein_path) or not os.path.exists(interface_path):
        _write_empty_alignment(output_path)
        return None
    try:
        with tempfile.TemporaryDirectory(prefix="tmalign-", dir=alignment_root) as scratch:
            matrix_path = os.path.join(scratch, "matrix.out")
            output_path_tm = os.path.join(scratch, "out.tm")
            command = [aligner, protein_path, interface_path, "-m", matrix_path]
            with open(output_path_tm, "w") as tm_output:
                result = subprocess.run(
                    command,
                    stdout=tm_output,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                    check=False,
                )
            if (
                result.returncode != 0
                or not os.path.isfile(matrix_path)
                or not os.path.isfile(output_path_tm)
                or not _valid_tmalign_outputs(matrix_path, output_path_tm)
            ):
                raise RuntimeError(
                    f"TM-align exit={result.returncode}: {(result.stderr or '').strip()[:500]}"
                )
            parse_tmalign(protein_path, interface_path, protein, template, chain,
                          matrix_path, output_path_tm, alignment_root)
    except (OSError, RuntimeError, ValueError, IndexError) as exc:
        logger.warning("TM-align failed for %s and %s_%s: %s",
                       protein, template, chain, exc)
        _write_empty_alignment(output_path)
    return None


def align(queries, templates):
    alignment_root = os.path.abspath("processed/alignment")
    os.makedirs(alignment_root, exist_ok=True)
    aligner = os.environ.get("PRISM_TMALIGN", "external_tools/TMalign")
    max_workers = int(os.environ.get("PRISM_TMALIGN_WORKERS", "8"))
    tasks = (
        (protein, template, chain, alignment_root, aligner)
        for protein in queries
        for template in templates
        for chain in template[4:]
    )
    total = len(queries) * sum(len(template[4:]) for template in templates)
    for i, (_, _) in enumerate(iter_bounded_results(tasks, _align_one, max_workers), 1):
        if i % 100 == 0:
            logger.info("TM-align progress: %d/%d", i, total)


def parse_tmalign(
    protein_path, interface_path, protein, template, chain,
    matrix_path=None, tm_path=None, output_dir="processed/alignment"
):
    matrix_path = matrix_path or "processed/alignment/matrix.out"
    tm_path = tm_path or "processed/alignment/out.tm"

    # --- Parse rotation matrix ---
    with open(matrix_path, "r") as matrix_file:
        translation = [0.0, 0.0, 0.0]
        rotation_mat = [[0.0, 0.0, 0.0] for _ in range(3)]
        for line in matrix_file:
            tokens = line.strip().split()
            if len(tokens) < 5:
                if line.startswith(" Code"):
                    break
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

    # --- Parse alignment ---
    with open(tm_path, "r") as tm_file:
        match_dict = {}
        tm_score_1 = 0.0
        tm_score_2 = 0.0
        match_count = 0
        seq1 = []
        match = []
        seq2 = []

        for line in tm_file:
            if line.startswith("Aligned length"):
                match_count = int(line.split("=")[1].split(",")[0].strip())
            elif line.startswith("TM-score"):
                tmscore = float(line.split()[1])
                if "Chain_1" in line:
                    tm_score_1 = tmscore
                elif "Chain_2" in line:
                    tm_score_2 = tmscore
            elif line.startswith('(":"'):
                seq1 = list(tm_file.readline().rstrip("\n"))
                match = list(tm_file.readline().rstrip("\n"))
                seq2 = list(tm_file.readline().rstrip("\n"))
                break

        seq1_res_ids, seq1_chain_ids = extract_chain_and_res_ids("protein", protein_path)
        seq2_res_ids, seq2_chain_ids = extract_chain_and_res_ids("interface", interface_path)
        index1 = 0
        index2 = 0

        usable = min(len(seq1), len(match), len(seq2))
        mapping_status = (
            MAPPING_SUCCESS_STATUS
            if len({len(seq1), len(match), len(seq2)}) == 1
            else MAPPING_TRUNCATED_STATUS
        )
        if mapping_status == MAPPING_TRUNCATED_STATUS:
            logger.warning(
                "Alignment mapping truncated for %s_%s_%s: "
                "seq1=%d match=%d seq2=%d (min=%d)",
                protein, template, chain,
                len(seq1), len(match), len(seq2), usable,
            )

        for i, s in enumerate(match[:usable]):
            if s in (":", "."):
                if (index1 < len(seq1_chain_ids) and index1 < len(seq1_res_ids)
                        and index2 < len(seq2_chain_ids) and index2 < len(seq2_res_ids)):
                    seq1_str = (
                        seq1_chain_ids[index1] + "." + seq1[i] + "."
                        + seq1_res_ids[index1]
                    )
                    seq2_str = (
                        seq2_chain_ids[index2] + "." + seq2[i] + "."
                        + seq2_res_ids[index2]
                    )
                    match_dict[seq2_str] = seq1_str
                else:
                    if mapping_status == MAPPING_SUCCESS_STATUS:
                        logger.warning(
                            "Index out of bounds during match_dict build for "
                            "%s_%s_%s (index1=%d/%d, index2=%d/%d)",
                            protein, template, chain,
                            index1, len(seq1_chain_ids),
                            index2, len(seq2_chain_ids),
                        )
                    mapping_status = MAPPING_TRUNCATED_STATUS
            if seq1[i] != "-":
                index1 += 1
            if seq2[i] != "-":
                index2 += 1

    multi_dict = {
        "match_count": match_count,
        "translation": translation,
        "rotation_mat": rotation_mat,
        "match_dict": match_dict,
        "tm_score": max(tm_score_1, tm_score_2),
        "status": mapping_status,
        "alignment_lengths": {
            "seq1": len(seq1),
            "match": len(match),
            "seq2": len(seq2),
        },
    }
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, f"{protein}_{template}_{chain}.json"), "w") as f:
        json.dump(multi_dict, f)


def extract_chain_and_res_ids(name, path):
    """Extract sequential residue IDs and chain IDs from a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(name, path)
    residue_ids = []
    chain_id = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    residue_ids.append(str(residue.id[1]))
                    chain_id.append(chain.id)
    return residue_ids, chain_id


def _has_ca_atoms(path):
    if not os.path.exists(path):
        return False
    with open(path, "r") as handle:
        return any(line.startswith("ATOM") and line[13:15].strip() == "CA" for line in handle)


def _valid_tmalign_outputs(matrix_path, tm_path):
    """Require the minimal records consumed by parse_tmalign."""
    try:
        with open(matrix_path, "r") as matrix_handle, open(tm_path, "r") as tm_handle:
            matrix_text = matrix_handle.read()
            tm_text = tm_handle.read()
    except OSError:
        return False
    matrix_rows = set()
    for line in matrix_text.splitlines():
        tokens = line.split()
        if len(tokens) >= 5 and tokens[0] in {"0", "1", "2"}:
            matrix_rows.add(tokens[0])
    alignment_marker = any(line.startswith('(":') for line in tm_text.splitlines())
    return (matrix_rows == {"0", "1", "2"}
            and "Aligned length" in tm_text
            and "TM-score" in tm_text
            and alignment_marker)


def _write_empty_alignment(path):
    with open(path, "w") as handle:
        json.dump({
            "match_count": 0,
            "translation": [0.0, 0.0, 0.0],
            "rotation_mat": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "match_dict": {},
            "tm_score": 0.0,
            "status": "alignment_unavailable",
        }, handle)
