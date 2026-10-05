import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.PDB import PDBParser
import json


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
        print(f"TM-align failed for {protein} and {template}_{chain}: {exc}")
        _write_empty_alignment(output_path)
    return None

def align(queries, templates, output_dir="processed/alignment_tmalign"):
    alignment_root = os.path.abspath(output_dir)
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
            print(f"TM-align progress: {i}/{total}")


def parse_tmalign(protein_path, interface_path, protein, template, chain, matrix_path=None, tm_path=None, output_dir="processed/alignment_tmalign"):
    matrix_path = matrix_path or "processed/alignment_tmalign/matrix.out"
    tm_path = tm_path or "processed/alignment_tmalign/out.tm"
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
            # Each matrix row is indicated by 0/1/2 for x/y/z, but original code uses 1/2/3
            # From matrix.out, line numbers are actually 0, 1, 2 in first column. Adjust accordingly.
            # tokens[0]: row (0,1,2); tokens[1]: t[m]; tokens[2-4]: u[m][0:2]
            if row_index in (0, 1, 2):
                translation[row_index] = float(tokens[1])
                rotation_mat[row_index][0] = float(tokens[2])
                rotation_mat[row_index][1] = float(tokens[3])
                rotation_mat[row_index][2] = float(tokens[4])

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
                # Next 3 lines: seq1 (Chain_1), match line, seq2 (Chain_2)
                seq1 = list(tm_file.readline().rstrip("\n"))
                match = list(tm_file.readline().rstrip("\n"))
                seq2 = list(tm_file.readline().rstrip("\n"))
                break

        seq1_res_ids, seq1_chain_ids = extract_chain_and_res_ids("protein", protein_path)
        seq2_res_ids, seq2_chain_ids = extract_chain_and_res_ids("interface", interface_path)
        index1 = 0
        index2 = 0
        # TM-align can emit an alignment whose sequence line is inconsistent
        # with residue records when the input PDB contains duplicate residue
        # numbers.  Preserve the valid prefix and record the truncation rather
        # than discarding the entire candidate with IndexError.
        usable = min(len(seq1), len(match), len(seq2))
        mapping_status = "success" if len({len(seq1), len(match), len(seq2)}) == 1 else "mapping_truncated"
        for i, s in enumerate(match[:usable]):
            if s == ":" or s == ".":
                if index1 < len(seq1_chain_ids) and index1 < len(seq1_res_ids) and index2 < len(seq2_chain_ids) and index2 < len(seq2_res_ids):
                    seq1_str = seq1_chain_ids[index1] + "." + seq1[i] + "." + seq1_res_ids[index1]
                    seq2_str = seq2_chain_ids[index2] + "." + seq2[i] + "." + seq2_res_ids[index2]
                    match_dict[seq2_str] = seq1_str
                else:
                    mapping_status = "mapping_truncated"
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
        "alignment_lengths": {"seq1": len(seq1), "match": len(match), "seq2": len(seq2)},
        "aligner": "TMalign",
    }
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, f"{protein}_{template}_{chain}.json"), "w") as f:
        json.dump(multi_dict, f)

def extract_chain_and_res_ids(name, path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(name, path)
    residue_ids = []
    chain_id = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    residue_ids.append(str(residue.id[1]))
                    chain_id.append(chain.id)
    return residue_ids, chain_id


def _has_ca_atoms(path):
    if not os.path.exists(path):
        return False
    with open(path, "r") as handle:
        return any(line.startswith("ATOM") and line[13:15].strip() == "CA" for line in handle)


def _valid_tmalign_outputs(matrix_path, tm_path):
    """Require the minimal records consumed by ``parse_tmalign``."""
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
    return matrix_rows == {"0", "1", "2"} and "Aligned length" in tm_text and "TM-score" in tm_text and alignment_marker


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

if __name__ == "__main__":
    protein = "1a28"
    template = "1a28AB"
    chains = ["A", "B"]
    align([protein], [template])
    for chain in chains:
        data = json.load(open(f"processed/alignment_tmalign/{protein}_{template}_{chain}.json", "r"))
        print(f"--- {protein}_{template}_{chain} ---")
        print(data)
        print()
