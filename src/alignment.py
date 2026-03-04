import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.PDB import PDBParser
import json
import csv
from tqdm import tqdm

os.makedirs("processed/alignment", exist_ok=True)

def _run_single_alignment(protein, template, chain):
    protein_path = f"processed/surface_extraction/{protein}.asa.pdb"
    interface_path = f"templates/interfaces/{template}_{chain}_int.pdb"
    matrix_path = f"processed/alignment/{protein}_{template}_{chain}_matrix.out"
    tm_path = f"processed/alignment/{protein}_{template}_{chain}_out.tm"
    try:
        with open(tm_path, "w") as f:
            subprocess.run(
                ["external_tools/TMalign", protein_path, interface_path, "-m", matrix_path],
                stdout=f,
                stderr=subprocess.DEVNULL,
                check=True,
            )
        result = parse_tmalign(protein_path, interface_path, protein, template, chain)
        return {
            "protein": protein,
            "template": template,
            "chain": chain,
            "match_count": result["match_count"],
            "tm_score": result["tm_score"],
            "translation": json.dumps(result["translation"]),
            "rotation_mat": json.dumps(result["rotation_mat"]),
            "match_dict": json.dumps(result["match_dict"]),
        }
    except Exception as e:
        print(f"TM-align did not run for {protein} and {template}_{chain}: {e}")
        return None
    finally:
        for p in (matrix_path, tm_path):
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass

def align(targets, templates):
    tasks = []
    for target in targets:
        for template in templates:
            for chain in template[4:]:
                interface_path = f"templates/interfaces/{template}_{chain}_int.pdb"
                try:
                    with open(interface_path, "r") as f:
                        if len(f.readlines()) <= 2:
                            print(f"Interface file {interface_path} is empty.")
                            continue
                except FileNotFoundError:
                    continue
                tasks.append((target, template, chain))

    if not tasks:
        print("No alignment tasks to run.")
        return

    workers = min(32, (os.cpu_count() or 4))
    rows = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_run_single_alignment, p, t, c): (p, t, c) for p, t, c in tasks}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Aligning"):
            row = future.result()
            if row is not None:
                rows.append(row)

    if rows:
        fieldnames = ["protein", "template", "chain", "match_count", "tm_score", "translation", "rotation_mat", "match_dict"]
        rows_by_target = {}
        for row in rows:
            target = row["protein"]
            rows_by_target.setdefault(target, []).append(row)
        for target, target_rows in rows_by_target.items():
            output_csv = f"processed/alignment/{target}.csv"
            with open(output_csv, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(target_rows)
            print(f"Wrote {len(target_rows)} alignments to {output_csv}")

def parse_tmalign(protein_path, interface_path, protein, template, chain):
    with open(f"processed/alignment/{protein}_{template}_{chain}_matrix.out", "r") as matrix_file:
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

    with open(f"processed/alignment/{protein}_{template}_{chain}_out.tm", "r") as tm_file:
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
        for i, s in enumerate(match):
            if s == ":" or s == ".":
                seq1_str = seq1_chain_ids[index1] + "." + seq1[i] + "." + seq1_res_ids[index1]
                seq2_str = seq2_chain_ids[index2] + "." + seq2[i] + "." + seq2_res_ids[index2]
                match_dict[seq2_str] = seq1_str
            if seq1[i] != "-":
                index1 += 1
            if seq2[i] != "-":
                index2 += 1

    return {
        "match_count": match_count,
        "translation": translation,
        "rotation_mat": rotation_mat,
        "match_dict": match_dict,
        "tm_score": max(tm_score_1, tm_score_2),
    }

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
