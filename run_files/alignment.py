import os
from Bio.PDB import PDBParser
import json

os.makedirs("alignment", exist_ok=True)

def align(queries, templates):
    for protein in queries:
        for template in templates:
            for chain in template[4:]:
                protein_path = f"surface_extraction/{protein}.asa.pdb"
                interface_path = f"template/interface/{template}_{chain}.int"
                try:
                    os.system(f"external_tools/TMalign {protein_path} {interface_path} -m alignment/matrix.out > alignment/out.tm")
                except:
                    print(f"TM-align did not run for {protein} and {template}_{chain}.")
                parse_tmalign(protein_path, interface_path, protein, template, chain)

    if os.path.exists("alignment/matrix.out"):
        os.remove("alignment/matrix.out")
    
    if os.path.exists("alignment/out.tm"):
        os.remove("alignment/out.tm")

def parse_tmalign(protein_path, interface_path, protein, template, chain):
    with open("alignment/matrix.out", "r") as matrix_file:
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

    with open("alignment/out.tm", "r") as tm_file:
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

    multi_dict = {
        "match_count": match_count,
        "translation": translation,
        "rotation_mat": rotation_mat,
        "match_dict": match_dict,
        "tm_score": max(tm_score_1, tm_score_2)
    }
    with open(f"alignment/{protein}_{template}_{chain}.json", "w") as f:
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

if __name__ == "__main__":
    protein = "1a28"
    template = "1a28AB"
    chains = ["A", "B"]
    align([protein], [template])
    for chain in chains:
        data = json.load(open(f"alignment/{protein}_{template}_{chain}.json", "r"))
        print(f"--- {protein}_{template}_{chain} ---")
        print(data)
        print()