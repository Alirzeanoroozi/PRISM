import json
import os
import pandas as pd
from .utils import read_ca_coordinates, distance_calculator

os.makedirs("processed/transformation", exist_ok=True)

MINIMUM_RESIDUE_MATCH_COUNT = 15
MINIMUM_RESIDUE_MATCH_PERCENTAGE = 50.0
MINIMUM_HOTSPOT_MATCH_NUMBER = 1
DIFF_PERCENTAGE = 20.0
CONTACT_COUNT = 5
CLASHING_DISTANCE = 3
MAX_CLASHING_COUNT = 5
TM_SCORE_THRESHOLD = 0.5
HOTSPOT_CRITERION = 2
HOTSPOT_COUNT = 1
TEMPLATE_RESIDUE_COUNT = 50
CONTACT_COUNT_THRESHOLD = 5

passed_pairs = []
template_size = {}

def transformer(templates):
    df = pd.read_csv("inputs.csv")

    for template in templates:
        chain1 = template[4]
        chain2 = template[5]

        with open(os.path.join("templates", "interfaces_lists", f"{template}.json"), "r") as f:
            data = json.load(f)

        template_size[f"{template}_{chain1}"] = len(data[chain1])
        template_size[f"{template}_{chain2}"] = len(data[chain2])

        for left_query, right_query in zip(df["Receptor"], df["Ligand"]):
            process_pair_for_template(template, chain1, chain2, left_query, right_query)

    return passed_pairs

def load_alignment(query_id, template, chain_id):
    path = os.path.join("processed/alignment", f"{query_id}_{template}_{chain_id}.json")
    with open(path, "r") as f:
        return json.load(f)

def hotspot_analysis(match_dict):
    return True

def alignment_passes_thresholds(template_key, alignment):
    match_count = alignment.get("match_count", 0)
    tm_score = alignment.get("tm_score", 0.0)
    match_dict = alignment.get("match_dict", {})

    protein_size = float(template_size.get(template_key, 0))
    if protein_size <= 0:
        # Without a size estimate we cannot compute match percentage; fall
        # back to simple count + TM-score checks.
        if match_count < MINIMUM_RESIDUE_MATCH_COUNT or tm_score < TM_SCORE_THRESHOLD:
            return False
        return True

    match_score = (match_count / protein_size) * 100.0

    if not hotspot_analysis(match_dict):
        return False

    if match_count < MINIMUM_RESIDUE_MATCH_COUNT or tm_score < TM_SCORE_THRESHOLD:
        return False

    if protein_size > TEMPLATE_RESIDUE_COUNT:
        return match_score > (MINIMUM_RESIDUE_MATCH_PERCENTAGE - DIFF_PERCENTAGE)
    else:
        return match_score > MINIMUM_RESIDUE_MATCH_PERCENTAGE

def process_pair_for_template(template, chain1, chain2, left_query, right_query):
    left_key_chain1 = f"{template}_{chain1}"
    left_key_chain2 = f"{template}_{chain2}"

    left_align_1 = load_alignment(left_query, template, chain1)
    right_align_1 = load_alignment(right_query, template, chain2)

    if alignment_passes_thresholds(left_key_chain1, left_align_1) and alignment_passes_thresholds(left_key_chain2, right_align_1):
        create_transformed_pair(template, left_query, right_query, left_align_1, right_align_1, passed_pairs, "o1")

    left_align_2 = load_alignment(left_query, template, chain2)
    right_align_2 = load_alignment(right_query, template, chain1)

    if alignment_passes_thresholds(left_key_chain2, left_align_2) and alignment_passes_thresholds(left_key_chain1, right_align_2):
        create_transformed_pair(template, left_query, right_query, left_align_2, right_align_2, passed_pairs, "o2")

def create_transformed_pair(template, left_query, right_query, left_alignment, right_alignment, passed_pairs, orientation_suffix):
    left_input = f"processed/pdbs/{left_query}.pdb"
    right_input = f"processed/pdbs/{right_query}.pdb"

    left_output = f"processed/transformation/{template}_{left_query}_{right_query}_{orientation_suffix}_L.pdb"
    right_output = f"processed/transformation/{template}_{left_query}_{right_query}_{orientation_suffix}_R.pdb"

    apply_tm_transform(left_input, left_output, left_alignment.get("translation", [0.0, 0.0, 0.0]), left_alignment.get("rotation_mat", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
    apply_tm_transform(right_input, right_output, right_alignment.get("translation", [0.0, 0.0, 0.0]), right_alignment.get("rotation_mat", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))

    if pair_has_acceptable_clashes(left_output, right_output):
        passed_pairs.append((left_output, right_output))

def apply_tm_transform(input_pdb, output_pdb, translation, rotation_mat):
    try:
        with open(input_pdb, "r") as in_f, open(output_pdb, "w") as out_f:
            for line in in_f:
                if line.startswith("ATOM"):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                    except ValueError:
                        out_f.write(line)
                        continue

                    new_x = (
                        x * rotation_mat[0][0]
                        + y * rotation_mat[0][1]
                        + z * rotation_mat[0][2]
                        + translation[0]
                    )
                    new_y = (
                        x * rotation_mat[1][0]
                        + y * rotation_mat[1][1]
                        + z * rotation_mat[1][2]
                        + translation[1]
                    )
                    new_z = (
                        x * rotation_mat[2][0]
                        + y * rotation_mat[2][1]
                        + z * rotation_mat[2][2]
                        + translation[2]
                    )

                    line = (
                        f"{line[:30]}"
                        f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}"
                        f"{line[54:]}"
                    )
                out_f.write(line)
    except Exception as exc:
        print(f"Error applying TM transform to {input_pdb}: {exc}")

def pair_has_acceptable_clashes(left_path, right_path):
    left_coords = read_ca_coordinates(left_path)
    right_coords = read_ca_coordinates(right_path)

    clash_count = 0
    for c1 in left_coords:
        for c2 in right_coords:
            if distance_calculator(c1, c2) < CLASHING_DISTANCE:
                clash_count += 1
                if clash_count >= MAX_CLASHING_COUNT:
                    return False
    return True
