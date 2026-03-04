import json
import os
import pandas as pd
from tqdm import tqdm

from .utils import read_ca_coordinates, distance_calculator

TRANSFORMATION_DIR = "processed/transformation"
os.makedirs(TRANSFORMATION_DIR, exist_ok=True)

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

_alignment_cache = {}

def transformer(receptor_targets, ligand_targets, templates):
    passed_pairs = []
    for template in tqdm(templates):
        for receptor, ligand in zip(receptor_targets, ligand_targets):
            pairs = process_pair_for_template(template, receptor, ligand)
            if pairs is not []:
                passed_pairs.extend(pairs)
    return passed_pairs

def load_alignment(protein, template, chain_id):
    if protein not in _alignment_cache:
        csv_path = f"processed/alignment/{protein}.csv"
        if not os.path.exists(csv_path):
            return None
        _alignment_cache[protein] = pd.read_csv(csv_path)
    df = _alignment_cache[protein]
    row = df[(df["template"] == template) & (df["chain"] == chain_id)]
    if row.empty:
        return None
    row = row.iloc[0]
    return {
        "match_count": int(row["match_count"]),
        "tm_score": float(row["tm_score"]),
        "translation": json.loads(row["translation"]),
        "rotation_mat": json.loads(row["rotation_mat"]),
    }

def process_pair_for_template(template, receptor, ligand):
    passed_pairs = []
    chain1 = template[4]
    chain2 = template[5]

    receptor_align_1 = load_alignment(receptor, template, chain1)
    ligand_align_2 = load_alignment(ligand, template, chain2)
    passed_pairs.append(create_transformed_pair(template, receptor, ligand, receptor_align_1, ligand_align_2, "o1"))

    receptor_align_2 = load_alignment(receptor, template, chain2)
    ligand_align_1 = load_alignment(ligand, template, chain1)
    passed_pairs.append(create_transformed_pair(template, receptor, ligand, receptor_align_2, ligand_align_1, "o2"))
    return passed_pairs

def create_transformed_pair(template, receptor, ligand, receptor_alignment, ligand_alignment, orientation_suffix):
    receptor_pdb = f"processed/pdbs/{receptor}.pdb"
    ligand_pdb = f"processed/pdbs/{ligand}.pdb"

    receptor_output = f"processed/transformation/{template}_{receptor}_{ligand}_{orientation_suffix}_R.pdb"
    ligand_output = f"processed/transformation/{template}_{receptor}_{ligand}_{orientation_suffix}_L.pdb"

    apply_tm_transform(receptor_pdb, receptor_output, receptor_alignment['translation'], receptor_alignment['rotation_mat'])
    apply_tm_transform(ligand_pdb, ligand_output, ligand_alignment['translation'], ligand_alignment['rotation_mat'])

    if pair_has_acceptable_clashes(receptor_output, ligand_output):
        return [receptor_output, ligand_output]
    return None

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
