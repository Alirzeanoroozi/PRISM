import os
import json
import pandas as pd

from .utils import read_ca_coordinates, distance_calculator

TRANSFORMATION_DIR = "processed/transformation"
os.makedirs(TRANSFORMATION_DIR, exist_ok=True)

MERGE_DIR = "processed/output"
os.makedirs(MERGE_DIR, exist_ok=True)

CLASHING_DISTANCE = 3

def transformer(receptor_targets, ligand_targets):
    all_passed_pairs = []
    for receptor, ligand in zip(receptor_targets, ligand_targets):
        passed_pair = process_pair_for_template(receptor, ligand)
        if passed_pair:
            all_passed_pairs.append(passed_pair)
    return all_passed_pairs

def process_pair_for_template(receptor, ligand):
    # ["protein", "template", "chain", "match_count", "tm_score", "len_target", "len_template", "translation", "rotation_mat"]
    receptor_df = pd.read_csv(f"processed/alignment/{receptor}.csv")
    ligand_df = pd.read_csv(f"processed/alignment/{ligand}.csv")

    receptor_templates = list(receptor_df["template"].unique())
    ligand_templates = list(ligand_df["template"].unique())

    for receptor_template in receptor_templates:
        for ligand_template in ligand_templates:
            receptor_alignments = receptor_df[receptor_df["template"] == receptor_template].iloc[0]
            ligand_alignments = ligand_df[ligand_df["template"] == ligand_template].iloc[0]

            if len(receptor_alignments) > 0 and len(ligand_alignments) > 0 and create_transformed_pair(receptor_template, receptor, ligand, receptor_alignments, ligand_alignments):
                return (receptor, ligand)

def create_transformed_pair(template, receptor, ligand, receptor_alignments, ligand_alignments):
    receptor_input = f"processed/pdbs/{receptor[:4].lower()}.pdb"
    ligand_input = f"processed/pdbs/{ligand[:4].lower()}.pdb"

    receptor_output = f"processed/transformation/{template}_{receptor}_{ligand}_R.pdb"
    ligand_output = f"processed/transformation/{template}_{receptor}_{ligand}_L.pdb"

    # translation / rotation_mat are stored as JSON strings in the CSV; decode them if needed.
    def _parse_vec(val, default):
        if isinstance(val, str):
            try:
                return json.loads(val)
            except Exception:
                return default
        return val if val is not None else default

    rec_translation = _parse_vec(receptor_alignments.get("translation"), [0.0, 0.0, 0.0])
    rec_rotation = _parse_vec(
        receptor_alignments.get("rotation_mat"),
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    )
    lig_translation = _parse_vec(ligand_alignments.get("translation"), [0.0, 0.0, 0.0])
    lig_rotation = _parse_vec(
        ligand_alignments.get("rotation_mat"),
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    )
    output_pdb = f"{MERGE_DIR}/{template}_{receptor}_{ligand}.pdb"
    apply_tm_transform(receptor_input, receptor_output, rec_translation, rec_rotation)
    apply_tm_transform(ligand_input, ligand_output, lig_translation, lig_rotation)
    merge_pdb_files(receptor_output, ligand_output, output_pdb)

    return pair_has_acceptable_clashes(receptor_output, ligand_output)

def merge_pdb_files(receptor_path, ligand_path, output_path):
    with open(output_path, "w") as out_f:
        with open(receptor_path, "r") as rec_f:
            for line in rec_f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    out_f.write(line)
                elif line.startswith("TER"):
                    out_f.write(line)
        with open(ligand_path, "r") as lig_f:
            for line in lig_f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    out_f.write(line)
                elif line.startswith("TER"):
                    out_f.write(line)

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

                    new_line = (
                        f"{line[:30]}"
                        f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}"
                        f"{line[54:]}"
                    )
                    out_f.write(new_line)
            out_f.write("ENDMDL\n")
        return True
    except Exception as exc:
        print(f"Error applying TM transform to {input_pdb}: {exc}")
        return False

def pair_has_acceptable_clashes(receptor_path, ligand_path):
    receptor_coords = read_ca_coordinates(receptor_path)
    ligand_coords = read_ca_coordinates(ligand_path)

    for receptor_coord in receptor_coords:
        for ligand_coord in ligand_coords:
            if distance_calculator(receptor_coord, ligand_coord) < CLASHING_DISTANCE:
                return False
    return True

if __name__ == "__main__":
    receptor_targets = ["1FGNH"]
    ligand_targets = ["1TFHA"]
    passed_pairs = transformer(receptor_targets, ligand_targets)
    print(passed_pairs)