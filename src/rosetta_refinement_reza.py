import os
import shutil

ROSETTA_PREPACK = os.environ.get(
    "PRISM_ROSETTA_PREPACK",
    "/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/bin/docking_prepack_protocol.static.linuxgccrelease",
)
ROSETTA_DOCK = os.environ.get(
    "PRISM_ROSETTA_DOCK",
    "/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/bin/docking_protocol.static.linuxgccrelease",
)
ROSETTA_DB = os.environ.get(
    "PRISM_ROSETTA_DB",
    "/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database/",
)
ROSETTA_INT_SCORE_THRESHOLD = -5.0

ROSETTA_DIR = "processed/rosetta_refinement"
ENERGY_DIR = f"{ROSETTA_DIR}/energies"
STRUCTURE_DIR = f"{ROSETTA_DIR}/structures"

os.makedirs(ROSETTA_DIR, exist_ok=True)
os.makedirs(ENERGY_DIR, exist_ok=True)
os.makedirs(STRUCTURE_DIR, exist_ok=True)


def refiner(passed_pairs):
    with open(f"{ROSETTA_DIR}/refinement_energies.txt", "w") as file_out:
        for passed0, passed1 in passed_pairs:
            totalscore, intscore, _ = calculate_energy(passed0, passed1)
            if intscore != "-":
                file_out.write(f"{passed0}\t{passed1}\t{intscore}\t{totalscore}\n")


def calculate_energy(passed0, passed1):
    try:
        combined_path = combine_pdb(passed0, passed1)
        if not combined_path:
            return "-", "-", "-"

        left_chains = _extract_chain_ids(passed0)
        right_chains = _extract_chain_ids(passed1)
        if not left_chains or not right_chains:
            print(f"Could not determine partner chains for {passed0} and {passed1}")
            return "-", "-", "-"
        partner_chains = f"{''.join(left_chains)}_{''.join(right_chains)}"

        out_name = os.path.splitext(os.path.basename(combined_path))[0]
        prepack_score = f"{ENERGY_DIR}/{out_name}_prepack_score.sc"
        prepack_cmd = (
            f"{ROSETTA_PREPACK} -database {ROSETTA_DB} -s {combined_path} -partners {partner_chains} "
            f"-ex1 -ex2aro -out:file:scorefile {prepack_score} -overwrite "
            f"-ignore_zero_occupancy false -detect_disulf false"
        )
        _run_cmd(prepack_cmd)

        prepacked_name = f"{out_name}_0001.pdb"
        prepacked_src = prepacked_name
        prepacked_path = f"{ROSETTA_DIR}/{prepacked_name}"
        if os.path.exists(prepacked_src):
            shutil.move(prepacked_src, prepacked_path)
        elif not os.path.exists(prepacked_path):
            print(f"Prepacked structure not found: {prepacked_name}")
            return "-", "-", "-"

        dock_cmd = (
            f"{ROSETTA_DOCK} -database {ROSETTA_DB} -s {prepacked_path} -docking_local_refine "
            f"-partners {partner_chains} -ex1 -ex2aro -overwrite "
            f"-ignore_zero_occupancy false -detect_disulf false -out:path:score {ENERGY_DIR}"
        )
        _run_cmd(dock_cmd)

        out_pdb = f"{os.path.splitext(prepacked_name)[0]}_0001.pdb"
        out_pdb_src = out_pdb
        rosetta_out_structure = f"{STRUCTURE_DIR}/{out_pdb}"
        if os.path.exists(out_pdb_src):
            shutil.move(out_pdb_src, rosetta_out_structure)

        totalscore = "-"
        interaction_score = "-"
        score_path = f"{ENERGY_DIR}/score.sc"
        renamed_score = f"{ENERGY_DIR}/{out_name}_score.sc"
        if os.path.exists(score_path):
            shutil.move(score_path, renamed_score)
            with open(renamed_score, "r") as scorefile:
                for i, line in enumerate(scorefile):
                    if i == 2:
                        temp = line.split()
                        totalscore = float(temp[1].strip())
                        interaction_score = float(temp[5].strip())
                        break
        else:
            print(f"Could not find score file for {out_pdb}")

        if os.path.exists(rosetta_out_structure) and interaction_score != "-" and interaction_score <= ROSETTA_INT_SCORE_THRESHOLD:
            final_copy = f"{ROSETTA_DIR}/{out_pdb}"
            shutil.copy2(rosetta_out_structure, final_copy)
            return totalscore, str(interaction_score), final_copy

        if not os.path.exists(rosetta_out_structure):
            print(f"structure file couldn't be found {rosetta_out_structure}!!")
        elif interaction_score == "-":
            print(f"interaction_score is '-' for structure {rosetta_out_structure}!!")
        else:
            print(f"interaction_score {interaction_score} exceeds threshold for structure {rosetta_out_structure}!!")
        return "-", "-", "-"
    except Exception as e:
        print(f"Exception during calculate_energy: {e}")
        return "-", "-", "-"


def combine_pdb(passed0, passed1):
    try:
        combined_path = (
            f"{ROSETTA_DIR}/"
            f"{os.path.splitext(os.path.basename(passed0))[0]}_"
            f"{os.path.splitext(os.path.basename(passed1))[0]}_rosetta.pdb"
        )

        with open(passed0, "r") as p0file, open(passed1, "r") as p1file, open(combined_path, "w") as combinedfile:
            for line in p0file:
                if line.startswith("ATOM"):
                    combinedfile.write(line)
            combinedfile.write("TER\n")
            for line in p1file:
                if line.startswith("ATOM"):
                    combinedfile.write(line)
            combinedfile.write("END\n")
        return combined_path
    except Exception as e:
        print(f"Exception during combine_pdb: {e}")
        return ""


def _extract_chain_ids(pdb_path):
    chain_ids = []
    seen = set()
    try:
        with open(pdb_path, "r") as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    chain = line[21]
                    if chain not in seen:
                        seen.add(chain)
                        chain_ids.append(chain)
    except Exception as exc:
        print(f"Could not read chain IDs from {pdb_path}: {exc}")
    return chain_ids


def _run_cmd(cmd):
    rc = os.system(cmd)
    if rc != 0:
        print(f"Command failed ({rc}): {cmd}")
    return rc
