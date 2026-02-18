import os
from .contact import get_contacts

ROSETTA_PREPACK = "docking_prepack_protocol.static.linuxgccrelease"
ROSETTA_DOCK = "docking_protocol.static.linuxgccrelease"
ROSETTA_DB = "/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database/"
ROSETTA_INT_SCORE_THRESHOLD = -5.0

os.makedirs("processed/rosetta_refinement", exist_ok=True)
os.makedirs("processed/rosetta_refinement/energies", exist_ok=True)
os.makedirs("processed/rosetta_refinement/structures", exist_ok=True)

def refiner(passed_pairs):
    with open("processed/rosetta_refinement/refinement_energies.txt", "w") as file_out:
        for passed0, passed1 in passed_pairs:
            #     left_output = f"processed/transformation/{template}_{left_query}_{right_query}_{orientation_suffix}_L.pdb"
            #     right_output = f"processed/transformation/{template}_{left_query}_{right_query}_{orientation_suffix}_R.pdb"
            totalscore, intscore, structure = calculate_energy(passed0, passed1)

            if intscore != "-":
                file_out.write(f"{passed0}\t{passed1}\t{intscore}\t{totalscore}\n")

def calculate_energy(passed0, passed1):
    try:
        combined_path = combine_pdb(passed0, passed1)
        partner_chains = f"{passed0[4]}_{passed1[4]}"
        os.system(f"{ROSETTA_PREPACK} -database {ROSETTA_DB} -s {combined_path} -partners {partner_chains} -ex1 -ex2aro \
                  -out:file:scorefile processed/rosetta_refinement/energies/{combined_path}_prepack_score.sc -overwrite -ignore_zero_occupancy false -detect_disulf false")
        prepacked_file = combined_path.split(".pdb")[0] + "_0001.pdb"
        os.system(f"mv {prepacked_file.split('/')[1]} processed/rosetta_refinement/")
        os.system(f"{ROSETTA_DOCK} -database {ROSETTA_DB} -s {prepacked_file} -docking_local_refine -partners {partner_chains} \
            -ex1 -ex2aro -overwrite -ignore_zero_occupancy false -detect_disulf false -out:path:score processed/rosetta_refinement/energies")
        out_name = (combined_path.split("/")[1]).split(".pdb")[0]
        out_pdb = out_name + "_0001_0001.pdb"

        totalscore = "-"
        interaction_score = "-"
        if os.path.exists("processed/rosetta_refinement/energies/score.sc"):
            os.system(f"mv processed/rosetta_refinement/energies/score.sc processed/rosetta_refinement/energies/{out_name}_score.sc")
            os.system(f"mv {out_pdb} processed/rosetta_refinement/structures/{out_pdb}")

            with open(f"processed/rosetta_refinement/energies/{out_name}_score.sc", "r") as scorefile:
                for i, line in enumerate(scorefile):
                    if i == 2:
                        temp = line.split()
                        totalscore = float(temp[1].strip())
                        interaction_score = float(temp[5].strip())
                        break
        else:
            print(f"Could not find a score file for {out_pdb}")

        structure_list = {0: [], 1: []}
        rosetta_out_structure = f"processed/rosetta_refinement/structures/{out_pdb}"
        if os.path.exists(rosetta_out_structure) and interaction_score != "-" and interaction_score <= ROSETTA_INT_SCORE_THRESHOLD:
            with open(rosetta_out_structure, "r") as fh:
                for l in fh:
                    if l[:3] == "TER":
                        break
                    if l[:4] == "ATOM" and l[21] == passed0[4]:
                        structure_list[0].append(l)
                    if l[:4] == "ATOM" and l[21] == passed1[4]:
                        structure_list[1].append(l)
            try:
                os.system(f"cp {rosetta_out_structure} processed/rosetta_refinement/{out_pdb}")
            except Exception as e:
                print(f"Exception during cp: {e}")

            int_res_path = f"processed/rosetta_refinement/{out_pdb}.intRes.txt"
            try:
                get_contacts(f"processed/rosetta_refinement/{out_pdb}", int_res_path, structure_list[0], structure_list[1])
            except Exception as e:
                print(f"Exception during get_contacts: {e}")

            return totalscore, str(interaction_score), f"processed/rosetta_refinement/{out_pdb}"
        else:
            if not os.path.exists(rosetta_out_structure):
                print(f"structure file couldn't be found {rosetta_out_structure}!!")
            elif interaction_score == "-":
                print(f"interaction_score is '-' for structure {rosetta_out_structure}!!")
            elif interaction_score > ROSETTA_INT_SCORE_THRESHOLD:
                print(f"interaction_score {interaction_score} exceeds threshold for structure {rosetta_out_structure}!!")
        return "-", "-", "-"
    except Exception as e:
        print(f"Exception during calculate_energy: {e}")
        return "-", "-", "-"

def combine_pdb(passed0, passed1):
    try:
        combined_path = f'processed/rosetta_refinement/{passed0.split("/")[-1].split(".")[0]}_{passed1.split("/")[-1].split(".")[0]}_rosetta.pdb'

        with open(passed0, "r") as p0file, open(passed1, "r") as p1file, open(combined_path, "w") as combinedfile:
            for line in p0file:
                if line[0:4] == "ATOM":
                    combinedfile.write(line)
            combinedfile.write("TER\n")
            for line in p1file:
                if line[0:4] == "ATOM":
                    combinedfile.write(line)
            combinedfile.write("END\n")
        return combined_path

    except Exception as e:
        print(f"Exception during combine_pdb: {e}")
        return ""

if __name__ == "__main__":
    combined_path = "templates/pdbs/2ai9.pdb"
    partner_chains = "A_B"
    os.system(f"{ROSETTA_PREPACK} \
        -database {ROSETTA_DB} \
            -s {combined_path} \
                -partners {partner_chains} \
                    -ex1 -ex2aro \
                -out:file:scorefile processed/rosetta_refinement/energies/{combined_path.split('/')[-1].split('.')[0]}_prepack_score.sc \
                    -overwrite -ignore_zero_occupancy false -detect_disulf false")
    
    # prepacked_file = "2ai9_0001.pdb"
    # os.system(f"{ROSETTA_DOCK} -database {ROSETTA_DB} -s {prepacked_file} -docking_local_refine -partners {partner_chains} \
    # -ex1 -ex2aro -overwrite -ignore_zero_occupancy false -detect_disulf false -out:path:score processed/rosetta_refinement/energies")