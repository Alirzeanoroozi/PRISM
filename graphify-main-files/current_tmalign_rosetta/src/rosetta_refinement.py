import os

from .contact import get_contacts_from_atom_lines

ROSETTA_PREPACK = os.environ.get("PRISM_ROSETTA_PREPACK", "docking_prepack_protocol.static.linuxgccrelease")
ROSETTA_DOCK = os.environ.get("PRISM_ROSETTA_DOCK", "docking_protocol.static.linuxgccrelease")
ROSETTA_DB = os.environ.get("PRISM_ROSETTA_DB", "/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database/")
ROSETTA_INT_SCORE_THRESHOLD = -5.0
ROSETTA_DIR = "processed/rosetta_refinement"
ENERGY_DIR = os.path.join(ROSETTA_DIR, "energies")
STRUCTURE_DIR = os.path.join(ROSETTA_DIR, "structures")
os.makedirs(ROSETTA_DIR, exist_ok=True)
os.makedirs(ENERGY_DIR, exist_ok=True)
os.makedirs(STRUCTURE_DIR, exist_ok=True)

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
        # Local change: guard against combine_pdb failures before launching Rosetta.
        if not combined_path:
            return "-", "-", "-"
        left_chains, right_chains = partner_chain_ids(passed0, passed1)
        partner_chains = f"{left_chains}_{right_chains}"
        os.system(f"{ROSETTA_PREPACK} -database {ROSETTA_DB} -s {combined_path} -partners {partner_chains} -ex1 -ex2aro \
                  -out:file:scorefile processed/rosetta_refinement/energies/{combined_path}_prepack_score.sc -overwrite -ignore_zero_occupancy false -detect_disulf false")
        prepacked_file = combined_path.split(".pdb")[0] + "_0001.pdb"
        os.system(f"mv {os.path.basename(prepacked_file)} processed/rosetta_refinement/")
        os.system(f"{ROSETTA_DOCK} -database {ROSETTA_DB} -s {prepacked_file} -docking_local_refine -partners {partner_chains} \
            -ex1 -ex2aro -overwrite -ignore_zero_occupancy false -detect_disulf false -out:path:score processed/rosetta_refinement/energies")
        out_name = os.path.basename(combined_path).split(".pdb")[0]
        out_pdb = out_name + "_0001_0001.pdb"

        totalscore = "-"
        interaction_score = "-"
        if os.path.exists("processed/rosetta_refinement/energies/score.sc"):
            os.system(f"mv processed/rosetta_refinement/energies/score.sc processed/rosetta_refinement/energies/{out_name}_score.sc")
            os.system(f"mv {out_pdb} processed/rosetta_refinement/structures/{out_pdb}")

            with open(f"processed/rosetta_refinement/energies/{out_name}_score.sc", "r") as scorefile:
                lines = [l for l in scorefile if l.strip() and not l.startswith("#")]
                # Score data is on the first non-header data line (Rosetta score.sc:
                # header comment, column names, then values).
                if lines:
                    data_line = lines[0] if lines[0].strip().startswith("SCORE:") else lines[-1]
                    temp = data_line.split()
                    totalscore = float(temp[1].strip())
                    interaction_score = float(temp[5].strip())
        else:
            print(f"Could not find a score file for {out_pdb}")

        structure_list = {0: [], 1: []}
        rosetta_out_structure = f"processed/rosetta_refinement/structures/{out_pdb}"
        if os.path.exists(rosetta_out_structure) and interaction_score != "-" and interaction_score <= ROSETTA_INT_SCORE_THRESHOLD:
            with open(rosetta_out_structure, "r") as fh:
                for l in fh:
                    if l[:3] == "TER":
                        break
                    if l.startswith("ATOM") and l[21].strip() in left_chains:
                        structure_list[0].append(l)
                    if l.startswith("ATOM") and l[21].strip() in right_chains:
                        structure_list[1].append(l)
            try:
                os.system(f"cp {rosetta_out_structure} processed/rosetta_refinement/{out_pdb}")
            except Exception as e:
                print(f"Exception during cp: {e}")

            int_res_path = f"processed/rosetta_refinement/{out_pdb}.intRes.txt"
            try:
                get_contacts_from_atom_lines(f"processed/rosetta_refinement/{out_pdb}", int_res_path, structure_list[0], structure_list[1])
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
        combined_path = os.path.join(ROSETTA_DIR, "{}_{}_rosetta.pdb".format(os.path.splitext(os.path.basename(passed0))[0], os.path.splitext(os.path.basename(passed1))[0]))
        left_chains, right_chains = partner_chain_ids(passed0, passed1)
        source_left_chains = source_chain_ids(passed0)
        source_right_chains = source_chain_ids(passed1)
        with open(passed0, "r") as p0file, open(passed1, "r") as p1file, open(combined_path, "w") as combinedfile:
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
    except Exception as e:
        print(f"Exception during combine_pdb: {e}")
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
    try:
        left = source_chain_ids(left_path)
        right = source_chain_ids(right_path)
    except ValueError:
        return "A", "B"
    if not set(left).intersection(right):
        return left, right
    available = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    used = set()
    new_left = ""
    for _ in left:
        candidate = next(candidate for candidate in available if candidate not in used)
        used.add(candidate)
        new_left += candidate
    new_right = ""
    for _ in right:
        candidate = next(candidate for candidate in available if candidate not in used)
        used.add(candidate)
        new_right += candidate
    return new_left, new_right


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
