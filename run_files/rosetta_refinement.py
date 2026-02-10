import os

from .fiberdockInterfaceExtractor import fiberdock_interface_extractor

ROSETTA_PREPACK = "/kuacc/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/bin/docking_prepack_protocol.static.linuxgccrelease"
ROSETTA_DOCK = "/kuacc/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/bin/docking_protocol.static.linuxgccrelease"
ROSETTA_DB = "/kuacc/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database/"
ROSETTA_INT_SCORE_THRESHOLD = -5.0

os.makedirs("processed/flexible_refinement", exist_ok=True)
os.makedirs("processed/flexible_refinement/energies", exist_ok=True)
os.makedirs("processed/flexible_refinement/structures", exist_ok=True)

class FlexibleRefinement:
    def __init__(self, passed_pairs):
        self.sizes = {}
        self.ca_num = {}
        self.last_residue = {}
        self.seen_structures = {}
        self.passed_pairs = passed_pairs

    def refiner(self):
        try:
            energy_structure = []
            energy_file = "flexibleRefinement/energies/refinement_energies"
            with open(energy_file, "w") as file_out:
                for passed0, passed1 in self.passed_pairs:
                    totalscore, intscore, structure = self.calculate_energy(passed0, passed1)

                    if intscore != "-":
                        energy_structure.append([f"{passed0}\t{passed1}\t{intscore}", structure])
                        file_out.write(f"{passed0}\t{passed1}\t{intscore}\n")

            pred_file = f"{self.global_ros}/refinement_energies"
            if energy_structure:
                os.system(f"cp {energy_file} {self.global_ros}/refinement_energies")
                os.system(f"chmod 775 {pred_file}")

            return energy_structure
        except Exception as e:
            print(f"Exception during refiner: {e}")
            return []

    def calculate_energy(self, passed0, passed1):
        """
        Prepare input PDBs, run Rosetta prepack and docking, parse scores,
        and generate interface visualization using FiberDock-style extractor.
        """
        try:
            target1_chain_name_mapping = {}
            target2_chain_name_mapping = {}
            original_chains_t1 = list(
                self.get_chains_from_pdb(os.path.join("transformation", passed0))
            )
            original_chains_t2 = list(
                self.get_chains_from_pdb(os.path.join("transformation", passed1))
            )

            all_characters = (
                [chr(ord("A") + x) for x in range(26)]
                + [chr(ord("a") + x) for x in range(26)]
                + [chr(ord("0") + x) for x in range(10)]
            )
            if len(original_chains_t1) + len(original_chains_t2) > 62:
                print(
                    "The number of chains in target1 and target2 exceeds 62, "
                    "can't assign unique IDs for chains. Rosetta run skipped"
                )
                return "-", "-", "-"

            new_chains_t1 = all_characters[: len(original_chains_t1)]
            new_chains_t2 = all_characters[
                len(original_chains_t1) : len(original_chains_t1) + len(original_chains_t2)
            ]
            chain_renames_mapping = {}
            for i, ch in enumerate(original_chains_t1):
                target1_chain_name_mapping[ch] = new_chains_t1[i]
                chain_renames_mapping[new_chains_t1[i]] = ch
            for i, ch in enumerate(original_chains_t2):
                target2_chain_name_mapping[ch] = new_chains_t2[i]
                chain_renames_mapping[new_chains_t2[i]] = ch

            self.change_chain_name(
                os.path.join("transformation", passed0), target1_chain_name_mapping
            )
            self.change_chain_name(
                os.path.join("transformation", passed1), target2_chain_name_mapping
            )

            combined_path = self.combine_pdb(passed0, passed1)
            if not combined_path:
                return "-", "-", "-"

            partner_chains = "".join(new_chains_t1) + "_" + "".join(new_chains_t2)
            cwd = os.getcwd()

            os.system(
                "%s -database %s -s %s -partners %s -ex1 -ex2aro "
                "-out:file:scorefile %s/%s_prepack_score.sc -overwrite "
                "-ignore_zero_occupancy false -detect_disulf false"
                % (
                    self.rosetta_prepack,
                    self.rosetta_db,
                    combined_path,
                    partner_chains,
                    os.path.join(cwd, "flexibleRefinement/energies"),
                    combined_path,
                )
            )
            prepacked_file = combined_path.split(".pdb")[0] + "_0001.pdb"
            os.system(f"mv {prepacked_file.split('/')[1]} flexibleRefinement/")
            os.system(
                "%s -database %s -s %s -docking_local_refine -partners %s "
                "-ex1 -ex2aro -overwrite -ignore_zero_occupancy false "
                "-detect_disulf false -out:path:score %s"
                % (
                    self.rosetta_dock,
                    self.rosetta_db,
                    prepacked_file,
                    partner_chains,
                    os.path.join(cwd, "flexibleRefinement/energies"),
                )
            )
            out_name = (combined_path.split("/")[1]).split(".pdb")[0]
            out_pdb = out_name + "_0001_0001.pdb"

            if not os.path.exists(self.global_ros):
                os.makedirs(self.global_ros, exist_ok=True)
                os.system(f"chmod 775 {self.global_ros}")

            totalscore = "-"
            interaction_score = "-"

            score_path = "flexibleRefinement/energies/score.sc"
            if os.path.exists(score_path):
                os.system(
                    "mv flexibleRefinement/energies/score.sc "
                    f"flexibleRefinement/energies/{out_name}_score.sc"
                )
                os.system(
                    f"mv {out_pdb} flexibleRefinement/structures/{out_pdb}"
                )

                self.change_chain_name(
                    os.path.join("flexibleRefinement/structures", out_pdb),
                    chain_renames_mapping,
                )
                with open(
                    f"flexibleRefinement/energies/{out_name}_score.sc", "r"
                ) as scorefile:
                    for i, line in enumerate(scorefile):
                        if i == 2:
                            temp = line.split()
                            totalscore = float(temp[1].strip())
                            interaction_score = float(temp[5].strip())
                            break
            else:
                print(f"Could not find a score file for {out_pdb}")

            structure_list = {0: [], 1: []}
            rosetta_out_structure = f"flexibleRefinement/structures/{out_pdb}"
            if (
                os.path.exists(rosetta_out_structure)
                and interaction_score != "-"
                and interaction_score <= self.rosetta_int_score_threshold
            ):
                with open(rosetta_out_structure, "r") as fh:
                    # First partner
                    for chain in original_chains_t1:
                        for l in fh:
                            if l[:3] == "TER":
                                break
                            if l[:4] == "ATOM" and l[21] == chain:
                                structure_list[0].append(l)

                    # Second partner
                    for chain in original_chains_t2:
                        for l in fh:
                            if l[:3] == "TER":
                                break
                            if l[:4] == "ATOM" and l[21] == chain:
                                structure_list[1].append(l)

                key = f"{passed0},{passed1},{interaction_score}"
                if key not in self.seen_structures:
                    self.seen_structures[key] = 1

                try:
                    os.system(
                        f"cp {rosetta_out_structure} "
                        f"{os.path.join(self.global_ros, out_pdb)}"
                    )
                except Exception as e:
                    print(e)

                int_res_path = os.path.join(
                    self.global_ros, f"{out_pdb}.intRes.txt"
                )
                try:
                    fiberdock_interface_extractor(
                        os.path.join(self.global_ros, out_pdb),
                        int_res_path,
                        structure_list[0],
                        structure_list[1],
                    )
                except Exception as e:
                    print(e)
                    print(
                        "Rosetta output visualization and interface residues "
                        "file were not created."
                    )

                return totalscore, str(interaction_score), os.path.join(
                    self.out_folder, out_pdb
                )
            else:
                if not os.path.exists(rosetta_out_structure):
                    print(f"structure file couldn't be found {rosetta_out_structure}!!")
                elif interaction_score == "-":
                    print(
                        f"interaction_score is '-' for structure "
                        f"{rosetta_out_structure}!!"
                    )
                elif interaction_score > self.rosetta_int_score_threshold:
                    print(
                        f"interaction_score {interaction_score} exceeds threshold "
                        f"for structure {rosetta_out_structure}!!"
                    )
            return "-", "-", "-"
        except Exception as e:
            print(f"Exception during calculate_energy: {e}")
            return "-", "-", "-"

    def get_chains_from_pdb(self, pdb_path):
        try:
            chains = set()
            with open(pdb_path, "r") as pdbfile:
                for line in pdbfile:
                    if line[0:4] == "ATOM":
                        chains.add(line[21])
            return chains
        except Exception as e:
            print(f"Exception during get_chains_from_pdb: {e}")
            return set()

    def change_chain_name(self, protein_path, chain_name_mapping):
        try:
            with open(protein_path, "r") as prevfile, open(
                protein_path + "_temp", "w"
            ) as newfile:
                for line in prevfile:
                    if len(line) >= 4 and line[0:4] == "ATOM":
                        if line[21] in chain_name_mapping:
                            line = (
                                line[0:21]
                                + chain_name_mapping[line[21]]
                                + line[22:]
                            )
                    newfile.write(line)
            os.system(f"cp {protein_path + '_temp'} {protein_path}")
        except Exception as e:
            print(f"Exception during change_chain_name: {e}")

    def combine_pdb(self, passed0, passed1):
        """
        Combine two transformed PDBs into a single complex PDB for Rosetta.

        The output filename is based only on the base names of the input
        transformed PDBs and no longer depends on the legacy underscore-based
        naming convention.
        """
        try:
            ppath0 = os.path.join("transformation", passed0)
            ppath1 = os.path.join("transformation", passed1)

            base0 = os.path.splitext(os.path.basename(passed0))[0]
            base1 = os.path.splitext(os.path.basename(passed1))[0]
            combined_name = f"{base0}__{base1}.rosetta.pdb"
            combined_path = os.path.join("flexibleRefinement", combined_name)

            with open(ppath0, "r") as p0file, open(ppath1, "r") as p1file, open(
                combined_path, "w"
            ) as combinedfile:
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
