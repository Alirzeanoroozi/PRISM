"""
Rosetta refinement stage using PyRosetta.
Supports multi-chain inputs: ligand and receptor can each have multiple chains,
or all/selected chains from the PDB.
"""
import os
import shutil

from .contact import get_contacts_from_atom_lines

ROSETTA_DB = os.environ.get(
    "ROSETTA_DB",
    "/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database/",
)
ROSETTA_INT_SCORE_THRESHOLD = -5.0

ROSETTA_DIR = "processed/rosetta_refinement"
ENERGY_DIR = os.path.join(ROSETTA_DIR, "energies")
STRUCTURE_DIR = os.path.join(ROSETTA_DIR, "structures")

os.makedirs(ROSETTA_DIR, exist_ok=True)
os.makedirs(ENERGY_DIR, exist_ok=True)
os.makedirs(STRUCTURE_DIR, exist_ok=True)


def refiner(passed_pairs):
    """Run Rosetta refinement on passed receptor-ligand pairs."""
    passed_pairs = [p for p in passed_pairs if p is not None]
    with open(os.path.join(ROSETTA_DIR, "refinement_energies.txt"), "w") as file_out:
        for passed0, passed1 in passed_pairs:
            totalscore, intscore, _ = calculate_energy(passed0, passed1)
            if intscore != "-":
                file_out.write(f"{passed0}\t{passed1}\t{intscore}\t{totalscore}\n")


def calculate_energy(passed0, passed1):
    """
    Run prepack + docking local refine via PyRosetta.
    Supports multi-chain receptor and ligand (e.g. AB_CD).
    """
    try:
        import pyrosetta
        from pyrosetta.rosetta.core.import_pose import pose_from_file
        from pyrosetta.rosetta.protocols.docking import setup_foldtree
        from pyrosetta.rosetta.utility import Vector1
    except ImportError as e:
        print(f"PyRosetta not available: {e}. Install with: pip install pyrosetta")
        return "-", "-", "-"

    try:
        combined_path = combine_pdb(passed0, passed1)
        if not combined_path:
            return "-", "-", "-"

        left_chains = _extract_chain_ids(passed0)
        right_chains = _extract_chain_ids(passed1)
        if not left_chains or not right_chains:
            print(f"Could not determine partner chains for {passed0} and {passed1}")
            return "-", "-", "-"
        # Multi-chain format: "AB_CD" for chains A,B vs C,D
        partner_chains = f"{''.join(left_chains)}_{''.join(right_chains)}"

        out_name = os.path.splitext(os.path.basename(combined_path))[0]

        # Initialize PyRosetta (no-op if already initialized)
        pyrosetta.init(
            extra_options="-ex1 -ex2aro -ignore_zero_occupancy false -detect_disulf false"
        )

        pose = pose_from_file(combined_path)
        jump_num = Vector1(1)
        setup_foldtree(pose, partner_chains, jump_num)

        # Prepack
        try:
            from pyrosetta.rosetta.protocols.docking import DockingPrepackProtocol
            prepack = DockingPrepackProtocol(jump_num[1])
            prepack.apply(pose)
        except (ImportError, AttributeError):
            _prepack_with_pack_rotamers(pose)

        prepacked_path = os.path.join(ROSETTA_DIR, f"{out_name}_0001.pdb")
        pose.dump_pdb(prepacked_path)

        # Docking local refine
        try:
            from pyrosetta.rosetta.protocols.docking import DockMCMProtocol
            from pyrosetta import create_score_function
            scorefxn = create_score_function("ref2015")
            dock = DockMCMProtocol(jump_num[1], scorefxn, scorefxn)
            try:
                dock.set_local_refine(True)
            except AttributeError:
                pass
            dock.apply(pose)
        except Exception as e:
            print(f"DockMCMProtocol failed: {e}; trying DockingProtocol fallback")
            _dock_local_refine_fallback(pose, jump_num[1])

        out_pdb_name = f"{out_name}_0001_0001.pdb"
        rosetta_out_path = os.path.join(STRUCTURE_DIR, out_pdb_name)
        final_out_path = os.path.join(ROSETTA_DIR, out_pdb_name)
        # DockMCM modifies pose in place; dump the refined pose
        pose.dump_pdb(rosetta_out_path)

        # Parse scores (try score.sc or from pose)
        totalscore = "-"
        interaction_score = "-"
        score_path = os.path.join(ENERGY_DIR, "score.sc")
        renamed_score = os.path.join(ENERGY_DIR, f"{out_name}_score.sc")
        if os.path.exists(score_path):
            shutil.move(score_path, renamed_score)
            with open(renamed_score, "r") as sf:
                for i, line in enumerate(sf):
                    if i == 2:
                        parts = line.split()
                        if len(parts) >= 6:
                            totalscore = float(parts[1].strip())
                            interaction_score = float(parts[5].strip())
                        break
        else:
            try:
                from pyrosetta.rosetta.core.scoring import Interface
                scorefxn = pyrosetta.create_score_function("ref2015")
                scorefxn(pose)
                totalscore = pose.energies().total_energy()
                interface = Interface(jump_num[1])
                interface.calculate(pose)
                interaction_score = interface.interface_energy(pose, scorefxn)
            except Exception:
                pass

        if (
            os.path.exists(rosetta_out_path)
            and interaction_score != "-"
            and float(interaction_score) <= ROSETTA_INT_SCORE_THRESHOLD
        ):
            shutil.copy2(rosetta_out_path, final_out_path)
            structure_list_0, structure_list_1 = _extract_atom_lines_by_partners(
                final_out_path, left_chains, right_chains
            )
            int_res_path = f"{final_out_path}.intRes.txt"
            try:
                get_contacts_from_atom_lines(
                    final_out_path, int_res_path, structure_list_0, structure_list_1
                )
            except Exception as e:
                print(f"Exception during get_contacts_from_atom_lines: {e}")
            return str(totalscore), str(interaction_score), final_out_path

        if not os.path.exists(rosetta_out_path):
            print(f"Structure file not found: {rosetta_out_path}")
        elif interaction_score == "-":
            print(f"Interaction score is '-' for {rosetta_out_path}")
        else:
            print(f"Interaction score {interaction_score} exceeds threshold for {rosetta_out_path}")
        return "-", "-", "-"

    except Exception as e:
        print(f"Exception during calculate_energy: {e}")
        import traceback
        traceback.print_exc()
        return "-", "-", "-"


def _prepack_with_pack_rotamers(pose):
    """Fallback prepack using PackRotamersMover."""
    from pyrosetta import create_score_function
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepackingRLT
    scorefxn = create_score_function("ref2015")
    tf = TaskFactory()
    tf.push_back(RestrictToRepackingRLT())
    packer = PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    packer.apply(pose)


def _dock_local_refine_fallback(pose, jump_num):
    """Fallback docking local refine using RigidBodyPerturbMover + MinMover."""
    from pyrosetta import create_score_function
    from pyrosetta.rosetta.protocols.rigid import RigidBodyPerturbMover
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover
    from pyrosetta.rosetta.core.kinematics import MoveMap
    scorefxn = create_score_function("ref2015")
    pert = RigidBodyPerturbMover(jump_num, 3, 1)
    pert.apply(pose)
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(True)
    mm.set_jump(True)
    min_mover = MinMover(mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.01, True)
    min_mover.apply(pose)


def combine_pdb(passed0, passed1):
    """Combine receptor and ligand PDBs into one file (multi-chain supported)."""
    try:
        base0 = os.path.splitext(os.path.basename(passed0))[0]
        base1 = os.path.splitext(os.path.basename(passed1))[0]
        combined_path = os.path.join(ROSETTA_DIR, f"{base0}_{base1}_rosetta.pdb")
        with (
            open(passed0, "r") as p0file,
            open(passed1, "r") as p1file,
            open(combined_path, "w") as combinedfile,
        ):
            for line in p0file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    combinedfile.write(line)
            combinedfile.write("TER\n")
            for line in p1file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    combinedfile.write(line)
            combinedfile.write("END\n")
        return combined_path
    except Exception as e:
        print(f"Exception during combine_pdb: {e}")
        return ""


def _extract_chain_ids(pdb_path):
    """Extract all chain IDs from a PDB (supports multi-chain)."""
    chain_ids = []
    seen = set()
    try:
        with open(pdb_path, "r") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if len(line) >= 22:
                        ch = line[21]
                        if ch not in seen:
                            seen.add(ch)
                            chain_ids.append(ch)
    except Exception as exc:
        print(f"Could not read chain IDs from {pdb_path}: {exc}")
    return chain_ids


def _extract_atom_lines_by_partners(pdb_path, left_chains, right_chains):
    """Extract ATOM lines for partner 0 (left) and partner 1 (right)."""
    left_set = set(left_chains)
    right_set = set(right_chains)
    lines_0 = []
    lines_1 = []
    try:
        with open(pdb_path, "r") as fh:
            for line in fh:
                if line.startswith("TER"):
                    continue
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if len(line) >= 22:
                        ch = line[21]
                        if ch in left_set:
                            lines_0.append(line)
                        elif ch in right_set:
                            lines_1.append(line)
    except Exception as exc:
        print(f"Could not read PDB {pdb_path}: {exc}")
    return lines_0, lines_1


if __name__ == "__main__":
    # Example: single structure with partners A_B
    combined_path = "templates/pdbs/2ai9.pdb"
    partner_chains = "A_B"
    print("Run refinement via: refiner([(receptor.pdb, ligand.pdb), ...])")
