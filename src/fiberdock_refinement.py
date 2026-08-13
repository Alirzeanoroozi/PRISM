#!/usr/bin/env python3
"""
FiberDock refinement module for the PRISM pipeline.

Integrates the FiberDock flexible refinement as a PRISM refiner option.
FiberDock adds hydrogens via Reduce, computes Normal Mode Analysis (NMA),
and performs energy-based docking refinement.

The FiberDock binary and helper tools must be staged under:
    external_tools/fiberdock/

Requires: FiberDock, NMA, Reduce, buildFiberDockParams.pl, addHydrogens.pl
"""

import os
import subprocess
import sys
from pathlib import Path
import shutil

REFINER_NAME = "fiberdock"
FIBERDOCK_DIR = os.environ.get(
    "PRISM_FIBERDOCK_DIR",
    os.path.abspath("external_tools/fiberdock"),
)
STRUCTURES_DIR = "processed/fiberdock_refinement/structures"
ENERGIES_DIR = "processed/fiberdock_refinement/energies"

FIBERDOCK_OUTPUT_PREFIX = "fiberdock_energies"


def _ensure_output_dirs():
    os.makedirs(STRUCTURES_DIR, exist_ok=True)
    os.makedirs(ENERGIES_DIR, exist_ok=True)


def _check_tools():
    """Verify that all required FiberDock tools are available."""
    required = {
        "FiberDock": os.path.join(FIBERDOCK_DIR, "FiberDock"),
        "buildFiberDockParams.pl": os.path.join(FIBERDOCK_DIR, "buildFiberDockParams.pl"),
        "addHydrogens.pl": os.path.join(FIBERDOCK_DIR, "addHydrogens.pl"),
        "nma": os.path.join(FIBERDOCK_DIR, "nma"),
    }
    missing = [name for name, path in required.items() if not os.path.exists(path)]
    if missing:
        print(f"FiberDock: missing tools: {missing}")
        return False
    return True


def _add_hydrogens(pdb_path, output_dir):
    """Add hydrogens to a PDB using Reduce."""
    base = os.path.basename(pdb_path).replace(".pdb", "")
    hb_path = os.path.join(output_dir, f"{base}.HB")

    # Try reduce.3 first (64-bit where available), fall back to reduce.2
    reduce_exe = os.path.join(FIBERDOCK_DIR, "reduce")
    if not os.path.exists(reduce_exe):
        for fn in ["reduce.3.23.130521", "reduce.2.21.030604", "addHydrogens.pl"]:
            cand = os.path.join(FIBERDOCK_DIR, fn)
            if os.path.exists(cand):
                reduce_exe = cand
                break

    reduce_dict = os.path.join(FIBERDOCK_DIR, "reduce_het_dict.txt")

    if not os.path.exists(hb_path):
        try:
            if reduce_exe.endswith(".pl"):
                subprocess.run(
                    ["perl", reduce_exe, pdb_path],
                    cwd=output_dir, timeout=120,
                )
            else:
                # Reduce outputs hydrogenated PDB to stdout
                # Save directly as .HB (matching legacy format expected by
                # buildFiberDockParams.pl)
                with open(hb_path, "w") as out:
                    subprocess.run(
                        [reduce_exe, "-OH", "-HIS", "-NOADJust", "-NOROTMET",
                         "-Quiet", "-DB", reduce_dict, pdb_path],
                        stdout=out, stderr=subprocess.DEVNULL,
                        timeout=120,
                    )
        except (OSError, subprocess.TimeoutExpired) as e:
            print(f"FiberDock: addHydrogens failed for {base}: {e}")
    return hb_path if os.path.exists(hb_path) and os.path.getsize(hb_path) > 0 else None


def _create_ca_pdb(pdb_path, output_dir):
    """Create a CA-only PDB for NMA input."""
    base = os.path.basename(pdb_path).replace(".pdb", "")
    ca_path = os.path.join(output_dir, f"{base}.ca.pdb")
    if not os.path.exists(ca_path):
        chain_sizes = {}
        with open(pdb_path) as f_in, open(ca_path, "w") as f_out:
            serial = 0
            for line in f_in:
                if line.startswith("ATOM") and line[13:15].strip() == "CA":
                    serial += 1
                    new_line = (
                        f"ATOM  {serial:5d}  CA  {line[17:20]} {line[21]}"
                        f"{int(line[22:26]):4d}    "
                        f"{float(line[30:38]):8.3f}{float(line[38:46]):8.3f}"
                        f"{float(line[46:54]):8.3f}  1.00  0.00           C  \n"
                    )
                    f_out.write(new_line)
                    chain = line[21]
                    chain_sizes[chain] = chain_sizes.get(chain, 0) + 1
            f_out.write("END\n")
    return ca_path, chain_sizes


def _run_nma(ca_path, output_dir, normal_modes=50):
    """Run Normal Mode Analysis."""
    base = os.path.basename(ca_path).replace(".ca.pdb", "")
    nma_path = os.path.join(output_dir, f"{base}.ca.nma")
    if not os.path.exists(nma_path):
        nma_exe = os.path.join(FIBERDOCK_DIR, "nma")
        try:
            subprocess.run(
                [nma_exe, ca_path, nma_path, str(normal_modes), "3", "10"],
                stdout=open(os.path.join(output_dir, f"{base}.NMA"), "w"),
                stderr=subprocess.PIPE,
                timeout=300,
            )
        except (OSError, subprocess.TimeoutExpired) as e:
            print(f"FiberDock: NMA failed for {base}: {e}")
            return None
    return nma_path if os.path.exists(nma_path) else None


def _build_fiberdock_params(receptor_hb, ligand_hb, output_dir, receptor, ligand):
    """Build FiberDock parameter file.

    receptor_hb/ligand_hb: Paths to hydrogenated PDBs (legacy .HB format).
    Must run from FIBERDOCK_DIR so FindBin resolves lib/ correctly.
    """
    params_script = os.path.join(FIBERDOCK_DIR, "buildFiberDockParams.pl")
    zero_trans = os.path.join(output_dir, "zero-transformation")
    if not os.path.exists(zero_trans):
        with open(zero_trans, "w") as f:
            f.write("1 0 0 0 0 0 0\n")

    params_file = os.path.join(output_dir, "fd_params.txt")
    pair_key = f"{receptor}_{ligand}"
    subprocess.run(
        [
            "perl", params_script,
            os.path.abspath(receptor_hb),
            os.path.abspath(ligand_hb),
            "U", "U", "Default",
            os.path.abspath(zero_trans),
            os.path.abspath(os.path.join(output_dir, FIBERDOCK_OUTPUT_PREFIX)),
            "0", "50", "0.80", "1", "glpk",
            os.path.abspath(params_file),
            "0.05",
            os.path.abspath(os.path.join(output_dir, f"{receptor}.ca.pdb")),
            os.path.abspath(os.path.join(output_dir, f"{receptor}.ca.nma")),
            os.path.abspath(os.path.join(output_dir, f"{ligand}.ca.pdb")),
            os.path.abspath(os.path.join(output_dir, f"{ligand}.ca.nma")),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        timeout=120,
        cwd=FIBERDOCK_DIR,  # Crucial: FindBin needs correct CWD
    )
    return params_file


def _run_fiberdock(params_file, fiberdock_dir, output_dir, receptor, ligand, pair_name=None):
    """Run FiberDock energy calculation.

    FiberDock must run from its own directory to find lib/ files.
    It creates output files (paramName.ref, paramName.pdb) in CWD.
    """
    _ensure_output_dirs()
    fib_out = os.path.join(output_dir, f"{receptor}_{ligand}.fib")
    expected_outputs = (
        os.path.join(output_dir, f"{FIBERDOCK_OUTPUT_PREFIX}.ref"),
        os.path.join(output_dir, f"{FIBERDOCK_OUTPUT_PREFIX}_1.ref.pdb"),
    )
    for output_path in expected_outputs:
        if os.path.isfile(output_path):
            os.unlink(output_path)
    # Run from fiberdock dir so FiberDock finds its lib/
    result = subprocess.run(
        [os.path.join(FIBERDOCK_DIR, "FiberDock"),
         os.path.abspath(params_file), "0"],
        capture_output=True, text=True, timeout=600,
        cwd=fiberdock_dir,
    )
    # Save stdout
    with open(fib_out, "w") as f:
        f.write(result.stdout)
    if result.returncode != 0:
        detail = result.stderr.strip() if result.stderr else "no stderr"
        raise RuntimeError(
            f"FiberDock failed with exit code {result.returncode}: {detail}"
        )

    # FiberDock names outputs after energiesOutFileName, not fd_params.txt.
    import glob as _glob
    for f in _glob.glob(os.path.join(fiberdock_dir, f"{FIBERDOCK_OUTPUT_PREFIX}*")):
        dst = os.path.join(output_dir, os.path.basename(f))
        if os.path.abspath(f) != os.path.abspath(dst):
            import shutil as _shutil
            _shutil.move(f, dst)

    # Parse energy
    energy = "-"
    ref_file = os.path.join(output_dir, f"{FIBERDOCK_OUTPUT_PREFIX}.ref")
    if os.path.exists(ref_file):
        with open(ref_file) as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) >= 3 and parts[0].strip().isdigit():
                    # Format: Sol # | glob | aVdW | rVdW | ...
                    energy = parts[1]  # glob = total energy
                    break

        aggregate_name = f"{pair_name}.ref" if pair_name else os.path.basename(ref_file)
        aggregate_path = os.path.join(ENERGIES_DIR, aggregate_name)
        if os.path.abspath(ref_file) != os.path.abspath(aggregate_path):
            shutil.copy2(ref_file, aggregate_path)

    refined = os.path.join(output_dir, f"{FIBERDOCK_OUTPUT_PREFIX}_1.ref.pdb")
    if os.path.isfile(refined):
        aggregate_name = (
            f"{pair_name}_fiberdock.ref.pdb"
            if pair_name else os.path.basename(refined)
        )
        aggregate_path = os.path.join(STRUCTURES_DIR, aggregate_name)
        if os.path.abspath(refined) != os.path.abspath(aggregate_path):
            shutil.copy2(refined, aggregate_path)

    return energy


def refine_pairs(passed_pairs):
    """FiberDock refinement entry point.

    Args:
        passed_pairs: List of (ligand_pdb, receptor_pdb) tuples from transformer.
    """
    _ensure_output_dirs()
    if not _check_tools():
        print("FiberDock: tools not available, skipping refinement")
        return []

    # Build a map: passed_pairs contains paths to R and L PDBs
    # Need to group them by (template_query_orientation)
    from collections import defaultdict
    pair_groups = defaultdict(lambda: {"R": None, "L": None})

    for pair in passed_pairs:
        for path in pair:
            if path.endswith("_R.pdb"):
                base = path.replace("_R.pdb", "")
                pair_groups[base]["R"] = path
            elif path.endswith("_L.pdb"):
                base = path.replace("_L.pdb", "")
                pair_groups[base]["L"] = path

    results = []
    for base, files in pair_groups.items():
        if not files["R"] or not files["L"]:
            continue

        r_path = files["R"]
        l_path = files["L"]
        pair_name = os.path.basename(base)

        print(f"FiberDock refining {pair_name}...")

        # Create working dir
        work_dir = os.path.join("processed/fiberdock_refinement", pair_name)
        os.makedirs(work_dir, exist_ok=True)

        # Step 1: Add hydrogens
        r_hb = _add_hydrogens(r_path, work_dir)
        l_hb = _add_hydrogens(l_path, work_dir)
        if not r_hb or not l_hb:
            print(f"  FiberDock: hydrogen addition failed for {pair_name}")
            continue

        # Step 2: Create CA PDBs
        r_ca, r_sizes = _create_ca_pdb(r_path, work_dir)
        l_ca, l_sizes = _create_ca_pdb(l_path, work_dir)

        # Step 3: Run NMA
        r_nma = _run_nma(r_ca, work_dir)
        l_nma = _run_nma(l_ca, work_dir)
        if not r_nma or not l_nma:
            print(f"  FiberDock: NMA failed for {pair_name}")
            continue

        # Step 4: Determine receptor/ligand (larger = receptor)
        r_total = sum(r_sizes.values())
        l_total = sum(l_sizes.values())

        if r_total >= l_total:
            receptor_base = os.path.basename(r_path).replace(".pdb", "")
            ligand_base = os.path.basename(l_path).replace(".pdb", "")
        else:
            receptor_base = os.path.basename(l_path).replace(".pdb", "")
            ligand_base = os.path.basename(r_path).replace(".pdb", "")
            r_hb, l_hb = l_hb, r_hb

        # Step 5: Build FiberDock params
        params_file = _build_fiberdock_params(
            r_hb, l_hb, work_dir, receptor_base, ligand_base
        )

        if not os.path.exists(params_file):
            print(f"  FiberDock: params file not created: {params_file}")
            continue

        # Step 6: Run FiberDock (from fiberdock dir so it finds lib/)
        fiberdock_cwd = FIBERDOCK_DIR
        energy = _run_fiberdock(
            params_file, fiberdock_cwd, work_dir, receptor_base, ligand_base,
            pair_name=pair_name,
        )

        results.append((pair_name, energy))
        print(f"  FiberDock energy: {energy}")

    return results


def _split_merged_candidate(input_pdb, receptor_chains, ligand_chains, left, right):
    """Split one merged PRISM model into FiberDock receptor/ligand inputs."""
    receptor_set = set(receptor_chains)
    ligand_set = set(ligand_chains)
    receptor_atoms = 0
    ligand_atoms = 0
    with open(input_pdb) as source, open(left, "w") as receptor_out, open(right, "w") as ligand_out:
        for line in source:
            if not line.startswith(("ATOM", "HETATM")) or len(line) <= 21:
                continue
            if line[21] in receptor_set:
                receptor_out.write(line)
                receptor_atoms += 1
            elif line[21] in ligand_set:
                ligand_out.write(line)
                ligand_atoms += 1
        receptor_out.write("END\n")
        ligand_out.write("END\n")
    if receptor_atoms == 0 or ligand_atoms == 0:
        raise ValueError(
            f"Cannot split {input_pdb}: receptor atoms={receptor_atoms}, ligand atoms={ligand_atoms}"
        )


def refine_merged_candidates(candidates):
    """Run FiberDock while preserving PRISM's merged candidate tuple model."""
    from .pdb_download import split_target_id

    global FIBERDOCK_DIR
    FIBERDOCK_DIR = os.path.abspath(
        os.environ.get("PRISM_FIBERDOCK_DIR", FIBERDOCK_DIR)
    )
    if not _check_tools():
        raise RuntimeError(
            f"FiberDock tools are unavailable under {FIBERDOCK_DIR}; set PRISM_FIBERDOCK_DIR"
        )
    accepted = []
    for receptor, ligand, template, input_pdb in candidates:
        _, receptor_chains = split_target_id(receptor)
        _, ligand_chains = split_target_id(ligand)
        key = Path(input_pdb).stem
        split_root = Path("processed/fiberdock_refinement/inputs")
        split_root.mkdir(parents=True, exist_ok=True)
        receptor_pdb = split_root / f"{key}_R.pdb"
        ligand_pdb = split_root / f"{key}_L.pdb"
        _split_merged_candidate(
            input_pdb, receptor_chains, ligand_chains, receptor_pdb, ligand_pdb,
        )
        results = refine_pairs([(str(receptor_pdb), str(ligand_pdb))])
        if not results:
            continue
        destination = Path(STRUCTURES_DIR) / f"{key}_fiberdock.ref.pdb"
        if not destination.is_file():
            continue
        accepted.append((receptor, ligand, template, str(destination)))
    return accepted
