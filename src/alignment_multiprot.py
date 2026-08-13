#!/usr/bin/env python3
"""
MultiProt structural alignment module for the PRISM pipeline.

Produces PRISM-compatible alignment JSONs (same format as TMalign/GTalign).
Current mode computes a Kabsch transform from matched coordinates; legacy
compatibility mode preserves MultiProt's native transforms and solutions.
"""

import csv
import json
import math
import os
import shutil
import subprocess
import tempfile
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

MULTIPROT = os.environ.get("PRISM_MULTIPROT", "external_tools/multiprot.Linux")


def _has_ca_atoms(path):
    try:
        with open(path) as handle:
            return any(line.startswith("ATOM") and " CA " in line[:20] for line in handle)
    except OSError:
        return False


def _write_empty_alignment(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        json.dump({
            "match_count": 0,
            "tm_score": 0.0,
            "match_dict": {},
            "status": "alignment_unavailable",
            "aligner": "MultiProt",
        }, handle)

# Resolve MultiProt path once at module level
_MP_ABS = None
if MULTIPROT.startswith("/"):
    _MP_ABS = MULTIPROT
else:
    # Try to find it relative to the repo root
    for _prefix in [os.getcwd(), os.path.dirname(os.path.dirname(os.path.abspath(__file__)))]:
        _candidate = os.path.join(_prefix, MULTIPROT)
        if os.path.exists(_candidate):
            _MP_ABS = _candidate
            break
    if _MP_ABS is None:
        _resolved = shutil.which(MULTIPROT)
        if _resolved:
            _MP_ABS = _resolved


def _check_seccomp():
    """Check if seccomp is enabled (blocks 32-bit binaries like MultiProt).

    Returns True if MultiProt can run, False if seccomp blocks it.
    """
    # Allow override for testing on compute nodes with known 32-bit support
    if os.environ.get("PRISM_MULTIPROT_FORCE", "0") == "1":
        return True
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("Seccomp:"):
                    seccomp_level = int(line.split(":")[1].strip())
                    if seccomp_level >= 2:
                        return False
                    return True
    except (IOError, ValueError):
        pass
    return True  # Assume OK if we can't check


def _parse_multiprot_solution(sol_res_path):
    """Parse MultiProt's 2_sol.res to extract the best alignment solution.

    Returns (match_dict, rmsd) where match_dict maps interface residues
    to query residues as {interface_chain.res_type.res_num: query_chain.res_type.res_num}
    """
    match_dict = {}
    best_rmsd = float('inf')

    with open(sol_res_path) as f:
        content = f.read()

    # Split into solutions
    solutions = re.split(r'Solution Num : \d+', content)

    # Skip header (first split before Solution Num : 0)
    if not solutions:
        return {}, float('inf')

    best_solution = None
    best_score = 0

    for sol_text in solutions:
        # Extract score
        score_match = re.search(r'Mult Corres Score : (\d+)', sol_text)
        if not score_match:
            continue
        score = int(score_match.group(1))

        if score > best_score:
            best_score = score
            best_solution = sol_text

    if not best_solution:
        return {}, float('inf')

    # Extract RMSD
    rmsd_match = re.search(r'RMSD : ([\d.e+\-]+)', best_solution)
    if rmsd_match:
        best_rmsd = float(rmsd_match.group(1))

    # Extract match list
    in_match_list = False
    for line in best_solution.split('\n'):
        line = line.strip()

        if line.startswith('Match List'):
            in_match_list = True
            continue
        if line.startswith('End of Match List'):
            break
        if not in_match_list:
            continue

        # Parse "A.N.103  A.N.103" format
        parts = line.split()
        if len(parts) >= 2:
            interface_res = parts[0].strip()  # e.g. "A.N.103"
            query_res = parts[1].strip()      # e.g. "A.N.103"
            if interface_res and query_res and '.' in interface_res:
                match_dict[interface_res] = query_res

    return match_dict, best_rmsd


def _parse_multiprot_solutions(sol_res_path, max_solutions=3):
    """Parse retained solutions using the historical MultiProt field contract."""
    with open(sol_res_path) as f:
        lines = f.readlines()

    solutions = []
    index = 0
    while index < len(lines) and len(solutions) < max_solutions:
        if not lines[index].startswith("Solution Num"):
            index += 1
            continue
        number = int(lines[index].split(":", 1)[1].strip())
        fields = {}
        index += 1
        while index < len(lines):
            line = lines[index].strip()
            if line.startswith("Solution Num"):
                break
            if line.startswith("Mult Corres Score"):
                fields["match_count"] = int(line.split(":", 1)[1].strip())
            elif line.startswith("Reference Molecule"):
                fields["reference_molecule"] = int(line.split(":", 1)[1].strip())
            elif line.startswith("Trans"):
                fields["trans"] = [float(value) for value in line.split(":", 1)[1].split()]
            elif line.startswith("RMSD"):
                fields["rmsd"] = float(line.split(":", 1)[1].strip())
            elif line.startswith("Match List"):
                match_dict = {}
                reference = fields.get("reference_molecule", 0)
                index += 1
                expected = fields.get("match_count", 0)
                while index < len(lines) and len(match_dict) < expected:
                    match_line = lines[index].strip()
                    if match_line.startswith("End of Match List"):
                        break
                    parts = match_line.split()
                    if len(parts) >= 2:
                        match_dict[parts[reference]] = parts[1 - reference]
                    index += 1
                fields["match_dict"] = match_dict
            index += 1

        fields.setdefault("match_count", len(fields.get("match_dict", {})))
        fields.setdefault("reference_molecule", 0)
        fields.setdefault("trans", [0.0] * 6)
        fields.setdefault("rmsd", float("inf"))
        fields.setdefault("match_dict", {})
        fields["solution_number"] = number
        solutions.append(fields)
    return solutions


def _legacy_solution_to_alignment(solution):
    """Convert legacy row-vector rotation fields to PRISM column transforms."""
    phi, theta, psi, x_translate, y_translate, z_translate = solution.get(
        "trans", [0.0] * 6
    )
    phi += math.pi
    theta += math.pi
    psi += math.pi
    cos_phi, sin_phi = math.cos(phi), math.sin(phi)
    cos_theta, sin_theta = math.cos(theta), math.sin(theta)
    cos_psi, sin_psi = math.cos(psi), math.sin(psi)
    legacy_rotation = np.array([
        [cos_theta * cos_psi, cos_theta * sin_psi, -sin_theta],
        [
            -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi,
            cos_phi * cos_psi + sin_phi * sin_theta * sin_psi,
            sin_phi * cos_theta,
        ],
        [
            sin_phi * sin_psi + cos_phi * sin_theta * cos_psi,
            -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi,
            cos_phi * cos_theta,
        ],
    ])
    trans = np.array([x_translate, y_translate, z_translate])
    if solution.get("reference_molecule", 0) == 0:
        rotation_mat = legacy_rotation.T
        translation = trans
    else:
        rotation_mat = legacy_rotation
        translation = -(legacy_rotation @ trans)
    return {
        "rotation_mat": rotation_mat.tolist(),
        "translation": translation.tolist(),
        "match_count": solution.get("match_count", 0),
        "reference_molecule": solution.get("reference_molecule", 0),
        "trans": solution.get("trans", [0.0] * 6),
        "rmsd": solution.get("rmsd", float("inf")),
        "match_dict": solution.get("match_dict", {}),
        "solution_number": solution.get("solution_number"),
    }


# 1-letter to 3-letter amino acid code mapping
_1TO3_AA = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}


def _parse_pdb_coords_by_resnum(pdb_path):
    """Parse CA coordinates keyed by resnum (int), ignoring chain ID.

    MultiProt labels molecules by molecule index (0, 1), not PDB chain IDs.
    We match by residue number plus amino acid to handle this.
    """
    coords_by_resnum = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[13:15].strip() == "CA":
                res_num = int(line[22:26].strip())
                res_name = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                # Use both resnum and resname as key to avoid ambiguity
                key = (res_num, res_name)
                if key not in coords_by_resnum:
                    coords_by_resnum[key] = np.array([x, y, z])
    return coords_by_resnum


def _compute_transform_from_matches(match_dict, query_path, interface_path):
    """Compute rotation and translation from the matched residue pairs using Kabsch algorithm.

    MultiProt labels residues as MolID.AA.ResNum where MolID is molecule index
    (0=first file=query, 1=second file=interface in MultiProt's invocation).
    We match by residue number + amino acid, ignoring the molecule ID prefix.

    Returns (translation, rotation_mat, rmsd).
    """
    q_coords = _parse_pdb_coords_by_resnum(query_path)
    i_coords = _parse_pdb_coords_by_resnum(interface_path)

    # Helper to extract (resnum, 3-letter AA) from a MultiProt key
    def _mp_key_to_res(mp_key):
        parts = mp_key.split('.')
        if len(parts) != 3:
            return None
        aa1 = parts[1].upper()
        try:
            res_num = int(parts[2])
        except ValueError:
            return None
        aa3 = _1TO3_AA.get(aa1)
        if aa3 is None:
            return None
        return (res_num, aa3)

    # Collect matched CA coordinates
    q_pts = []
    i_pts = []

    for mp_interface_key, mp_query_key in match_dict.items():
        q_res = _mp_key_to_res(mp_query_key)
        i_res = _mp_key_to_res(mp_interface_key)
        if q_res is None or i_res is None:
            continue
        if i_res in i_coords and q_res in q_coords:
            i_pts.append(i_coords[i_res])
            q_pts.append(q_coords[q_res])

    if len(q_pts) < 3:
        return [0.0, 0.0, 0.0], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], float('inf')

    q_pts = np.array(q_pts)
    i_pts = np.array(i_pts)

    # Centroids
    q_centroid = np.mean(q_pts, axis=0)
    i_centroid = np.mean(i_pts, axis=0)

    # Cross-covariance matrix
    H = (q_pts - q_centroid).T @ (i_pts - i_centroid)

    # SVD
    U, _, Vt = np.linalg.svd(H)

    # Ensure proper rotation (no reflection)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Rotation matrix (as nested list)
    rotation_mat = R.tolist()

    # Translation from query coordinates into the template/interface frame.
    t_vec = i_centroid - R @ q_centroid
    translation = t_vec.tolist()

    # RMSD
    transformed_q = (R @ q_pts.T).T + t_vec
    rmsd = np.sqrt(np.mean(np.sum((i_pts - transformed_q)**2, axis=1)))

    return translation, rotation_mat, rmsd


def _align_one(task):
    """Align one (query, template, chain) pair using only MultiProt.

    Runs MultiProt for structural alignment, parses the residue correspondence
    and computes the rotation/translation matrix directly from the match list
    using the Kabsch algorithm. No TMalign fallback.
    """
    (
        query, template, chain, alignment_root, multiprot_mode,
        multiprot_params, multiprot_solutions,
    ) = task
    query_path = f"processed/surface_extraction/{query}.asa.pdb"
    interface_path = f"templates/interfaces/{template}_{chain}_int.pdb"
    output_path = os.path.join(alignment_root, f"{query}_{template}_{chain}.json")

    if not _has_ca_atoms(query_path) or not os.path.exists(interface_path):
        _write_empty_alignment(output_path)
        return None

    # Check seccomp before running 32-bit MultiProt
    if not _check_seccomp():
        print(f"MultiProt skipped for {query} and {template}_{chain}: "
              "seccomp blocks 32-bit binaries on this node")
        _write_empty_alignment(output_path)
        return None

    try:
        with tempfile.TemporaryDirectory(prefix="multiprot-", dir=alignment_root) as scratch:
            scratch = Path(scratch)

            # Step 1: Copy PDBs with short names for MultiProt
            q_pdb = scratch / "query.pdb"
            i_pdb = scratch / "interface.pdb"
            shutil.copy2(query_path, str(q_pdb))
            shutil.copy2(interface_path, str(i_pdb))
            if multiprot_params:
                shutil.copy2(multiprot_params, str(scratch / "params.txt"))

            # Step 2: Use pre-resolved absolute path to MultiProt
            # (ThreadPoolExecutor threads inherit CWD but the binary path
            # must be absolute since subprocess cwd is set to a temp dir)
            global _MP_ABS
            mp_exe = _MP_ABS or shutil.which(MULTIPROT) or os.path.abspath(MULTIPROT)
            if not os.path.exists(mp_exe):
                raise RuntimeError(f"MultiProt not found at {mp_exe}")

            mp_inputs = (
                [str(i_pdb), str(q_pdb)]
                if multiprot_mode == "legacy_compatible"
                else [str(q_pdb), str(i_pdb)]
            )
            mp_result = subprocess.run(
                [mp_exe, *mp_inputs],
                capture_output=True, text=True, timeout=300,
                cwd=str(scratch)
            )

            if mp_result.returncode != 0 and not mp_result.stdout:
                raise RuntimeError(f"MultiProt failed (exit={mp_result.returncode}): "
                                   f"{(mp_result.stderr or '').strip()[:300]}")

            # Parse aligned residue count from stdout
            largest_solution = 0
            for line in mp_result.stdout.split('\n'):
                if "Largest Solution" in line:
                    try:
                        largest_solution = int(line.split(":")[-1].strip())
                    except (ValueError, IndexError):
                        pass

            if multiprot_mode == "current" and largest_solution < 5:
                _write_empty_alignment(output_path)
                return None

            # Step 3: Parse 2_sol.res for the actual alignment
            sol_res = scratch / "2_sol.res"
            if not sol_res.exists():
                raise RuntimeError("MultiProt did not create 2_sol.res")

            if multiprot_mode == "legacy_compatible":
                parsed_solutions = _parse_multiprot_solutions(
                    str(sol_res), max_solutions=multiprot_solutions
                )
                if not parsed_solutions:
                    _write_empty_alignment(output_path)
                    return None
                converted = [
                    _legacy_solution_to_alignment(solution)
                    for solution in parsed_solutions
                ]
                primary = converted[0]
                multi_dict = {
                    "match_count": primary["match_count"],
                    "translation": primary["translation"],
                    "rotation_mat": primary["rotation_mat"],
                    "match_dict": primary["match_dict"],
                    "tm_score": 0.0,
                    "tm_score_contract": "multiprot_legacy_native",
                    "score_gate_contract": "native_match_count_and_coverage",
                    "rmsd": primary["rmsd"],
                    "status": "success",
                    "aligner": "MultiProt",
                    "multiprot_mode": "legacy_compatible",
                    "multiprot_solution": primary["solution_number"],
                    "reference_molecule": primary["reference_molecule"],
                    "multiprot_trans": primary["trans"],
                    "multiprot_solutions": converted,
                }
                with open(output_path, "w") as f:
                    json.dump(multi_dict, f)
                return _summary_row(query, template, chain, query_path, interface_path, multi_dict)

            match_dict, rmsd = _parse_multiprot_solution(str(sol_res))

            if len(match_dict) < 3:
                _write_empty_alignment(output_path)
                return None

            # Step 4: Compute transform from matched residues via Kabsch
            translation, rotation_mat, kabsch_rmsd = _compute_transform_from_matches(
                match_dict, query_path, interface_path
            )

            if kabsch_rmsd == float('inf'):
                _write_empty_alignment(output_path)
                return None

            # Step 5: Write alignment JSON
            # Keep the compatibility field for existing consumers, but mark it
            # explicitly: this RMSD-derived value is not a TMalign TM-score and
            # must not be used by the shared TMalign score gate.
            tm_score = max(0.0, 1.0 - (kabsch_rmsd / 10.0))
            multi_dict = {
                "match_count": len(match_dict),
                "translation": translation,
                "rotation_mat": rotation_mat,
                "match_dict": match_dict,
                "tm_score": tm_score,
                "tm_score_contract": "multiprot_kabsch_rmsd_proxy",
                "score_gate_contract": "native_match_count_and_coverage",
                "rmsd": kabsch_rmsd,
                "status": "success",
                "aligner": "MultiProt",
                "multiprot_solution": largest_solution,
            }
            os.makedirs(alignment_root, exist_ok=True)
            with open(output_path, "w") as f:
                json.dump(multi_dict, f)
            return _summary_row(query, template, chain, query_path, interface_path, multi_dict)

    except (OSError, RuntimeError, ValueError, IndexError, subprocess.TimeoutExpired, KeyError, AttributeError) as exc:
        import traceback
        print(f"MultiProt failed for {query} and {template}_{chain}: {exc}")
        traceback.print_exc()
        _write_empty_alignment(output_path)

    return None


def _count_ca(path):
    with open(path) as handle:
        return sum(1 for line in handle if line.startswith("ATOM") and " CA " in line[:20])


def _summary_row(query, template, chain, query_path, interface_path, payload):
    return {
        "protein": query,
        "template": template,
        "chain": chain,
        "match_count": int(payload["match_count"]),
        "tm_score": float(payload.get("tm_score", 0.0)),
        "len_target": _count_ca(query_path),
        "len_template": _count_ca(interface_path),
        "translation": json.dumps(payload["translation"]),
        "rotation_mat": json.dumps(payload["rotation_mat"]),
    }


def _passes_native_gate(row):
    template_length = row["len_template"]
    coverage = row["match_count"] / template_length if template_length else 0.0
    return row["match_count"] >= 10 and coverage >= 0.30


def align_multiprot(
    queries,
    templates,
    output_dir="processed/alignment",
    max_workers=8,
    multiprot_path=None,
    multiprot_mode="current",
    multiprot_params=None,
    multiprot_solutions=3,
):
    """MultiProt alignment entry point for the PRISM pipeline.

    Aligns query proteins against template interfaces and writes the JSON and
    CSV contracts consumed by the transformation stage.

    Args:
        queries: List of query protein IDs (e.g. ['1fgnHL', '1tfhA'])
        templates: List of template IDs (e.g. ['1cl7HL', '1a0cCD'])
        output_dir: Directory for alignment JSONs
        max_workers: Number of parallel worker threads
    """
    global MULTIPROT, _MP_ABS
    if multiprot_mode not in {"current", "legacy_compatible"}:
        raise ValueError("multiprot_mode must be current or legacy_compatible")
    if multiprot_solutions < 1:
        raise ValueError("multiprot_solutions must be at least 1")
    if multiprot_params:
        multiprot_params = os.path.abspath(str(multiprot_params))
        if not os.path.isfile(multiprot_params):
            raise FileNotFoundError(f"MultiProt params file not found: {multiprot_params}")
    if multiprot_path:
        MULTIPROT = str(multiprot_path)
        _MP_ABS = os.path.abspath(MULTIPROT) if os.path.sep in MULTIPROT else shutil.which(MULTIPROT)

    alignment_root = os.path.abspath(output_dir)
    os.makedirs(alignment_root, exist_ok=True)

    tasks = []
    for query in queries:
        for template in templates:
            for chain in template[4:]:
                tasks.append((
                    query, template, chain, alignment_root, multiprot_mode,
                    multiprot_params, multiprot_solutions,
                ))

    total = len(tasks)
    if total == 0:
        return

    print(f"MultiProt alignment: {total} pairs to process")

    rows = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_align_one, t): t for t in tasks}
        for i, future in enumerate(as_completed(futures), 1):
            row = future.result()
            if row is not None and _passes_native_gate(row):
                rows.append(row)
            if i % 100 == 0:
                print(f"MultiProt progress: {i}/{total}")

    fieldnames = [
        "protein", "template", "chain", "match_count", "tm_score",
        "len_target", "len_template", "translation", "rotation_mat",
    ]
    for target in queries:
        target_rows = [row for row in rows if row["protein"] == target]
        if not target_rows:
            continue
        target_rows.sort(key=lambda row: (row["match_count"], row["tm_score"]), reverse=True)
        with open(os.path.join(alignment_root, f"{target}.csv"), "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(target_rows)
    print(f"MultiProt alignment finished: {total} pairs processed; {len(rows)} passed native gates")
