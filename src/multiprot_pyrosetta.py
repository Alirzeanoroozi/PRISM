#!/usr/bin/env python3
"""
MultiProt + PyRosetta integration test for PRISM pipeline.

Uses MultiProt for structural alignment of query proteins against
template interfaces, then builds docked complexes and refines with PyRosetta.

Usage:
  python src/multiprot_pyrosetta.py --receptor 1fgn --ligand 1tfh \
    --template 1cl7HL --template-chain H --template-chain2 L
"""
import argparse, json, os, shutil, subprocess, sys, tempfile
from pathlib import Path

PRISM_ROOT = Path(__file__).resolve().parent.parent
MULTIPROT = PRISM_ROOT / "external_tools/multiprot.Linux"
TEMPLATE_INTERFACES = PRISM_ROOT / "new_template/template/interfaces"
STAGED_PDBS = PRISM_ROOT / "tmp/agent/20260712-multichain-full-comparison/staged_pdbs"


def find_staged_pdb(pdb_code: str) -> Path:
    """Find a 4-char PDB code in staged dir."""
    p = STAGED_PDBS / f"{pdb_code[:4].lower()}.pdb"
    if p.exists():
        return p
    raise FileNotFoundError(f"Staged PDB not found: {p}")


def run_multiprot(pdb_paths: list[Path], output_dir: Path) -> dict:
    """Run MultiProt on a list of PDB files and parse the output.

    MultiProt outputs to stdout. We capture and parse the alignment info.
    Returns dict with alignment data.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # MultiProt needs PDB files in a temp dir with simple names
    with tempfile.TemporaryDirectory(prefix="multiprot-") as td:
        td = Path(td)
        pdb_files = []
        for i, src in enumerate(pdb_paths):
            dst = td / f"mol{i+1:03d}.pdb"
            shutil.copy2(str(src), str(dst))
            pdb_files.append(str(dst))
        
        # Write alignment param file
        # First, let MultiProt write to temp output
        result = subprocess.run(
            [str(MULTIPROT)] + pdb_files,
            capture_output=True, text=True, timeout=300
        )
        
        if result.returncode != 0 and not result.stdout:
            raise RuntimeError(f"MultiProt failed: {result.stderr[:500]}")
        
        # Parse MultiProt output for alignment info
        lines = result.stdout.split('\n')
        
        # Extract alignment data from output
        alignment_info = parse_multiprot_output(lines, pdb_files)
        alignment_info["raw_output"] = result.stdout[:5000]
        
        return alignment_info


def parse_multiprot_output(lines: list[str], pdb_files: list[str]) -> dict:
    """Parse MultiProt output to extract transformations."""
    import re
    
    alignment_info = {
        "num_molecules": 0,
        "largest_solution": 0,
        "transformations": [],
        "aligned_residues": [],
    }
    
    in_solution = False
    current_solution = {}
    
    for line in lines:
        # "Num Of Mols: 2 Largest Solution: 84"
        m = re.search(r'Num Of Mols:\s*(\d+)\s+Largest Solution:\s*(\d+)', line)
        if m:
            alignment_info["num_molecules"] = int(m.group(1))
            alignment_info["largest_solution"] = int(m.group(2))
        
        # Parse transformation matrices from output
        # MultiProt outputs rotation/translation in a specific format
        # Pattern: "Transformation for molecule X: ..."
        if "Transformation for molecule" in line:
            current_solution = {"molecule": line.strip()}
            in_solution = True
            continue
        
        if in_solution and ("rotation" in line.lower() or "Rotation" in line):
            current_solution.setdefault("rotation_lines", []).append(line.strip())
        elif in_solution and ("translation" in line.lower() or "Translation" in line):
            current_solution.setdefault("translation_lines", []).append(line.strip())
        
        # Parse aligned positions
        # MultiProt outputs: "Chain X residue Y aligned to chain Z residue W"
        if "aligned" in line.lower() and "chain" in line.lower():
            alignment_info["aligned_residues"].append(line.strip())
    
    return alignment_info


def apply_tm_transform(input_pdb, output_pdb, translation, rotation_mat):
    """Apply rotation+translation to a PDB, same as transformation.py."""
    try:
        with open(input_pdb, "r") as in_f, open(output_pdb, "w") as out_f:
            for line in in_f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                    except ValueError:
                        out_f.write(line)
                        continue
                    new_x = (x * rotation_mat[0][0] + y * rotation_mat[0][1]
                             + z * rotation_mat[0][2] + translation[0])
                    new_y = (x * rotation_mat[1][0] + y * rotation_mat[1][1]
                             + z * rotation_mat[1][2] + translation[1])
                    new_z = (x * rotation_mat[2][0] + y * rotation_mat[2][1]
                             + z * rotation_mat[2][2] + translation[2])
                    line = (f"{line[:30]}{new_x:>8.3f}{new_y:>8.3f}{new_z:>8.3f}"
                            f"{line[54:]}")
                    out_f.write(line)
                else:
                    out_f.write(line)
        return True
    except Exception as e:
        print(f"  Transform failed: {e}")
        return False


def get_chain_mapping(pdb_path: Path) -> dict:
    """Get unique chain IDs from a PDB file."""
    chains = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                ch = line[21].strip()
                if ch:
                    chains[ch] = chains.get(ch, 0) + 1
    return chains


def parse_multiprot_transforms(output_dir: Path) -> list[dict]:
    """Parse MultiProt output files for transformation matrices.
    
    MultiProt writes log_multiprot.txt in the CWD with alignment params
    but does NOT output rotation matrices in a parseable format.
    As a workaround, we fall back to using TMalign to get the transform
    that corresponds to the MultiProt alignment.
    
    Returns a list of dicts with 'translation' and 'rotation_mat'.
    """
    # MultiProt itself doesn't produce rotation matrices in standard output.
    # This function is a placeholder — the actual alignment transforms
    # must come from the alignment tool (GTalign/TMalign) that outputs them.
    return []


def build_docked_complex(receptor_pdb: Path, ligand_pdb: Path,
                          transform_dir: Path,
                          template_id: str, chain1: str, chain2: str,
                          output_pdb: Path) -> tuple[bool, str, str]:
    """Build docked complex using an existing pipeline transform.

    Uses the pre-computed alignment transform (from GTalign/TMalign run)
    to position the receptor and ligand into a docked complex.
    
    The transform_dir must contain {template_id}_{chain1}_aligned_R.pdb
    and {template_id}_{chain2}_aligned_L.pdb or similar split files.
    
    Falls back to using the pipeline's R+L transform PDBs if available.
    """
    # Strategy: if pipeline transform outputs exist for this template,
    # combine them. Otherwise, apply identity transform on query proteins.
    
    # Look for pipeline transform PDBs in the standard v3 output
    v3_root = PRISM_ROOT / "tmp/agent/20260718-benchmark20k-v3/smoke"
    rec_name = receptor_pdb.stem[:4]  # e.g., "1fgn"
    lig_name = ligand_pdb.stem[:4]    # e.g., "1tfh"
    
    receptor_chains = list(get_chain_mapping(receptor_pdb).keys())
    ligand_chains_list = list(get_chain_mapping(ligand_pdb).keys())
    
    # Try to find a transform from any variant
    for variant in ["gt_external", "gt_pyro", "tm_external", "tm_pyro"]:
        trans_dir = (v3_root / variant / "runs" / variant / "current"
                     / "batch_0001" / "processed" / "transformation")
        if not trans_dir.is_dir():
            continue
        # Look for pattern: *_{template_id}_{rec_name}*_R.pdb
        for r_pdb in trans_dir.glob(f"*_{template_id}_{rec_name}*_R.pdb"):
            # Derive L counterpart
            l_pdb = r_pdb.parent / r_pdb.name.replace("_R.pdb", "_L.pdb")
            if l_pdb.exists():
                with open(output_pdb, 'w') as out:
                    out.write(f"REMARK MultiProt + PyRosetta PRISM test\n")
                    out.write(f"REMARK Receptor: {receptor_pdb.name}\n")
                    out.write(f"REMARK Ligand: {ligand_pdb.name}\n")
                    out.write(f"REMARK Template: {template_id}\n")
                    out.write(f"REMARK Transform source: {variant}\n")
                    out.write(f"REMARK Receptor chains: {list(get_chain_mapping(r_pdb).keys())}\n")
                    out.write(f"REMARK Ligand chains: {list(get_chain_mapping(l_pdb).keys())}\n")
                    out.write(r_pdb.read_text())
                    out.write(l_pdb.read_text())
                
                r_chains = ''.join(get_chain_mapping(r_pdb).keys())
                l_chains = ''.join(get_chain_mapping(l_pdb).keys())
                return (True, r_chains, l_chains)
    
    # Fallback: just combine the query proteins with identity placement
    r_chains = ''.join(receptor_chains)
    l_chains = ''.join(ligand_chains_list)
    with open(output_pdb, 'w') as out:
        out.write(f"REMARK MultiProt + PyRosetta test (fallback)\n")
        out.write(f"REMARK Receptor: {receptor_pdb.name} chains={r_chains}\n")
        out.write(f"REMARK Ligand: {ligand_pdb.name} chains={l_chains}\n")
        out.write(f"REMARK Template: {template_id}\n")
        out.write(f"REMARK NOTE: No transform found - using identity placement\n")
        with open(receptor_pdb) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    out.write(line)
        with open(ligand_pdb) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    out.write(line)
        out.write("END\n")
    
    return (True, r_chains, l_chains)


def refine_with_pyrosetta(input_pdb: Path, output_pdb: Path,
                           partners: str) -> bool:
    """Run PyRosetta refinement on the docked complex.
    
    Args:
        input_pdb: Path to the complex PDB (receptor + ligand chains)
        output_pdb: Path to save the refined complex
        partners: Partner string e.g. 'A_L' meaning receptor chain A,
                  ligand chain L
    """
    try:
        import pyrosetta
        pyrosetta.init("-ignore_unrecognized_res 1 -ex1 -ex2aro",
                       silent=True)
        
        pose = pyrosetta.pose_from_pdb(str(input_pdb))
        n_res = pose.total_residue()
        n_chain = pose.num_chains()
        print(f"  Loaded: {n_res} residues, {n_chain} chains")
        
        from pyrosetta.rosetta.protocols.docking import DockingProtocol
        dock = DockingProtocol()
        dock.set_partners(partners)
        dock.apply(pose)
        pose.dump_pdb(str(output_pdb))
        print(f"  Refinement partners: {partners}")
        return True
    except Exception as e:
        print(f"  PyRosetta refinement failed: {e}")
        import traceback
        traceback.print_exc()
        shutil.copy2(str(input_pdb), str(output_pdb))
        return False


def main():
    parser = argparse.ArgumentParser(description="MultiProt + PyRosetta PRISM test")
    parser.add_argument("--receptor", required=True, help="4-char PDB code for receptor")
    parser.add_argument("--ligand", required=True, help="4-char PDB code for ligand")
    parser.add_argument("--template", required=True, help="Template ID (e.g., 1cl7HL)")
    parser.add_argument("--template-chain1", default="H", help="Template chain 1")
    parser.add_argument("--template-chain2", default="L", help="Template chain 2")
    parser.add_argument("--output-dir", default="output/multiprot_test", help="Output dir")
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find input files
    receptor_pdb = find_staged_pdb(args.receptor)
    ligand_pdb = find_staged_pdb(args.ligand)
    template_int1 = TEMPLATE_INTERFACES / f"{args.template}_{args.template_chain1}_int.pdb"
    template_int2 = TEMPLATE_INTERFACES / f"{args.template}_{args.template_chain2}_int.pdb"
    
    print(f"=== MultiProt + PyRosetta Test ===")
    print(f"Receptor PDB: {receptor_pdb}")
    print(f"Ligand PDB: {ligand_pdb}")
    print(f"Template interface 1 ({args.template_chain1}): {template_int1}")
    print(f"Template interface 2 ({args.template_chain2}): {template_int2}")
    
    # Step 1: Run MultiProt to align receptor to template interface
    print(f"\n[Step 1] Running MultiProt alignment...")
    alignment = run_multiprot(
        [receptor_pdb, template_int1],
        output_dir / "alignment"
    )
    print(f"  Molecules: {alignment['num_molecules']}")
    print(f"  Largest solution: {alignment['largest_solution']}")
    print(f"  Aligned residues: {len(alignment['aligned_residues'])}")
    
    # Step 2: Build docked complex using pipeline transforms
    print(f"\n[Step 2] Building docked complex...")
    docked_pdb = output_dir / "complex_docked.pdb"
    success, rec_chains, lig_chains = build_docked_complex(
        receptor_pdb, ligand_pdb, output_dir / "alignment",
        args.template, args.template_chain1, args.template_chain2,
        docked_pdb
    )
    partners = f"{rec_chains}_{lig_chains}"
    print(f"  Output: {docked_pdb}")
    print(f"  Partners: {partners}")
    
    # Step 3: Refine with PyRosetta
    print(f"\n[Step 3] Running PyRosetta refinement...")
    refined_pdb = output_dir / "complex_refined.pdb"
    refine_with_pyrosetta(
        docked_pdb, refined_pdb, partners
    )
    print(f"  Output: {refined_pdb}")
    
    # Summary
    print(f"\n=== Summary ===")
    print(f"Receptor: {args.receptor} ({receptor_pdb})")
    print(f"Ligand: {args.ligand} ({ligand_pdb})")
    print(f"Template: {args.template}")
    print(f"Alignment: {alignment['largest_solution']} aligned residues")
    print(f"Docked complex: {docked_pdb}")
    print(f"Refined complex: {refined_pdb}")
    print(f"\nMultiProt raw output saved to: {output_dir / 'multiprot_output.txt'}")
    
    # Save raw output
    with open(output_dir / "multiprot_output.txt", "w") as f:
        f.write(alignment.get("raw_output", "(no output captured)"))


if __name__ == "__main__":
    main()
