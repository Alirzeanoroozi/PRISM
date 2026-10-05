#!/usr/bin/env python3
"""
Compute true TM-score from MultiProt match_dict using the standard formula.

TM-score = (1/L_target) * sum_i 1/(1 + (d_i/d0)^2)
where d0 = 1.24 * (L_target - 15)^(1/3) - 1.8

This enables fair comparison with TMalign/GTalign which report native TM-scores.
"""

import json
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, Selection

_1TO3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}

def parse_mp_key(key: str):
    """Parse MultiProt key like 'C.T.110' -> (chain, resnum, resname3)"""
    parts = key.split('.')
    if len(parts) != 3:
        return None
    chain = parts[0]
    aa1 = parts[1].upper()
    try:
        resnum = int(parts[2])
    except ValueError:
        return None
    aa3 = _1TO3.get(aa1)
    if aa3 is None:
        return None
    return (chain, resnum, aa3)

def get_ca_coords(pdb_path):
    """Return dict: (chain, resnum, resname) -> CA coord array"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    coords = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ':  # skip hetero
                    continue
                resnum = residue.id[1]
                resname = residue.resname
                if 'CA' in residue:
                    ca = residue['CA']
                    coords[(chain.id, resnum, resname)] = ca.get_coord()
    return coords

def compute_tm_score(match_dict, query_path, interface_path):
    """Compute true TM-score from MultiProt match_dict."""
    q_coords = get_ca_coords(query_path)
    i_coords = get_ca_coords(interface_path)
    
    # Collect matched CA pairs
    q_pts = []
    i_pts = []
    
    for mp_i_key, mp_q_key in match_dict.items():
        i_parsed = parse_mp_key(mp_i_key)
        q_parsed = parse_mp_key(mp_q_key)
        if i_parsed is None or q_parsed is None:
            continue
        if i_parsed in i_coords and q_parsed in q_coords:
            i_pts.append(i_coords[i_parsed])
            q_pts.append(q_coords[q_parsed])
    
    if len(q_pts) < 3:
        return 0.0
    
    q_pts = np.array(q_pts)
    i_pts = np.array(i_pts)
    
    # Kabsch alignment
    q_centroid = q_pts.mean(axis=0)
    i_centroid = i_pts.mean(axis=0)
    H = (q_pts - q_centroid).T @ (i_pts - i_centroid)
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Align interface to query
    i_aligned = (R @ i_pts.T).T + q_centroid - R @ i_centroid
    
    # Compute distances
    d = np.sqrt(np.sum((q_pts - i_aligned)**2, axis=1))
    
    # TM-score formula
    L_target = len(q_pts)  # length of aligned region (could use full query length)
    if L_target > 15:
        d0 = 1.24 * ((L_target - 15)**(1/3)) - 1.8
    else:
        d0 = 0.5
    if d0 <= 0:
        d0 = 0.5
    
    tm = np.mean(1.0 / (1.0 + (d / d0)**2))
    return float(tm)

def process_alignment_file(json_path):
    """Add true_tm_score to an alignment JSON."""
    with open(json_path) as f:
        data = json.load(f)
    
    if data.get('status') != 'success':
        return
    
    # Extract paths from filename pattern: query_template_chain.json
    # e.g., 5zngA_1a0cCD_C.json
    fname = Path(json_path).stem
    parts = fname.split('_')
    if len(parts) < 3:
        return
    
    query = parts[0]
    template = parts[1]
    chain = parts[2]
    
    query_path = f"processed/surface_extraction/{query}.asa.pdb"
    interface_path = f"templates/interfaces/{template}_{chain}_int.pdb"
    
    if not Path(query_path).exists() or not Path(interface_path).exists():
        return
    
    match_dict = data.get('match_dict', {})
    if not match_dict:
        return
    
    true_tm = compute_tm_score(match_dict, query_path, interface_path)
    data['true_tm_score'] = true_tm
    data['tm_score_contract'] = 'standard_length_normalized'
    
    with open(json_path, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"{fname}: true_tm={true_tm:.4f}, proxy_tm={data.get('tm_score', 0):.4f}, matches={data.get('match_count', 0)}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        for f in sys.argv[1:]:
            process_alignment_file(f)
    else:
        # Process all
        for f in Path("processed/alignment").glob("*.json"):
            process_alignment_file(f)
