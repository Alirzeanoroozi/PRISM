#!/usr/bin/env python3
"""Find contacting residues between chains in a PDB structure."""

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa


def get_contacting_residues(
    pdb_path: str,
    chain_id1: str,
    chain_id2: str,
    distance_threshold: float = 5.5,
    model_id: int = 0,
    use_backbone_only: bool = False,
):
    """Get residues that contact between two chains.

    Two residues are in contact if any atom pair from the residues is within
    the distance threshold. By default uses all heavy atoms; set use_backbone_only
    for N, CA, C, O only.

    Args:
        pdb_path: Path to PDB file.
        chain_id1: First chain ID.
        chain_id2: Second chain ID.
        distance_threshold: Maximum distance (Å) for atoms to be in contact.
            Default 5.5 Å is common for protein-protein interfaces.
        model_id: Model index (default 0 for single-model structures).
        use_backbone_only: If True, use only N, CA, C, O atoms.

    Returns:
        List of tuples (res1, res2) where each is a Bio.PDB Residue object
        that forms a contact. Also returns (set of chain1 res_ids, set of chain2 res_ids)
        for convenience.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    model = structure[model_id]

    if chain_id1 not in model or chain_id2 not in model:
        raise ValueError(f"Chains {chain_id1} or {chain_id2} not found in structure")

    chain1 = model[chain_id1]
    chain2 = model[chain_id2]

    if use_backbone_only:
        atom_names = {"N", "CA", "C", "O"}
    else:
        atom_names = None  # all heavy atoms

    def get_atoms(residue):
        for atom in residue:
            if atom.element == "H":
                continue
            if atom_names is not None and atom.get_name() not in atom_names:
                continue
            yield atom

    contacts = []
    seen_res1, seen_res2 = set(), set()

    for res1 in chain1:
        if not is_aa(res1):
            continue
        coords1 = np.array([a.get_coord() for a in get_atoms(res1)])
        if len(coords1) == 0:
            continue

        for res2 in chain2:
            if not is_aa(res2):
                continue
            coords2 = np.array([a.get_coord() for a in get_atoms(res2)])
            if len(coords2) == 0:
                continue

            # Minimum distance between any atom pair
            dists = np.linalg.norm(coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :], axis=2)
            min_dist = np.min(dists)

            if min_dist <= distance_threshold:
                contacts.append((res1, res2))
                seen_res1.add(res1.get_id())
                seen_res2.add(res2.get_id())

    return contacts, seen_res1, seen_res2


def get_contacting_residue_ids(
    pdb_path: str,
    chain_id1: str,
    chain_id2: str,
    distance_threshold: float = 5.5,
    model_id: int = 0,
    use_backbone_only: bool = False,
) -> tuple[list[tuple], set, set]:
    """Get contacting residue IDs (convenience wrapper).

    Returns:
        (contact_pairs, chain1_res_ids, chain2_res_ids)
        where contact_pairs is [(res_id1, res_id2), ...] and res_ids are
        (hetflag, resseq, icode) tuples.
    """
    contacts, res_ids1, res_ids2 = get_contacting_residues(
        pdb_path, chain_id1, chain_id2, distance_threshold, model_id, use_backbone_only
    )
    pairs = [(r1.get_id(), r2.get_id()) for r1, r2 in contacts]
    return pairs, res_ids1, res_ids2
