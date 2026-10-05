#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Isolated iRMSD alignment testcase that compares the current index-based
alignment logic with a paired-residue alignment variant.

This script is intentionally standalone and does NOT modify any production
pipeline code. It is safe to run ad hoc.
"""

import argparse
import os
import sys
from typing import List, Tuple

from Bio.Align import PairwiseAligner, substitution_matrices

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPTS_DIR)

import irmsd  # noqa: E402


def _parse_chain_list(chains: str) -> List[str]:
    chains = chains.strip()
    if not chains:
        return []
    return list(chains)


def _get_aligner() -> PairwiseAligner:
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    return aligner


def get_common_residues_paired(
    ref_file: str,
    ref_chains: str,
    model_file: str,
    model_chains: str,
) -> Tuple[list, list]:
    """
    Returns paired residue lists where each index maps ref<->model residues.
    This avoids index drift after gap handling.
    """
    aligner = _get_aligner()
    ref_component = []
    model_component = []

    for refch, mdlch in zip(_parse_chain_list(ref_chains), _parse_chain_list(model_chains)):
        model_res_list, model_seq = irmsd.parseChainResiduesFromStructure(model_file, mdlch)
        ref_res_list, ref_seq = irmsd.parseChainResiduesFromStructure(ref_file, refch)

        alignment = aligner.align(model_seq, ref_seq)[0]
        model_aligned = alignment[0]
        ref_aligned = alignment[1]

        m_idx = 0
        r_idx = 0
        for i in range(len(model_aligned)):
            m_char = model_aligned[i]
            r_char = ref_aligned[i]
            if m_char != '-' and r_char != '-':
                model_component.append(model_res_list[m_idx])
                ref_component.append(ref_res_list[r_idx])
                m_idx += 1
                r_idx += 1
            elif m_char == '-' and r_char != '-':
                r_idx += 1
            elif m_char != '-' and r_char == '-':
                m_idx += 1

    return ref_component, model_component


def compute_irmsd(
    model_file: str,
    model_rch: str,
    model_lch: str,
    ref_file: str,
    ref_rch: str,
    ref_lch: str,
    use_paired: bool,
) -> Tuple[float, int]:
    if len(model_rch) > 1:
        mdl_receptor_chain_combinations = irmsd.getSymmetricChainsList(model_file, model_rch)
    else:
        mdl_receptor_chain_combinations = [model_rch]

    if len(model_lch) > 1:
        mdl_ligand_chain_combinations = irmsd.getSymmetricChainsList(model_file, model_lch)
    else:
        mdl_ligand_chain_combinations = [model_lch]

    least_irmsd = 100000.0
    least_interface_len = 0

    for mdl_rchain in mdl_receptor_chain_combinations:
        for mdl_lchain in mdl_ligand_chain_combinations:
            if use_paired:
                ref_receptor, model_receptor = get_common_residues_paired(
                    ref_file, ref_rch, model_file, mdl_rchain
                )
                ref_ligand, model_ligand = get_common_residues_paired(
                    ref_file, ref_lch, model_file, mdl_lchain
                )
            else:
                ref_receptor, model_receptor = irmsd.getCommonResiduesList(
                    ref_file, ref_rch, model_file, mdl_rchain
                )
                ref_ligand, model_ligand = irmsd.getCommonResiduesList(
                    ref_file, ref_lch, model_file, mdl_lchain
                )

            ref_coords, model_coords = irmsd.getInterfaceResidues(
                ref_receptor, ref_ligand, model_receptor, model_ligand
            )
            interface_len = len(ref_coords)
            if interface_len == 0 or len(model_coords) == 0:
                curr_irmsd = 100.0
            else:
                curr_irmsd = irmsd.getRMSD(model_coords, ref_coords)

            if curr_irmsd < least_irmsd:
                least_irmsd = curr_irmsd
                least_interface_len = interface_len

    return least_irmsd, least_interface_len


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare current iRMSD alignment vs paired-residue variant."
    )
    parser.add_argument("model_file")
    parser.add_argument("model_receptor_chains")
    parser.add_argument("model_ligand_chains")
    parser.add_argument("ref_file")
    parser.add_argument("ref_receptor_chains")
    parser.add_argument("ref_ligand_chains")
    args = parser.parse_args()

    legacy_irmsd, legacy_if_len = compute_irmsd(
        args.model_file,
        args.model_receptor_chains,
        args.model_ligand_chains,
        args.ref_file,
        args.ref_receptor_chains,
        args.ref_ligand_chains,
        use_paired=False,
    )
    paired_irmsd, paired_if_len = compute_irmsd(
        args.model_file,
        args.model_receptor_chains,
        args.model_ligand_chains,
        args.ref_file,
        args.ref_receptor_chains,
        args.ref_ligand_chains,
        use_paired=True,
    )

    print("legacy_irmsd=%.3f interface_len=%d" % (legacy_irmsd, legacy_if_len))
    print("paired_irmsd=%.3f interface_len=%d" % (paired_irmsd, paired_if_len))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
