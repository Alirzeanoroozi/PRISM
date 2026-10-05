#!/usr/bin/env python3
"""Compute grouped benchmark iRMSD with safe symmetric-chain expansion.

This compatibility entry point deliberately reuses the established residue
pairing, interface, and RMSD implementation from ``irmsd.py``.  It replaces
only that script's mutable symmetric-chain grouping, which fails for three or
more near-identical chains and does not form the Cartesian product when more
than one symmetry group is present.
"""

from __future__ import annotations

import itertools
import sys
import warnings
from functools import lru_cache
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import irmsd as legacy  # noqa: E402


def symmetric_chain_orders(model_file: str | Path, model_chains: str) -> list[str]:
    """Return every independent symmetric-chain ordering exactly once."""

    chains = list(model_chains)
    if len(chains) < 2:
        return [model_chains]

    sequences = {
        chain: legacy.parseChainResiduesFromStructure(str(model_file), chain)[1]
        for chain in chains
    }
    parent = {chain: chain for chain in chains}

    def find(chain: str) -> str:
        while parent[chain] != chain:
            parent[chain] = parent[parent[chain]]
            chain = parent[chain]
        return chain

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for index, left in enumerate(chains):
        for right in chains[index + 1 :]:
            if legacy.almostIdentical(sequences[left], sequences[right]):
                union(left, right)

    components: list[list[str]] = []
    component_by_root: dict[str, list[str]] = {}
    for chain in chains:
        component = component_by_root.setdefault(find(chain), [])
        if not component:
            components.append(component)
        component.append(chain)

    permutation_groups = [list(itertools.permutations(group)) for group in components]
    orders: list[str] = []
    for selected in itertools.product(*permutation_groups):
        replacement = {
            original: replacement
            for group, permutation in zip(components, selected)
            for original, replacement in zip(group, permutation)
        }
        order = "".join(replacement[chain] for chain in chains)
        if order not in orders:
            orders.append(order)
    return orders


@lru_cache(maxsize=None)
def _aligned_chain_pair(
    reference_file: str,
    reference_chain: str,
    model_file: str,
    model_chain: str,
):
    """Cache the established single-chain residue pairing.

    ``legacy.getCommonResiduesList`` is chain-wise even when called with
    multi-chain strings. Caching its single-chain result removes repeated PDB
    parsing and sequence alignment without changing gap-removal semantics.
    """

    return legacy.getCommonResiduesList(
        reference_file, reference_chain, model_file, model_chain
    )


def _pair_interface_indices(reference_receptor_residues, reference_ligand_residues):
    """Return interface residue indices for one reference chain pair."""

    inter_dis = 10.0
    safe_dis = 20.0
    receptor_interface = set()
    ligand_interface = set()

    for receptor_index, reference_receptor in enumerate(reference_receptor_residues):
        for ligand_index, reference_ligand in enumerate(reference_ligand_residues):
            distance = reference_ligand["CA"] - reference_receptor["CA"]
            if distance > safe_dis:
                continue
            if distance < inter_dis:
                receptor_interface.add(receptor_index)
                ligand_interface.add(ligand_index)
                continue

            found_interface = False
            for receptor_atom in reference_receptor:
                if found_interface:
                    break
                for ligand_atom in reference_ligand:
                    if receptor_atom - ligand_atom < inter_dis:
                        receptor_interface.add(receptor_index)
                        ligand_interface.add(ligand_index)
                        found_interface = True
                        break

    return tuple(sorted(receptor_interface)), tuple(sorted(ligand_interface))


@lru_cache(maxsize=None)
def _paired_interface_indices(
    reference_file: str,
    reference_receptor_chain: str,
    model_file: str,
    model_receptor_chain: str,
    reference_ligand_chain: str,
    model_ligand_chain: str,
):
    """Cache interface masks for one reference/model chain-pair mapping."""

    reference_receptor, _ = _aligned_chain_pair(
        reference_file, reference_receptor_chain, model_file, model_receptor_chain
    )
    reference_ligand, _ = _aligned_chain_pair(
        reference_file, reference_ligand_chain, model_file, model_ligand_chain
    )
    return _pair_interface_indices(reference_receptor, reference_ligand)


def _interface_for_orders(
    reference_file: str,
    reference_receptor: str,
    model_file: str,
    model_receptor: str,
    reference_ligand: str,
    model_ligand: str,
):
    """Build interface coordinates from cached chain-pair masks."""

    receptor_blocks = [
        (
            reference_chain,
            model_chain,
            *_aligned_chain_pair(reference_file, reference_chain, model_file, model_chain),
        )
        for reference_chain, model_chain in zip(reference_receptor, model_receptor)
    ]
    ligand_blocks = [
        (
            reference_chain,
            model_chain,
            *_aligned_chain_pair(reference_file, reference_chain, model_file, model_chain),
        )
        for reference_chain, model_chain in zip(reference_ligand, model_ligand)
    ]
    receptor_indices = [set() for _ in receptor_blocks]
    ligand_indices = [set() for _ in ligand_blocks]

    for receptor_block_index, (reference_receptor_chain, model_receptor_chain, _, _) in enumerate(receptor_blocks):
        for ligand_block_index, (reference_ligand_chain, model_ligand_chain, _, _) in enumerate(ligand_blocks):
            pair_receptor_indices, pair_ligand_indices = _paired_interface_indices(
                reference_file,
                reference_receptor_chain,
                model_file,
                model_receptor_chain,
                reference_ligand_chain,
                model_ligand_chain,
            )
            receptor_indices[receptor_block_index].update(pair_receptor_indices)
            ligand_indices[ligand_block_index].update(pair_ligand_indices)

    reference_receptor_interface = {}
    reference_ligand_interface = {}
    model_receptor_interface = {}
    model_ligand_interface = {}
    offset = 0
    for (_, _, reference_residues, model_residues), indices in zip(receptor_blocks, receptor_indices):
        for index in indices:
            reference_receptor_interface[offset + index] = reference_residues[index]["CA"].get_coord()
            model_receptor_interface[offset + index] = model_residues[index]["CA"].get_coord()
        offset += len(reference_residues)
    offset = 0
    for (_, _, reference_residues, model_residues), indices in zip(ligand_blocks, ligand_indices):
        for index in indices:
            reference_ligand_interface[offset + index] = reference_residues[index]["CA"].get_coord()
            model_ligand_interface[offset + index] = model_residues[index]["CA"].get_coord()
        offset += len(reference_residues)

    return (
        legacy.combineCoords(reference_receptor_interface, reference_ligand_interface),
        legacy.combineCoords(model_receptor_interface, model_ligand_interface),
    )


def grouped_irmsd(
    model_file: str | Path,
    model_receptor: str,
    model_ligand: str,
    reference_file: str | Path,
    reference_receptor: str,
    reference_ligand: str,
) -> float:
    """Return the least legacy iRMSD over valid partner symmetry orders."""

    receptor_orders = symmetric_chain_orders(model_file, model_receptor)
    ligand_orders = symmetric_chain_orders(model_file, model_ligand)
    least = 100000.0
    for receptor_order in receptor_orders:
        for ligand_order in ligand_orders:
            ref_coords, model_coords = _interface_for_orders(
                str(reference_file),
                reference_receptor,
                str(model_file),
                receptor_order,
                reference_ligand,
                ligand_order,
            )
            current = 100.0 if len(ref_coords) == 0 or len(model_coords) == 0 else legacy.getRMSD(model_coords, ref_coords)
            least = min(least, current)
    return least


def main() -> int:
    if len(sys.argv) != 7:
        raise SystemExit(
            "usage: irmsd_grouped_safe.py MODEL MODEL_R MODEL_L NATIVE NATIVE_R NATIVE_L"
        )
    warnings.filterwarnings("ignore")
    value = grouped_irmsd(*sys.argv[1:])
    print(f"{value:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
