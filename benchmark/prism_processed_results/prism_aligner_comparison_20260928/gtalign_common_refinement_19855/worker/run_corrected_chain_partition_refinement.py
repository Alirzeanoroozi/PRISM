#!/usr/bin/env python3
"""Run the common-refinement worker with source/native side partitions decoupled.

The original worker used the native receptor/ligand chain counts to normalize
each transformed source half.  Some valid candidates partition the same total
chains differently between source halves and native interface groups.  This
wrapper preserves source-side chain counts for refinement and restores the
true native mapping for DockQ/iRMSD.  It is intentionally separate from the
active worker so the running array retains its original provenance.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


_WORKER_PATH = Path(__file__).with_name("run_small_refinement_comparison_one.py")
_SPEC = importlib.util.spec_from_file_location("prism_common_refinement_worker", _WORKER_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError(f"cannot load worker module: {_WORKER_PATH}")
_WORKER = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_WORKER)
_ORIGINAL_CANONICALIZE_MODEL = _WORKER.canonicalize_model


def _clean(value: str) -> str:
    return "".join(str(value).split())


def prepare_candidate_inputs(row: dict[str, str], candidate_root: Path) -> dict[str, object]:
    """Normalize each source half using its own chain count.

    The native groups still define the final DockQ/iRMSD scope.  Only the
    source-side split used by FiberDock/Rosetta is preserved here.
    """
    left = Path(row["left"]).resolve()
    right = Path(row["right"]).resolve()
    left_source = _WORKER.chain_ids(left)
    right_source = _WORKER.chain_ids(right)
    native_left = _clean(row["native_receptor_chains"])
    native_right = _clean(row["native_ligand_chains"])
    total = len(left_source) + len(right_source)
    if total != len(native_left) + len(native_right):
        raise ValueError(
            "transformed/native total chain-count mismatch: "
            f"source {left_source}+{right_source}, native {native_left}+{native_right}"
        )
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    if total < 2 or total > len(alphabet):
        raise ValueError(f"unsupported total source chain count: {total}")
    left_chains = alphabet[: len(left_source)]
    right_chains = alphabet[len(left_source) : total]
    input_dir = candidate_root / "inputs"
    basename = f"normalized_srcx{left_chains}_srcy{right_chains}_o1"
    left_out = input_dir / f"{basename}_L.pdb"
    right_out = input_dir / f"{basename}_R.pdb"
    left_info = _WORKER.normalize_partner(left, left_out, left_chains)
    right_info = _WORKER.normalize_partner(right, right_out, right_chains)
    return {
        "prepared_left": str(left_out.resolve()),
        "prepared_right": str(right_out.resolve()),
        "prepared_left_chains": left_chains,
        "prepared_right_chains": right_chains,
        "normalization_left": left_info,
        "normalization_right": right_info,
        "source_chain_groups": {"left": "".join(left_source), "right": "".join(right_source)},
        "native_chain_groups": {"left": native_left, "right": native_right},
        "chain_partition_mode": "source_partition_native_scope",
    }


def canonicalize_model(
    model: Path,
    left: Path,
    right: Path,
    native_receptor: str,
    native_ligand: str,
    output: Path,
) -> dict[str, object]:
    """Canonicalize with source-side model groups and true native mapping."""
    left_count = len(_WORKER.input_chain_specs(left))
    right_count = len(_WORKER.input_chain_specs(right))
    native_receptor = _clean(native_receptor)
    native_ligand = _clean(native_ligand)
    total = left_count + right_count
    if total != len(native_receptor) + len(native_ligand):
        raise ValueError(
            "transformed/native total chain-count mismatch: "
            f"source {left_count}+{right_count}, native {native_receptor}+{native_ligand}"
        )
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    fake_receptor = alphabet[:left_count]
    fake_ligand = alphabet[left_count:total]
    result = _ORIGINAL_CANONICALIZE_MODEL(
        model,
        left,
        right,
        fake_receptor,
        fake_ligand,
        output,
    )
    model_receptor = str(result["model_receptor_chains"])
    model_ligand = str(result["model_ligand_chains"])
    result.update(
        native_receptor_chains=native_receptor,
        native_ligand_chains=native_ligand,
        mapping=f"{model_receptor}{model_ligand}:{native_receptor}{native_ligand}",
        chain_partition_mode="source_partition_native_scope",
    )
    return result


_WORKER.prepare_candidate_inputs = prepare_candidate_inputs
_WORKER.canonicalize_model = canonicalize_model


if __name__ == "__main__":
    raise SystemExit(_WORKER.main())
