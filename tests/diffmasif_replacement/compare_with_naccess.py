import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, Iterable, Optional, Tuple

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from export_residue_accessibility import STANDARD_ASA, build_proxy_accessibility


def parse_rsa_file(rsa_path: str, allowed_chains: Iterable[str]) -> Dict[str, float]:
    allowed = set(allowed_chains)
    relative_asa: Dict[str, float] = {}
    with open(rsa_path, "r") as handle:
        for line in handle:
            if not line.startswith("RES"):
                continue
            items = line.split()
            if len(items) < 5:
                continue
            resname = items[1].upper()
            chain = items[2]
            if chain not in allowed:
                continue
            resseq = items[3]
            absolute_asa = float(items[4])
            standard = STANDARD_ASA.get(resname)
            if not standard:
                continue
            relative_asa[f"{resname}_{resseq}_{chain}"] = absolute_asa * 100.0 / standard
    return relative_asa


def run_naccess_for_structure(pdb_path: str, structure_id: str, allowed_chains: Iterable[str]) -> Dict[str, float]:
    naccess_exec = os.path.abspath("external_tools/naccess/naccess")
    if not os.path.exists(naccess_exec):
        raise FileNotFoundError(f"Naccess executable not found: {naccess_exec}")

    with tempfile.TemporaryDirectory(prefix="diffmasif_naccess_", dir=".") as tmpdir:
        local_pdb = os.path.join(tmpdir, os.path.basename(pdb_path))
        shutil.copy2(pdb_path, local_pdb)
        result = subprocess.run(
            [naccess_exec, os.path.basename(local_pdb)],
            cwd=tmpdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            check=False,
        )
        rsa_path = os.path.join(tmpdir, f"{structure_id[:4].lower()}.rsa")
        if result.returncode != 0 or not os.path.exists(rsa_path):
            raise RuntimeError(
                f"Naccess failed for {structure_id} (exit code {result.returncode}).\n{result.stdout}"
            )
        return parse_rsa_file(rsa_path, allowed_chains=allowed_chains)


def resolve_structure(structure_kind: str, structure_id: str) -> Tuple[str, Iterable[str]]:
    if structure_kind == "target":
        if len(structure_id) < 5:
            raise ValueError("Target ids must include a chain, for example `1fgnA`.")
        return f"processed/pdbs/{structure_id[:4].lower()}.pdb", [structure_id[4]]
    if len(structure_id) < 6:
        raise ValueError("Template ids must include two chains, for example `1a28AB`.")
    return f"templates/pdbs/{structure_id[:4].lower()}.pdb", [structure_id[4], structure_id[5]]


def summarize_comparison(
    naccess: Dict[str, float],
    proxy_stats: Dict[str, Dict[str, float]],
) -> Dict[str, object]:
    proxy = {
        residue: stats["relative_accessibility_proxy"]
        for residue, stats in proxy_stats.items()
        if residue != "_meta"
    }
    common = sorted(set(naccess) & set(proxy))
    naccess_only = sorted(set(naccess) - set(proxy))
    proxy_only = sorted(set(proxy) - set(naccess))
    absolute_errors = [abs(naccess[key] - proxy[key]) for key in common]
    mae = sum(absolute_errors) / len(absolute_errors) if absolute_errors else None

    buried_match_total = 0
    buried_match_hits = 0
    surface_match_total = 0
    surface_match_hits = 0
    for key in common:
        buried_match_total += 1
        if (naccess[key] <= 20.0) == (proxy[key] <= 20.0):
            buried_match_hits += 1
        surface_match_total += 1
        if (naccess[key] > 15.0) == (proxy[key] > 15.0):
            surface_match_hits += 1

    return {
        "common_residue_count": len(common),
        "naccess_only_residue_count": len(naccess_only),
        "proxy_only_residue_count": len(proxy_only),
        "mean_absolute_error": mae,
        "buried_threshold_agreement_at_20": (
            buried_match_hits / buried_match_total if buried_match_total else None
        ),
        "surface_threshold_agreement_at_15": (
            surface_match_hits / surface_match_total if surface_match_total else None
        ),
        "naccess_only_residues": naccess_only[:20],
        "proxy_only_residues": proxy_only[:20],
    }


def build_per_residue_table(
    naccess: Dict[str, float],
    proxy_stats: Dict[str, Dict[str, float]],
) -> Dict[str, Dict[str, Optional[float]]]:
    residue_keys = sorted((set(naccess) | set(proxy_stats)) - {"_meta"})
    table: Dict[str, Dict[str, Optional[float]]] = {}
    for residue in residue_keys:
        proxy_entry = proxy_stats.get(residue, {})
        proxy_value = proxy_entry.get("relative_accessibility_proxy")
        naccess_value = naccess.get(residue)
        table[residue] = {
            "naccess_relative_asa": naccess_value,
            "diffmasif_proxy_relative_asa": proxy_value,
            "absolute_difference": (
                abs(naccess_value - proxy_value)
                if naccess_value is not None and proxy_value is not None
                else None
            ),
            "point_count": proxy_entry.get("point_count"),
            "point_density_per_standard_asa": proxy_entry.get("point_density_per_standard_asa"),
            "mean_score": proxy_entry.get("mean_score"),
        }
    return table


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure-kind", required=True, choices=["target", "template"])
    parser.add_argument("--structure-id", required=True)
    parser.add_argument("--surface-points", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--assignment-cutoff", type=float, default=4.0)
    parser.add_argument("--score-column", default="score")
    args = parser.parse_args()

    pdb_path, chains = resolve_structure(args.structure_kind, args.structure_id)
    naccess = run_naccess_for_structure(pdb_path, args.structure_id, allowed_chains=chains)
    proxy_stats = build_proxy_accessibility(
        pdb_path=pdb_path,
        surface_points_path=args.surface_points,
        allowed_chains=chains,
        assignment_cutoff=args.assignment_cutoff,
        score_column=args.score_column,
    )
    proxy_relative_asa = {
        residue: stats["relative_accessibility_proxy"]
        for residue, stats in proxy_stats.items()
        if residue != "_meta"
    }

    result = {
        "structure_kind": args.structure_kind,
        "structure_id": args.structure_id,
        "pdb_path": pdb_path,
        "chains": list(chains),
        "surface_points_path": args.surface_points,
        "naccess_relative_asa": naccess,
        "diffmasif_proxy_relative_asa": proxy_relative_asa,
        "per_residue": build_per_residue_table(naccess, proxy_stats),
        "summary": summarize_comparison(naccess, proxy_stats),
        "proxy_meta": proxy_stats["_meta"],
    }
    ensure_parent_dir(args.output)
    with open(args.output, "w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
