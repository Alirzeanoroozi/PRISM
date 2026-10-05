#!/usr/bin/env python3
"""Build a deterministic target/template leakage manifest.

This runner performs only sequence/asset provenance checks.  It never selects
templates using native scores and it never authorizes a confirmatory run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import re
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

from benchmark.scripts.build_template_source_gate import (
    SENSITIVITY_THRESHOLDS,
    classify_template,
    confirmatory_row_eligibility,
)


SIMILARITY_FIELDS = (
    "dataset_row_id", "pair_id", "target_role", "target_selector", "target_pdb",
    "target_chain", "template_id", "template_chain", "orientation", "asset_status",
    "identity_percent", "query_coverage_percent", "template_coverage_percent",
    "shorter_coverage_percent", "exclusion_reason", "similarity_eligible",
    "source_gate_status", "source_gate_reason",
) + tuple(f"excluded_identity_gt_{threshold}" for threshold in SENSITIVITY_THRESHOLDS)

ELIGIBLE_FIELDS = (
    "dataset_row_id", "pair_id", "source_gate_status", "source_gate_reason",
    "target_asset_status", "template_list_sha256", "eligible_template_count",
    "eligible_templates", "excluded_template_count", "missing_template_count",
    "similarity_tsv_sha256",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_selector(value: str) -> tuple[str, str]:
    token = (value or "").strip()
    pdb_id, separator, chains = token.partition("_")
    if not separator:
        chains = ""
    chains = re.sub(r"[^A-Za-z0-9]", "", chains)
    if len(pdb_id) < 4:
        raise ValueError(f"invalid selector: {value!r}")
    return pdb_id[:4].lower(), chains


def chain_sequences(path: Path, wanted: str | None = None) -> dict[str, str]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("source_gate", str(path))
    allowed = set(wanted or "")
    result: dict[str, str] = {}
    model = next(structure.get_models())
    for chain in model:
        if allowed and chain.id not in allowed:
            continue
        sequence: list[str] = []
        for residue in chain:
            if "CA" not in residue:
                continue
            sequence.append(protein_letters_3to1.get(residue.resname.upper(), "X"))
        if sequence:
            result[chain.id] = "".join(sequence)
    return result


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_templates(path: Path) -> list[str]:
    values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()
              if line.strip() and not line.lstrip().startswith("#")]
    if len(values) != len(set(values)):
        raise ValueError("template list contains duplicate IDs")
    return values


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def best_assignment(target: dict[str, str], template: dict[str, str]) -> list[tuple[str, str]]:
    """Match target/template chains by maximum total identity deterministically."""

    pairs = list(itertools.product(sorted(target), sorted(template)))
    if not pairs:
        return []
    target_ids, template_ids = sorted(target), sorted(template)
    assignments: list[list[tuple[str, str]]] = []
    for selected_targets in itertools.permutations(target_ids, min(len(target_ids), len(template_ids))):
        for selected_templates in itertools.permutations(template_ids, len(selected_targets)):
            assignments.append(list(zip(selected_targets, selected_templates)))
    def score(assignment: list[tuple[str, str]]) -> tuple[float, tuple[tuple[str, str], ...]]:
        total = 0.0
        for left, right in assignment:
            total += float(classify_template(
                target_pdb="", target_chains=left, template_id="0000" + right,
                target_sequence=target[left], template_sequence=template[right]
            )["identity_percent"])
        return total, tuple(assignment)
    return max(assignments, key=lambda item: (score(item)[0], tuple(reversed(score(item)[1]))))


def build_manifest(
    *, inputs: Path, staged_pdb_dir: Path, template_list: Path, interface_root: Path,
    source_policy: Path, similarity_output: Path, eligible_output: Path,
) -> dict[str, object]:
    decision = json.loads(source_policy.read_text(encoding="utf-8"))["decision"]
    templates = read_templates(template_list)
    template_hash = sha256_file(template_list)
    rows = read_rows(inputs)
    alignment_engine = PairwiseAligner()
    alignment_engine.mode = "global"
    similarity_cache: dict[tuple[str, str], dict[str, float | int]] = {}
    # Parse each immutable template interface exactly once.  Without this
    # cache, the 19,855-template matrix would repeatedly reparse tens of
    # thousands of PDB files for every benchmark row.
    template_cache: dict[tuple[str, str], tuple[str, str]] = {}
    for template_id in templates:
        for chain in template_id[4:].replace("_", ""):
            path = interface_root / f"{template_id}_{chain}_int.pdb"
            try:
                seqs = chain_sequences(path)
                template_cache[(template_id, chain)] = ("complete", next(iter(seqs.values())))
            except (OSError, ValueError, StopIteration):
                template_cache[(template_id, chain)] = ("missing_or_invalid", "")
    similarity_rows: list[dict[str, object]] = []
    eligible_rows: list[dict[str, object]] = []
    for row in rows:
        pair_id = row.get("pair_id", "")
        dataset_row_id = row.get("dataset_row_id", pair_id)
        source_gate_ok, source_gate_reason = confirmatory_row_eligibility(dataset_row_id, decision)
        selectors = [("receptor", row.get("Receptor", "")), ("ligand", row.get("Ligand", ""))]
        target_assets: dict[str, dict[str, str]] = {}
        target_asset_status = "complete"
        for role, selector in selectors:
            try:
                pdb_id, chains = parse_selector(selector)
                target_assets[role] = chain_sequences(staged_pdb_dir / f"{pdb_id}.pdb", chains)
                if set(chains) - set(target_assets[role]):
                    target_asset_status = "missing_chain"
            except (OSError, ValueError, StopIteration):
                target_assets[role] = {}
                target_asset_status = "missing_or_invalid"
        eligible: list[str] = []
        excluded: list[str] = []
        missing: list[str] = []
        for template_id in templates:
            template_chains = template_id[4:].replace("_", "")
            template_assets: dict[str, str] = {}
            asset_status = "complete"
            for chain in template_chains:
                status, sequence = template_cache.get((template_id, chain), ("missing_or_invalid", ""))
                if status != "complete":
                    asset_status = status
                else:
                    template_assets[chain] = sequence
            if asset_status != "complete":
                missing.append(template_id)
                for role, _selector in selectors:
                    similarity_rows.append({"dataset_row_id": dataset_row_id, "pair_id": pair_id, "target_role": role, "target_selector": row.get("Receptor" if role == "receptor" else "Ligand", ""), "target_pdb": "", "target_chain": "", "template_id": template_id, "template_chain": "", "orientation": "", "asset_status": asset_status, "exclusion_reason": "template_asset_missing", "similarity_eligible": False, "source_gate_status": decision.get("status", ""), "source_gate_reason": source_gate_reason})
                continue
            for orientation, mapping in (("o1", list(zip(("receptor", "ligand"), template_chains))), ("o2", list(zip(("receptor", "ligand"), reversed(template_chains))))):
                orientation_reasons: list[str] = []
                orientation_eligible = True
                for role, template_chain in mapping:
                    target_sequences = target_assets.get(role, {})
                    if not target_sequences or template_chain not in template_assets:
                        orientation_eligible = False
                        continue
                    for target_chain, target_sequence in sorted(target_sequences.items()):
                        cache_key = (target_sequence, template_assets[template_chain])
                        if cache_key not in similarity_cache:
                            from benchmark.scripts.build_template_source_gate import sequence_similarity
                            similarity_cache[cache_key] = sequence_similarity(*cache_key, aligner=alignment_engine)
                        decision_row = classify_template(target_pdb=parse_selector(row.get("Receptor" if role == "receptor" else "Ligand", ""))[0], target_chains=target_chain, template_id=template_id, target_sequence=target_sequence, template_sequence=template_assets[template_chain], aligner=alignment_engine, similarity=similarity_cache[cache_key])
                        if decision_row["exclusion_reason"] != "eligible":
                            orientation_eligible = False
                            orientation_reasons.append(str(decision_row["exclusion_reason"]))
                        similarity_rows.append({"dataset_row_id": dataset_row_id, "pair_id": pair_id, "target_role": role, "target_selector": row.get("Receptor" if role == "receptor" else "Ligand", ""), "target_pdb": parse_selector(row.get("Receptor" if role == "receptor" else "Ligand", ""))[0], "target_chain": target_chain, "template_id": template_id, "template_chain": template_chain, "orientation": orientation, "asset_status": asset_status, "similarity_eligible": decision_row["confirmatory_eligible"], "source_gate_status": decision.get("status", ""), "source_gate_reason": source_gate_reason, **{key: decision_row.get(key, "") for key in SIMILARITY_FIELDS if key in decision_row}, **{f"excluded_identity_gt_{threshold}": decision_row.get(f"excluded_identity_gt_{threshold}", "") for threshold in SENSITIVITY_THRESHOLDS}})
                if source_gate_ok and orientation_eligible and not orientation_reasons:
                    eligible.append(template_id)
            if template_id not in eligible and template_id not in missing:
                excluded.append(template_id)
        eligible = sorted(set(eligible))
        eligible_rows.append({"dataset_row_id": dataset_row_id, "pair_id": pair_id, "source_gate_status": decision.get("status", ""), "source_gate_reason": source_gate_reason, "target_asset_status": target_asset_status, "template_list_sha256": template_hash, "eligible_template_count": len(eligible), "eligible_templates": ",".join(eligible), "excluded_template_count": len(excluded), "missing_template_count": len(missing), "similarity_tsv_sha256": ""})
    write_tsv(similarity_output, SIMILARITY_FIELDS, similarity_rows)
    similarity_hash = sha256_file(similarity_output)
    for row in eligible_rows:
        row["similarity_tsv_sha256"] = similarity_hash
    write_tsv(eligible_output, ELIGIBLE_FIELDS, eligible_rows)
    return {"row_count": len(rows), "template_count": len(templates), "similarity_row_count": len(similarity_rows), "similarity_sha256": similarity_hash, "template_list_sha256": template_hash, "source_gate_status": decision.get("status", ""), "confirmatory_run_authorized": bool(decision.get("confirmatory_run_authorized"))}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--staged-pdb-dir", type=Path, required=True)
    parser.add_argument("--template-list", type=Path, required=True)
    parser.add_argument("--interface-root", type=Path, required=True)
    parser.add_argument("--source-policy", type=Path, required=True)
    parser.add_argument("--similarity-output", type=Path, required=True)
    parser.add_argument("--eligible-output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_manifest(inputs=args.inputs, staged_pdb_dir=args.staged_pdb_dir, template_list=args.template_list, interface_root=args.interface_root, source_policy=args.source_policy, similarity_output=args.similarity_output, eligible_output=args.eligible_output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
