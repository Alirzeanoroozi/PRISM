#!/usr/bin/env python3
"""
Score extension protein-DNA predictions against a native complex using chain-qualified contacts.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa


DNA_RESIDUES = {"DA", "DC", "DG", "DT", "DI", "A", "C", "G", "T", "I"}
CONTACT_DISTANCE = 4.5


def is_dna_residue(residue):
    return residue.get_resname().strip().upper() in DNA_RESIDUES


def parse_chain_list(value):
    if value in ("", None, "NA"):
        return []
    if isinstance(value, list):
        return value
    return [item.strip() for item in str(value).split(",") if item.strip()]


def infer_chain_groups(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("infer_groups_ext", pdb_path)
    protein = []
    dna = []
    for model in structure:
        for chain in model:
            aa_count = 0
            dna_count = 0
            for residue in chain:
                if is_aa(residue, standard=True):
                    aa_count += 1
                elif is_dna_residue(residue):
                    dna_count += 1
            if aa_count:
                protein.append(chain.id)
            if dna_count:
                dna.append(chain.id)
        break
    return protein, dna


def collect_atoms(pdb_path, chain_ids, predicate):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("collect_atoms_ext", pdb_path)
    grouped = {}
    for model in structure:
        for chain in model:
            if chain_ids and chain.id not in chain_ids:
                continue
            for residue in chain:
                if not predicate(residue):
                    continue
                key = (chain.id, residue.id[1], residue.get_resname().strip())
                grouped.setdefault(key, [])
                for atom in residue:
                    if atom.element == "H":
                        continue
                    grouped[key].append(np.array(atom.get_coord()))
        break
    return grouped


def collect_ordered_residue_ranks(pdb_path, chain_ids, predicate):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("collect_ranks_ext", pdb_path)
    ranks = {}
    for model in structure:
        for chain in model:
            if chain_ids and chain.id not in chain_ids:
                continue
            rank = 0
            for residue in chain:
                if not predicate(residue):
                    continue
                rank += 1
                key = (chain.id, residue.id[1], residue.get_resname().strip())
                ranks[key] = rank
        break
    return ranks


def compute_contacts(pdb_path, protein_chains, dna_chains):
    protein_atoms = collect_atoms(pdb_path, protein_chains, lambda residue: is_aa(residue, standard=True))
    dna_atoms = collect_atoms(pdb_path, dna_chains, is_dna_residue)
    dna_ranks = collect_ordered_residue_ranks(pdb_path, dna_chains, is_dna_residue)

    contacts = []
    protein_residues = set()
    dna_residues = set()
    for (p_chain, p_resnum, p_resname), p_coords in protein_atoms.items():
        arr1 = np.array(p_coords)
        for (d_chain, d_resnum, d_resname), d_coords in dna_atoms.items():
            arr2 = np.array(d_coords)
            dists = np.linalg.norm(arr1[:, np.newaxis, :] - arr2[np.newaxis, :, :], axis=2)
            min_dist = float(np.min(dists))
            if min_dist <= CONTACT_DISTANCE:
                protein_key = f"{p_chain}.{p_resnum}.{p_resname}"
                dna_key = f"{d_chain}.{d_resnum}.{d_resname}"
                protein_residues.add(protein_key)
                dna_residues.add(dna_key)
                contacts.append(
                    {
                        "protein_key": protein_key,
                        "dna_key": dna_key,
                        "dna_rank": dna_ranks.get((d_chain, d_resnum, d_resname)),
                        "protein_chain": p_chain,
                        "protein_residue_number": p_resnum,
                        "protein_residue_name": p_resname,
                        "dna_chain": d_chain,
                        "dna_residue_number": d_resnum,
                        "dna_residue_name": d_resname,
                        "min_distance": min_dist,
                    }
                )
    return contacts, protein_residues, dna_residues


def compute_protein_alignment(model_pdb, native_pdb, model_protein_chains, native_protein_chains):
    parser = PDBParser(QUIET=True)
    model_structure = parser.get_structure("model_align_ext", model_pdb)
    native_structure = parser.get_structure("native_align_ext", native_pdb)

    def ordered_ca_atoms(structure, chain_ids):
        atoms = []
        for model in structure:
            for chain in model:
                if chain_ids and chain.id not in chain_ids:
                    continue
                for residue in chain:
                    if is_aa(residue, standard=True) and "CA" in residue:
                        atoms.append(residue["CA"])
            break
        return atoms

    model_atoms = ordered_ca_atoms(model_structure, model_protein_chains)
    native_atoms = ordered_ca_atoms(native_structure, native_protein_chains)
    n = min(len(model_atoms), len(native_atoms))
    if n == 0:
        return None, None

    sup = Superimposer()
    sup.set_atoms(native_atoms[:n], model_atoms[:n])
    rmsd = float(sup.rms)
    tm_like = 1.0 / (1.0 + (rmsd / 5.0) ** 2)
    return tm_like, rmsd


def load_rosetta_score(sidecar_path):
    if not sidecar_path:
        return None
    path = Path(sidecar_path)
    if not path.exists():
        return None
    with open(path, "r") as handle:
        payload = json.load(handle)
    return payload.get("rosetta_dna_interface_score")


def precision_recall_f1(overlap, predicted_size, native_size):
    precision = overlap / float(predicted_size) if predicted_size else 0.0
    recall = overlap / float(native_size) if native_size else 0.0
    f1 = 0.0
    if precision + recall:
        f1 = 2.0 * precision * recall / (precision + recall)
    return precision, recall, f1


def score_one_ext(
    model_pdb,
    native_pdb,
    model_protein_chains=None,
    model_dna_chains=None,
    native_protein_chains=None,
    native_dna_chains=None,
    score_json=None,
):
    if model_protein_chains is None or model_dna_chains is None:
        protein, dna = infer_chain_groups(model_pdb)
        model_protein_chains = model_protein_chains or protein
        model_dna_chains = model_dna_chains or dna
    if native_protein_chains is None or native_dna_chains is None:
        protein, dna = infer_chain_groups(native_pdb)
        native_protein_chains = native_protein_chains or protein
        native_dna_chains = native_dna_chains or dna

    model_contacts, model_protein_if, model_dna_if = compute_contacts(model_pdb, model_protein_chains, model_dna_chains)
    native_contacts, native_protein_if, native_dna_if = compute_contacts(native_pdb, native_protein_chains, native_dna_chains)

    model_pairs = {(contact["protein_key"], contact["dna_key"]) for contact in model_contacts}
    native_pairs = {(contact["protein_key"], contact["dna_key"]) for contact in native_contacts}
    model_register_pairs = {
        (contact["protein_key"], contact["dna_chain"], contact["dna_rank"])
        for contact in model_contacts
        if contact["dna_rank"] is not None
    }
    native_register_pairs = {
        (contact["protein_key"], contact["dna_chain"], contact["dna_rank"])
        for contact in native_contacts
        if contact["dna_rank"] is not None
    }
    overlap = len(model_pairs.intersection(native_pairs))
    register_overlap = len(model_register_pairs.intersection(native_register_pairs))

    precision, recall, contact_f1 = precision_recall_f1(overlap, len(model_pairs), len(native_pairs))
    register_precision, register_recall, register_f1 = precision_recall_f1(
        register_overlap, len(model_register_pairs), len(native_register_pairs)
    )
    protein_overlap = len(model_protein_if.intersection(native_protein_if))
    dna_overlap = len(model_dna_if.intersection(native_dna_if))
    protein_precision, protein_recall, protein_f1 = precision_recall_f1(protein_overlap, len(model_protein_if), len(native_protein_if))
    dna_precision, dna_recall, dna_f1 = precision_recall_f1(dna_overlap, len(model_dna_if), len(native_dna_if))
    tm_score, rmsd = compute_protein_alignment(model_pdb, native_pdb, model_protein_chains, native_protein_chains)

    return {
        "model_pdb": str(Path(model_pdb).resolve()),
        "native_pdb": str(Path(native_pdb).resolve()),
        "model_protein_chains": ",".join(model_protein_chains),
        "model_dna_chains": ",".join(model_dna_chains),
        "native_protein_chains": ",".join(native_protein_chains),
        "native_dna_chains": ",".join(native_dna_chains),
        "contact_precision": precision,
        "contact_recall": recall,
        "contact_f1": contact_f1,
        "dna_register_contact_precision": register_precision,
        "dna_register_contact_recall": register_recall,
        "dna_register_contact_f1": register_f1,
        "protein_interface_precision": protein_precision,
        "protein_interface_recall": protein_recall,
        "protein_interface_f1": protein_f1,
        "nucleotide_contact_precision": dna_precision,
        "nucleotide_contact_recall": dna_recall,
        "nucleotide_contact_f1": dna_f1,
        "alignment_tm_score": tm_score,
        "alignment_rmsd": rmsd,
        "rosetta_dna_interface_score": load_rosetta_score(score_json),
        "status": "passed",
        "reason": "",
    }


def write_csv(rows, output_path):
    fieldnames = list(rows[0].keys()) if rows else []
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Score one extension protein-DNA model or directory against a native complex")
    parser.add_argument("model_input", help="Model PDB file or directory")
    parser.add_argument("native_pdb", help="Native/reference PDB")
    parser.add_argument("--model-protein-chains", default=None)
    parser.add_argument("--model-dna-chains", default=None)
    parser.add_argument("--native-protein-chains", default=None)
    parser.add_argument("--native-dna-chains", default=None)
    parser.add_argument("--score-json", default=None)
    parser.add_argument("--out-csv", default=None)
    args = parser.parse_args()

    model_path = Path(args.model_input)
    models = [model_path] if model_path.is_file() else sorted(model_path.glob("*.pdb"))
    rows = []
    for model_pdb in models:
        rows.append(
            score_one_ext(
                str(model_pdb),
                args.native_pdb,
                model_protein_chains=parse_chain_list(args.model_protein_chains) or None,
                model_dna_chains=parse_chain_list(args.model_dna_chains) or None,
                native_protein_chains=parse_chain_list(args.native_protein_chains) or None,
                native_dna_chains=parse_chain_list(args.native_dna_chains) or None,
                score_json=args.score_json,
            )
        )
    if args.out_csv:
        write_csv(rows, args.out_csv)
    if len(rows) == 1:
        print(json.dumps(rows[0], indent=2))
    else:
        print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
