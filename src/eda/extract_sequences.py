"""Extract sequences from template PDBs using Biopython (ATOM records).

Requires the conda env ``new_bg`` (or any env with Biopython installed)::

    conda activate new_bg
    python extract_sequences.py

Uses ``Bio.SeqIO`` with format ``pdb-atom`` (coordinates, not SEQRES). For each
PDB, every protein chain becomes one entry ``{chain_id: sequence}``.

Reads ``templates/templates_analysis_sorted_unique.csv`` for template ids,
loads ``templates/pdbs/<pdb_id>.pdb`` once per structure (cached), and writes
``templates/template_sequences.json``::

    {
        "104lAB": {"A": "...", "B": "..."},
        ...
    }

Each value dict contains **all** chains present in that PDB file. Parsing many
large structures can take tens of minutes; run under ``conda activate new_bg``.
"""

import csv
import json
import os
import sys
import warnings

try:
    from Bio import BiopythonWarning
    from Bio import SeqIO
except ImportError as e:
    raise SystemExit(
        "Biopython is required. Activate your env first, e.g.\n"
        "  conda activate new_bg\n"
        f"Import error: {e}"
    ) from e

CSV_PATH = "templates/templates_analysis_sorted_unique.csv"
PDB_DIR = "templates/pdbs"
OUTPUT_PATH = "templates/template_sequences.json"


def parse_template(template):
    """Split a template id like ``1a0dAB`` into (pdb_id, chain1, chain2)."""
    if len(template) < 6:
        raise ValueError(f"template id too short: {template!r}")
    return template[:4], template[4], template[5]


def read_all_chains_from_pdb(pdb_path):
    """Return {chain_id: one-letter sequence} for all chains (model 0, pdb-atom)."""
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        for record in SeqIO.parse(pdb_path, "pdb-atom"):
            cid = record.annotations.get("chain")
            if cid is None:
                continue
            seq = str(record.seq).replace("*", "X")
            if seq:
                out[cid] = seq
    return out


def build_sequence_dict(csv_path, pdb_dir):
    sequences = {}
    missing_pdb = []
    missing_pair = []
    pdb_cache = {}

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            template = row["template"].strip()
            try:
                pdb_id, c1, c2 = parse_template(template)
            except ValueError as err:
                print(f"  skip: {err}", file=sys.stderr)
                continue

            pdb_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
            if not os.path.exists(pdb_path):
                missing_pdb.append(template)
                continue

            if pdb_id not in pdb_cache:
                pdb_cache[pdb_id] = read_all_chains_from_pdb(pdb_path)

            chain_seqs = pdb_cache[pdb_id]
            if c1 not in chain_seqs or c2 not in chain_seqs:
                missing_pair.append(template)
            sequences[template] = dict(chain_seqs)

    return sequences, missing_pdb, missing_pair, len(pdb_cache)


def main():
    base = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(base, CSV_PATH)
    pdb_dir = os.path.join(base, PDB_DIR)
    output_path = os.path.join(base, OUTPUT_PATH)

    sequences, missing_pdb, missing_pair, n_pdbs = build_sequence_dict(
        csv_path, pdb_dir
    )

    with open(output_path, "w") as f:
        json.dump(sequences, f, indent=2, sort_keys=True)

    print(f"templates written:        {len(sequences)}", flush=True)
    print(f"unique PDBs parsed:       {n_pdbs}", flush=True)
    print(f"missing pdb files:        {len(missing_pdb)}", flush=True)
    print(f"template pair not in pdb: {len(missing_pair)}", flush=True)
    print(f"saved -> {output_path}", flush=True)

    if missing_pdb[:5]:
        print(f"  sample missing pdb: {missing_pdb[:5]}")
    if missing_pair[:5]:
        print(f"  sample pair missing: {missing_pair[:5]}")


if __name__ == "__main__":
    main()
