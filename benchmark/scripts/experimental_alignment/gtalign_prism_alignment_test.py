"""
Standalone GTalign-based prototype for PRISM's alignment stage.

This script does NOT modify the PRISM pipeline. It reads the same PRISM inputs:
- queries from processed/surface_extraction/*.asa.pdb
- template interfaces from templates/interfaces/*_int.pdb

It runs GTalign in batch mode and writes PRISM-compatible alignment JSON files
to a separate output directory (default: processed/alignment_gtalign_test).
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio.PDB import PDBParser


def load_queries_from_inputs_csv(inputs_csv="inputs.csv"):
    # Matches PRISM's expectation of 5-char IDs (e.g., 1abcA) but only reads file.
    import csv

    ids = []
    with open(inputs_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for col in ("Receptor", "Ligand"):
                val = (row.get(col) or "").strip()
                if len(val) == 5:
                    ids.append(val)
    return sorted(set(ids))


def load_templates(generate_templates=False):
    path = "templates/calculated_templates.txt"
    if generate_templates:
        # The script is test-like; template generation is expected to be run through PRISM.
        raise RuntimeError("Run PRISM template generation first, then use this script on existing artifacts.")
    with open(path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def build_match_dict_from_aligned_sequences(query_seq, ref_seq, protein_path, interface_path):
    if not query_seq or not ref_seq:
        return {}

    seq1_res_ids, seq1_chain_ids = extract_chain_and_res_ids("protein", protein_path)
    seq2_res_ids, seq2_chain_ids = extract_chain_and_res_ids("interface", interface_path)

    match_dict = {}
    index1 = 0
    index2 = 0
    for q_char, r_char in zip(query_seq, ref_seq):
        q_gap = q_char == "-"
        r_gap = r_char == "-"
        if not q_gap and not r_gap and index1 < len(seq1_res_ids) and index2 < len(seq2_res_ids):
            seq1_str = seq1_chain_ids[index1] + "." + q_char + "." + seq1_res_ids[index1]
            seq2_str = seq2_chain_ids[index2] + "." + r_char + "." + seq2_res_ids[index2]
            match_dict[seq2_str] = seq1_str
        if not q_gap:
            index1 += 1
        if not r_gap:
            index2 += 1
    return match_dict


def extract_chain_and_res_ids(name, path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(name, path)
    residue_ids = []
    chain_ids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    residue_ids.append(str(residue.id[1]))
                    chain_ids.append(chain.id)
    return residue_ids, chain_ids


def write_alignment_json(out_dir, protein, template, chain, match_count, translation, rotation_mat, match_dict, tm_score):
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "match_count": int(match_count) if match_count is not None else 0,
        "translation": translation or [0.0, 0.0, 0.0],
        "rotation_mat": rotation_mat or [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "match_dict": match_dict or {},
        "tm_score": float(tm_score) if tm_score is not None else 0.0,
    }
    with open(out_dir / f"{protein}_{template}_{chain}.json", "w") as f:
        json.dump(payload, f)


def _symlink_or_copy(src, dst):
    try:
        os.symlink(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def _extract_floats(line):
    return [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", line)]


def _extract_gtalign_alignment_seq(line, label):
    pattern = r"^\s*" + re.escape(label) + r":\s+\d+\s+([A-Za-z\-]+)\s+\d+\s*$"
    m = re.match(pattern, line)
    if m:
        return m.group(1)
    toks = line.strip().split()
    if len(toks) >= 4:
        return toks[-2]
    return None


def _extract_gtalign_query_path(lines):
    for i, line in enumerate(lines):
        if line.startswith(" Query ("):
            j = i + 1
            while j < len(lines) and not lines[j].startswith(" Searched:"):
                candidate = lines[j].strip()
                if candidate and not candidate.startswith("Chn:"):
                    return candidate.split(" Chn:", 1)[0].strip()
                j += 1
    return None


def _parse_gtalign_hit_block(lines, start_idx):
    i = start_idx
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines) or not lines[i].lstrip().startswith(">"):
        return None, i + 1

    ref_path = lines[i].strip()[1:].split(" Chn:", 1)[0].strip()
    i += 1

    hit = {
        "ref_path": ref_path,
        "tm_ref": None,
        "tm_query": None,
        "rmsd": None,
        "aligned_length": None,
        "rotation_mat": None,
        "translation": None,
        "query_aln": "",
        "ref_aln": "",
    }
    query_chunks = []
    ref_chunks = []

    while i < len(lines):
        line = lines[i]
        if re.match(r"^\s*\d+\.\s*$", line) or line.startswith("Query length:"):
            break

        if "TM-score (Refn./Query)" in line:
            m = re.search(
                r"TM-score \(Refn\./Query\)\s*=\s*([0-9.]+)\s*/\s*([0-9.]+).*RMSD\s*=\s*([0-9.]+)",
                line,
            )
            if m:
                hit["tm_ref"] = float(m.group(1))
                hit["tm_query"] = float(m.group(2))
                hit["rmsd"] = float(m.group(3))

        if "Matched =" in line:
            m = re.search(r"Matched\s*=\s*(\d+)/", line)
            if m:
                hit["aligned_length"] = int(m.group(1))

        if line.lstrip().startswith("Query:") and i + 2 < len(lines) and lines[i + 2].lstrip().startswith("Refn.:"):
            q = _extract_gtalign_alignment_seq(line, "Query")
            r = _extract_gtalign_alignment_seq(lines[i + 2], "Refn.")
            if q is not None and r is not None:
                query_chunks.append(q)
                ref_chunks.append(r)
                i += 3
                continue

        if line.strip().startswith("Rotation [3,3] and translation [3,1] for Query:"):
            rot = []
            trans = []
            ok = True
            for j in range(1, 4):
                if i + j >= len(lines):
                    ok = False
                    break
                vals = _extract_floats(lines[i + j])
                if len(vals) < 4:
                    ok = False
                    break
                rot.append(vals[:3])
                trans.append(vals[3])
            if ok:
                hit["rotation_mat"] = rot
                hit["translation"] = trans
            i += 4
            continue

        i += 1

    hit["query_aln"] = "".join(query_chunks)
    hit["ref_aln"] = "".join(ref_chunks)
    if hit["aligned_length"] is None and hit["query_aln"] and hit["ref_aln"]:
        hit["aligned_length"] = sum(1 for a, b in zip(hit["query_aln"], hit["ref_aln"]) if a != "-" and b != "-")
    return hit, i


def parse_gtalign_output_text(text):
    lines = text.splitlines()
    query_path = _extract_gtalign_query_path(lines)
    if query_path is None:
        raise ValueError("Could not parse GTalign query path")
    hits = []
    i = 0
    while i < len(lines):
        if re.match(r"^\s*\d+\.\s*$", lines[i]):
            hit, i = _parse_gtalign_hit_block(lines, i + 1)
            if hit:
                hits.append(hit)
            continue
        i += 1
    return query_path, hits


def run_gtalign_stage(
    gtalign_path,
    queries,
    templates,
    output_dir,
    dev_min_length=3,
    pre_score=0.0,
    speed=None,
    refinement=None,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    selected_query_paths = {}
    selected_ref_paths = {}
    for protein in queries:
        qpath = Path(f"processed/surface_extraction/{protein}.asa.pdb")
        if qpath.exists():
            selected_query_paths[protein] = qpath.resolve()
        else:
            print(f"Missing query surface file: {qpath}")
    for template in templates:
        for chain in template[4:]:
            rpath = Path(f"templates/interfaces/{template}_{chain}_int.pdb")
            if rpath.exists():
                selected_ref_paths[(template, chain)] = rpath.resolve()
            else:
                print(f"Missing template interface file: {rpath}")

    if not selected_query_paths or not selected_ref_paths:
        raise RuntimeError("No valid query/reference files found for GTalign test stage.")

    parsed_pairs = set()
    with tempfile.TemporaryDirectory(prefix="gtalign_prism_test_", dir=str(output_dir.parent)) as td:
        td = Path(td)
        qdir = td / "queries"
        rdir = td / "refs"
        outdir = td / "out"
        qdir.mkdir()
        rdir.mkdir()
        outdir.mkdir()

        query_basename_to_real = {}
        ref_basename_to_real = {}
        ref_basename_to_key = {}

        for _, src in selected_query_paths.items():
            dst = qdir / src.name
            _symlink_or_copy(str(src), str(dst))
            query_basename_to_real[dst.name] = str(src)

        for key, src in selected_ref_paths.items():
            dst = rdir / src.name
            _symlink_or_copy(str(src), str(dst))
            ref_basename_to_real[dst.name] = str(src)
            ref_basename_to_key[dst.name] = key

        cmd = [
            gtalign_path,
            f"--qrs={qdir}",
            f"--rfs={rdir}",
            "-o",
            str(outdir),
            "-s",
            "0",
            f"--nhits={len(selected_ref_paths)}",
            f"--nalns={len(selected_ref_paths)}",
            f"--dev-min-length={int(dev_min_length)}",
            f"--pre-score={float(pre_score)}",
        ]
        if speed is not None:
            cmd.append(f"--speed={int(speed)}")
        if refinement is not None:
            cmd.append(f"--refinement={int(refinement)}")
        print("Running:", " ".join(cmd))
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode != 0:
            raise RuntimeError(f"GTalign failed ({result.returncode}): {(result.stderr or result.stdout)[:2000]}")

        for out_file in sorted(outdir.glob("*.out")):
            query_path, hits = parse_gtalign_output_text(out_file.read_text())
            qbase = os.path.basename(query_path)
            if not qbase.endswith(".asa.pdb"):
                raise ValueError(f"Unexpected query filename: {qbase}")
            protein = qbase[:-8]
            protein_path = query_basename_to_real[qbase]

            for hit in hits:
                ref_base = os.path.basename(hit["ref_path"])
                if ref_base not in ref_basename_to_key:
                    continue
                template, chain = ref_basename_to_key[ref_base]
                interface_path = ref_basename_to_real[ref_base]
                match_dict = build_match_dict_from_aligned_sequences(
                    hit["query_aln"],
                    hit["ref_aln"],
                    protein_path,
                    interface_path,
                )
                tm_candidates = [v for v in (hit["tm_ref"], hit["tm_query"]) if isinstance(v, (float, int))]
                tm_score = max(tm_candidates) if tm_candidates else 0.0
                match_count = hit["aligned_length"] or len(match_dict)
                write_alignment_json(
                    output_dir,
                    protein,
                    template,
                    chain,
                    match_count,
                    hit["translation"],
                    hit["rotation_mat"],
                    match_dict,
                    tm_score,
                )
                parsed_pairs.add((protein, template, chain))

    # Write empty JSONs for missing pairs to mimic PRISM stage behavior
    for protein in queries:
        for template in templates:
            for chain in template[4:]:
                if (protein, template, chain) not in parsed_pairs:
                    write_alignment_json(
                        output_dir,
                        protein,
                        template,
                        chain,
                        0,
                        [0.0, 0.0, 0.0],
                        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                        {},
                        0.0,
                    )

    print(f"Wrote {len(list(output_dir.glob('*.json')))} JSON files to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Standalone GTalign prototype for PRISM alignment stage")
    parser.add_argument("--gtalign_path", required=True, help="Path to GTalign binary")
    parser.add_argument(
        "--output_dir",
        default="processed/alignment_gtalign_test",
        help="Output directory for PRISM-format alignment JSON files",
    )
    parser.add_argument(
        "--inputs_csv",
        default="inputs.csv",
        help="CSV file with Receptor/Ligand columns (default: inputs.csv)",
    )
    parser.add_argument(
        "--queries",
        nargs="*",
        default=None,
        help="Optional explicit list of query IDs (default: derive from inputs.csv)",
    )
    parser.add_argument(
        "--templates",
        nargs="*",
        default=None,
        help="Optional explicit list of template IDs (default: templates/calculated_templates.txt)",
    )
    parser.add_argument(
        "--dev_min_length",
        type=int,
        default=3,
        help="GTalign --dev-min-length value (default: 3; useful for short PRISM interface fragments)",
    )
    parser.add_argument(
        "--pre_score",
        type=float,
        default=0.0,
        help="GTalign --pre-score value (default: 0.0 to reduce pruning of weak alignments)",
    )
    parser.add_argument(
        "--speed",
        type=int,
        default=None,
        help="Optional GTalign --speed override (lower is slower but potentially more sensitive)",
    )
    parser.add_argument(
        "--refinement",
        type=int,
        default=None,
        help="Optional GTalign --refinement override (0-3)",
    )
    args = parser.parse_args()

    queries = args.queries if args.queries else load_queries_from_inputs_csv(args.inputs_csv)
    templates = args.templates if args.templates else load_templates()
    run_gtalign_stage(
        args.gtalign_path,
        queries,
        templates,
        args.output_dir,
        dev_min_length=args.dev_min_length,
        pre_score=args.pre_score,
        speed=args.speed,
        refinement=args.refinement,
    )


if __name__ == "__main__":
    main()
