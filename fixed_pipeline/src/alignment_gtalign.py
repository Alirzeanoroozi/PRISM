"""
Fixed version of src/alignment_gtalign.py

Key fixes:
1. Missing query surface files AND missing template interface files are
   counted and reported as a warning summary at the end.
2. Added a pre-flight report so operators know coverage was reduced.
3. Added logging module integration.
"""

import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio.PDB import PDBParser

import alignment as _align
extract_chain_and_res_ids = _align.extract_chain_and_res_ids

logger = logging.getLogger(__name__)


def align_gtalign(
    queries,
    templates,
    gtalign_path="gtalign",
    output_dir="processed/alignment_gtalign",
    dev_min_length=3,
    pre_score=0.0,
    speed=0,
    refinement=3,
):
    """
    GTalign-backed alignment stage that writes PRISM-compatible JSONs.
    """
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise RuntimeError(
            f"refusing to reuse non-empty GTalign output directory: {output_dir}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    selected_query_paths = {}
    selected_ref_paths = {}
    missing_queries = []
    missing_templates = []

    for protein in queries:
        qpath = Path(f"processed/surface_extraction/{protein}.asa.pdb")
        if qpath.exists():
            selected_query_paths[protein] = qpath.resolve()
        else:
            missing_queries.append(protein)
            logger.warning("Missing query surface file: %s", qpath)

    for template in templates:
        for chain in template[4:]:
            rpath = Path(f"templates/interfaces/{template}_{chain}_int.pdb")
            if rpath.exists():
                selected_ref_paths[(template, chain)] = rpath.resolve()
            else:
                missing_templates.append(f"{template}_{chain}")
                logger.warning("Missing template interface file: %s", rpath)

    # --- Coverage summary ---
    total_expected = len(queries) * sum(len(t[4:]) for t in templates)
    total_actual = len(selected_query_paths) * len(selected_ref_paths)
    if missing_queries or missing_templates:
        logger.warning(
            "GTalign coverage reduced: %d/%d queries missing, %d/%d template-chains missing "
            "(expected %d pairs, have %d pairs)",
            len(missing_queries), len(queries),
            len(missing_templates), sum(len(t[4:]) for t in templates),
            total_expected, total_actual,
        )
        if missing_queries:
            logger.info("Missing queries (%d): %s", len(missing_queries), missing_queries[:10])
        if missing_templates:
            logger.info("Missing templates (%d): %s", len(missing_templates), missing_templates[:10])

    if not selected_query_paths or not selected_ref_paths:
        raise RuntimeError(
            f"GTalign alignment stage: no valid query/reference files found. "
            f"({len(selected_query_paths)} queries, {len(selected_ref_paths)} refs)"
        )

    parsed_pairs = set()
    raw_hash_by_protein = {}
    with tempfile.TemporaryDirectory(prefix="gtalign_prism_stage_", dir="processed") as td:
        td = Path(td)
        qdir = td / "queries"
        rdir = td / "refs"
        outdir = td / "out"
        qdir.mkdir(parents=True, exist_ok=True)
        rdir.mkdir(parents=True, exist_ok=True)
        outdir.mkdir(parents=True, exist_ok=True)

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

        pre_score = float(os.environ.get("PRISM_GTALIGN_PRE_SCORE") or pre_score)
        cmd = [
            str(gtalign_path),
            f"--qrs={qdir}",
            f"--rfs={rdir}",
            "-o", str(outdir),
            "-s", "0",
            f"--nhits=2000",
            f"--nalns=2000",
            f"--dev-min-length={int(dev_min_length)}",
            f"--pre-score={pre_score}",
        ]
        if speed is not None:
            cmd.append(f"--speed={int(speed)}")
        if refinement is not None:
            cmd.append(f"--refinement={int(refinement)}")

        logger.info("Running GTalign with %d queries, %d refs", len(query_basename_to_real), len(ref_basename_to_real))
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"GTalign failed ({result.returncode}): {(result.stderr or result.stdout)[:2000]}"
            )

        total_hits = 0
        skipped_no_match = 0
        skipped_low_tm = 0
        skipped_low_matches = 0
        written_jsons = 0

        for out_file in sorted(outdir.glob("*.out")):
            raw_output_sha256 = hashlib.sha256(out_file.read_bytes()).hexdigest()
            try:
                query_path, hits = parse_gtalign_output_text(out_file.read_text(errors="replace"))
            except (ValueError, IndexError, KeyError, TypeError, AttributeError) as exc:
                logger.warning("Skipping malformed GTalign output %s: %s", out_file, exc)
                continue
            qbase = os.path.basename(query_path)
            if not qbase.endswith(".asa.pdb"):
                logger.warning("Skipping unexpected GTalign query filename: %s", qbase)
                continue
            protein = qbase[:-8]
            protein_path = query_basename_to_real.get(qbase)
            if not protein_path:
                logger.warning("GTalign query not mapped to PRISM input: %s", qbase)
                continue
            raw_hash_by_protein[protein] = raw_output_sha256

            for hit in hits:
                total_hits += 1
                ref_base = os.path.basename(hit["ref_path"])
                if ref_base not in ref_basename_to_key:
                    skipped_no_match += 1
                    continue
                template, chain = ref_basename_to_key[ref_base]
                interface_path = ref_basename_to_real[ref_base]
                match_dict = build_match_dict_from_aligned_sequences(
                    hit["query_aln"],
                    hit["ref_aln"],
                    protein_path,
                    interface_path,
                )
                tm_score = (
                    hit.get("tm_ref", 0.0)
                    if isinstance(hit.get("tm_ref"), (float, int))
                    else 0.0
                )
                match_count = hit["aligned_length"] or len(match_dict)

                min_tm = pre_score
                min_matches = int(os.environ.get("PRISM_MINIMUM_RESIDUE_MATCH_COUNT", "15"))
                tm_candidates = [
                    v for v in (hit["tm_ref"], hit["tm_query"])
                    if isinstance(v, (float, int))
                ]

                if len(tm_candidates) != 2 or tm_score < min_tm or match_count < min_matches:
                    if tm_score < min_tm:
                        skipped_low_tm += 1
                    elif match_count < min_matches:
                        skipped_low_matches += 1
                    else:
                        skipped_no_match += 1
                    parsed_pairs.add((protein, template, chain))
                    continue

                write_alignment_json(
                    output_dir, protein, template, chain,
                    match_count, hit["translation"], hit["rotation_mat"],
                    match_dict, tm_score,
                    tm_score_ref=hit.get("tm_ref"),
                    tm_score_query=hit.get("tm_query"),
                    raw_output_sha256=raw_output_sha256,
                )
                written_jsons += 1
                parsed_pairs.add((protein, template, chain))

    logger.info(
        "GTalign done: %d total hits, %d JSONs written, "
        "%d skipped (no ref match), %d skipped (low TM), %d skipped (low matches)",
        total_hits, written_jsons, skipped_no_match, skipped_low_tm, skipped_low_matches,
    )


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


def write_alignment_json(
    out_dir, protein, template, chain,
    match_count, translation, rotation_mat, match_dict, tm_score,
    *, tm_score_ref=None, tm_score_query=None,
    raw_output_sha256=None, status="success",
):
    payload = {
        "match_count": int(match_count) if match_count is not None else 0,
        "translation": translation or [0.0, 0.0, 0.0],
        "rotation_mat": rotation_mat or [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "match_dict": match_dict or {},
        "tm_score": float(tm_score) if tm_score is not None else 0.0,
        "tm_score_ref": tm_score_ref,
        "tm_score_query": tm_score_query,
        "raw_output_sha256": raw_output_sha256,
        "aligner": "GTalign",
        "status": status,
    }
    with open(Path(out_dir) / f"{protein}_{template}_{chain}.json", "w") as f:
        json.dump(payload, f)


def parse_gtalign_hits(raw_output, *, min_tm_score=0.4, min_match_count=15):
    """Parse and filter GTAlign hit records without fabricating missing hits."""
    if isinstance(raw_output, str):
        _, hits = parse_gtalign_output_text(raw_output)
    else:
        hits = list(raw_output)
    accepted = []
    for hit in hits:
        if not isinstance(hit, dict):
            continue
        tm_candidates = [
            value for value in (hit.get("tm_ref"), hit.get("tm_query"))
            if isinstance(value, (int, float))
        ]
        tm_score = min(tm_candidates, default=0.0)
        match_count = hit.get("aligned_length") or 0
        if len(tm_candidates) == 2 and tm_score >= float(min_tm_score) and int(match_count) >= int(min_match_count):
            accepted.append({**hit, "tm_score": tm_score, "status": "accepted"})
    return accepted


def _symlink_or_copy(src, dst):
    try:
        os.symlink(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def _extract_floats(line):
    return [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", line)]


def _extract_gtalign_alignment_seq(line, label):
    pattern = r"^\s*" + re.escape(label) + r":\s+\d+\s+([A-Za-z\\-]+)\s+\d+\s*$"
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
                r"TM-score \(Refn\./Query\)\s*=\s*([0-9.]+)\s*/\s*([0-9.]+)", line
            )
            if m:
                hit["tm_ref"] = float(m.group(1))
                hit["tm_query"] = float(m.group(2))

        if "Matched =" in line:
            m = re.search(r"Matched\s*=\s*(\d+)/", line)
            if m:
                hit["aligned_length"] = int(m.group(1))

        if line.lstrip().startswith("Query:") and i + 2 < len(lines) and lines[i + 2].lstrip().startswith("Refn."):
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
        hit["aligned_length"] = sum(
            1 for a, b in zip(hit["query_aln"], hit["ref_aln"])
            if a != "-" and b != "-"
        )
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
