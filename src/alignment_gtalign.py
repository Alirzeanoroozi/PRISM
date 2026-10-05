import json
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio.PDB import PDBParser


def align_gtalign(
    queries,
    templates,
    gtalign_path="gtalign",
    output_dir="processed/alignment_gtalign",
    dev_min_length=3,
    pre_score=0.0,
    speed=0,
    refinement=3,
    min_match_count=None,
):
    """
    GTalign-backed alignment stage that writes PRISM-compatible JSONs.

    This is a separate implementation to keep the original TMalign pipeline
    intact. Downstream stages can read the generated JSONs from output_dir.
    """
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise RuntimeError(f"refusing to reuse non-empty GTalign output directory: {output_dir}")
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
        raise RuntimeError("GTalign alignment stage: no valid query/reference files found.")

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

        # Honor the CLI-provided pre-filter unless the environment explicitly
        # overrides it.  GTalign reports ALL hits when --pre-score=0.0,
        # producing 50+ MB output files per query that take hours to parse in
        # Python.  A moderate threshold keeps raw output manageable; the
        # downstream transformer() stage independently enforces the full
        # production thresholds, so no double-filtering leak occurs.
        pre_score = float(os.environ.get("PRISM_GTALIGN_PRE_SCORE", str(pre_score)))
        cmd = [
            str(gtalign_path),
            f"--qrs={qdir}",
            f"--rfs={rdir}",
            "-o",
            str(outdir),
            "-s",
            "0",
            # Limit output to GTalign's default of 2000 hits per query.
            # Using len(selected_ref_paths)=39710 would dump ALL hits into
            # 50+ MB output files that the Python parser takes hours to
            # process.  With --pre-score=0.2, the top 2000 hits per query
            # capture all potentially useful candidates.
            f"--nhits=2000",
            f"--nalns=2000",
            f"--dev-min-length={int(dev_min_length)}",
            f"--pre-score={pre_score}",
        ]
        if speed is not None:
            cmd.append(f"--speed={int(speed)}")
        if refinement is not None:
            cmd.append(f"--refinement={int(refinement)}")

        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode != 0:
            raise RuntimeError(f"GTalign failed ({result.returncode}): {(result.stderr or result.stdout)[:2000]}")

        for out_file in sorted(outdir.glob("*.out")):
            raw_output_sha256 = hashlib.sha256(out_file.read_bytes()).hexdigest()
            try:
                query_path, hits = parse_gtalign_output_text(out_file.read_text(errors="replace"))
            except (ValueError, IndexError, KeyError, TypeError, AttributeError) as exc:
                print(f"Skipping malformed GTalign output {out_file}: {exc}")
                continue
            qbase = os.path.basename(query_path)
            if not qbase.endswith(".asa.pdb"):
                print(f"Skipping unexpected GTalign query filename: {qbase}")
                continue
            protein = qbase[:-8]
            protein_path = query_basename_to_real.get(qbase)
            if not protein_path:
                print(f"GTalign query not mapped to PRISM input: {qbase}")
                continue
            raw_hash_by_protein[protein] = raw_output_sha256

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
                tm_score = max(tm_candidates) if len(tm_candidates) == 2 else 0.0
                match_count = hit["aligned_length"] or len(match_dict)

                # Pre-filter: skip writing JSON for low-quality hits.
                # Honor the CLI pre-score unless the environment overrides it.
                # The downstream transformer() stage independently enforces
                # the full production thresholds.
                min_tm = float(os.environ.get("PRISM_GTALIGN_PRE_SCORE", str(pre_score)))
                min_matches = (
                    int(min_match_count)
                    if min_match_count is not None
                    else int(os.environ.get("PRISM_MINIMUM_RESIDUE_MATCH_COUNT", "15"))
                )
                if (
                    len(tm_candidates) != 2
                    or min(tm_candidates) < min_tm
                    or match_count < min_matches
                ):
                    parsed_pairs.add((protein, template, chain))
                    continue

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
                    tm_score_ref=hit.get("tm_ref"),
                    tm_score_query=hit.get("tm_query"),
                    raw_output_sha256=raw_output_sha256,
                    return_code=result.returncode,
                )
                parsed_pairs.add((protein, template, chain))

    # Only write JSONs for actual hits.  transformer.py handles
    # missing files gracefully via try/except around load_alignment().


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


def write_alignment_json(
    out_dir,
    protein,
    template,
    chain,
    match_count,
    translation,
    rotation_mat,
    match_dict,
    tm_score,
    *,
    tm_score_ref=None,
    tm_score_query=None,
    raw_output_sha256=None,
    return_code=None,
    status="success",
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
        "return_code": return_code,
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
        tm_candidates = [value for value in (hit.get("tm_ref"), hit.get("tm_query")) if isinstance(value, (int, float))]
        tm_score = max(tm_candidates, default=0.0)
        match_count = hit.get("aligned_length") or 0
        if (
            len(tm_candidates) == 2
            and min(tm_candidates) >= float(min_tm_score)
            and int(match_count) >= int(min_match_count)
        ):
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
    """Parse a single GTalign hit block from plain text output.

    GTalign plain text hit format:
        [spaces]N ...path/to/file.pdb Chn:X tm_score_query tm_score_ref rmsd n_aligned ...
    """
    i = start_idx
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines):
        return None, i

    line = lines[i]
    # Check if this is a hit summary line (starts with spaces, then number, then dots)
    if not re.match(r"^\s*\d+\s*\.\s*$", line):
        return None, i + 1

    # GTalign 0.19 emits the reference path on the line immediately after
    # the numbered hit marker (prefixed with ``>``).  Older output placed it
    # on the marker line itself, so accept both forms.
    ref_path = ""
    inline_path = line.strip().split("...", 1)
    if len(inline_path) == 2:
        ref_path = inline_path[1].split(" Chn:", 1)[0].strip()

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
        if not ref_path and line.lstrip().startswith(">"):
            ref_path = line.lstrip()[1:].split(" Chn:", 1)[0].strip()
            hit["ref_path"] = ref_path
            i += 1
            continue
        if re.match(r"^\s*\d+\.\s*$", line) or line.startswith("Query length:"):
            break

        if "TM-score (Refn./Query)" in line:
            m = re.search(r"TM-score \(Refn\./Query\)\s*=\s*([0-9.]+)\s*/\s*([0-9.]+)", line)
            if m:
                hit["tm_ref"] = float(m.group(1))
                hit["tm_query"] = float(m.group(2))

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
        # Match lines like "     1 ...tes/interfaces/..." (number, spaces, then dots/content)
        # GTalign output format: spaces + number + spaces + dots + path
        if re.match(r"^\s*\d+\s*\.\s*$", lines[i]):
            hit, i = _parse_gtalign_hit_block(lines, i)
            if hit:
                hits.append(hit)
            continue
        i += 1
    return query_path, hits
