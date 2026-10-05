"""
Fixed version of src/transformation.py

Key fixes:
1. apply_tm_transform() returns a TransformResult with diagnostics
   (file missing, parse error, NaN coords, success).
2. create_transformed_pair() propagates diagnostics.
3. pair_has_acceptable_clashes() logs clash counts.
4. Added transform diagnostics to audit records.
5. alignment_passes_thresholds() checks for mapping_truncated status.
"""

import json
import logging
import math
import os

import pandas as pd

from .alignment import MAPPING_TRUNCATED_STATUS
import candidate_audit as _ca
import pdb_download as _pd
import template_filtering as _tf
import utils as _ut

logger = logging.getLogger(__name__)

os.makedirs("processed/transformation", exist_ok=True)

MINIMUM_RESIDUE_MATCH_COUNT = int(os.environ.get("PRISM_MINIMUM_RESIDUE_MATCH_COUNT", "15"))
MINIMUM_RESIDUE_MATCH_PERCENTAGE = float(os.environ.get("PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE", "50.0"))
MINIMUM_HOTSPOT_MATCH_NUMBER = 1
DIFF_PERCENTAGE = float(os.environ.get("PRISM_DIFF_PERCENTAGE", "20.0"))
CONTACT_COUNT = 5
CLASHING_DISTANCE = float(os.environ.get("PRISM_CLASHING_DISTANCE", "3"))
MAX_CLASHING_COUNT = int(os.environ.get("PRISM_MAX_CLASHING_COUNT", "5"))
TM_SCORE_THRESHOLD = float(os.environ.get("PRISM_TM_SCORE_THRESHOLD", "0.5"))
HOTSPOT_CRITERION = 2
HOTSPOT_COUNT = 1
TEMPLATE_RESIDUE_COUNT = 50
CONTACT_COUNT_THRESHOLD = 5
AUDIT_PATH = os.environ.get("PRISM_CANDIDATE_AUDIT_PATH")

passed_pairs = []
template_size = {}
INPUTS_CSV = os.environ.get("PRISM_INPUTS_CSV", "inputs.csv")
FILTER_MODE = os.environ.get("PRISM_FILTER_MODE", "geometry_only_experimental")
FILTER_ASSET_ROOT = os.environ.get("PRISM_FILTER_ASSET_ROOT", "")
if FILTER_MODE not in {"published_protocol", "geometry_only_experimental"}:
    raise ValueError("PRISM_FILTER_MODE must be published_protocol or geometry_only_experimental")


class TransformResult:
    """Structured result from a transform operation."""

    def __init__(self, success=False, error_code=None, error_detail=None, nan_coords=0):
        self.success = success
        self.error_code = error_code
        self.error_detail = error_detail
        self.nan_coords = nan_coords

    @classmethod
    def ok(cls):
        return cls(success=True)

    @classmethod
    def failure(cls, code, detail=None):
        return cls(success=False, error_code=code, error_detail=detail)

    def __bool__(self):
        return self.success


def transformer(templates, alignment_dir="processed/alignment"):
    df = pd.read_csv(INPUTS_CSV)

    for template in templates:
        chain1 = template[4]
        chain2 = template[5]

        with open(os.path.join("templates", "interfaces_lists", f"{template}.json"), "r") as f:
            data = json.load(f)

        template_size[f"{template}_{chain1}"] = len(data[chain1])
        template_size[f"{template}_{chain2}"] = len(data[chain2])

        for left_query, right_query in zip(df["Receptor"], df["Ligand"]):
            process_pair_for_template(
                template, chain1, chain2, left_query, right_query,
                alignment_dir=alignment_dir,
            )

    return passed_pairs


def load_alignment(query_id, template, chain_id, alignment_dir="processed/alignment"):
    candidates = [_pd.normalize_target_id(query_id)]
    if str(query_id) not in candidates:
        candidates.append(str(query_id))
    for candidate in candidates:
        path = os.path.join(alignment_dir, f"{candidate}_{template}_{chain_id}.json")
        if os.path.isfile(path):
            with open(path, "r") as f:
                return json.load(f)
    raise FileNotFoundError(
        f"alignment not found for query {query_id} (tried {', '.join(candidates)})"
    )


def hotspot_analysis(match_dict, alignment=None, hotspots=None):
    if FILTER_MODE == "geometry_only_experimental":
        return True
    alignment = alignment or {}
    if hotspots is None:
        hotspots = alignment.get("hotspots")
    return _tf.evaluate_hotspots(match_dict, hotspots, minimum=MINIMUM_HOTSPOT_MATCH_NUMBER).passed


def alignment_passes_thresholds(template_key, alignment, protocol_hotspots=None):
    match_count = alignment.get("match_count", 0)
    tm_score = alignment.get("tm_score", 0.0)
    match_dict = alignment.get("match_dict", {})
    align_status = alignment.get("status", "unknown")

    # --- FIX: Reject truncated mappings ---
    if align_status == MAPPING_TRUNCATED_STATUS:
        logger.warning(
            "Rejecting alignment with truncated mapping for %s: "
            "match_count=%d tm_score=%.3f",
            template_key, match_count, tm_score,
        )
        return False

    protein_size = float(template_size.get(template_key, 0))
    if protein_size <= 0:
        if match_count < MINIMUM_RESIDUE_MATCH_COUNT or tm_score < TM_SCORE_THRESHOLD:
            return False
        return True

    match_score = (match_count / protein_size) * 100.0

    if not hotspot_analysis(match_dict, alignment, hotspots=protocol_hotspots):
        return False

    if match_count < MINIMUM_RESIDUE_MATCH_COUNT or tm_score < TM_SCORE_THRESHOLD:
        return False

    if protein_size > TEMPLATE_RESIDUE_COUNT:
        return match_score > (MINIMUM_RESIDUE_MATCH_PERCENTAGE - DIFF_PERCENTAGE)
    else:
        return match_score > MINIMUM_RESIDUE_MATCH_PERCENTAGE


def _protocol_hotspots(filter_assets, chain_id):
    if not filter_assets:
        return None
    by_chain = filter_assets.get("hotspots_by_chain")
    if isinstance(by_chain, dict):
        return by_chain.get(str(chain_id), [])
    return filter_assets.get("hotspots", [])


def _protocol_contacts(filter_assets, left_chain, right_chain):
    if not filter_assets:
        return []
    contacts = filter_assets.get("contacts", [])
    if not contacts:
        return []
    oriented = []
    for contact in contacts:
        if isinstance(contact, dict):
            left = contact.get("left", contact.get("a"))
            right = contact.get("right", contact.get("b"))
        elif isinstance(contact, (list, tuple)) and len(contact) >= 2:
            left, right = contact[0], contact[1]
        else:
            continue
        left_text, right_text = str(left), str(right)
        if left_text.startswith(f"{left_chain}.") and right_text.startswith(f"{right_chain}."):
            oriented.append([left, right])
        elif left_text.startswith(f"{right_chain}.") and right_text.startswith(f"{left_chain}."):
            oriented.append([right, left])
    return oriented


def process_pair_for_template(
    template, chain1, chain2, left_query, right_query,
    alignment_dir="processed/alignment",
):
    left_key_chain1 = f"{template}_{chain1}"
    left_key_chain2 = f"{template}_{chain2}"
    empty_align = {
        "match_count": 0, "tm_score": 0.0, "match_dict": {},
        "translation": [0.0, 0.0, 0.0],
        "rotation_mat": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    }

    try:
        left_align_1 = load_alignment(left_query, template, chain1, alignment_dir=alignment_dir)
    except FileNotFoundError:
        left_align_1 = dict(empty_align)
    try:
        right_align_1 = load_alignment(right_query, template, chain2, alignment_dir=alignment_dir)
    except FileNotFoundError:
        right_align_1 = dict(empty_align)

    filter_assets = None
    if FILTER_MODE == "published_protocol":
        try:
            filter_assets = _tf.load_filter_assets(template, FILTER_ASSET_ROOT)
        except (FileNotFoundError, ValueError) as exc:
            write_audit_record(
                template, left_query, right_query, chain1, chain2, "o1",
                left_align_1, right_align_1,
                f"protocol_assets_failed:{exc}",
            )
            filter_assets = {"hotspots": [], "hotspots_by_chain": {}, "contacts": []}

    left1_hotspots = _protocol_hotspots(filter_assets, chain1)
    right1_hotspots = _protocol_hotspots(filter_assets, chain2)
    contacts_o1 = _protocol_contacts(filter_assets, chain1, chain2)
    published_pair = (
        FILTER_MODE != "published_protocol"
        or _tf.evaluate_protocol_candidate(
            left_align_1.get("match_dict", {}),
            right_align_1.get("match_dict", {}),
            left1_hotspots if filter_assets else left_align_1.get("hotspots"),
            right1_hotspots if filter_assets else right_align_1.get("hotspots"),
            contacts_o1 if filter_assets else left_align_1.get("contacts", []),
            minimum_contacts=CONTACT_COUNT_THRESHOLD,
        ).passed
    )

    if (published_pair
            and alignment_passes_thresholds(left_key_chain1, left_align_1, left1_hotspots)
            and alignment_passes_thresholds(left_key_chain2, right_align_1, right1_hotspots)):
        status = create_transformed_pair(
            template, left_query, right_query, left_align_1, right_align_1,
            passed_pairs, "o1",
        )
    else:
        status = "alignment_failed"
    write_audit_record(
        template, left_query, right_query, chain1, chain2, "o1",
        left_align_1, right_align_1, status,
    )

    try:
        left_align_2 = load_alignment(left_query, template, chain2, alignment_dir=alignment_dir)
    except FileNotFoundError:
        left_align_2 = dict(empty_align)
    try:
        right_align_2 = load_alignment(right_query, template, chain1, alignment_dir=alignment_dir)
    except FileNotFoundError:
        right_align_2 = dict(empty_align)

    left2_hotspots = _protocol_hotspots(filter_assets, chain2)
    right2_hotspots = _protocol_hotspots(filter_assets, chain1)
    contacts_o2 = _protocol_contacts(filter_assets, chain2, chain1)
    published_pair = (
        FILTER_MODE != "published_protocol"
        or _tf.evaluate_protocol_candidate(
            left_align_2.get("match_dict", {}),
            right_align_2.get("match_dict", {}),
            left2_hotspots if filter_assets else left_align_2.get("hotspots"),
            right2_hotspots if filter_assets else right_align_2.get("hotspots"),
            contacts_o2 if filter_assets else left_align_2.get("contacts", []),
            minimum_contacts=CONTACT_COUNT_THRESHOLD,
        ).passed
    )

    if (published_pair
            and alignment_passes_thresholds(left_key_chain2, left_align_2, left2_hotspots)
            and alignment_passes_thresholds(left_key_chain1, right_align_2, right2_hotspots)):
        status = create_transformed_pair(
            template, left_query, right_query, left_align_2, right_align_2,
            passed_pairs, "o2",
        )
    else:
        status = "alignment_failed"
    write_audit_record(
        template, left_query, right_query, chain2, chain1, "o2",
        left_align_2, right_align_2, status,
    )


def write_audit_record(
    template, left_query, right_query, chain_left, chain_right,
    orientation, left_alignment, right_alignment, status,
):
    if not AUDIT_PATH:
        return
    audit = _ca.CandidateAudit(AUDIT_PATH)
    _ca.record_alignment_pair(
        audit,
        query_left=left_query,
        query_right=right_query,
        template=template,
        chain_left=chain_left,
        chain_right=chain_right,
        orientation=orientation,
        left_alignment=left_alignment,
        right_alignment=right_alignment,
        metadata={"native_complex_id": os.environ.get("PRISM_NATIVE_COMPLEX_ID")},
        status=status,
    )


def create_transformed_pair(
    template, left_query, right_query, left_alignment, right_alignment,
    passed_pairs, orientation_suffix,
):
    left_id = _pd.normalize_target_id(left_query)
    right_id = _pd.normalize_target_id(right_query)
    left_input = f"processed/pdbs/{left_id}.pdb"
    right_input = f"processed/pdbs/{right_id}.pdb"
    if not os.path.exists(left_input):
        left_input = f"processed/pdbs/{left_id[:4]}.pdb"
    if not os.path.exists(right_input):
        right_input = f"processed/pdbs/{right_id[:4]}.pdb"

    left_output = (
        f"processed/transformation/{template}_{left_query}_{right_query}"
        f"_{orientation_suffix}_L.pdb"
    )
    right_output = (
        f"processed/transformation/{template}_{left_query}_{right_query}"
        f"_{orientation_suffix}_R.pdb"
    )

    left_result = apply_tm_transform(
        left_input, left_output,
        left_alignment.get("translation", [0.0, 0.0, 0.0]),
        left_alignment.get("rotation_mat", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    right_result = apply_tm_transform(
        right_input, right_output,
        right_alignment.get("translation", [0.0, 0.0, 0.0]),
        right_alignment.get("rotation_mat", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    )

    if not left_result:
        return f"transformation_failed:left:{left_result.error_code}"
    if not right_result:
        return f"transformation_failed:right:{right_result.error_code}"
    if not os.path.exists(left_output):
        return "transformation_failed:left:missing_output"
    if not os.path.exists(right_output):
        return "transformation_failed:right:missing_output"

    clash_ok = pair_has_acceptable_clashes(left_output, right_output)
    if not clash_ok:
        return "clash_rejected"

    if os.path.exists(left_output) and os.path.exists(right_output):
        passed_pairs.append((left_output, right_output))
        return "generated"

    return "transformation_failed:unknown"


def apply_tm_transform(input_pdb, output_pdb, translation, rotation_mat):
    """
    Apply a TM-align rotation+translation to a PDB file.

    Returns a TransformResult with diagnostics on failure.
    """
    if not os.path.exists(input_pdb):
        return TransformResult.failure(
            "INPUT_MISSING",
            f"Input PDB does not exist: {input_pdb}",
        )

    nan_count = 0
    try:
        with open(input_pdb, "r") as in_f, open(output_pdb, "w") as out_f:
            for line in in_f:
                if line.startswith("ATOM"):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                    except ValueError as exc:
                        logger.warning(
                            "Coordinate parse error in %s at line: %s",
                            input_pdb, line[:60],
                        )
                        out_f.write(line)
                        continue

                    new_x = (
                        x * rotation_mat[0][0]
                        + y * rotation_mat[0][1]
                        + z * rotation_mat[0][2]
                        + translation[0]
                    )
                    new_y = (
                        x * rotation_mat[1][0]
                        + y * rotation_mat[1][1]
                        + z * rotation_mat[1][2]
                        + translation[1]
                    )
                    new_z = (
                        x * rotation_mat[2][0]
                        + y * rotation_mat[2][1]
                        + z * rotation_mat[2][2]
                        + translation[2]
                    )

                    # Check for NaN in output coordinates
                    if math.isnan(new_x) or math.isnan(new_y) or math.isnan(new_z):
                        nan_count += 1
                        out_f.write(line)
                        continue

                    line = (
                        f"{line[:30]}"
                        f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}"
                        f"{line[54:]}"
                    )
                out_f.write(line)

        if nan_count > 0:
            logger.warning(
                "Transform produced %d NaN coordinates in %s",
                nan_count, input_pdb,
            )
            return TransformResult(
                success=False,
                error_code="NAN_COORDS",
                error_detail=f"{nan_count} atoms had NaN coordinates after transform",
                nan_coords=nan_count,
            )

        return TransformResult.ok()

    except Exception as exc:
        logger.error("Error applying TM transform to %s: %s", input_pdb, exc)
        return TransformResult.failure(
            "TRANSFORM_EXCEPTION",
            f"{type(exc).__name__}: {exc}",
        )


def pair_has_acceptable_clashes(left_path, right_path):
    left_coords = _ut.read_ca_coordinates(left_path)
    right_coords = _ut.read_ca_coordinates(right_path)

    clash_count = 0
    for c1 in left_coords:
        for c2 in right_coords:
            if _ut.distance_calculator(c1, c2) < CLASHING_DISTANCE:
                clash_count += 1
                if clash_count >= MAX_CLASHING_COUNT:
                    logger.info(
                        "Clash threshold reached: %d clashes (max %d) between %s and %s",
                        clash_count, MAX_CLASHING_COUNT, left_path, right_path,
                    )
                    return False

    if clash_count > 0:
        logger.debug("Pair has %d clashes (threshold: %d)", clash_count, MAX_CLASHING_COUNT)

    return True
