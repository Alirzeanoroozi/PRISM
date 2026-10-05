import json
import os
import pandas as pd

from .utils import read_ca_coordinates, distance_calculator
from .candidate_audit import CandidateAudit, record_alignment_pair
from .pdb_download import normalize_target_id
from .template_filtering import evaluate_protocol_candidate, evaluate_hotspots, load_filter_assets

os.makedirs("processed/transformation", exist_ok=True)

MINIMUM_RESIDUE_MATCH_COUNT = int(os.environ.get("PRISM_MINIMUM_RESIDUE_MATCH_COUNT", "15"))
MINIMUM_RESIDUE_MATCH_PERCENTAGE = float(os.environ.get("PRISM_MINIMUM_RESIDUE_MATCH_PERCENTAGE", "50.0"))
MINIMUM_HOTSPOT_MATCH_NUMBER = 1
DIFF_PERCENTAGE = float(os.environ.get("PRISM_DIFF_PERCENTAGE", "20.0"))
CONTACT_COUNT = 5
CLASHING_DISTANCE = float(os.environ.get("PRISM_CLASHING_DISTANCE", "3"))
MAX_CLASHING_COUNT = int(os.environ.get("PRISM_MAX_CLASHING_COUNT", "5"))
TM_SCORE_THRESHOLD = float(os.environ.get("PRISM_TM_SCORE_THRESHOLD", "0.5"))
# MultiProt does not emit a TMalign-compatible TM-score. Its alignment JSON
# contains a solver correspondence count and a Kabsch RMSD instead. Keep the
# TMalign threshold for TMalign/GTalign records and use the native MultiProt
# match/coverage gates below.
MULTIPROT_SCORE_CONTRACT = "native_match_count_and_coverage"
HOTSPOT_CRITERION = 2
HOTSPOT_COUNT = 1
TEMPLATE_RESIDUE_COUNT = 50
CONTACT_COUNT_THRESHOLD = 5
# Kept as an override for compatibility with callers/tests that set the module
# attribute.  Environment lookup happens at write time so a CLI option parsed
# after this module is imported can still configure the audit destination.
AUDIT_PATH = None

passed_pairs = []
template_size = {}
# Local change: allow test runs to override the pair list without editing inputs.csv.
INPUTS_CSV = os.environ.get("PRISM_INPUTS_CSV", "inputs.csv")
FILTER_MODE = os.environ.get("PRISM_FILTER_MODE", "geometry_only_experimental")
FILTER_ASSET_ROOT = os.environ.get("PRISM_FILTER_ASSET_ROOT", "")
if FILTER_MODE not in {"published_protocol", "geometry_only_experimental"}:
    raise ValueError("PRISM_FILTER_MODE must be published_protocol or geometry_only_experimental")

def transformer(templates, alignment_dir="processed/alignment", audit_path=None):
    # Local change: read from env-configurable CSV path.
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
                alignment_dir=alignment_dir, audit_path=audit_path,
            )

    return passed_pairs

def load_alignment(query_id, template, chain_id, alignment_dir="processed/alignment"):
    # Alignment emits canonical PDB+chain tokens (for example ``1fgnH``),
    # while the benchmark CSV intentionally preserves raw selector spelling
    # (for example ``1FGNH``).  Resolve the canonical spelling first and keep
    # the raw spelling as a compatibility fallback for historical outputs.
    candidates = [normalize_target_id(query_id)]
    if str(query_id) not in candidates:
        candidates.append(str(query_id))
    for candidate in candidates:
        path = os.path.join(alignment_dir, f"{candidate}_{template}_{chain_id}.json")
        if os.path.isfile(path):
            with open(path, "r") as f:
                return json.load(f)
    raise FileNotFoundError(
        "alignment not found for query {} (tried {})".format(query_id, ", ".join(candidates))
    )

def hotspot_analysis(match_dict, alignment=None, hotspots=None):
    if FILTER_MODE == "geometry_only_experimental":
        return True
    alignment = alignment or {}
    if hotspots is None:
        hotspots = alignment.get("hotspots")
    return evaluate_hotspots(match_dict, hotspots, minimum=MINIMUM_HOTSPOT_MATCH_NUMBER).passed

def alignment_score_passes(alignment):
    """Apply the aligner-specific alignment score contract.

    TMalign and GTalign records expose a TM-score on the shared tm_score
    field. MultiProt's field is now a true TM-score (computed from match_dict),
    but distributions differ from TMalign. Lowered thresholds based on
    diagnostic run: true TM ~0.05-0.6, matches 4-68.
    """

    aligner = str(alignment.get("aligner", "")).strip().lower()
    tm_score = float(alignment.get("tm_score", 0.0) or 0.0)
    match_count = alignment.get("match_count", 0)
    match_dict = alignment.get("match_dict", {})
    template_size = len(match_dict) if match_dict else 0
    match_pct = (match_count / template_size * 100.0) if template_size > 0 else 0.0

    if aligner == "multiprot":
        # MultiProt's wrapper does NOT emit a TMalign-compatible TM-score. It
        # writes an RMSD-proxy value (``1 - kabsch_rmsd/10``) on the ts_score
        # field and declares ``score_gate_contract: native_match_count_and_coverage``.
        # Gate on the native match count + interface coverage, NOT the proxy
        # score, so real MultiProt solutions are not discarded. When no
        # match_dict is available (coverage not computable), fall back to the
        # match-count gate alone.
        if template_size <= 0:
            return match_count >= 10
        return match_count >= 10 and match_pct >= 30.0
    return tm_score >= TM_SCORE_THRESHOLD


def alignment_passes_thresholds(template_key, alignment, protocol_hotspots=None):
    match_count = alignment.get("match_count", 0)
    tm_score = alignment.get("tm_score", 0.0)
    match_dict = alignment.get("match_dict", {})

    protein_size = float(template_size.get(template_key, 0))
    if protein_size <= 0:
        # Without a size estimate we cannot compute match percentage; fall
        # back to simple count + TM-score checks.
        if match_count < MINIMUM_RESIDUE_MATCH_COUNT or not alignment_score_passes(alignment):
            return False
        return True

    match_score = (match_count / protein_size) * 100.0

    if not hotspot_analysis(match_dict, alignment, hotspots=protocol_hotspots):
        return False

    if match_count < MINIMUM_RESIDUE_MATCH_COUNT or not alignment_score_passes(alignment):
        return False

    if protein_size > TEMPLATE_RESIDUE_COUNT:
        return match_score > (MINIMUM_RESIDUE_MATCH_PERCENTAGE - DIFF_PERCENTAGE)
    else:
        return match_score > MINIMUM_RESIDUE_MATCH_PERCENTAGE


def _protocol_hotspots(filter_assets, chain_id):
    """Return hotspots for one template chain, preserving legacy fallback."""

    if not filter_assets:
        return None
    by_chain = filter_assets.get("hotspots_by_chain")
    if isinstance(by_chain, dict):
        return by_chain.get(str(chain_id), [])
    # Older callers may provide only a flattened list.  Keep that API
    # compatible, but never use it when chain-specific modern assets exist.
    return filter_assets.get("hotspots", [])


def _protocol_contacts(filter_assets, left_chain, right_chain):
    """Orient template contact pairs to match the two alignment sides."""

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


def _alignment_variants(alignment):
    """Return primary plus legacy MultiProt solution variants.

    Current and GTalign records have one transform.  The opt-in MultiProt
    compatibility record carries the legacy solver's retained solutions; each
    solution is promoted to the same alignment schema so existing threshold
    and clash checks can evaluate it without changing the default path.
    """
    solutions = alignment.get("multiprot_solutions")
    if not isinstance(solutions, list) or not solutions:
        return [alignment]

    variants = []
    for solution in solutions:
        variant = dict(alignment)
        variant.update(solution)
        variant["aligner"] = alignment.get("aligner", "MultiProt")
        variant["tm_score"] = alignment.get("tm_score", 0.0)
        variant["tm_score_contract"] = alignment.get(
            "tm_score_contract", "multiprot_legacy_native"
        )
        variant["score_gate_contract"] = alignment.get(
            "score_gate_contract", MULTIPROT_SCORE_CONTRACT
        )
        variant["multiprot_mode"] = alignment.get("multiprot_mode")
        variants.append(variant)
    return variants

def process_pair_for_template(
    template,
    chain1,
    chain2,
    left_query,
    right_query,
    alignment_dir="processed/alignment",
    audit_path=None,
):
    left_key_chain1 = f"{template}_{chain1}"
    left_key_chain2 = f"{template}_{chain2}"

    try:
        left_align_1 = load_alignment(left_query, template, chain1, alignment_dir=alignment_dir)
    except FileNotFoundError:
        left_align_1 = {"match_count": 0, "tm_score": 0.0, "match_dict": {}, "translation": [0.0, 0.0, 0.0], "rotation_mat": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]}
    try:
        right_align_1 = load_alignment(right_query, template, chain2, alignment_dir=alignment_dir)
    except FileNotFoundError:
        right_align_1 = {"match_count": 0, "tm_score": 0.0, "match_dict": {}, "translation": [0.0, 0.0, 0.0], "rotation_mat": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]}

    filter_assets = None
    if FILTER_MODE == "published_protocol":
        try:
            filter_assets = load_filter_assets(template, FILTER_ASSET_ROOT)
        except (FileNotFoundError, ValueError) as exc:
            write_audit_record(
                template, left_query, right_query, chain1, chain2, "o1",
                left_align_1, right_align_1, f"protocol_assets_failed:{exc}",
                audit_path=audit_path,
            )
            filter_assets = {"hotspots": [], "hotspots_by_chain": {}, "contacts": []}
    left1_hotspots = _protocol_hotspots(filter_assets, chain1)
    right1_hotspots = _protocol_hotspots(filter_assets, chain2)
    contacts_o1 = _protocol_contacts(filter_assets, chain1, chain2)
    left_variants_1 = _alignment_variants(left_align_1)
    right_variants_1 = _alignment_variants(right_align_1)
    status = "alignment_failed"
    for left_index, left_variant in enumerate(left_variants_1):
        for right_index, right_variant in enumerate(right_variants_1):
            published_pair = FILTER_MODE != "published_protocol" or evaluate_protocol_candidate(
                left_variant.get("match_dict", {}), right_variant.get("match_dict", {}),
                left1_hotspots if filter_assets else left_variant.get("hotspots"),
                right1_hotspots if filter_assets else right_variant.get("hotspots"),
                contacts_o1 if filter_assets else left_variant.get("contacts", []),
                minimum_contacts=CONTACT_COUNT_THRESHOLD,
            ).passed
            if not published_pair:
                continue
            if not (
                alignment_passes_thresholds(left_key_chain1, left_variant, left1_hotspots)
                and alignment_passes_thresholds(left_key_chain2, right_variant, right1_hotspots)
            ):
                continue
            suffix = "o1" if len(left_variants_1) == len(right_variants_1) == 1 else f"o1_s{left_index}_{right_index}"
            candidate_status = create_transformed_pair(
                template, left_query, right_query, left_variant, right_variant,
                passed_pairs, suffix,
            )
            if candidate_status == "generated":
                status = "generated"
            elif status != "generated":
                status = candidate_status
    write_audit_record(
        template, left_query, right_query, chain1, chain2, "o1",
        left_align_1, right_align_1, status, audit_path=audit_path,
    )

    try:
        left_align_2 = load_alignment(left_query, template, chain2, alignment_dir=alignment_dir)
    except FileNotFoundError:
        left_align_2 = {"match_count": 0, "tm_score": 0.0, "match_dict": {}, "translation": [0.0, 0.0, 0.0], "rotation_mat": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]}
    try:
        right_align_2 = load_alignment(right_query, template, chain1, alignment_dir=alignment_dir)
    except FileNotFoundError:
        right_align_2 = {"match_count": 0, "tm_score": 0.0, "match_dict": {}, "translation": [0.0, 0.0, 0.0], "rotation_mat": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]}

    left2_hotspots = _protocol_hotspots(filter_assets, chain2)
    right2_hotspots = _protocol_hotspots(filter_assets, chain1)
    contacts_o2 = _protocol_contacts(filter_assets, chain2, chain1)
    left_variants_2 = _alignment_variants(left_align_2)
    right_variants_2 = _alignment_variants(right_align_2)
    status = "alignment_failed"
    for left_index, left_variant in enumerate(left_variants_2):
        for right_index, right_variant in enumerate(right_variants_2):
            published_pair = FILTER_MODE != "published_protocol" or evaluate_protocol_candidate(
                left_variant.get("match_dict", {}), right_variant.get("match_dict", {}),
                left2_hotspots if filter_assets else left_variant.get("hotspots"),
                right2_hotspots if filter_assets else right_variant.get("hotspots"),
                contacts_o2 if filter_assets else left_variant.get("contacts", []),
                minimum_contacts=CONTACT_COUNT_THRESHOLD,
            ).passed
            if not published_pair:
                continue
            if not (
                alignment_passes_thresholds(left_key_chain2, left_variant, left2_hotspots)
                and alignment_passes_thresholds(left_key_chain1, right_variant, right2_hotspots)
            ):
                continue
            suffix = "o2" if len(left_variants_2) == len(right_variants_2) == 1 else f"o2_s{left_index}_{right_index}"
            candidate_status = create_transformed_pair(
                template, left_query, right_query, left_variant, right_variant,
                passed_pairs, suffix,
            )
            if candidate_status == "generated":
                status = "generated"
            elif status != "generated":
                status = candidate_status
    write_audit_record(
        template, left_query, right_query, chain2, chain1, "o2",
        left_align_2, right_align_2, status, audit_path=audit_path,
    )


def write_audit_record(
    template,
    left_query,
    right_query,
    chain_left,
    chain_right,
    orientation,
    left_alignment,
    right_alignment,
    status,
    audit_path=None,
):
    resolved_audit_path = audit_path or AUDIT_PATH or os.environ.get("PRISM_CANDIDATE_AUDIT_PATH")
    if not resolved_audit_path:
        return
    audit = CandidateAudit(resolved_audit_path)
    record_alignment_pair(
        audit,
        query_left=left_query,
        query_right=right_query,
        template=template,
        chain_left=chain_left,
        chain_right=chain_right,
        orientation=orientation,
        left_alignment=left_alignment,
        right_alignment=right_alignment,
        template_size_left=template_size.get(f"{template}_{chain_left}"),
        template_size_right=template_size.get(f"{template}_{chain_right}"),
        metadata={"native_complex_id": os.environ.get("PRISM_NATIVE_COMPLEX_ID")},
        status=status,
    )

def create_transformed_pair(template, left_query, right_query, left_alignment, right_alignment, passed_pairs, orientation_suffix):
    # Local change: input IDs include chain suffix (e.g., 1abcA), but PDB files are stored by 4-char PDB ID.
    left_id = normalize_target_id(left_query)
    right_id = normalize_target_id(right_query)
    left_input = f"processed/pdbs/{left_id}.pdb"
    right_input = f"processed/pdbs/{right_id}.pdb"
    if not os.path.exists(left_input):
        left_input = f"processed/pdbs/{left_id[:4]}.pdb"
    if not os.path.exists(right_input):
        right_input = f"processed/pdbs/{right_id[:4]}.pdb"

    left_output = f"processed/transformation/{template}_{left_query}_{right_query}_{orientation_suffix}_L.pdb"
    right_output = f"processed/transformation/{template}_{left_query}_{right_query}_{orientation_suffix}_R.pdb"

    # Local change: check transform success before downstream clash filtering.
    left_ok = apply_tm_transform(left_input, left_output, left_alignment.get("translation", [0.0, 0.0, 0.0]), left_alignment.get("rotation_mat", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
    right_ok = apply_tm_transform(right_input, right_output, right_alignment.get("translation", [0.0, 0.0, 0.0]), right_alignment.get("rotation_mat", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))

    if not left_ok or not right_ok or not os.path.exists(left_output) or not os.path.exists(right_output):
        return "transformation_failed"
    if not pair_has_acceptable_clashes(left_output, right_output):
        return "clash_rejected"
    if left_ok and right_ok and os.path.exists(left_output) and os.path.exists(right_output):
        passed_pairs.append((left_output, right_output))
        return "generated"

def apply_tm_transform(input_pdb, output_pdb, translation, rotation_mat):
    try:
        with open(input_pdb, "r") as in_f, open(output_pdb, "w") as out_f:
            for line in in_f:
                if line.startswith("ATOM"):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                    except ValueError:
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

                    line = (
                        f"{line[:30]}"
                        f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}"
                        f"{line[54:]}"
                    )
                out_f.write(line)
        # Local change: return success/failure so caller can skip invalid outputs.
        return True
    except Exception as exc:
        print(f"Error applying TM transform to {input_pdb}: {exc}")
        # Local change: signal transform failure to caller.
        return False

def pair_has_acceptable_clashes(left_path, right_path):
    left_coords = read_ca_coordinates(left_path)
    right_coords = read_ca_coordinates(right_path)

    clash_count = 0
    for c1 in left_coords:
        for c2 in right_coords:
            if distance_calculator(c1, c2) < CLASHING_DISTANCE:
                clash_count += 1
                if clash_count >= MAX_CLASHING_COUNT:
                    return False
    return True
