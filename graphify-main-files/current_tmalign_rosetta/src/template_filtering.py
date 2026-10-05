"""Pure PRISM published-protocol hotspot/contact decisions."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

from Bio.PDB.Polypeptide import protein_letters_3to1


@dataclass(frozen=True)
class FilterDecision:
    passed: bool
    reason: str
    hotspot_matches: int = 0
    complementary_contacts: int = 0


def load_filter_assets(template_id: str, root: str | Path):
    """Load and normalize one protocol asset pair without mixing sources.

    Current JSON assets store hotspots as ``{chain: [[number, resname]]}``
    and contacts as numeric ``[left_number, right_number]`` pairs.  Derived
    legacy assets store chain-qualified, one-letter records.  Both forms are
    normalized here, while ``asset_format`` and raw hashes preserve which
    source was actually used.
    """

    root = Path(root)
    hotspot_path = root / "hotspots" / f"{template_id}.json"
    contact_path = root / "contacts" / f"{template_id}.json"
    if not hotspot_path.is_file() or not contact_path.is_file():
        raise FileNotFoundError(f"missing protocol filter assets for {template_id}")
    try:
        hotspots = json.loads(hotspot_path.read_text(encoding="utf-8"))
        contacts = json.loads(contact_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid protocol filter assets for {template_id}: {exc}") from exc
    if not isinstance(hotspots, (dict, list)) or not isinstance(contacts, list):
        raise ValueError(f"invalid protocol filter asset schema for {template_id}")
    chains = template_id[4:].replace("_", "")
    hotspots_by_chain: dict[str, list] = {}
    asset_format = "derived_legacy_json"
    if isinstance(hotspots, dict):
        asset_format = "modern_json"
        for chain, entries in hotspots.items():
            if not isinstance(entries, list):
                raise ValueError(f"invalid hotspot records for {template_id}_{chain}")
            normalized = []
            for entry in entries:
                if not isinstance(entry, (list, tuple)) or len(entry) < 2:
                    continue
                number, residue = str(entry[0]), str(entry[1]).upper()
                residue = protein_letters_3to1.get(residue, residue[:1])
                normalized.append([str(chain), residue, number])
            hotspots_by_chain[str(chain)] = normalized
    else:
        for entry in hotspots:
            if not isinstance(entry, (list, tuple)):
                continue
            if len(entry) >= 3:
                chain, residue, number = str(entry[0]), str(entry[1]).upper(), str(entry[2])
            elif len(entry) >= 2:
                chain, residue, number = "", str(entry[1]).upper(), str(entry[0])
            else:
                continue
            hotspots_by_chain.setdefault(chain, []).append([chain, residue[:1], number])

    normalized_contacts = []
    for pair in contacts:
        if not isinstance(pair, (list, tuple)) or len(pair) < 2:
            continue
        left, right = pair[0], pair[1]
        if all(isinstance(value, (int, float)) or str(value).strip().isdigit() for value in (left, right)) and len(chains) >= 2:
            asset_format = "modern_json"
            left, right = f"{chains[0]}..{int(left)}", f"{chains[1]}..{int(right)}"
        normalized_contacts.append([left, right])
    flattened_hotspots = [record for records in hotspots_by_chain.values() for record in records]
    return {
        "template_id": template_id,
        "asset_format": asset_format,
        "hotspots": flattened_hotspots,
        "hotspots_by_chain": hotspots_by_chain,
        "contacts": normalized_contacts,
        "hotspots_sha256": hashlib.sha256(hotspot_path.read_bytes()).hexdigest(),
        "contacts_sha256": hashlib.sha256(contact_path.read_bytes()).hexdigest(),
    }


def _residue_key(value):
    if isinstance(value, (tuple, list)) and len(value) >= 3:
        return str(value[2]).strip(), str(value[1]).upper().strip()
    if isinstance(value, (tuple, list)) and len(value) >= 2:
        return str(value[0]).strip(), str(value[1]).upper().strip()
    if isinstance(value, dict):
        return str(value.get("residue", value.get("number", ""))).strip(), str(value.get("name", value.get("resname", ""))).upper().strip()
    text = str(value).strip()
    parts = text.replace(".", " ").split()
    return (parts[-1], parts[-2].upper()) if len(parts) >= 2 else (text, "")


def _match_residue_keys(match_dict):
    return {_residue_key(key): _residue_key(value) for key, value in (match_dict or {}).items()}


def evaluate_hotspots(match_dict, hotspots, criterion=2, minimum=1):
    """Count matching template hotspots using residue number and type."""

    if not isinstance(hotspots, (list, tuple, set)) or not hotspots:
        return FilterDecision(False, "hotspot_assets_missing")
    mapped = _match_residue_keys(match_dict)
    matches = sum(1 for hotspot in hotspots if _residue_key(hotspot) in mapped and mapped[_residue_key(hotspot)] != ("", ""))
    if matches < minimum:
        return FilterDecision(False, "hotspot_threshold_failed", hotspot_matches=matches)
    return FilterDecision(True, "hotspot_threshold_passed", hotspot_matches=matches)


def count_matched_complementary_contacts(left_match, right_match, contacts):
    left_keys = {_template_residue_key(key) for key in (left_match or {})}
    right_keys = {_template_residue_key(key) for key in (right_match or {})}
    count = 0
    for contact in contacts or []:
        if isinstance(contact, dict):
            pair = (contact.get("left", contact.get("a")), contact.get("right", contact.get("b")))
        elif isinstance(contact, (list, tuple)) and len(contact) >= 2:
            pair = (contact[0], contact[1])
        else:
            continue
        if _template_residue_key(pair[0]) in left_keys and _template_residue_key(pair[1]) in right_keys:
            count += 1
    return count


def _template_residue_key(value):
    """Return chain and residue number for a template-side match/contact."""

    if isinstance(value, (tuple, list)) and len(value) >= 3:
        return str(value[0]).strip(), str(value[2]).strip()
    text = str(value).strip()
    parts = text.split(".")
    if len(parts) >= 3:
        return parts[0].strip(), parts[-1].strip()
    if len(parts) == 2:
        return parts[0].strip(), parts[1].strip()
    return "", text


def evaluate_protocol_candidate(left_match, right_match, left_hotspots, right_hotspots, contacts, *, minimum_contacts=5):
    left_hotspot = evaluate_hotspots(left_match, left_hotspots)
    right_hotspot = evaluate_hotspots(right_match, right_hotspots)
    contact_count = count_matched_complementary_contacts(left_match, right_match, contacts)
    if not left_hotspot.passed or not right_hotspot.passed:
        return FilterDecision(False, "hotspot_threshold_failed", left_hotspot.hotspot_matches + right_hotspot.hotspot_matches, contact_count)
    if contact_count < minimum_contacts:
        return FilterDecision(False, "complementary_contact_threshold_failed", left_hotspot.hotspot_matches + right_hotspot.hotspot_matches, contact_count)
    return FilterDecision(True, "published_protocol_passed", left_hotspot.hotspot_matches + right_hotspot.hotspot_matches, contact_count)
