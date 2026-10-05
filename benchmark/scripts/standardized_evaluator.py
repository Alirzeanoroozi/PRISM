"""Dependency-light normalization of DockQ 2.1.3 JSON output.

This module deliberately does not execute DockQ.  It turns an already-produced
JSON document into stable global and interface records, retaining provenance
needed to compare historical and current scoring paths.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path


DOCKQ_METRIC_FIELDS = ("DockQ", "iRMSD", "LRMSD", "fnat", "F1", "clashes")
METRIC_RANGES = {
    "GlobalDockQ": (0.0, 1.0),
    "DockQ": (0.0, 1.0),
    "iRMSD": (0.0, None),
    "LRMSD": (0.0, None),
    "fnat": (0.0, 1.0),
    "F1": (0.0, 1.0),
    "clashes": (0.0, None),
    "grouped_iRMSD": (0.0, None),
}
GLOBAL_TSV_FIELDS = (
    "record_type",
    "interface",
    "GlobalDockQ",
    *DOCKQ_METRIC_FIELDS,
    "grouped_iRMSD",
    "raw_json_sha256",
)
INTERFACE_TSV_FIELDS = GLOBAL_TSV_FIELDS


class DockQJSONError(ValueError):
    """Raised when a document is not a supported DockQ JSON record."""


class EvaluationRecords(list):
    """Ordered records with the retained source document attached."""

    def __init__(self, records, *, raw_json, raw_json_sha256):
        super().__init__(records)
        self.raw_json = raw_json
        self.raw_json_sha256 = raw_json_sha256


@dataclass(frozen=True)
class MappingValidation:
    """Result of validating a frozen model-to-native residue mapping."""

    valid: bool
    errors: tuple[str, ...]
    normalized_chain_mapping: tuple[tuple[str, str], ...]
    symmetric_chain_equivalence_declared: bool = False

    def __bool__(self):
        return self.valid

    def as_dict(self):
        return {
            "valid": self.valid,
            "errors": list(self.errors),
            "normalized_chain_mapping": dict(self.normalized_chain_mapping),
            "symmetric_chain_equivalence_declared": self.symmetric_chain_equivalence_declared,
        }

    def __getitem__(self, key):
        return self.as_dict()[key]

    def get(self, key, default=None):
        return self.as_dict().get(key, default)


def validate_raw_pdb_chain_contract(model_pdb, model_receptor, model_ligand):
    """Reject raw model files that cannot represent the declared partners.

    Biopython may collapse repeated chain/residue identifiers while parsing a
    PDB.  This guard intentionally inspects the fixed-width records first, so
    a legacy file containing two partner segments under one chain ID is
    rejected instead of being scored as a one-chain model.
    """

    def declared_chains(value):
        if isinstance(value, str):
            tokens = value.replace(",", " ").split()
            if len(tokens) == 1:
                return list(tokens[0])
            return [token for token in tokens if token]
        return [str(item) for item in value or []]

    receptor = declared_chains(model_receptor)
    ligand = declared_chains(model_ligand)
    errors = []
    if not receptor or not ligand:
        errors.append("model receptor and ligand chains must be declared")
    overlap = sorted(set(receptor) & set(ligand))
    if overlap:
        errors.append("model receptor/ligand chains overlap: " + ",".join(overlap))

    observed = set()
    residue_order = {}
    try:
        with open(model_pdb, encoding="ascii", errors="replace") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.startswith("ATOM  ") or len(line) < 27:
                    continue
                chain = line[21].strip() or "_"
                observed.add(chain)
                try:
                    residue_number = int(line[22:26])
                except ValueError:
                    continue
                insertion_code = line[26].strip()
                key = (residue_number, insertion_code)
                order = residue_order.setdefault(chain, [])
                if not order or order[-1] != key:
                    order.append(key)
    except OSError as exc:
        return MappingValidation(False, (f"cannot read model PDB: {exc}",), tuple())

    required = set(receptor + ligand)
    missing = sorted(required - observed)
    if missing:
        errors.append("declared model chains absent from raw PDB: " + ",".join(missing))
    if len(observed) < 2:
        errors.append(f"raw PDB contains fewer than two chains: {sorted(observed)}")

    for chain, order in residue_order.items():
        numbers = [item[0] for item in order]
        if any(later < earlier for earlier, later in zip(numbers, numbers[1:])):
            errors.append(f"raw PDB residue numbering resets within chain {chain}")

    return MappingValidation(
        valid=not errors,
        errors=tuple(dict.fromkeys(errors)),
        normalized_chain_mapping=tuple(),
    )


def _canonical_json_bytes(data):
    return json.dumps(
        data, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _read_json(source):
    if isinstance(source, Mapping):
        data = copy.deepcopy(dict(source))
        raw = _canonical_json_bytes(data)
    else:
        path = Path(source)
        raw = path.read_bytes()
        try:
            data = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise DockQJSONError(f"invalid DockQ JSON file {path}: {exc}") from exc
    if not isinstance(data, Mapping):
        raise DockQJSONError("DockQ JSON root must be an object")
    return dict(data), raw, hashlib.sha256(raw).hexdigest()


def _numeric(value, field, *, required=False):
    if value is None:
        if required:
            raise DockQJSONError(f"missing required DockQ field: {field}")
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise DockQJSONError(f"DockQ field {field} must be numeric")
    if not math.isfinite(value):
        raise DockQJSONError(f"DockQ field {field} must be finite")
    lower, upper = METRIC_RANGES.get(field, (None, None))
    if (lower is not None and value < lower) or (upper is not None and value > upper):
        if upper is None:
            expected = f">={lower:g}"
        elif lower is None:
            expected = f"<={upper:g}"
        else:
            expected = f"between {lower:g} and {upper:g}"
        raise DockQJSONError(f"DockQ field {field} is out of range; expected {expected}")
    return value


def parse_metric_output(output, field):
    """Parse a successful scalar-metric subprocess result fail-closed."""

    if not isinstance(output, str):
        raise DockQJSONError(f"{field} subprocess output must be text")
    value = output.strip()
    if not value:
        raise DockQJSONError(f"{field} subprocess output is empty")
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise DockQJSONError(f"{field} subprocess output must be numeric") from exc
    return _numeric(numeric, field, required=True)


def _interface_name(item, key, index):
    if isinstance(item, Mapping):
        for field in ("interface", "interface_id", "chains", "mapping"):
            if item.get(field) is not None:
                return str(item[field])
        model = item.get("model_chains", item.get("model"))
        native = item.get("native_chains", item.get("native"))
        if model is not None and native is not None:
            return f"{model}:{native}"
    if key is not None:
        return str(key)
    return f"interface_{index + 1}"


def _interface_items(data):
    if "best_result" in data:
        raw_interfaces = data["best_result"]
        field = "best_result"
    elif "interfaces" in data:
        raw_interfaces = data["interfaces"]
        field = "interfaces"
    else:
        raise DockQJSONError("missing DockQ interface field: best_result")

    if isinstance(raw_interfaces, Mapping):
        return [(_interface_name(value, key, index), value) for index, (key, value) in enumerate(raw_interfaces.items())]
    if isinstance(raw_interfaces, list):
        if not all(isinstance(item, Mapping) for item in raw_interfaces):
            raise DockQJSONError(f"{field} must contain interface objects")
        return [(_interface_name(item, None, index), item) for index, item in enumerate(raw_interfaces)]
    raise DockQJSONError(f"{field} must be an object or list")


def evaluate_dockq_json(source, *, grouped_irmsd=None):
    """Return one global record followed by one record for every interface.

    ``source`` is a mapping or a path to a JSON file.  Mapping input is hashed
    from deterministic JSON bytes; file input is hashed from its exact bytes.
    The returned list also exposes ``raw_json`` and ``raw_json_sha256``.
    """

    data, raw, raw_hash = _read_json(source)
    pairwise_cross_fallback = data.get("evaluation_mode") == "pairwise_cross_fallback"
    global_dockq = _numeric(
        data.get("GlobalDockQ"),
        "GlobalDockQ",
        required=not pairwise_cross_fallback,
    )
    interfaces = _interface_items(data)
    if grouped_irmsd is None:
        grouped_irmsd = data.get("grouped_iRMSD")
    grouped_irmsd = _numeric(grouped_irmsd, "grouped_iRMSD")
    interface_metrics = {
        field: _numeric(data.get(field), field)
        for field in ("iRMSD", "LRMSD", "fnat", "F1", "clashes")
    }
    if len(interfaces) > 1:
        interface_metrics = {field: None for field in interface_metrics}

    global_record = {
        "record_type": "global",
        "interface": "",
        "GlobalDockQ": global_dockq,
        "DockQ": global_dockq,
        **interface_metrics,
        "grouped_iRMSD": grouped_irmsd,
        "raw_json_sha256": raw_hash,
    }
    records = [global_record]
    for interface_name, item in interfaces:
        record = {
            "record_type": "interface",
            "interface": interface_name,
            "GlobalDockQ": global_dockq,
            "grouped_iRMSD": None,
            "raw_json_sha256": raw_hash,
        }
        for field in DOCKQ_METRIC_FIELDS:
            record[field] = _numeric(item.get(field), field)
        records.append(record)
    return EvaluationRecords(records, raw_json=data, raw_json_sha256=raw_hash)


def parse_dockq_json(source, *, grouped_irmsd=None):
    """Compatibility name for :func:`evaluate_dockq_json`."""

    return evaluate_dockq_json(source, grouped_irmsd=grouped_irmsd)


def standardize_dockq_json(source, *, grouped_irmsd=None):
    """Explicitly named alias for the standardized evaluator contract."""

    return evaluate_dockq_json(source, grouped_irmsd=grouped_irmsd)


def _chains(value):
    if isinstance(value, str):
        value = value.replace(",", " ").split()
        if len(value) == 1 and len(value[0]) > 1:
            value = list(value[0])
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        value = list(value)
    else:
        return []
    return [str(chain) for chain in value if str(chain)]


def _chain_groups(mapping):
    groups = mapping.get("chain_groups", []) if isinstance(mapping, Mapping) else []
    if not groups and isinstance(mapping, Mapping):
        model = mapping.get("model_chains", mapping.get("model"))
        native = mapping.get("native_chains", mapping.get("native"))
        if model is not None or native is not None:
            groups = [{"model": model, "native": native}]
    result = []
    for group in groups:
        if isinstance(group, Mapping):
            model = group.get("model", group.get("model_chains"))
            native = group.get("native", group.get("native_chains"))
        elif isinstance(group, Sequence) and len(group) == 2:
            model, native = group
        else:
            result.append(([], []))
            continue
        result.append((_chains(model), _chains(native)))
    return result


def _explicit_chain_mapping(mapping, groups):
    declared = mapping.get("chain_mapping", mapping.get("mapping")) if isinstance(mapping, Mapping) else None
    if declared is None and isinstance(mapping, Mapping):
        reserved = {"chain_groups", "model", "native", "model_chains", "native_chains", "residue_correspondence", "residue_mapping", "residues", "residue_correspondence_complete", "symmetric_chain_equivalence", "symmetric_chain_equivalence_declared"}
        candidate = {key: value for key, value in mapping.items() if key not in reserved}
        if candidate and all(isinstance(key, str) and isinstance(value, str) for key, value in candidate.items()):
            declared = candidate
    if isinstance(declared, str):
        if "," not in declared and ":" in declared:
            left, right = declared.split(":", 1)
            if len(left) == len(right) and len(left) > 1:
                return dict(zip(left, right))
        pairs = [part.split(":", 1) for part in declared.split(",") if ":" in part]
        return {left: right for left, right in pairs}
    if isinstance(declared, Mapping):
        return {str(key): str(value) for key, value in declared.items()}
    if isinstance(declared, Sequence) and not isinstance(declared, (bytes, bytearray, str)):
        return {str(pair[0]): str(pair[1]) for pair in declared if isinstance(pair, Sequence) and len(pair) == 2}
    result = {}
    for model, native in groups:
        if len(model) == len(native):
            result.update(zip(model, native))
    return result


def _residue_parts(item):
    if not isinstance(item, Mapping):
        return None, None
    model = item.get("model") if isinstance(item.get("model"), Mapping) else item
    native = item.get("native") if isinstance(item.get("native"), Mapping) else item

    def part(side, names):
        for name in names:
            if f"{side}_{name}" in item:
                return item[f"{side}_{name}"]
            if name in item and side == "model":
                return item[name]
        for name in names:
            if name in model and side == "model":
                return model[name]
            if name in native and side == "native":
                return native[name]
        return None

    return (
        {"chain": part("model", ("chain",)), "number": part("model", ("number", "resnum", "residue_number")), "name": part("model", ("name", "resname", "residue_name"))},
        {"chain": part("native", ("chain",)), "number": part("native", ("number", "resnum", "residue_number")), "name": part("native", ("name", "resname", "residue_name"))},
    )


def _same_number(left, right):
    return str(left).strip() == str(right).strip()


def validate_mapping(
    mapping,
    residue_correspondence=None,
    *,
    symmetric_chain_equivalence=False,
    symmetric_chain_equivalence_declared=None,
    residue_correspondence_complete=None,
):
    """Validate chain groups and residue identity/numbering for ``--no_align``.

    The validator is intentionally structural: it reads only the supplied
    mapping declaration and residue pairs.  It never opens PDB files or runs a
    scorer.  A valid result is the sole input accepted by
    :func:`no_align_is_safe`.
    """

    if isinstance(mapping, str):
        groups = _chain_groups({"model_chains": mapping.split(":", 1)[0], "native_chains": mapping.split(":", 1)[1] if ":" in mapping else None})
        mapping_object = {}
    elif isinstance(mapping, Mapping):
        mapping_object = mapping
        groups = _chain_groups(mapping)
    else:
        mapping_object = {"chain_mapping": mapping}
        groups = []
    if residue_correspondence_complete is None and isinstance(mapping_object, Mapping):
        residue_correspondence_complete = mapping_object.get("residue_correspondence_complete")
    if symmetric_chain_equivalence_declared is None:
        symmetric_chain_equivalence_declared = bool(
            mapping_object.get("symmetric_chain_equivalence_declared", mapping_object.get("symmetric_chain_equivalence", symmetric_chain_equivalence))
        )
    chain_map = _explicit_chain_mapping(mapping_object, groups)
    errors = []
    if residue_correspondence_complete is not True:
        errors.append("complete residue correspondence must be explicitly declared")

    if not chain_map:
        errors.append("missing chain mapping")
    if len(chain_map) != len(set(chain_map)) or len(set(chain_map.values())) != len(chain_map):
        errors.append("chain mapping must be one-to-one")
    for model, native in chain_map.items():
        if not model or not native:
            errors.append("chain mapping contains an empty chain")

    for model_group, native_group in groups:
        if not model_group or not native_group or len(model_group) != len(native_group):
            errors.append("chain groups must be non-empty and have equal size")
            continue
        if len(set(model_group)) != len(model_group) or len(set(native_group)) != len(native_group):
            errors.append("chain groups must not repeat chains")
        mapped = [chain_map.get(chain) for chain in model_group]
        if any(value is None for value in mapped):
            errors.append("chain group contains a chain missing from chain_mapping")
        elif symmetric_chain_equivalence_declared:
            if set(mapped) != set(native_group):
                errors.append("symmetric chain group does not have equivalent members")
        elif mapped != native_group:
            errors.append("chain order differs; declare symmetric chain equivalence")

    if residue_correspondence is None and isinstance(mapping_object, Mapping):
        residue_correspondence = mapping_object.get(
            "residue_correspondence",
            mapping_object.get("residue_mapping", mapping_object.get("residues")),
        )
    if not isinstance(residue_correspondence, Sequence) or isinstance(residue_correspondence, (str, bytes, bytearray)) or not residue_correspondence:
        errors.append("missing residue correspondence")
    else:
        seen_model_residues = set()
        seen_native_residues = set()
        for index, item in enumerate(residue_correspondence):
            model, native = _residue_parts(item)
            if model is None or native is None:
                errors.append(f"residue correspondence {index} is malformed")
                continue
            if model["chain"] not in chain_map:
                errors.append(f"residue correspondence {index} uses unmapped model chain")
            elif chain_map[model["chain"]] != native["chain"]:
                errors.append(f"residue correspondence {index} has inconsistent chain mapping")
            if model["name"] is None or native["name"] is None or str(model["name"]).upper() != str(native["name"]).upper():
                errors.append(f"residue correspondence {index} has different residue identity")
            if model["number"] is None or native["number"] is None or not _same_number(model["number"], native["number"]):
                errors.append(f"residue correspondence {index} has different residue numbering")
            if model["chain"] is not None and model["number"] is not None:
                model_key = (model["chain"], str(model["number"]).strip())
                if model_key in seen_model_residues:
                    errors.append(f"residue correspondence {index} duplicates a model residue")
                seen_model_residues.add(model_key)
            if native["chain"] is not None and native["number"] is not None:
                native_key = (native["chain"], str(native["number"]).strip())
                if native_key in seen_native_residues:
                    errors.append(f"residue correspondence {index} duplicates a native residue")
                seen_native_residues.add(native_key)

        covered_model_chains = set()
        covered_native_chains = set()
        for item in residue_correspondence:
            model, native = _residue_parts(item)
            if model and native:
                covered_model_chains.add(model.get("chain"))
                covered_native_chains.add(native.get("chain"))
        for model_chain, native_chain in chain_map.items():
            if model_chain not in covered_model_chains:
                errors.append(f"chain {model_chain} has no residue correspondence")
            if native_chain not in covered_native_chains:
                errors.append(f"native chain {native_chain} has no residue correspondence")

    return MappingValidation(
        valid=not errors,
        errors=tuple(dict.fromkeys(errors)),
        normalized_chain_mapping=tuple(sorted(chain_map.items())),
        symmetric_chain_equivalence_declared=bool(symmetric_chain_equivalence_declared),
    )


def validate_chain_mapping(*args, **kwargs):
    """Compatibility alias for :func:`validate_mapping`."""

    return validate_mapping(*args, **kwargs)


def validate_pdb_mapping(
    model_pdb,
    native_pdb,
    chain_mapping,
    *,
    symmetric_chain_equivalence=False,
):
    """Validate chain, residue identity, and numbering correspondence in PDBs.

    This is the concrete precondition for using DockQ ``--no_align``.  The
    comparison is intentionally strict: standard residues must exist on both
    mapped chains at the same sequence number and have the same residue name.
    Gaps, insertions, or identity changes therefore remain visible as mapping
    failures rather than being silently aligned by DockQ.
    """

    try:
        from Bio.PDB import PDBParser
    except ImportError as exc:  # pragma: no cover - exercised by env checks
        raise RuntimeError("validate_pdb_mapping requires Biopython") from exc

    parser = PDBParser(QUIET=True)
    model_structure = parser.get_structure("model", str(model_pdb))[0]
    native_structure = parser.get_structure("native", str(native_pdb))[0]
    mapping = _explicit_chain_mapping({}, []) if not isinstance(chain_mapping, Mapping) else {
        str(key): str(value) for key, value in chain_mapping.items()
    }
    if not mapping and isinstance(chain_mapping, str):
        mapping = _explicit_chain_mapping({"mapping": chain_mapping}, [])

    correspondence = []
    errors = []
    for model_chain, native_chain in sorted(mapping.items()):
        if model_chain not in model_structure or native_chain not in native_structure:
            errors.append(f"missing mapped chain {model_chain}:{native_chain}")
            continue

        def standard_residues(chain):
            return {
                (residue.id[1], str(residue.id[2]).strip()): residue
                for residue in chain
                if residue.id[0] == " "
            }

        model_residues = standard_residues(model_structure[model_chain])
        native_residues = standard_residues(native_structure[native_chain])
        for key in sorted(set(model_residues) | set(native_residues)):
            model_residue = model_residues.get(key)
            native_residue = native_residues.get(key)
            if model_residue is None or native_residue is None:
                errors.append(f"missing residue correspondence {model_chain}:{native_chain}:{key}")
                continue
            correspondence.append(
                {
                    "model_chain": model_chain,
                    "native_chain": native_chain,
                    "model_number": f"{key[0]}{key[1]}",
                    "native_number": f"{key[0]}{key[1]}",
                    "model_name": model_residue.resname,
                    "native_name": native_residue.resname,
                }
            )

    covered_model_chains = {item["model_chain"] for item in correspondence}
    covered_native_chains = {item["native_chain"] for item in correspondence}
    for model_chain, native_chain in mapping.items():
        if model_chain not in covered_model_chains:
            errors.append(f"chain {model_chain} has no residue correspondence")
        if native_chain not in covered_native_chains:
            errors.append(f"native chain {native_chain} has no residue correspondence")

    result = validate_mapping(
        {"chain_mapping": mapping},
        correspondence,
        symmetric_chain_equivalence=symmetric_chain_equivalence,
        residue_correspondence_complete=True,
    )
    if errors:
        return MappingValidation(
            valid=False,
            errors=tuple(dict.fromkeys((*result.errors, *errors))),
            normalized_chain_mapping=result.normalized_chain_mapping,
            symmetric_chain_equivalence_declared=result.symmetric_chain_equivalence_declared,
        )
    return result


def no_align_is_safe(mapping_validation):
    """Return whether the validated declaration permits DockQ ``--no_align``."""

    if isinstance(mapping_validation, MappingValidation):
        return mapping_validation.valid
    return False


def is_no_align_safe(mapping_validation):
    """Compatibility alias for :func:`no_align_is_safe`."""

    return no_align_is_safe(mapping_validation)


def can_use_no_align(mapping_validation):
    """Compatibility alias for :func:`no_align_is_safe`."""

    return no_align_is_safe(mapping_validation)
