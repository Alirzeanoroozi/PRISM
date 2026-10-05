import hashlib
import json

import pytest

from benchmark.scripts.standardized_evaluator import (
    GLOBAL_TSV_FIELDS,
    INTERFACE_TSV_FIELDS,
    evaluate_dockq_json,
    no_align_is_safe,
    parse_metric_output,
    validate_mapping,
    validate_pdb_mapping,
)


def _write_pdb(path, chain, residues, *, offset=(0.0, 0.0, 0.0), reverse_x=False):
    lines = []
    for serial, (number, name) in enumerate(residues, 1):
        x = float(serial) * (-1.0 if reverse_x else 1.0) + offset[0]
        y = offset[1]
        z = offset[2]
        lines.append(
            f"ATOM  {serial:5d}  CA  {name:>3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  "
        )
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def _interface(interface, dockq, irmsd):
    return {
        "interface": interface,
        "DockQ": dockq,
        "iRMSD": irmsd,
        "LRMSD": irmsd + 1.0,
        "fnat": 0.5,
        "F1": 0.6,
        "clashes": 2,
    }


def test_evaluator_keeps_global_score_and_repeated_interfaces_separate():
    payload = {
        "GlobalDockQ": 0.66,
        "best_result": [
            _interface("H,L:H,L", 0.41, 1.2),
            _interface("H,L:H,L", 0.22, 3.4),
        ],
    }

    records = evaluate_dockq_json(payload, grouped_irmsd=1.679)

    assert [record["record_type"] for record in records] == ["global", "interface", "interface"]
    assert records[0]["GlobalDockQ"] == 0.66
    assert records[0]["DockQ"] == 0.66
    assert [record["DockQ"] for record in records[1:]] == [0.41, 0.22]
    assert [record["iRMSD"] for record in records[1:]] == [1.2, 3.4]
    assert records[0]["grouped_iRMSD"] == 1.679
    assert all(record["raw_json_sha256"] == records[0]["raw_json_sha256"] for record in records)


def test_evaluator_hashes_exact_json_file_bytes(tmp_path):
    raw = b'{\n  "GlobalDockQ": 0.5,\n  "best_result": []\n}\n'
    source = tmp_path / "dockq.json"
    source.write_bytes(raw)

    records = evaluate_dockq_json(source)

    assert records[0]["raw_json_sha256"] == hashlib.sha256(raw).hexdigest()


def test_missing_interface_structural_metric_remains_null_not_zero():
    records = evaluate_dockq_json(
        {
            "GlobalDockQ": 0.4,
            "best_result": [{"interface": "A:B", "DockQ": 0.4, "LRMSD": 2.0, "fnat": 0.2}],
        }
    )
    assert records[0]["iRMSD"] is None
    assert records[0]["clashes"] is None


def test_multi_interface_metrics_stay_per_interface_and_global_fields_are_ambiguous():
    payload = {
        "GlobalDockQ": 0.66,
        "iRMSD": 1.1,
        "LRMSD": 2.2,
        "fnat": 0.3,
        "F1": 0.4,
        "clashes": 2,
        "best_result": [
            _interface("A:B", 0.41, 1.2),
            _interface("C:D", 0.22, 3.4),
        ],
    }

    records = evaluate_dockq_json(payload)

    assert records[0]["record_type"] == "global"
    assert all(records[0][field] is None for field in ("iRMSD", "LRMSD", "fnat", "F1", "clashes"))
    assert [(record["interface"], record["iRMSD"]) for record in records[1:]] == [
        ("A:B", 1.2),
        ("C:D", 3.4),
    ]


def test_metric_ranges_are_validated():
    invalid_payloads = [
        {"GlobalDockQ": -0.01, "best_result": []},
        {"GlobalDockQ": 1.01, "best_result": []},
        {"GlobalDockQ": 0.5, "best_result": [_interface("A:B", 0.5, -1.0)]},
        {"GlobalDockQ": 0.5, "best_result": [_interface("A:B", 0.5, 1.0) | {"fnat": 1.1}]},
        {"GlobalDockQ": 0.5, "best_result": [_interface("A:B", 0.5, 1.0) | {"clashes": -1}]},
    ]

    for payload in invalid_payloads:
        with pytest.raises(ValueError, match="range"):
            evaluate_dockq_json(payload)


@pytest.mark.parametrize("output", ["", "not-a-number", "1.0 Angstroms"])
def test_invalid_metric_subprocess_output_fails_closed(output):
    with pytest.raises(ValueError, match="iRMSD"):
        parse_metric_output(output, "iRMSD")


def test_metric_subprocess_output_accepts_one_finite_number():
    assert parse_metric_output(" 1.25\n", "iRMSD") == 1.25


@pytest.mark.parametrize(
    "payload, message",
    [
        ({"best_result": []}, "GlobalDockQ"),
        ({"GlobalDockQ": None, "best_result": []}, "GlobalDockQ"),
        ({"GlobalDockQ": 0.5, "best_result": "not a result list"}, "best_result"),
    ],
)
def test_evaluator_rejects_malformed_or_missing_global_fields(payload, message):
    with pytest.raises(ValueError, match=message):
        evaluate_dockq_json(payload)


def test_mapping_validation_checks_chain_groups_residue_identity_and_numbering():
    mapping = {
        "chain_groups": [{"model": "AB", "native": "HL"}],
        "chain_mapping": {"A": "H", "B": "L"},
        "residue_correspondence_complete": True,
        "residue_correspondence": [
            {
                "model_chain": "A",
                "model_number": 10,
                "model_name": "ALA",
                "native_chain": "H",
                "native_number": 10,
                "native_name": "ALA",
            },
            {
                "model_chain": "B",
                "model_number": 10,
                "model_name": "GLY",
                "native_chain": "L",
                "native_number": 10,
                "native_name": "GLY",
            },
        ],
    }

    result = validate_mapping(mapping)

    assert result.valid
    assert result.errors == ()
    assert no_align_is_safe(result)


def test_symmetric_chain_equivalence_requires_explicit_declaration():
    mapping = {
        "chain_groups": [{"model": "AB", "native": "HL"}],
        "chain_mapping": {"A": "L", "B": "H"},
        "residue_correspondence_complete": True,
        "residue_correspondence": [
            {
                "model_chain": "A",
                "model_number": 10,
                "model_name": "ALA",
                "native_chain": "L",
                "native_number": 10,
                "native_name": "ALA",
            },
            {
                "model_chain": "B",
                "model_number": 10,
                "model_name": "GLY",
                "native_chain": "H",
                "native_number": 10,
                "native_name": "GLY",
            },
        ],
    }

    undeclared = validate_mapping(mapping)
    declared = validate_mapping(mapping, symmetric_chain_equivalence=True)

    assert not undeclared.valid
    assert not no_align_is_safe(undeclared)
    assert declared.valid
    assert declared.symmetric_chain_equivalence_declared
    assert no_align_is_safe(declared)


def test_sparse_residue_correspondence_does_not_authorize_no_align():
    mapping = {
        "chain_mapping": {"A": "H"},
        "residue_correspondence": [
            {
                "model_chain": "A",
                "model_number": 1,
                "model_name": "ALA",
                "native_chain": "H",
                "native_number": 1,
                "native_name": "ALA",
            }
        ],
    }

    result = validate_mapping(mapping)

    assert not result.valid
    assert any("complete" in error for error in result.errors)
    assert not no_align_is_safe(result)

    mapping["residue_correspondence_complete"] = True
    assert no_align_is_safe(validate_mapping(mapping))


def test_no_align_gate_rejects_identity_or_numbering_mismatch():
    mapping = {
        "chain_mapping": {"A": "H"},
        "residue_correspondence": [
            {
                "model_chain": "A",
                "model_number": 10,
                "model_name": "ALA",
                "native_chain": "H",
                "native_number": 11,
                "native_name": "GLY",
            }
        ],
    }

    result = validate_mapping(mapping)

    assert not result.valid
    assert any("identity" in error or "number" in error for error in result.errors)
    assert not no_align_is_safe(result)
    assert not no_align_is_safe({"valid": True})


def test_pdb_mapping_validator_checks_residue_identity_and_numbering(tmp_path):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native.pdb"
    _write_pdb(model, "A", [(1, "ALA"), (2, "GLY")])
    _write_pdb(native, "B", [(1, "ALA"), (2, "GLY")])
    valid = validate_pdb_mapping(model, native, {"A": "B"})
    assert valid.valid
    assert no_align_is_safe(valid)

    _write_pdb(native, "B", [(1, "ALA"), (3, "GLY")])
    invalid = validate_pdb_mapping(model, native, {"A": "B"})
    assert not invalid.valid
    assert any("missing residue correspondence" in error for error in invalid.errors)


def test_pdb_mapping_contract_is_invariant_to_rigid_coordinate_transform(tmp_path):
    model = tmp_path / "model.pdb"
    native = tmp_path / "native-transformed.pdb"
    residues = [(1, "ALA"), (2, "GLY"), (3, "SER")]
    _write_pdb(model, "A", residues)
    _write_pdb(native, "B", residues, offset=(37.0, -11.0, 4.5), reverse_x=True)

    result = validate_pdb_mapping(model, native, {"A": "B"})

    assert result.valid
    assert no_align_is_safe(result)


def test_tsv_field_lists_are_stable_and_cover_required_metrics():
    assert GLOBAL_TSV_FIELDS == (
        "record_type",
        "interface",
        "GlobalDockQ",
        "DockQ",
        "iRMSD",
        "LRMSD",
        "fnat",
        "F1",
        "clashes",
        "grouped_iRMSD",
        "raw_json_sha256",
    )
    assert INTERFACE_TSV_FIELDS == (
        "record_type",
        "interface",
        "GlobalDockQ",
        "DockQ",
        "iRMSD",
        "LRMSD",
        "fnat",
        "F1",
        "clashes",
        "grouped_iRMSD",
        "raw_json_sha256",
    )


def test_pairwise_cross_fallback_allows_explicitly_unavailable_global_score(tmp_path):
    path = tmp_path / "fallback.json"
    path.write_text(json.dumps({
        "evaluation_mode": "pairwise_cross_fallback",
        "GlobalDockQ": None,
        "best_result": {
            "BA": {"DockQ": 0.4, "iRMSD": 2.0, "LRMSD": 3.0,
                   "fnat": 0.3, "F1": 0.4, "clashes": 0},
        },
    }))

    records = evaluate_dockq_json(path)

    assert records[0]["record_type"] == "global"
    assert records[0]["GlobalDockQ"] is None
    assert records[1]["DockQ"] == 0.4
