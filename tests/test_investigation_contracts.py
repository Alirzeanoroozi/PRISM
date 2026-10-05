import csv
import json

import pytest

from benchmark.scripts.investigation_contracts import (
    CONTRACT_SCORE_FIELDS,
    FOREIGN_KEY_FIELDS,
    STRUCTURAL_METRIC_FIELDS,
    ContractViolation,
    evaluate_dockq_with_frozen_mapping,
    require_foreign_keys,
    select_top_k_candidates,
    summarize_top20,
    validate_native_independent_ranking_keys,
    write_contract_scores_tsv,
)
from benchmark.scripts.investigation_artifacts import write_contract_pair_summary_tsv


def _foreign_keys(**overrides):
    values = {
        "cohort": "cohort-1",
        "dataset_row_id": "row-1",
        "arm": "current",
        "experiment": "exp-1",
        "attempt": "attempt-1",
        "pose": "pose-1",
        "mapping": "caller-mapping-is-replaced",
        "native": "native-1",
        "evaluator": "dockq-2.1.3",
    }
    values.update(overrides)
    return values


def _mapping():
    return {
        "chain_mapping": {"B": "L", "A": "H"},
        "residue_correspondence_complete": True,
        "residue_correspondence": [
            {
                "model_chain": "A",
                "model_number": 1,
                "model_name": "ALA",
                "native_chain": "H",
                "native_number": 1,
                "native_name": "ALA",
            },
            {
                "model_chain": "B",
                "model_number": 1,
                "model_name": "GLY",
                "native_chain": "L",
                "native_number": 1,
                "native_name": "GLY",
            },
        ],
    }


def _dockq_payload(global_dockq=0.7):
    return {
        "GlobalDockQ": global_dockq,
        "best_result": [
            {"interface": "A:B", "DockQ": global_dockq, "iRMSD": 1.2, "LRMSD": 2.3, "fnat": 0.4, "F1": 0.5, "clashes": 0}
        ],
    }


def _candidate(pose, score, *, model_valid=True, global_dockq=None):
    row = _foreign_keys(pose=pose, dataset_row_id=f"row-{pose}")
    row.update({
        "model_valid": model_valid,
        "tm_score_left": score,
        "tm_score_right": score,
        "match_count_left": 20,
        "match_count_right": 20,
    })
    if global_dockq is not None:
        row["GlobalDockQ"] = global_dockq
    return row


def test_contract_requires_all_foreign_key_columns():
    assert tuple(require_foreign_keys(_foreign_keys()))[0:9] == FOREIGN_KEY_FIELDS
    with pytest.raises(ContractViolation, match="mapping"):
        require_foreign_keys({key: value for key, value in _foreign_keys().items() if key != "mapping"})


def test_mapping_is_canonical_and_frozen_before_score_values_are_read():
    valid_keys = {key: value for key, value in _foreign_keys().items() if key != "mapping"}
    records = evaluate_dockq_with_frozen_mapping(
        _dockq_payload(), mapping=_mapping(), foreign_keys=valid_keys
    )
    assert records.frozen_mapping.serialized == "A:H;B:L"
    assert records[0]["mapping"] == records.frozen_mapping.mapping_id
    assert records[0]["chain_mapping"] == "A:H;B:L"
    assert all(set(record) >= set(FOREIGN_KEY_FIELDS) for record in records)

    with pytest.raises(ContractViolation, match="valid before score evaluation"):
        evaluate_dockq_with_frozen_mapping(
            {"malformed": True},
            mapping={"chain_mapping": {"A": "H"}},
            foreign_keys=valid_keys,
        )


def test_no_valid_model_nulls_structural_metrics_and_only_utility_is_zero():
    summary = summarize_top20(
        [_candidate("pose-failed", 0.99, model_valid=False, global_dockq=0.99)],
        foreign_keys=_foreign_keys(dataset_row_id="row-pose-failed", pose="pose-failed"),
    )
    assert summary["best_GlobalDockQ_at_20"] == 0.0
    assert all(summary[field] is None for field in STRUCTURAL_METRIC_FIELDS)
    assert summary["model_valid_count"] == 0
    assert summary["selected_count"] == 0


def test_pair_summary_uses_native_independent_top20_and_nulls_unscored_metrics(tmp_path):
    rows = []
    for pose, ranking, dockq in (("pose-a", 0.8, 0.1), ("pose-b", 0.9, 0.7)):
        row = _candidate(pose, ranking, global_dockq=dockq)
        row["dataset_row_id"] = "row-1"
        rows.append(row)
    summary = summarize_top20(rows, foreign_keys=_foreign_keys(dataset_row_id="row-1"))
    assert summary["best_GlobalDockQ_at_20"] == 0.7
    assert summary["model_valid_count"] == 2
    assert summary["scoreable_model_count"] == 2
    assert summary["GlobalDockQ"] == 0.7

    path = write_contract_pair_summary_tsv(
        tmp_path / "pair_summary.tsv",
        rows,
        foreign_keys=_foreign_keys(dataset_row_id="row-1"),
    )
    with path.open(newline="") as handle:
        output = next(csv.DictReader(handle, delimiter="\t"))
    assert output["best_GlobalDockQ_at_20"] == "0.7"
    assert "pose" not in output


def test_summary_rejects_mixed_foreign_key_context():
    first = _candidate("pose-a", 0.9)
    second = _candidate("pose-b", 0.8)
    with pytest.raises(ContractViolation, match="summary context"):
        summarize_top20([first, second], foreign_keys={**_foreign_keys(), "dataset_row_id": "row-pose-a"})


def test_top20_selection_is_native_independent_and_deterministic():
    rows = [_candidate(f"pose-{index:02d}", 1.0 - index / 100.0, global_dockq=index / 100.0) for index in range(21)]
    reordered = list(reversed(rows))
    reordered[0]["GlobalDockQ"] = 1.0
    reordered[-1]["GlobalDockQ"] = 0.0

    first = select_top_k_candidates(rows)
    second = select_top_k_candidates(reordered)
    assert [row["pose"] for row in first] == [row["pose"] for row in second]
    assert len(first) == 20
    with pytest.raises(ContractViolation, match="native-derived"):
        validate_native_independent_ranking_keys(["tm_score_left", "GlobalDockQ"])


def test_contract_score_writer_populates_foreign_keys(tmp_path):
    global_path, interface_path = write_contract_scores_tsv(
        _dockq_payload(),
        tmp_path / "global.tsv",
        tmp_path / "interfaces.tsv",
        mapping=_mapping(),
        foreign_keys={key: value for key, value in _foreign_keys().items() if key != "mapping"},
    )
    with global_path.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert tuple(row[field] for field in FOREIGN_KEY_FIELDS)[:8] == tuple(
        _foreign_keys()[field] if field != "mapping" else row[field] for field in FOREIGN_KEY_FIELDS[:8]
    )
    assert row["mapping"].startswith("mapping-")
    assert tuple(row) == CONTRACT_SCORE_FIELDS
    assert "A:H;B:L" in interface_path.read_text()


def test_contract_score_writer_verifies_raw_json_hash(tmp_path):
    payload_path = tmp_path / "dockq.json"
    payload_path.write_text(json.dumps(_dockq_payload()) + "\n", encoding="utf-8")
    global_path, _ = write_contract_scores_tsv(
        payload_path,
        tmp_path / "global.tsv",
        tmp_path / "interfaces.tsv",
        mapping=_mapping(),
        foreign_keys={key: value for key, value in _foreign_keys().items() if key != "mapping"},
    )
    with global_path.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["raw_json_path"] == str(payload_path)
    assert len(row["raw_json_sha256"]) == 64
