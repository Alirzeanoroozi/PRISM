from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from scripts.build_gtalign_clusters import (
    ClusterAccountingError,
    PanelHashMismatch,
    build_command,
    build_index,
    load_panel_members,
    parse_gtalign_clusters,
    select_representatives,
    validate_panel_hash,
)


def _member(member_id: str, *, template_id: str | None = None, chain_id: str = "A", path: str | None = None):
    return {
        "member_id": member_id,
        "template_id": template_id or member_id,
        "chain_id": chain_id,
        "path": path or f"/panel/{member_id}.pdb",
    }


def test_build_command_captures_required_gtalign_arguments():
    command = build_command(
        input_dir="/panel",
        raw_output="/run/raw",
        cache_dir="/run/cache",
        threshold=0.62,
        coverage=0.81,
        algorithm=1,
        speed=13,
        gtalign_path="/opt/bin/gtalign",
    )
    assert command == [
        "/opt/bin/gtalign",
        "--cls=/panel",
        "-o",
        "/run/raw",
        "-c",
        "/run/cache",
        "--cls-threshold=0.62",
        "--cls-coverage=0.81",
        "--cls-algorithm=1",
        "--speed=13",
    ]


def test_build_command_captures_sensitivity_and_input_eligibility_controls():
    command = build_command(
        input_dir="/panel",
        raw_output="/run/raw",
        cache_dir="/run/cache",
        threshold=0.6,
        coverage=0.8,
        algorithm=0,
        speed=9,
        pre_score=0.0,
        min_length=3,
        cpu_threads_reading=8,
        sort=3,
        gtalign_path="/opt/bin/gtalign",
    )
    assert command[-4:] == [
        "--pre-score=0",
        "--dev-min-length=3",
        "--cpu-threads-reading=8",
        "--sort=3",
    ]


def test_parse_actual_gtalign_cluster_list_preserves_singletons():
    text = """gtalign 0.19.00\n\n===============================\n\n1wte_A_1 2pkh_A\n3iji_A\n"""
    declarations = parse_gtalign_clusters(text)
    assert [list(d.members) for d in declarations] == [["1wte_A_1", "2pkh_A"], ["3iji_A"]]
    assert declarations[1].declared_status == "SINGLETON"


def test_parse_explicit_empty_and_json_cluster_states():
    declarations = parse_gtalign_clusters("cluster-empty:\ncluster-two:\tfoo bar\n")
    assert list(declarations[0].members) == []
    assert declarations[0].declared_status == "EMPTY_DECLARED"

    payload = {"clusters": [{"cluster_id": "c1", "members": []}, {"cluster_id": "c2", "members": ["x"]}]}
    declarations = parse_gtalign_clusters(json.dumps(payload))
    assert [(d.cluster_id, list(d.members)) for d in declarations] == [("c1", []), ("c2", ["x"])]


def test_panel_hash_mismatch_fails_closed(tmp_path):
    panel = tmp_path / "panel.tsv"
    panel.write_text("a\tA\t/panel/a.pdb\n")
    actual = hashlib.sha256(panel.read_bytes()).hexdigest()
    assert validate_panel_hash(panel, actual) == actual
    with pytest.raises(PanelHashMismatch):
        validate_panel_hash(panel, "0" * 64)


def test_json_panel_manifest_accepts_gtalign_repeated_chain_suffix(tmp_path: Path):
    structure = tmp_path / "1lyqAB_A.pdb"
    structure.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C\nEND\n"
    )
    panel = tmp_path / "panel.json"
    panel.write_text(
        json.dumps(
            [{
                "member_id": "1lyqAB:A",
                "template_id": "1lyqAB",
                "chain_id": "A",
                "path": str(structure),
            }]
        )
    )
    members, _ = load_panel_members(tmp_path, panel)
    assert "1lyqAB_A_A" in members[0]["aliases"]


def test_json_panel_manifest_accepts_staging_int_repeated_chain_suffix(tmp_path: Path):
    structure = tmp_path / "2axtAI_I_int.pdb"
    structure.write_text(
        "ATOM      1  CA  ALA I   1       0.000   0.000   0.000  1.00 10.00           C\nEND\n"
    )
    panel = tmp_path / "panel.json"
    panel.write_text(
        json.dumps(
            [{
                "member_id": "2axtAI:I",
                "template_id": "2axtAI",
                "chain_id": "I",
                "path": str(structure),
            }]
        )
    )
    members, _ = load_panel_members(tmp_path, panel)
    assert "2axtAI_I_I" in members[0]["aliases"]


def test_duplicate_and_missing_members_fail_closed():
    panel = [_member("a"), _member("b")]
    with pytest.raises(ClusterAccountingError, match="duplicate"):
        build_index(
            panel_members=panel,
            declarations=parse_gtalign_clusters("a b\na\n"),
            panel_sha256="h",
            tool_version="0.19.00",
            parameters={},
        )
    with pytest.raises(ClusterAccountingError, match="missing"):
        build_index(
            panel_members=panel,
            declarations=parse_gtalign_clusters("a\n"),
            panel_sha256="h",
            tool_version="0.19.00",
            parameters={},
        )


def test_missing_member_singleton_policy_preserves_complete_accounting():
    panel = [_member("a"), _member("b")]
    index = build_index(
        panel_members=panel,
        declarations=parse_gtalign_clusters("a\n"),
        panel_sha256="h",
        tool_version="0.19.00",
        parameters={},
        missing_member_policy="singleton",
    )
    fallback = next(c for c in index["clusters"] if c["cluster_id"] == "fallback-b")
    assert fallback["status"] == "SINGLETON"
    assert fallback["metadata"]["index_membership_source"] == "singleton_fallback"
    assert index["parameters"]["missing_member_fallback_count"] == 1


def test_representatives_are_deterministic_and_native_first():
    members = [_member("a"), _member("b"), _member("c")]
    metadata = {
        "a": {"family": "f1", "species": "s1", "role": "native", "is_native": True},
        "b": {"family": "f2", "species": "s1", "role": "decoy"},
        "c": {"family": "f1", "species": "s2", "role": "decoy"},
    }
    similarity = {
        ("a", "b"): 0.2,
        ("a", "c"): 0.4,
        ("b", "c"): 0.3,
    }
    first = select_representatives(members, similarity=similarity, metadata=metadata, coverage_threshold=0.8)
    second = select_representatives(list(reversed(members)), similarity=similarity, metadata=metadata, coverage_threshold=0.8)
    assert first == second
    assert [item["member_id"] for item in first.representatives] == ["a", "b", "c"]
    assert first.metadata["member_metadata"]["a"]["family"] == "f1"
    assert "inferred_family" not in first.metadata


def test_singleton_and_empty_cluster_states_are_preserved():
    singleton = build_index(
        panel_members=[_member("a")],
        declarations=parse_gtalign_clusters("a\n"),
        panel_sha256="h",
        tool_version="0.19.00",
        parameters={},
    )
    assert singleton["clusters"][0]["status"] == "SINGLETON"
    assert len(singleton["clusters"][0]["representatives"]) == 1

    empty = build_index(
        panel_members=[],
        declarations=parse_gtalign_clusters("cluster-empty:\n"),
        panel_sha256="h",
        tool_version="0.19.00",
        parameters={},
    )
    assert empty["clusters"][0]["status"] == "EMPTY_DECLARED"
    assert empty["clusters"][0]["members"] == []
    assert empty["clusters"][0]["representatives"] == []


def test_adaptive_representatives_mark_split_required_after_three():
    members = [_member(name) for name in ["a", "b", "c", "d"]]
    similarity = {(left, right): 0.1 for left in "abcd" for right in "abcd" if left < right}
    metadata = {name: {"family": name} for name in "abcd"}
    result = select_representatives(members, similarity=similarity, metadata=metadata, coverage_threshold=0.8)
    assert len(result.representatives) == 3
    assert result.status == "SPLIT_REQUIRED"


def test_shared_index_contract_has_required_top_level_fields():
    index = build_index(
        panel_members=[_member("a")],
        declarations=parse_gtalign_clusters("a\n"),
        panel_sha256="panel-hash",
        tool_version="0.19.00",
        parameters={"threshold": 0.5},
        provenance={"command": ["gtalign"]},
    )
    assert set(["schema_version", "tool", "tool_version", "panel_sha256", "parameters", "clusters"]).issubset(index)
    assert index["schema_version"] == "prism-template-cluster-index/v1"
    assert index["tool"] == "GTalign"
    assert index["clusters"][0]["members"][0].keys() >= {"template_id", "chain_id", "path"}
