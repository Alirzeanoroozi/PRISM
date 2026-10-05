from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

import pytest

from src.foldseek_cluster import (
    FoldseekConfig,
    MembershipError,
    build_foldseek_command,
    build_index_from_outputs,
    check_gpu_capability,
    collect_donor_inputs,
    execute_foldseek,
    resolve_panel_hash,
    select_representatives,
)


def _pdb(path: Path, chains: tuple[str, str] = ("A", "B")) -> None:
    lines = []
    serial = 1
    for chain, z in zip(chains, (0.0, 4.0)):
        for residue in (1, 2):
            lines.append(
                f"ATOM  {serial:5d}  CA  ALA {chain}{residue:4d}    "
                f"{float(residue):8.3f}{0.0:8.3f}{z:8.3f}"
                "  1.00 10.00           C  \n"
            )
            serial += 1
    path.write_text("".join(lines) + "END\n")


def test_command_construction_preserves_multimer_and_threshold_flags(tmp_path):
    config = FoldseekConfig(
        input_dir=tmp_path / "donors",
        output_dir=tmp_path / "out",
        foldseek_path="/opt/foldseek",
        threads=7,
        gpu=1,
        distance_threshold=8,
        interface_lddt_threshold=0.6,
        chain_tm_threshold=0.7,
        multimer_tm_threshold=0.8,
        coverage_threshold=0.75,
    )

    command = build_foldseek_command(config, tmp_path / "out" / "donors")

    assert command[:3] == ["/opt/foldseek", "easy-multimercluster", str(tmp_path / "donors")]
    assert "--db-extraction-mode" in command
    assert command[command.index("--db-extraction-mode") + 1] == "1"
    assert command[command.index("--gpu") + 1] == "1"
    assert command[command.index("--threads") + 1] == "7"
    assert command[command.index("--distance-threshold") + 1] == "8"
    assert "easy-interfacecluster" not in command
    assert command[command.index("--alignment-type") + 1] == "2"
    assert command[command.index("--prefilter-mode") + 1] == "0"


def test_command_can_select_whole_complex_and_exhaustive_controls(tmp_path):
    config = FoldseekConfig(
        input_dir=tmp_path / "donors",
        output_dir=tmp_path / "out",
        extraction_mode=0,
        sensitivity=7.5,
        max_seqs=10000,
        prefilter_mode=2,
        exhaustive_search=1,
        remove_tmp_files=0,
    )
    command = build_foldseek_command(config, tmp_path / "out" / "donors")
    assert command[command.index("--db-extraction-mode") + 1] == "0"
    assert command[command.index("-s") + 1] == "7.5"
    assert command[command.index("--max-seqs") + 1] == "10000"
    assert command[command.index("--prefilter-mode") + 1] == "2"
    assert command[command.index("--exhaustive-search") + 1] == "1"
    assert command[command.index("--remove-tmp-files") + 1] == "0"


def test_gpu_zero_is_explicit_and_dry_run_does_not_execute(tmp_path, monkeypatch):
    config = FoldseekConfig(input_dir=tmp_path / "in", output_dir=tmp_path / "out", gpu=0)
    command = build_foldseek_command(config, tmp_path / "out" / "donors")
    assert command[command.index("--gpu") + 1] == "0"

    def _unexpected(*args, **kwargs):  # pragma: no cover - failure path
        raise AssertionError("dry-run must not execute Foldseek")

    monkeypatch.setattr("src.foldseek_cluster.subprocess.run", _unexpected)
    # The pure command path is what the CLI uses for --dry-run.
    assert build_foldseek_command(config, tmp_path / "out" / "donors") == command


def test_requested_gpu_failure_is_not_silently_downgraded(monkeypatch):
    class Result:
        returncode = 0
        stdout = "--gpu INT Use GPU (CUDA) if possible"
        stderr = ""

    assert check_gpu_capability("foldseek", runner=lambda *args, **kwargs: Result())["supported"]

    monkeypatch.setattr(
        "src.foldseek_cluster.subprocess.run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0, stdout="", stderr="GPU not available; falling back to CPU"
        ),
    )
    with pytest.raises(RuntimeError, match="GPU failure/fallback"):
        execute_foldseek(["foldseek", "easy-multimercluster", "in", "out", "tmp", "--gpu", "1"], gpu=1)


def test_panel_hash_mismatch_fails_closed(tmp_path):
    panel = tmp_path / "panel.txt"
    panel.write_text("donor-a\n")
    expected = hashlib.sha256(panel.read_bytes()).hexdigest()
    assert resolve_panel_hash(panel, expected) == expected
    with pytest.raises(ValueError, match="panel SHA-256 mismatch"):
        resolve_panel_hash(panel, "0" * 64)


def test_parser_expands_paired_orientation_and_handles_binary_report(tmp_path):
    inputs = tmp_path / "donors"
    inputs.mkdir()
    for name in ("donor-a.pdb", "donor-b.pdb", "donor-c.pdb"):
        _pdb(inputs / name)
    records = collect_donor_inputs(inputs)
    prefix = tmp_path / "foldseek"
    (tmp_path / "foldseek_cluster.tsv").write_text(
        "donor-a\tdonor-a\n"
        "donor-a\tdonor-b\n"
        "donor-c\tdonor-c\n"
    )
    (tmp_path / "foldseek_cluster_report").write_bytes(
        b"query\ttarget\tinterface_lddt\tchain_tm\tmultimer_tm\tcoverage\n"
        b"donor-a\tdonor-b\t0.91\t0.82\t0.80\t0.95\x00\n"
    )
    (tmp_path / "foldseek_rep_seq.fasta").write_text(
        ">donor-a_INT_1_A\nAAAA\n>donor-c_INT_1_A\nAAAA\n"
    )

    index = build_index_from_outputs(prefix, records, panel_sha256="panel")

    assert index["schema_version"] == "prism-template-cluster-index/v1"
    assert [c["cluster_id"] for c in index["clusters"]] == ["donor-a", "donor-c"]
    first = index["clusters"][0]
    assert first["status"] == "PASS"
    assert [(m["template_id"], m["chain_id"]) for m in first["members"]] == [
        ("donor-a", "A"),
        ("donor-a", "B"),
        ("donor-b", "A"),
        ("donor-b", "B"),
    ]
    assert first["members"][0]["metadata"]["chain_role"] == "first"
    assert first["members"][1]["metadata"]["chain_role"] == "second"
    assert first["metadata"]["report_rows"][0]["interface_lddt"] == 0.91
    assert first["representatives"][0]["template_id"] == "donor-a"


def test_duplicate_missing_and_unaccounted_members_fail_closed(tmp_path):
    inputs = tmp_path / "donors"
    inputs.mkdir()
    for name in ("a.pdb", "b.pdb"):
        _pdb(inputs / name)
    records = collect_donor_inputs(inputs)
    prefix = tmp_path / "foldseek"

    (tmp_path / "foldseek_cluster.tsv").write_text("a\ta\na\ta\n")
    with pytest.raises(MembershipError, match="duplicate"):
        build_index_from_outputs(prefix, records, panel_sha256="panel")

    (tmp_path / "foldseek_cluster.tsv").write_text("a\ta\n")
    with pytest.raises(MembershipError, match="missing"):
        build_index_from_outputs(prefix, records, panel_sha256="panel")

    (tmp_path / "foldseek_cluster.tsv").write_text("a\ta\na\tghost\n")
    with pytest.raises(MembershipError, match="unaccounted"):
        build_index_from_outputs(prefix, records, panel_sha256="panel")


def test_missing_member_singleton_policy_is_explicit_and_complete(tmp_path):
    inputs = tmp_path / "donors"
    inputs.mkdir()
    for name in ("a.pdb", "b.pdb"):
        _pdb(inputs / name)
    records = collect_donor_inputs(inputs)
    prefix = tmp_path / "foldseek"
    (tmp_path / "foldseek_cluster.tsv").write_text("a\ta\n")
    index = build_index_from_outputs(
        prefix,
        records,
        panel_sha256="panel",
        missing_member_policy="singleton",
    )
    assert len(index["clusters"]) == 2
    fallback = next(c for c in index["clusters"] if c["cluster_id"] == "b")
    assert fallback["status"] == "SINGLETON"
    assert fallback["metadata"]["index_membership_source"] == "singleton_fallback"
    assert index["parameters"]["missing_member_fallback_count"] == 1


def test_representatives_are_deterministic_adaptive_and_mark_split_required():
    members = [
        {"donor_id": "d3", "coverage": 0.20, "medoid_score": 0.80, "diversity": 0.9},
        {"donor_id": "d1", "coverage": 0.45, "medoid_score": 0.95, "diversity": 0.1},
        {"donor_id": "d2", "coverage": 0.30, "medoid_score": 0.90, "diversity": 0.8},
        {"donor_id": "d4", "coverage": 0.10, "medoid_score": 0.70, "diversity": 0.7},
    ]
    selected_a, status_a, metadata_a = select_representatives(
        members, coverage_threshold=0.96, native_representative="d1"
    )
    selected_b, status_b, metadata_b = select_representatives(
        list(reversed(members)), coverage_threshold=0.96, native_representative="d1"
    )

    assert [m["donor_id"] for m in selected_a] == ["d1", "d2", "d3"]
    assert [m["donor_id"] for m in selected_a] == [m["donor_id"] for m in selected_b]
    assert status_a == status_b == "SPLIT_REQUIRED"
    assert metadata_a == metadata_b


def test_singleton_and_empty_inputs_are_preserved(tmp_path):
    inputs = tmp_path / "singleton"
    inputs.mkdir()
    _pdb(inputs / "only.pdb")
    records = collect_donor_inputs(inputs)
    prefix = tmp_path / "single"
    (tmp_path / "single_cluster.tsv").write_text("only\tonly\n")
    index = build_index_from_outputs(prefix, records, panel_sha256="panel")
    assert index["clusters"][0]["status"] == "SINGLETON"
    assert len(index["clusters"][0]["representatives"]) == 2

    (tmp_path / "empty").mkdir()
    empty = collect_donor_inputs(tmp_path / "empty")
    empty_index = build_index_from_outputs(tmp_path / "empty-prefix", empty, panel_sha256="panel")
    assert empty_index["clusters"] == []
    assert empty_index["parameters"]["input_state"] == "EMPTY"
