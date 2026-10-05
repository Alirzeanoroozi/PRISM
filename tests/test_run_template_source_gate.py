import csv
import json
from pathlib import Path

from benchmark.scripts.run_template_source_gate import build_manifest


def _atom(serial: int, residue: int, chain: str, resname: str) -> str:
    return (
        f"ATOM  {serial:5d}  CA  {resname:>3s} {chain}{residue:4d}    "
        "   0.000   0.000   0.000  1.00 20.00           C\n"
    )


def _write_pdb(path: Path, chain: str, residues: list[str]) -> None:
    path.write_text("".join(_atom(i, i, chain, r) for i, r in enumerate(residues, 1)), encoding="ascii")


def test_source_gate_outputs_similarity_rows_and_respects_policy(tmp_path):
    staged = tmp_path / "staged"
    interfaces = tmp_path / "interfaces"
    (interfaces).mkdir(parents=True)
    staged.mkdir()
    _write_pdb(staged / "1abc.pdb", "A", ["ALA", "ALA", "ALA", "ALA"])
    _write_pdb(staged / "1abc.pdb", "B", ["ALA", "ALA", "ALA", "ALA"])
    # Keep the fixture single-chain so the source-gate mapping is explicit.
    (staged / "1abc.pdb").write_text("".join(_atom(i, i, "A", "ALA") for i in range(1, 5)), encoding="ascii")
    _write_pdb(interfaces / "2defAB_A_int.pdb", "A", ["CYS", "CYS", "CYS", "CYS"])
    _write_pdb(interfaces / "2defAB_B_int.pdb", "B", ["CYS", "CYS", "CYS", "CYS"])

    inputs = tmp_path / "inputs.csv"
    inputs.write_text("pair_id,dataset_row_id,Receptor,Ligand\np1,rigid:1,1ABC_A,1ABC_A\n", encoding="utf-8")
    templates = tmp_path / "templates.txt"
    templates.write_text("2defAB\n", encoding="utf-8")
    policy = tmp_path / "policy.json"
    policy.write_text(json.dumps({"decision": {"status": "blocked_source_authority", "confirmatory_run_authorized": False, "excluded_dataset_row_ids": []}}), encoding="utf-8")
    similarity = tmp_path / "similarity.tsv"
    eligible = tmp_path / "eligible.tsv"

    summary = build_manifest(
        inputs=inputs,
        staged_pdb_dir=staged,
        template_list=templates,
        interface_root=interfaces,
        source_policy=policy,
        similarity_output=similarity,
        eligible_output=eligible,
    )

    assert summary["row_count"] == 1
    assert summary["template_count"] == 1
    with eligible.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["source_gate_status"] == "blocked_source_authority"
    assert row["eligible_template_count"] == "0"
    with similarity.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert {row["exclusion_reason"] for row in rows} == {"eligible"}
