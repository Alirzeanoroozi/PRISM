import csv
import hashlib
import json
from pathlib import Path

from benchmark.scripts.assemble_bm55_native_complexes import assemble_native_complexes


def _atom(chain: str, serial: int) -> str:
    return (
        f"ATOM  {serial:5d}  CA  ALA {chain}   1    "
        "   0.000   0.000   0.000  1.00 20.00           C\n"
    )


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_table(path: Path, rows: list[dict[str, str]], delimiter: str) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def test_assemble_curated_roles_and_enrich_stage_manifest(tmp_path):
    receptor = tmp_path / "1ABC_r_b.pdb"
    ligand = tmp_path / "1ABC_l_b.pdb"
    receptor.write_text(_atom("A", 1) + _atom("B", 2) + "END\n", encoding="ascii")
    ligand.write_text(
        _atom("C", 1)
        + "HETATM    2  O   HOH     2       0.000   0.000   0.000  1.00 20.00           O\n"
        + "END\n",
        encoding="ascii",
    )

    stage_manifest = tmp_path / "stage.csv"
    _write_table(
        stage_manifest,
        [{"dataset_row_id": "rigid:000001", "status": "staged_symlink", "complex": "1ABC_AB:C"}],
        ",",
    )
    source_manifest = tmp_path / "source.tsv"
    identity = {
        "dataset_row_id": "rigid:000001",
        "native_complex": "1ABC_AB:C",
        "native_receptor_chains": "AB",
        "native_ligand_chains": "C",
    }
    _write_table(
        source_manifest,
        [
            identity | {"source_role": "native_receptor", "sha256": _sha(receptor)},
            identity | {"source_role": "native_ligand", "sha256": _sha(ligand)},
        ],
        "\t",
    )
    staged_sources = tmp_path / "staged.tsv"
    _write_table(
        staged_sources,
        [
            {
                "dataset_row_id": "rigid:000001", "source_role": "native_receptor", "status": "staged",
                "staged_path": str(receptor), "staged_sha256": _sha(receptor),
                "source_payload_sha256": _sha(receptor), "archive_member": "benchmark5.5/structures/1ABC_r_b.pdb",
            },
            {
                "dataset_row_id": "rigid:000001", "source_role": "native_ligand", "status": "staged",
                "staged_path": str(ligand), "staged_sha256": _sha(ligand),
                "source_payload_sha256": _sha(ligand), "archive_member": "benchmark5.5/structures/1ABC_l_b.pdb",
            },
        ],
        "\t",
    )
    policy = tmp_path / "policy.json"
    policy.write_text(json.dumps({"decision": {"excluded_dataset_row_ids": []}}), encoding="utf-8")

    assemblies, enriched = assemble_native_complexes(
        stage_manifest,
        source_manifest,
        staged_sources,
        policy,
        tmp_path / "native",
        tmp_path / "assemblies.tsv",
        tmp_path / "enriched.csv",
    )

    assert assemblies[0]["assembly_status"] == "assembled"
    native = Path(assemblies[0]["native_pdb_path"])
    assert native.is_file()
    assert " A   1" in native.read_text(encoding="ascii")
    assert " C   1" in native.read_text(encoding="ascii")
    assert "HETATM" not in native.read_text(encoding="ascii")
    assert enriched[0]["source_gate_status"] == "strict_clean"
    assert enriched[0]["native_pdb_path"] == str(native)
    assert enriched[0]["native_receptor_archive_member"].endswith("1ABC_r_b.pdb")
