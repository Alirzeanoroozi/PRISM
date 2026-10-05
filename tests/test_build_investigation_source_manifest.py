import csv
import hashlib
import io
import json
import tarfile
from pathlib import Path

from benchmark.scripts.build_investigation_source_manifest import ARCHIVE_PREFIX_ALIASES, _local_candidates, main
from benchmark.scripts.stage_curated_sources import stage_sources


HEADER = [
    "Complex",
    "Cat.",
    "PDB ID 1",
    "Protein 1",
    "PDB ID 2",
    "Protein 2",
    "I-RMSD (Å)",
    "ΔASA(Å2)",
    "BM version introduced",
]


def _pdb(chain: str, residue: int = 1) -> bytes:
    return (
        f"ATOM      1  CA  ALA {chain}{residue:4d}      1.000   2.000   3.000  1.00 20.00           C  \n"
        "TER\nEND\n"
    ).encode("ascii")


def _write_dataset(root: Path, name: str, row: list[str]) -> None:
    path = root / "benchmark" / "data" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        writer.writerow(row)


def _write_archive(path: Path, members: dict[str, bytes]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as archive:
        for name, payload in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _make_repo(tmp_path: Path, *, include_archive: bool = True) -> Path:
    root = tmp_path / "repo"
    row = [
        "1NAT_A:C",
        "AA",
        "1ABC_A",
        "raw receptor",
        "2DEF_B",
        "raw ligand",
        "1.23",
        "456",
        "5.5",
    ]
    for filename in ("T_Rigid.csv", "T_medium.csv", "T_difficult.csv"):
        _write_dataset(root, filename, row)

    # Two exact local selector candidates intentionally collide across source sets.
    for set_name in ("rigid", "medium"):
        selector = root / "benchmark" / "data" / "pdbs" / set_name / "chainwise" / "1ABC_A" / "1abc_A.pdb"
        selector.parent.mkdir(parents=True, exist_ok=True)
        selector.write_bytes(_pdb("A"))
    ligand = root / "benchmark" / "data" / "pdbs" / "rigid" / "chainwise" / "2DEF_B" / "2def_B.pdb"
    ligand.parent.mkdir(parents=True, exist_ok=True)
    ligand.write_bytes(_pdb("B"))

    if include_archive:
        _write_archive(
            root / "benchmark" / "originals" / "benchmark5.5.tgz",
            {
                "benchmark5.5/structures/1NAT_r_b.pdb": _pdb("A", 1),
                "benchmark5.5/structures/1NAT_l_b.pdb": _pdb("C", 2),
                "benchmark5.5/structures/1NAT_r_u.pdb": _pdb("A", 3),
                "benchmark5.5/structures/1NAT_l_u.pdb": _pdb("C", 4),
            },
        )
    return root


def test_cli_preserves_rows_hashes_parser_metadata_and_collisions(tmp_path):
    root = _make_repo(tmp_path)
    output = tmp_path / "out"

    assert main(["--repo-root", str(root), "--output-dir", str(output), "--limit", "1"]) == 0
    source_rows = _read_tsv(output / "source_manifest.tsv")
    validation_rows = _read_tsv(output / "structure_validation.tsv")

    assert len(source_rows) == 7
    assert {row["dataset_row_id"] for row in source_rows} == {"rigid:000001"}
    assert {row["difficulty"] for row in source_rows} == {"rigid"}
    assert {row["raw_receptor_selector"] for row in source_rows} == {"1ABC_A"}
    assert {row["raw_ligand_selector"] for row in source_rows} == {"2DEF_B"}
    assert {row["native_complex"] for row in source_rows} == {"1NAT_A:C"}
    assert {row["native_receptor_chains"] for row in source_rows} == {"A"}
    assert {row["native_ligand_chains"] for row in source_rows} == {"C"}

    source_row = json.loads(source_rows[0]["source_row_json"])
    assert source_row["PDB ID 1"] == "1ABC_A"
    assert source_row["PDB ID 2"] == "2DEF_B"
    local_rows = [row for row in source_rows if row["source_kind"] == "benchmark_data_pdbs"]
    assert len(local_rows) == 3
    assert {row["candidate_status"] for row in local_rows if row["source_role"] == "receptor_selector"} == {"collision"}

    archive_rows = [row for row in source_rows if row["source_kind"] == "benchmark5.5_archive_member"]
    assert len(archive_rows) == 4
    assert {row["archive_prefix"] for row in archive_rows} == {"1NAT"}
    assert all(len(row["sha256"]) == 64 for row in source_rows)
    assert all(row["resolution_status"] == "resolved" for row in source_rows)

    assert len(validation_rows) == len(source_rows)
    assert all(row["parse_status"] == "ok" for row in validation_rows)
    assert {row["parser_name"] for row in validation_rows} == {"Bio.PDB.PDBParser"}
    assert {row["biopython_version"] for row in validation_rows} == {"1.84"}
    assert {row["chain_ids"] for row in validation_rows} == {"A", "B", "C"}
    assert all(int(row["atom_count"]) == 1 for row in validation_rows)
    assert all(row["residue_ids_sha256"] for row in validation_rows)
    assert all(row["sequence_hashes"] for row in validation_rows)
    assert all(row["source_scope"] == "pipeline" for row in validation_rows if row["source_role"].startswith(("pipeline_", "native_")))

    first = (output / "source_manifest.tsv").read_bytes(), (output / "structure_validation.tsv").read_bytes()
    output_again = tmp_path / "out-again"
    assert main(["--repo-root", str(root), "--output-dir", str(output_again), "--limit", "1"]) == 0
    assert first == (
        (output_again / "source_manifest.tsv").read_bytes(),
        (output_again / "structure_validation.tsv").read_bytes(),
    )


def test_strict_mode_writes_explicit_unresolved_rows_and_fails_closed(tmp_path):
    root = _make_repo(tmp_path, include_archive=False)
    output = tmp_path / "out"

    assert main(["--repo-root", str(root), "--output-dir", str(output), "--limit", "1", "--strict"]) == 2
    source_rows = _read_tsv(output / "source_manifest.tsv")
    validation_rows = _read_tsv(output / "structure_validation.tsv")
    native_rows = [row for row in source_rows if row["source_role"].startswith("native_")]
    assert len(native_rows) == 2
    assert all(row["resolution_status"] == "unresolved" for row in native_rows)
    assert all(row["archive_prefix"] == "1NAT" for row in native_rows)
    native_validation = [row for row in validation_rows if row["source_role"].startswith("native_")]
    assert all(row["parse_status"] == "unresolved" for row in native_validation)
    assert any("archive" in row["error"] for row in native_validation)


def test_hash_in_fixture_matches_manifest(tmp_path):
    root = _make_repo(tmp_path)
    output = tmp_path / "out"
    assert main(["--repo-root", str(root), "--output-dir", str(output), "--limit", "1"]) == 0
    rows = _read_tsv(output / "source_manifest.tsv")
    local = next(row for row in rows if row["source_kind"] == "benchmark_data_pdbs")
    assert local["sha256"] == hashlib.sha256(Path(root / local["source_path"]).read_bytes()).hexdigest()


def test_exact_dataset_row_selection_is_deterministic(tmp_path):
    root = _make_repo(tmp_path)
    output = tmp_path / "selected"
    assert main([
        "--repo-root", str(root),
        "--output-dir", str(output),
        "--dataset-row-id", "rigid:000001",
    ]) == 0
    source_rows = _read_tsv(output / "source_manifest.tsv")
    assert source_rows
    assert {row["dataset_row_id"] for row in source_rows} == {"rigid:000001"}

    assert main([
        "--repo-root", str(root),
        "--output-dir", str(tmp_path / "unknown"),
        "--dataset-row-id", "rigid:999999",
    ]) == 2


def test_curated_roles_stage_to_hash_named_files_only(tmp_path):
    root = _make_repo(tmp_path)
    output = tmp_path / "out"
    assert main(["--repo-root", str(root), "--output-dir", str(output), "--limit", "1"]) == 0
    staged, failures = stage_sources(
        output / "source_manifest.tsv",
        tmp_path / "staged",
        repo_root=root,
        strict=True,
    )
    assert failures == 0
    rows = _read_tsv(staged)
    assert {row["source_role"] for row in rows} == {
        "pipeline_receptor",
        "pipeline_ligand",
        "native_receptor",
        "native_ligand",
    }
    assert all(row["status"] == "staged" for row in rows)
    assert all(Path(row["staged_path"]).is_file() for row in rows)
    assert all(row["staged_sha256"] == row["source_payload_sha256"] for row in rows)


def test_qualified_selector_does_not_accept_unqualified_full_pdb(tmp_path):
    pdb_root = tmp_path / "pdbs" / "rigid"
    pdb_root.mkdir(parents=True)
    (pdb_root / "1abc.pdb").write_bytes(_pdb("A"))
    assert _local_candidates(pdb_root.parent, "1ABC_A") == []


def test_archive_alias_table_covers_readme_aliases():
    assert ARCHIVE_PREFIX_ALIASES["1QFW"] == ("1QFW", "9QFW")
    assert ARCHIVE_PREFIX_ALIASES["3AAD"] == ("3AAD", "BAAD")
    assert ARCHIVE_PREFIX_ALIASES["1OYV"] == ("1OYV", "BOYV")
    assert ARCHIVE_PREFIX_ALIASES["3P57"] == ("3P57", "BP57", "CP57")
