from pathlib import Path

from benchmark.scripts.stage_current_models_for_main_benchmark import parse_current_model_name, stage_models


def _pdb_line(chain: str, residue: int) -> str:
    return (
        f"ATOM      1  CA  ALA {chain}{residue:4d}    "
        "   0.000   0.000   0.000  1.00 20.00           C\n"
    )


def _write_inputs(batch_root: Path, receptor: str = "1abcAB", ligand: str = "3ghiB") -> None:
    batch = batch_root / "batch_0001"
    batch.mkdir(parents=True)
    (batch / "inputs.csv").write_text(
        "pair_id,complex,benchmark_set,source_row,pdb_id_1_raw,pdb_id_2_raw,Receptor,Ligand\n"
        f"unit_0002,1abc_A:B,unit,2,1ABC_AB,3GHI_B,{receptor},{ligand}\n",
        encoding="utf-8",
    )


def test_parse_external_rosetta_suffix_variants():
    prefix = "1bjaAB_1abcAB_3ghiB_o2_L_1bjaAB_1abcAB_3ghiB_o2_R_rosetta"
    for suffix in (".pdb", "_0001.pdb", "_0001_0001.pdb"):
        parsed = parse_current_model_name(prefix + suffix)
        assert parsed is not None
        assert parsed["template_1"] == "1bjaA"
        assert parsed["template_2"] == "1bjaB"
        assert parsed["target_left"] == "1abcAB"
        assert parsed["target_right"] == "3ghiB"


def test_stage_rejects_l_half_even_when_filename_matches_input(tmp_path):
    batch_root = tmp_path / "batches"
    current_root = tmp_path / "current"
    output_root = tmp_path / "staged"
    _write_inputs(batch_root)
    model = current_root / "batch_0001/processed/rosetta_refinement"
    model.mkdir(parents=True)
    (model / "1abcAB_2defAB_3ghiAB_o1_L_rosetta_0001.pdb").write_text(
        _pdb_line("A", 1) + _pdb_line("B", 1), encoding="ascii"
    )

    records = stage_models(batch_root, current_root, output_root)

    assert records[0]["status"] == "transformation_intermediate"
    assert records[0]["integrity_reason"] == "incomplete_transformation_half"
    assert not output_root.exists()


def test_stage_accepts_pyrosetta_model_with_compact_template_tokens(tmp_path):
    batch_root = tmp_path / "batches"
    current_root = tmp_path / "current"
    output_root = tmp_path / "staged"
    _write_inputs(batch_root)
    model_dir = current_root / "batch_0001/processed/pyrosetta_refinement/structures"
    model_dir.mkdir(parents=True)
    model = model_dir / "1abcAB_2defA_3ghiB_o1_L_1abcAB_2defA_3ghiB_o1_R_rosetta.pdb"
    model.write_text(_pdb_line("A", 1) + _pdb_line("B", 1) + _pdb_line("C", 1), encoding="ascii")

    records = stage_models(batch_root, current_root, output_root)

    assert records[0]["status"] == "staged_symlink"
    assert records[0]["template_1"] == "2defA"
    assert records[0]["template_2"] == "3ghiB"
    assert records[0]["observed_chain_order"] == "A,B,C"
    assert records[0]["model_receptor_chains"] == "AB"
    assert records[0]["model_ligand_chains"] == "C"
    assert records[0]["ca_counts"] == "A:1;B:1;C:1"
    assert len(records[0]["source_model_sha256"]) == 64
    assert Path(records[0]["staged_model_path"]).is_symlink()


def test_stage_accepts_query_first_external_rosetta_name(tmp_path):
    batch_root = tmp_path / "batches"
    current_root = tmp_path / "current"
    output_root = tmp_path / "staged"
    _write_inputs(batch_root)
    model_dir = current_root / "batch_0001/processed/rosetta_refinement"
    model_dir.mkdir(parents=True)
    model = model_dir / "1bjaAB_1abcAB_3ghiB_o2_L_1bjaAB_1abcAB_3ghiB_o2_R_rosetta_0001_0001.pdb"
    model.write_text(_pdb_line("A", 1) + _pdb_line("B", 1) + _pdb_line("C", 1), encoding="ascii")

    records = stage_models(batch_root, current_root, output_root)

    assert records[0]["status"] == "staged_symlink"
    assert records[0]["template_1"] == "1bjaA"
    assert records[0]["template_2"] == "1bjaB"
    assert records[0]["target_left"] == "1abcAB"
    assert records[0]["target_right"] == "3ghiB"
    assert records[0]["model_receptor_chains"] == "AB"
    assert records[0]["model_ligand_chains"] == "C"
    assert records[0]["dataset_row_id"] == "unit:000001"
    assert records[0]["raw_receptor_selector"] == "1ABC_AB"


def test_stage_selects_only_most_refined_external_artifact(tmp_path):
    batch_root = tmp_path / "batches"
    current_root = tmp_path / "current"
    output_root = tmp_path / "staged"
    _write_inputs(batch_root)
    model_dir = current_root / "batch_0001/processed/rosetta_refinement"
    model_dir.mkdir(parents=True)
    prefix = "1bjaAB_1abcAB_3ghiB_o2_L_1bjaAB_1abcAB_3ghiB_o2_R_rosetta"
    pdb = _pdb_line("A", 1) + _pdb_line("B", 1) + _pdb_line("C", 1)
    for suffix in (".pdb", "_0001.pdb", "_0001_0001.pdb"):
        (model_dir / (prefix + suffix)).write_text(pdb, encoding="ascii")

    records = stage_models(batch_root, current_root, output_root)

    assert len(records) == 1
    assert records[0]["status"] == "staged_symlink"
    assert records[0]["source_model_path"].endswith("_rosetta_0001_0001.pdb")


def test_stage_rejects_external_model_with_template_sized_chain_groups(tmp_path):
    batch_root = tmp_path / "batches"
    current_root = tmp_path / "current"
    output_root = tmp_path / "staged"
    _write_inputs(batch_root)
    model_dir = current_root / "batch_0001/processed/rosetta_refinement"
    model_dir.mkdir(parents=True)
    model = model_dir / "1bjaAB_1abcAB_3ghiB_o2_L_1bjaAB_1abcAB_3ghiB_o2_R_rosetta_0001_0001.pdb"
    model.write_text(_pdb_line("A", 1) + _pdb_line("B", 1), encoding="ascii")

    records = stage_models(batch_root, current_root, output_root)

    assert records[0]["status"] == "invalid_model_chain_contract"
    assert records[0]["integrity_reason"].startswith("unexpected_model_chain_count")
    assert "expected=3 receptor=2 ligand=1" in records[0]["integrity_reason"]
    assert not output_root.exists()


def test_stage_rejects_models_that_lack_a_declared_partner_chain(tmp_path):
    batch_root = tmp_path / "batches"
    current_root = tmp_path / "current"
    output_root = tmp_path / "staged"
    _write_inputs(batch_root)
    model_dir = current_root / "batch_0001/processed/rosetta_refinement"
    model_dir.mkdir(parents=True)
    model = model_dir / "1abcAB_2defA_3ghiB_o1_L_1abcAB_2defA_3ghiB_o1_R_rosetta_0001.pdb"
    model.write_text(_pdb_line("A", 1), encoding="ascii")

    records = stage_models(batch_root, current_root, output_root)

    assert records[0]["status"] == "invalid_model_chain_contract"
    assert "unexpected_model_chain_count" in records[0]["chain_contract_errors"]
    assert not output_root.exists()
