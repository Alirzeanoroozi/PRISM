#!/usr/bin/env python3
import tempfile
import unittest
from pathlib import Path
import os

import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.contact import get_contacts_from_atom_lines
import src.alignment as alignment
import src.pdb_download as pdb_download
import src.rosetta_refinement as rosetta_refinement
import src.surface_extract as surface_extract
from src.transformation import create_transformed_pair
from src.transformation import load_alignment
from unittest.mock import patch


class PrismPipelineHelperTests(unittest.TestCase):
    def test_transformation_resolves_raw_benchmark_selector_to_alignment_token(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            alignment_dir = Path(tmpdir)
            (alignment_dir / "1fgnH_1kcaCH_C.json").write_text('{"match_count": 15}\n')
            self.assertEqual(
                load_alignment("1FGNH", "1kcaCH", "C", alignment_dir=str(alignment_dir))["match_count"],
                15,
            )

    def test_alignment_writes_unavailable_records_for_empty_surface_inputs(self):
        original_cwd = os.getcwd()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                os.chdir(tmpdir)
                (Path(tmpdir) / "processed" / "alignment").mkdir(parents=True)
                (Path(tmpdir) / "processed" / "surface_extraction").mkdir(parents=True)
                (Path(tmpdir) / "processed" / "surface_extraction" / "1ABCA.asa.pdb").write_text("END\n")
                alignment.align(["1ABCA"], ["1tmplAB"])
                for chain in "AB":
                    output = Path(tmpdir) / "processed" / "alignment" / f"1ABCA_1tmplAB_{chain}.json"
                    self.assertEqual(__import__("json").loads(output.read_text())["match_count"], 0)
        finally:
            os.chdir(original_cwd)

    def test_target_id_accepts_multiple_chains_and_canonicalizes_chain_order(self):
        self.assertEqual(pdb_download.normalize_target_id("1abcBAA"), "1abcAB")
        self.assertEqual(pdb_download.target_chain_ids("1abcAB"), ("A", "B"))
        self.assertEqual(pdb_download.normalize_target_id("1ABCB"), "1abcB")
        self.assertEqual(pdb_download.normalize_target_id("1IK0_A(10)"), "1ik0A")
        self.assertEqual(pdb_download.normalize_target_id("3LZT_"), "3lzt")

    def test_downloader_materializes_all_requested_chains_for_multichain_target(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            source = tmp_path / "1abc.pdb"
            source.write_text(
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
                "ATOM      2  CA  GLY B   2       1.000   0.000   0.000  1.00  0.00           C  \n"
                "ATOM      3  CA  SER C   3       2.000   0.000   0.000  1.00  0.00           C  \nEND\n"
            )

            output = Path(pdb_download.materialize_target_pdb("1abcBA", source, tmp_path))

            self.assertEqual(output.name, "1abcAB.pdb")
            atom_lines = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]
            self.assertEqual([line[21] for line in atom_lines], ["A", "B"])
            self.assertNotIn(" SER ", output.read_text())

    def test_surface_scaffold_threshold_matches_working_pipeline_default(self):
        self.assertEqual(surface_extract.SCFFTHRESHOLD, 5.0)

    def test_surface_extraction_uses_only_the_requested_chain_pdb(self):
        original_cwd = os.getcwd()
        original_surface_dir = surface_extract.SURFACE_EXTRACTION_DIR
        original_scaffold_threshold = surface_extract.SCFFTHRESHOLD

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir)
                os.chdir(tmp_path)
                (tmp_path / "processed" / "pdbs").mkdir(parents=True)
                (tmp_path / "processed" / "surface_extraction").mkdir(parents=True)
                (tmp_path / "processed" / "pdbs" / "1abc.pdb").write_text(
                    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
                    "ATOM      2  CA  GLY B   2       0.000   0.000   0.000  1.00  0.00           C  \nEND\n"
                )
                (tmp_path / "processed" / "pdbs" / "1ABCB.pdb").write_text(
                    "ATOM      2  CA  GLY B   2       0.000   0.000   0.000  1.00  0.00           C  \nEND\n"
                )
                surface_extract.SURFACE_EXTRACTION_DIR = "processed/surface_extraction"
                surface_extract.SCFFTHRESHOLD = 5.0

                with patch.object(
                    surface_extract,
                    "get_asa_complex_target",
                    return_value={"GLY_2_B": 20.0},
                ):
                    surface_extract.extract_surface("1ABCB")

                atom_lines = [
                    line
                    for line in (tmp_path / "processed" / "surface_extraction" / "1ABCB.asa.pdb").read_text().splitlines()
                    if line.startswith("ATOM")
                ]
                self.assertEqual({line[21] for line in atom_lines}, {"B"})
        finally:
            os.chdir(original_cwd)
            surface_extract.SURFACE_EXTRACTION_DIR = original_surface_dir
            surface_extract.SCFFTHRESHOLD = original_scaffold_threshold

    def test_surface_extraction_writes_empty_pdb_when_no_residues_pass_rsa(self):
        original_surface_dir = surface_extract.SURFACE_EXTRACTION_DIR
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                surface_extract.SURFACE_EXTRACTION_DIR = str(Path(tmpdir) / "surface")
                with patch.object(surface_extract, "get_asa_complex_target", return_value={}):
                    surface_extract.extract_surface("1ABCA")
                self.assertEqual(
                    (Path(tmpdir) / "surface" / "1ABCA.asa.pdb").read_text(),
                    "END\n",
                )
        finally:
            surface_extract.SURFACE_EXTRACTION_DIR = original_surface_dir

    def test_downloader_materializes_chain_pdbs_consumed_by_transformer(self):
        original_cwd = os.getcwd()
        original_target_dir = pdb_download.TARGET_DIR
        original_inputs_csv = os.environ.get("PRISM_INPUTS_CSV")

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir)
                os.chdir(tmp_path)
                (tmp_path / "processed" / "pdbs").mkdir(parents=True)
                (tmp_path / "processed" / "transformation").mkdir(parents=True)
                (tmp_path / "inputs.csv").write_text("Receptor,Ligand\n1ABCB,2DEFB\n")
                (tmp_path / "processed" / "pdbs" / "1abc.pdb").write_text(
                    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
                    "ATOM      2  CA  GLY B   2       0.000   0.000   0.000  1.00  0.00           C  \nEND\n"
                )
                (tmp_path / "processed" / "pdbs" / "2def.pdb").write_text(
                    "ATOM      1  CA  SER B   3      20.000   0.000   0.000  1.00  0.00           C  \nEND\n"
                )

                pdb_download.TARGET_DIR = "processed/pdbs"
                os.environ["PRISM_INPUTS_CSV"] = "inputs.csv"
                pdb_download.pdb_downloader()

                identity = {
                    "translation": [0.0, 0.0, 0.0],
                    "rotation_mat": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                }
                passed_pairs = []
                create_transformed_pair(
                    "1tmplAB", "1ABCB", "2DEFB", identity, identity, passed_pairs, "test"
                )

                self.assertEqual(len(passed_pairs), 1)
                receptor_text = (tmp_path / "processed" / "pdbs" / "1ABCB.pdb").read_text()
                self.assertIn("GLY B", receptor_text)
                self.assertNotIn("ALA A", receptor_text)
        finally:
            os.chdir(original_cwd)
            pdb_download.TARGET_DIR = original_target_dir
            if original_inputs_csv is None:
                os.environ.pop("PRISM_INPUTS_CSV", None)
            else:
                os.environ["PRISM_INPUTS_CSV"] = original_inputs_csv

    def test_get_contacts_from_atom_lines_writes_expected_pair(self):
        atom_lines_0 = [
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n",
        ]
        atom_lines_1 = [
            "ATOM      2  CA  GLY B   2       0.000   0.000   3.000  1.00  0.00           C  \n",
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "contacts.txt"
            get_contacts_from_atom_lines("dummy.pdb", str(output_path), atom_lines_0, atom_lines_1)
            self.assertEqual(output_path.read_text().strip(), "1\t2")

    def test_combine_pdb_reassigns_partner_chains_to_a_and_b(self):
        original_rosetta_dir = rosetta_refinement.ROSETTA_DIR
        original_energy_dir = rosetta_refinement.ENERGY_DIR
        original_structure_dir = rosetta_refinement.STRUCTURE_DIR

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir)
                rosetta_refinement.ROSETTA_DIR = str(tmp_path / "rosetta_refinement")
                rosetta_refinement.ENERGY_DIR = str(tmp_path / "rosetta_refinement" / "energies")
                rosetta_refinement.STRUCTURE_DIR = str(tmp_path / "rosetta_refinement" / "structures")
                Path(rosetta_refinement.ROSETTA_DIR).mkdir(parents=True, exist_ok=True)
                Path(rosetta_refinement.ENERGY_DIR).mkdir(parents=True, exist_ok=True)
                Path(rosetta_refinement.STRUCTURE_DIR).mkdir(parents=True, exist_ok=True)

                left = tmp_path / "left.pdb"
                right = tmp_path / "right.pdb"
                left.write_text(
                    "ATOM      1  CA  ALA H   1       1.000   2.000   3.000  1.00  0.00           C  \nEND\n"
                )
                right.write_text(
                    "ATOM      1  CA  GLY P   2       4.000   5.000   6.000  1.00  0.00           C  \nEND\n"
                )

                combined_path = rosetta_refinement.combine_pdb(str(left), str(right))
                combined_text = Path(combined_path).read_text().splitlines()

                self.assertEqual(combined_text[0][21], "A")
                self.assertEqual(combined_text[2][21], "B")
                self.assertEqual(combined_text[1], "TER")
                self.assertEqual(combined_text[-1], "END")
        finally:
            rosetta_refinement.ROSETTA_DIR = original_rosetta_dir
            rosetta_refinement.ENERGY_DIR = original_energy_dir
            rosetta_refinement.STRUCTURE_DIR = original_structure_dir

    def test_combine_pdb_preserves_all_chains_in_multichain_partners(self):
        original_rosetta_dir = rosetta_refinement.ROSETTA_DIR
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir)
                rosetta_refinement.ROSETTA_DIR = str(tmp_path / "rosetta_refinement")
                Path(rosetta_refinement.ROSETTA_DIR).mkdir(parents=True)
                left = tmp_path / "1tmpl_1abcAB_2defC_o1_L.pdb"
                right = tmp_path / "1tmpl_1abcAB_2defC_o1_R.pdb"
                left.write_text(
                    "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C  \n"
                    "ATOM      2  CA  GLY B   2       2.000   2.000   3.000  1.00  0.00           C  \nEND\n"
                )
                right.write_text(
                    "ATOM      1  CA  SER C   3       4.000   5.000   6.000  1.00  0.00           C  \nEND\n"
                )

                combined_path = rosetta_refinement.combine_pdb(str(left), str(right))
                atom_lines = [line for line in Path(combined_path).read_text().splitlines() if line.startswith("ATOM")]

                self.assertEqual([line[21] for line in atom_lines], ["A", "B", "C"])
                self.assertEqual(rosetta_refinement.partner_chain_ids(str(left), str(right)), ("AB", "C"))
        finally:
            rosetta_refinement.ROSETTA_DIR = original_rosetta_dir


if __name__ == "__main__":
    unittest.main()
