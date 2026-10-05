import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from Bio.PDB import PDBParser

from benchmark.scripts.protein_dna_output.frontier_model_adapters import (
    load_registry,
    normalize_prediction_file,
    probe_tool,
    alphafold3_database_root,
    alphafold3_database_available,
    alphafold3_model_root,
    write_alphafold3_input,
    write_boltz_input,
    write_chai_input,
)


PROTEIN_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.400   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   1.200   0.000  1.00  0.00           C
ATOM      4  N   GLY A   2       3.000   1.400   0.000  1.00  0.00           N
ATOM      5  CA  GLY A   2       4.300   1.400   0.000  1.00  0.00           C
ATOM      6  C   GLY A   2       5.000   2.500   0.000  1.00  0.00           C
END
"""


DNA_PDB = """\
ATOM      1  P   DA  C   1       8.000   0.000   0.000  1.00  0.00           P
ATOM      2  O1P DA  C   1       8.900   0.000   0.000  1.00  0.00           O
ATOM      3  P   DT  C   2      10.000   0.000   0.000  1.00  0.00           P
ATOM      4  O1P DT  C   2      10.900   0.000   0.000  1.00  0.00           O
END
"""


MODEL_PDB = """\
ATOM      1  N   ALA X   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA X   1       1.400   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA X   1       2.000   1.200   0.000  1.00  0.00           C
ATOM      4  P   DA  Y   1       8.000   0.000   0.000  1.00  0.00           P
ATOM      5  O1P DA  Y   1       8.900   0.000   0.000  1.00  0.00           O
END
"""


class FrontierModelAdapterTests(unittest.TestCase):
    def test_registry_contains_frontier_tools(self):
        registry = load_registry("benchmark/data/protein_dna_frontier_tools.json")
        tool_ids = {row["tool_id"] for row in registry}
        self.assertIn("chai1", tool_ids)
        self.assertIn("boltz2", tool_ids)
        self.assertIn("alphafold3", tool_ids)
        self.assertIn("rosettafoldna", tool_ids)
        self.assertIn("rosettafold_all_atom", tool_ids)
        af3 = next(row for row in registry if row["tool_id"] == "alphafold3")
        self.assertEqual(af3["database_root_env"], "PRISM_AF3_DB_DIR")
        self.assertEqual(af3["database_root_default"], "/datasets/alphafold3")

    def test_prepare_frontier_inputs_and_normalize_chains(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            protein_path = tmpdir / "protein.pdb"
            dna_path = tmpdir / "dna.pdb"
            protein_path.write_text(PROTEIN_PDB)
            dna_path.write_text(DNA_PDB)
            manifest_row = {
                "pair_id": "case1",
                "protein_unbound_pdb": str(protein_path),
                "dna_unbound_pdb": str(dna_path),
                "protein_chain_ids": "A",
                "dna_chain_ids": "C",
            }

            chai_artifacts = write_chai_input(tmpdir, manifest_row)
            boltz_artifacts = write_boltz_input(tmpdir, manifest_row)
            af3_artifacts = write_alphafold3_input(tmpdir, manifest_row)

            self.assertTrue(Path(chai_artifacts["input_path"]).exists())
            self.assertTrue(Path(boltz_artifacts["input_path"]).exists())
            self.assertTrue(Path(af3_artifacts["input_path"]).exists())
            self.assertIn(">protein|name=A", Path(chai_artifacts["input_path"]).read_text())
            self.assertIn("version: 1", Path(boltz_artifacts["input_path"]).read_text())
            self.assertIn('"dialect": "alphafold3"', Path(af3_artifacts["input_path"]).read_text())

            model_path = tmpdir / "model.pdb"
            model_path.write_text(MODEL_PDB)
            normalized = normalize_prediction_file(model_path, tmpdir / "normalized", manifest_row)
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("normalized", normalized)
            chains = [chain.id for chain in next(structure.get_models()).get_chains()]
            self.assertEqual(chains, ["A", "C"])

    def test_probe_tool_marks_missing_binaries_skipped(self):
        availability = probe_tool(
            {
                "tool_id": "chai1",
                "integration_kind": "local_cli",
                "binary_candidates": ["__definitely_missing_binary__"],
                "python_module": "__definitely_missing_module__",
            },
            python_executable="python3",
        )
        self.assertFalse(availability.available)
        self.assertTrue(availability.reason)

    @mock.patch("benchmark.scripts.protein_dna_output.frontier_model_adapters.alphafold3_container_available", return_value=True)
    @mock.patch("benchmark.scripts.protein_dna_output.frontier_model_adapters.alphafold3_supported_gpu", return_value=(True, "A100"))
    @mock.patch("benchmark.scripts.protein_dna_output.frontier_model_adapters.alphafold3_database_available", return_value=(False, "missing_af3_database_root:/home/rshadi25/public_databases"))
    def test_probe_tool_skips_alphafold3_without_database(self, *_):
        availability = probe_tool(
            {
                "tool_id": "alphafold3",
                "integration_kind": "reference_local",
                "binary_candidates": ["run_alphafold.py"],
                "python_module": "alphafold3",
            },
            python_executable="python3",
        )
        self.assertFalse(availability.available)
        self.assertEqual(availability.reason, "missing_af3_database_root:/home/rshadi25/public_databases")

    @mock.patch("benchmark.scripts.protein_dna_output.frontier_model_adapters._python_module_available", return_value=False)
    def test_probe_tool_skips_rosettafold_all_atom_without_module(self, *_):
        availability = probe_tool(
            {
                "tool_id": "rosettafold_all_atom",
                "integration_kind": "reference_local",
                "binary_candidates": ["python"],
                "python_module": "rf2aa",
            },
            python_executable="python3",
        )
        self.assertFalse(availability.available)
        self.assertEqual(availability.reason, "tool_not_installed_or_missing_runner")

    def test_alphafold3_database_root_prefers_env_override(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            db_root = tmpdir / "af3_db"
            db_root.mkdir(parents=True, exist_ok=True)
            (db_root / "bfd-first_non_consensus_sequences.fasta").write_text(">test\nACGT\n")
            with mock.patch.dict(os.environ, {"PRISM_AF3_DB_DIR": str(db_root)}, clear=False):
                self.assertEqual(alphafold3_database_root(), db_root)
                available, reason = alphafold3_database_available()
                self.assertTrue(available)
                self.assertEqual(reason, str(db_root))

    def test_alphafold3_model_root_prefers_env_override(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            model_root = tmpdir / "af3_models"
            model_root.mkdir(parents=True, exist_ok=True)
            (model_root / "params.bin").write_text("model")
            with mock.patch.dict(os.environ, {"PRISM_AF3_MODEL_DIR": str(model_root)}, clear=False):
                self.assertEqual(alphafold3_model_root(), model_root)


if __name__ == "__main__":
    unittest.main()
