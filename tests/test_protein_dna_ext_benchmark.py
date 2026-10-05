import csv
import tempfile
import unittest
from pathlib import Path

from benchmark.scripts.protein_dna_output.run_protein_dna_dockground_ext import case_to_target_ids, stage_workspace
from benchmark.scripts.protein_dna_output.score_single_protein_dna_pair_ext import score_one_ext


MODEL_CHAIN_B = """\
ATOM      1  N   LYS B   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  LYS B   1       1.200   0.000   0.000  1.00  0.00           C
ATOM      3  NZ  LYS B   1       2.300   0.000   0.000  1.00  0.00           N
ATOM      4  P   DA  C   1       4.000   0.000   0.000  1.00  0.00           P
ATOM      5  O1P DA  C   1       5.000   0.000   0.000  1.00  0.00           O
END
"""

MODEL_SHIFTED_DNA = """\
ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  LYS A   1       1.200   0.000   0.000  1.00  0.00           C
ATOM      3  NZ  LYS A   1       2.300   0.000   0.000  1.00  0.00           N
ATOM      4  P   DA  C  10       4.000   0.000   0.000  1.00  0.00           P
ATOM      5  O1P DA  C  10       5.000   0.000   0.000  1.00  0.00           O
END
"""

NATIVE_SHIFTED_DNA = """\
ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  LYS A   1       1.200   0.000   0.000  1.00  0.00           C
ATOM      3  NZ  LYS A   1       2.300   0.000   0.000  1.00  0.00           N
ATOM      4  P   DA  C   1       4.000   0.000   0.000  1.00  0.00           P
ATOM      5  O1P DA  C   1       5.000   0.000   0.000  1.00  0.00           O
END
"""

NATIVE_CHAIN_A = """\
ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  LYS A   1       1.200   0.000   0.000  1.00  0.00           C
ATOM      3  NZ  LYS A   1       2.300   0.000   0.000  1.00  0.00           N
ATOM      4  P   DA  C   1       4.000   0.000   0.000  1.00  0.00           P
ATOM      5  O1P DA  C   1       5.000   0.000   0.000  1.00  0.00           O
END
"""


class ProteinDnaExtBenchmarkTests(unittest.TestCase):
    def test_score_one_ext_uses_chain_qualified_contacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model.pdb"
            native = Path(tmpdir) / "native.pdb"
            model.write_text(MODEL_CHAIN_B)
            native.write_text(NATIVE_CHAIN_A)
            result = score_one_ext(
                str(model),
                str(native),
                model_protein_chains=["B"],
                model_dna_chains=["C"],
                native_protein_chains=["A"],
                native_dna_chains=["C"],
            )
        self.assertAlmostEqual(result["contact_precision"], 0.0)
        self.assertAlmostEqual(result["contact_recall"], 0.0)
        self.assertIn("protein_interface_precision", result)
        self.assertIn("protein_interface_f1", result)
        self.assertIn("dna_register_contact_precision", result)
        self.assertIn("nucleotide_contact_f1", result)

    def test_score_one_ext_reports_register_aware_dna_overlap(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model_shifted.pdb"
            native = Path(tmpdir) / "native_shifted.pdb"
            model.write_text(MODEL_SHIFTED_DNA)
            native.write_text(NATIVE_SHIFTED_DNA)
            result = score_one_ext(
                str(model),
                str(native),
                model_protein_chains=["A"],
                model_dna_chains=["C"],
                native_protein_chains=["A"],
                native_dna_chains=["C"],
            )
        self.assertAlmostEqual(result["contact_precision"], 0.0)
        self.assertAlmostEqual(result["contact_recall"], 0.0)
        self.assertAlmostEqual(result["dna_register_contact_precision"], 1.0)
        self.assertAlmostEqual(result["dna_register_contact_recall"], 1.0)

    def test_stage_workspace_writes_inputs_and_hints(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            protein = tmpdir / "protein.pdb"
            dna = tmpdir / "dna.pdb"
            complex_pdb = tmpdir / "complex.pdb"
            protein.write_text(MODEL_CHAIN_B)
            dna.write_text("ATOM      1  P   DA  C   1       4.000   0.000   0.000  1.00  0.00           P\nEND\n")
            complex_pdb.write_text(NATIVE_CHAIN_A)
            manifest = tmpdir / "manifest.csv"
            manifest.write_text(
                "pair_id,case_id,template_id,protein_unbound_pdb,dna_unbound_pdb,native_complex_pdb,template_complex_pdb,protein_chain_ids,dna_chain_ids,template_protein_chain_ids,template_dna_chain_ids\n"
                f"case1,1R4O,r4oaAC,{protein},{dna},{complex_pdb},{complex_pdb},A,C,A,C\n"
            )
            workspace = tmpdir / "workspace"
            rows = stage_workspace(manifest, workspace)
            self.assertEqual(len(rows), 1)
            self.assertEqual(case_to_target_ids("1R4O"), ("r4op", "r4od"))
            self.assertTrue((workspace / "processed/pdbs/r4op.pdb").exists())
            self.assertTrue((workspace / "processed/pdbs/r4od.pdb").exists())
            self.assertTrue((workspace / "templates/pdbs/r4oa.pdb").exists())
            with open(workspace / "inputs.csv", newline="") as handle:
                inputs_rows = list(csv.DictReader(handle))
            self.assertEqual(inputs_rows[0]["Receptor"], "r4op")
            self.assertEqual(inputs_rows[0]["Ligand"], "r4od")
            self.assertTrue((workspace / "templates/dna_ext_template_hints.json").exists())


if __name__ == "__main__":
    unittest.main()
