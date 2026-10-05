import json
import tempfile
import unittest
from pathlib import Path

from benchmark.scripts.protein_dna_output.validate_protein_dna_frontier_pipeline import (
    parse_stage_list,
    run_manifest_stage,
    run_workspace_stage,
)


PROTEIN_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.400   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   1.200   0.000  1.00  0.00           C
END
"""


DNA_PDB = """\
ATOM      1  P   DA  C   1       8.000   0.000   0.000  1.00  0.00           P
ATOM      2  O1P DA  C   1       8.900   0.000   0.000  1.00  0.00           O
END
"""


class FrontierPipelineValidationTests(unittest.TestCase):
    def test_parse_stage_list_normalizes_aliases(self):
        self.assertEqual(parse_stage_list("manifest,stage,probe,dry_run"), ["manifest", "workspace", "probe", "dry-run"])

    def test_manifest_and_workspace_stages_write_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            protein = tmpdir / "protein.pdb"
            dna = tmpdir / "dna.pdb"
            template = tmpdir / "template.pdb"
            protein.write_text(PROTEIN_PDB)
            dna.write_text(DNA_PDB)
            template.write_text(PROTEIN_PDB.replace("A   1", "A   2"))

            manifest = tmpdir / "manifest.csv"
            manifest.write_text(
                "pair_id,case_id,template_id,protein_unbound_pdb,dna_unbound_pdb,native_complex_pdb,template_complex_pdb,protein_chain_ids,dna_chain_ids,template_protein_chain_ids,template_dna_chain_ids\n"
                f"case1,1R4O,r4oaAC,{protein},{dna},{protein},{template},A,C,A,C\n"
            )

            registry = tmpdir / "registry.json"
            registry.write_text(
                json.dumps(
                    [
                        {
                            "tool_id": "chai1",
                            "name": "Chai-1",
                            "integration_kind": "local_cli",
                        }
                    ]
                )
            )

            output_root = tmpdir / "output"
            work_root = tmpdir / "work"

            manifest_summary = run_manifest_stage(manifest, registry, output_root)
            workspace_summary = run_workspace_stage(manifest, output_root, work_root)

            self.assertEqual(manifest_summary["manifest_row_count"], 1)
            self.assertEqual(manifest_summary["registry_row_count"], 1)
            self.assertFalse(manifest_summary["missing_manifest_fields"])
            self.assertFalse(manifest_summary["missing_registry_fields"])

            self.assertEqual(workspace_summary["pair_count"], 1)
            self.assertTrue(Path(workspace_summary["inputs_csv"]).exists())
            self.assertTrue(Path(workspace_summary["checked_templates"]).exists())
            self.assertTrue(Path(workspace_summary["template_hints"]).exists())


if __name__ == "__main__":
    unittest.main()
