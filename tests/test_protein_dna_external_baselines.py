import csv
import json
import tempfile
import unittest
from pathlib import Path

from benchmark.scripts.protein_dna_output.run_protein_dna_external_baselines import run_probe


class ProteinDnaExternalBaselineTests(unittest.TestCase):
    def test_probe_emits_plan_and_skips_unavailable_tools(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            manifest = tmpdir / "manifest.csv"
            registry = tmpdir / "registry.json"
            output_root = tmpdir / "out"
            manifest.write_text(
                "pair_id,case_id,template_id,protein_unbound_pdb,dna_unbound_pdb,native_complex_pdb,template_complex_pdb,protein_chain_ids,dna_chain_ids,template_protein_chain_ids,template_dna_chain_ids\n"
                "case1,1R4O,r4oaAC,/tmp/p.pdb,/tmp/d.pdb,/tmp/n.pdb,/tmp/t.pdb,A,C,A,C\n"
            )
            registry.write_text(
                json.dumps(
                    [
                        {
                            "tool_id": "lightdock",
                            "name": "LightDock",
                            "integration_kind": "local_cli",
                            "supports": ["protein-dna"],
                            "binary_candidates": ["lightdock"],
                            "official_url": "https://lightdock.org/",
                            "source_url": "https://github.com/lightdock/lightdock",
                            "notes": "test",
                        },
                        {
                            "tool_id": "pydockdna",
                            "name": "pyDockDNA",
                            "integration_kind": "web_server",
                            "supports": ["protein-dna"],
                            "binary_candidates": [],
                            "official_url": "https://example.org",
                            "source_url": "https://example.org",
                            "notes": "test",
                        },
                    ],
                    indent=2,
                )
            )

            summary = run_probe(manifest, registry, output_root, tools=["lightdock", "pydockdna"])

            self.assertEqual(summary["tool_count"], 2)
            self.assertEqual(summary["pair_count"], 1)
            self.assertIn("lightdock", summary["skipped_tools"])
            self.assertIn("pydockdna", summary["skipped_tools"])
            self.assertTrue((output_root / "tool_capability_report.csv").exists())
            self.assertTrue((output_root / "tool_baseline_plan.csv").exists())
            with open(output_root / "tool_baseline_plan.csv", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["pair_id"], "case1")


if __name__ == "__main__":
    unittest.main()
