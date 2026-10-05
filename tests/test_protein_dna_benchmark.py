import tempfile
import unittest
from pathlib import Path

from benchmark.scripts.protein_dna_output.score_single_protein_dna_pair import score_one


class ProteinDnaBenchmarkTests(unittest.TestCase):
    def test_score_one_fixture_pair(self):
        repo_root = Path(__file__).resolve().parents[1]
        native = repo_root / "benchmark/data/protein_dna_fixtures/fixture_native.pdb"
        positive = repo_root / "benchmark/data/protein_dna_fixtures/fixture_positive_model.pdb"
        result = score_one(
            str(positive),
            str(native),
            model_protein_chains=["A"],
            model_dna_chains=["B"],
            native_protein_chains=["A"],
            native_dna_chains=["B"],
        )
        self.assertAlmostEqual(result["contact_precision"], 1.0)
        self.assertAlmostEqual(result["contact_recall"], 1.0)
        self.assertGreater(result["alignment_tm_score"], 0.9)


if __name__ == "__main__":
    unittest.main()
