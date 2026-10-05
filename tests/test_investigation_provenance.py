import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from benchmark.scripts.investigation_provenance import (
    REDACTED,
    build_provenance_manifest,
    hash_files,
    preflight_template_assets,
    sha256_file,
    write_template_preflight_tsv,
)


class InvestigationProvenanceTests(unittest.TestCase):
    def test_sha256_and_hash_records_are_stable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            first = root / "first.txt"
            second = root / "second.txt"
            first.write_text("PRISM\n", encoding="utf-8")
            second.write_text("current-vs-legacy\n", encoding="utf-8")

            expected = hashlib.sha256(first.read_bytes()).hexdigest()
            self.assertEqual(sha256_file(first), expected)
            records_a = hash_files([second, first], root=root)
            records_b = hash_files([first, second], root=root)
            self.assertEqual(records_a, records_b)
            self.assertEqual(records_a[0]["path"], "first.txt")
            self.assertEqual(records_a[0]["sha256"], expected)

    def test_manifest_redacts_secrets_and_keeps_effective_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tracked = root / "input.csv"
            tracked.write_text("template_id\n1abcAB\n", encoding="utf-8")
            environ = {
                "SAFE_SETTING": "visible",
                "API_TOKEN": "do-not-record",
                "RANDOM_SEED_TOKEN": "do-not-record-either",
                "SLURM_JOB_ID": "1234",
                "SLURM_CPUS_PER_TASK": "8",
                "PYTHONHASHSEED": "7",
            }
            manifest = build_provenance_manifest(
                root,
                files=[tracked],
                environment=["SAFE_SETTING", "API_TOKEN"],
                packages=["definitely-not-installed-prism-package"],
                command=["python", "run.py", "--api-token", "secret-value"],
                seeds={"numpy": 7},
                config={"threshold": 0.4, "database_password": "hidden"},
                environ=environ,
            )

            self.assertEqual(manifest["environment"]["SAFE_SETTING"], "visible")
            self.assertEqual(manifest["environment"]["API_TOKEN"], REDACTED)
            self.assertEqual(manifest["command"], ["python", "run.py", "--api-token", REDACTED])
            self.assertEqual(manifest["config"]["database_password"], REDACTED)
            self.assertEqual(manifest["seeds"]["numpy"], 7)
            self.assertEqual(manifest["slurm"]["SLURM_CPUS_PER_TASK"], "8")
            self.assertIsNone(manifest["packages"]["definitely-not-installed-prism-package"])
            serialized = json.dumps(manifest, sort_keys=True)
            self.assertNotIn("do-not-record", serialized)
            self.assertNotIn("secret-value", serialized)
            self.assertNotIn("hidden", serialized)

            seeded = build_provenance_manifest(
                root,
                seeds=["RANDOM_SEED_TOKEN"],
                environ=environ,
            )
            assert seeded["seeds"]["RANDOM_SEED_TOKEN"] == REDACTED

    def test_template_preflight_counts_duplicates_and_missing_assets(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            assets = root / "templates"
            (assets / "contacts").mkdir(parents=True)
            (assets / "interfaces_lists").mkdir()
            (assets / "interfaces").mkdir()
            (assets / "contacts" / "oldAB.txt").write_text("1 2\n", encoding="utf-8")
            (assets / "contacts" / "newCD.json").write_text("[[1, 2]]\n", encoding="utf-8")
            (assets / "interfaces_lists" / "newCD.json").write_text('{"C": [1]}\n', encoding="utf-8")
            (assets / "interfaces" / "newCD_C_int.pdb").write_text("ATOM\n", encoding="utf-8")
            (assets / "contacts" / "badGH.json").write_text("not json\n", encoding="utf-8")
            (assets / "interfaces_lists" / "badGH.json").write_text('{"G": [1]}\n', encoding="utf-8")
            (assets / "interfaces" / "badGH_G_int.pdb").write_text("ATOM\n", encoding="utf-8")
            manifest = root / "templates.csv"
            manifest.write_text("template_id\nnewCD\noldAB\nnewCD\nmissingEF\nbadGH\n", encoding="utf-8")

            report = preflight_template_assets(manifest, assets)

            self.assertEqual(report["listed"], 5)
            self.assertEqual(report["unique"], 4)
            self.assertEqual(report["valid"], 2)
            self.assertEqual(report["fully_resolvable"], 2)
            self.assertEqual(report["missing"], 2)
            new_template = next(item for item in report["templates"] if item["template_id"] == "newCD")
            self.assertEqual(new_template["listed_count"], 2)
            missing_template = next(item for item in report["templates"] if item["template_id"] == "missingEF")
            self.assertFalse(missing_template["fully_resolvable"])
            self.assertIn("contact_txt_or_json", missing_template["missing"])
            bad_template = next(item for item in report["templates"] if item["template_id"] == "badGH")
            self.assertFalse(bad_template["fully_resolvable"])
            self.assertIn("contact_json_invalid", bad_template["missing"])
            self.assertTrue(any(row["sha256"] for row in report["rows"] if row["template_id"] == "newCD"))

    def test_template_tsv_is_stable_and_has_hashes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "contacts").mkdir()
            (root / "contacts" / "1abcAB.txt").write_text("contact\n", encoding="utf-8")
            manifest = root / "list.txt"
            manifest.write_text("1abcAB\n", encoding="utf-8")
            report = preflight_template_assets(manifest, root)
            first = root / "first.tsv"
            second = root / "second.tsv"
            write_template_preflight_tsv(first, report)
            write_template_preflight_tsv(second, report)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with first.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["template_id"], "1abcAB")
            self.assertEqual(rows[0]["asset_type"], "contact_txt")
            self.assertEqual(rows[0]["exists"], "True")
            self.assertEqual(len(rows[0]["sha256"]), 64)

    def test_template_preflight_checks_every_interface_pdb(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "contacts").mkdir()
            (root / "interfaces").mkdir()
            (root / "contacts" / "1abcAB.json").write_text("{}\n", encoding="utf-8")
            (root / "interfaces" / "1abcAB_A_int.pdb").write_text("ATOM\n", encoding="utf-8")
            (root / "interfaces" / "1abcAB_B_int.pdb").write_text("", encoding="utf-8")
            manifest = root / "list.txt"
            manifest.write_text("1abcAB\n", encoding="utf-8")

            report = preflight_template_assets(manifest, root)

            assert report["valid"] == 0
            assert report["fully_resolvable"] == 0
            template = report["templates"][0]
            assert any("interface_pdb" in item for item in template["missing"])


if __name__ == "__main__":
    unittest.main()
