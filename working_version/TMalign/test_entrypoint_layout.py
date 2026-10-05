#!/usr/bin/env python2.7
import os
import subprocess
import sys
import unittest


ROOT = os.path.abspath(os.path.dirname(__file__))
RUN_FILES = os.path.join(ROOT, "run_files")
if RUN_FILES not in sys.path:
    sys.path.insert(0, RUN_FILES)

from checkTemplate import TemplateChecker


class WorkingVersionEntrypointTests(unittest.TestCase):
    def test_slurm_entrypoint_resolves_project_root_from_script_location(self):
        script_path = os.path.join(RUN_FILES, "tmalignscript.sh")
        with open(script_path, "r") as script_file:
            script = script_file.read()

        self.assertNotIn("/scratch/users/rshadi25/hpc_run/fatma/prism-fixed-vFatma", script)
        self.assertNotIn("conda init bash", script)
        self.assertIn("project_root=", script)

    def test_root_entrypoint_prints_usage_without_arguments(self):
        process = subprocess.Popen(
            [sys.executable, os.path.join(ROOT, "prism.py")],
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        self.assertEqual(process.returncode, 0, stderr)
        self.assertIn("usage: python prism.py", stdout)

    def test_template_checker_reads_bundled_default_independent_of_cwd(self):
        original_cwd = os.getcwd()
        try:
            os.chdir("/")
            status, templates = TemplateChecker("unused", ["1b27AD"]).checker()
        finally:
            os.chdir(original_cwd)

        self.assertEqual(status, 2)
        self.assertEqual(templates, ["1b27AD"])


if __name__ == "__main__":
    unittest.main()
