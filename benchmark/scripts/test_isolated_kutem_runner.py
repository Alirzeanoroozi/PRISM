#!/usr/bin/env python3
"""Focused tests for the isolated KUTEM execution design."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import isolated_kutem_runner as runner
import prepare_investigation_tasks


class IsolatedKutemRunnerTests(unittest.TestCase):
    def make_manifest(self, root: Path, command: str = "printf ok > result.txt", array_size: int = 10) -> Path:
        config = root / "config.ini"
        input_file = root / "input.dat"
        config.write_text("setting=value\n", encoding="utf-8")
        input_file.write_text("input\n", encoding="utf-8")
        manifest = root / "tasks.csv"
        fields = [
            "array_index",
            "task_id",
            "command",
            "config_paths",
            "input_paths",
            "output_paths",
            "scientific_retry_id",
        ]
        with manifest.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for index in range(1, array_size + 1):
                writer.writerow(
                    {
                        "array_index": index,
                        "task_id": f"pair-{index:02d}",
                        "command": command,
                        "config_paths": "config.ini",
                        "input_paths": "input.dat",
                        "output_paths": "result.txt",
                        "scientific_retry_id": f"sci-{index}",
                    }
                )
        return manifest

    def test_task_directories_are_isolated_and_outputs_stay_inside(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.make_manifest(root)
            run_root = root / "run"
            first = runner.task_directory(run_root, 1)
            second = runner.task_directory(run_root, 2)
            self.assertNotEqual(first, second)
            self.assertEqual(first.name, "task_0001")
            self.assertEqual(second.name, "task_0002")
            output = runner.resolve_output_path("nested/result.txt", first)
            self.assertEqual(output, first / "nested" / "result.txt")
            with self.assertRaises(ValueError):
                runner.resolve_output_path("../escape.txt", first)
            self.assertEqual(manifest.name, "tasks.csv")

    def test_task_execution_writes_exit_schema_and_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = self.make_manifest(root)
            run_root = root / "run"
            env = {
                **os.environ,
                "SLURM_JOB_ID": "12345",
                "SLURM_ARRAY_JOB_ID": "12345",
                "SLURM_ARRAY_TASK_ID": "1",
                "SLURM_PARTITION": "kutem",
                "SLURM_ACCOUNT": "kutem",
                "SLURM_QOS": "kutem",
                "SLURM_CPUS_PER_TASK": "2",
                "SLURM_MEM_PER_NODE": "2G",
                "SLURM_NNODES": "1",
                "SLURM_NTASKS": "1",
                "SLURM_RESTART_COUNT": "3",
            }
            status, exit_path = runner.execute_task(
                manifest=manifest,
                run_root=run_root,
                array_index=1,
                env=env,
            )
            self.assertEqual(status, 0)
            record = json.loads(exit_path.read_text(encoding="utf-8"))
            self.assertEqual(record["schema_version"], "isolated-kutem-exit/v1")
            self.assertEqual(record["termination"]["return_code"], 0)
            self.assertIsNone(record["termination"]["signal"])
            self.assertIsNone(record["scientific_result"]["pair_success"])
            self.assertIn("batch_completion", record["scientific_result"]["reason"])
            self.assertEqual(record["retry_ids"], {"scientific_retry_id": "sci-1", "scheduler_retry_id": "12345:restart:3"})
            self.assertEqual(record["slurm"]["requested"], runner.PROFILE)
            self.assertEqual(record["slurm"]["observed"]["ids"]["SLURM_JOB_ID"], "12345")
            self.assertEqual(record["command"]["sha256"], runner.sha256_bytes(b"printf ok > result.txt"))
            self.assertEqual(record["hashes"]["command_sha256"], record["command"]["sha256"])
            self.assertEqual(record["hashes"]["config_sha256"], record["artifacts"]["config"]["aggregate_sha256"])
            self.assertEqual(record["hashes"]["input_sha256"], record["artifacts"]["inputs"]["aggregate_sha256"])
            self.assertEqual(record["artifacts"]["inputs"]["files"][0]["sha256"], runner.sha256_file(root / "input.dat"))
            output_record = record["artifacts"]["outputs"]["files"][0]
            self.assertTrue(output_record["exists"])
            self.assertEqual(output_record["sha256"], runner.sha256_file(Path(output_record["resolved_path"])))
            self.assertEqual(record["hashes"]["output_sha256"], record["artifacts"]["outputs"]["aggregate_sha256"])
            self.assertTrue(record["timestamps"]["started_at"].endswith("Z"))
            self.assertTrue((exit_path.parent / "stdout.log").is_file())
            self.assertTrue((exit_path.parent / "stderr.log").is_file())

    def test_dry_run_is_non_submitting_and_does_not_execute_commands(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            marker = root / "should-not-exist"
            manifest = self.make_manifest(root, command=f"touch {marker}")
            run_root = root / "run"
            with mock.patch.object(runner.subprocess, "run", side_effect=AssertionError("sbatch must not run")):
                plan = runner.dry_run_plan(manifest.resolve(), run_root.resolve(), Path("template.sbatch").resolve())
            self.assertFalse(plan["submits_job"])
            self.assertEqual(len(plan["tasks"]), 10)
            self.assertFalse(marker.exists())
            self.assertFalse((run_root / "tasks").exists())

    def test_partial_array_manifest_is_exactly_sized_and_submits_partial_array(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = prepare_investigation_tasks.build_task_rows(
                [f"rigid:{index:06d}" for index in range(1, 8)],
                "printf ok > result.txt",
            )
            self.assertEqual(len(rows), 7)
            manifest = root / "tasks.csv"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=prepare_investigation_tasks.FIELDS)
                writer.writeheader()
                writer.writerows(rows)
            loaded = runner.load_manifest(manifest, array_size=7)
            self.assertEqual(len(loaded), 7)
            plan = runner.dry_run_plan(manifest, root / "run", root / "template.sbatch", array_size=7)
            self.assertEqual(plan["array_size"], 7)
            self.assertIn("--array=1-7", plan["submission_command"])
            with self.assertRaises(ValueError):
                prepare_investigation_tasks.build_task_rows(["rigid:000001"] * 11, "printf ok")


if __name__ == "__main__":
    unittest.main()
