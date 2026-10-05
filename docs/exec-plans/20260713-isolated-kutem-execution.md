# Add isolated KUTEM array execution design

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` aligned with the implementation.

## Purpose / Big Picture

Add a standalone, manifest-driven Slurm execution surface for ten isolated KUTEM tasks. The design will be usable locally in dry-run mode, and the batch template will request exactly one node, one task, two CPUs, 2G memory, five minutes, and array indices 1–10 with partition/account/QOS `kutem`. Each task will have independent paths and a machine-readable `exit.json`; the runner will never turn batch completion into a scientific pair-success claim.

## Progress

- [x] Route to `PRISM-prescript` and read project guidance/memory without changing memory files.
- [x] Add the standard-library runner and exact KUTEM array template.
- [x] Add focused tests for path isolation and the exit schema.
- [x] Run shell syntax, dry-run, and focused tests without submitting Slurm jobs.
- [x] Review the critical changed regions and report caveats.

## Surprises & Discoveries

- Observation: The repository has extensive pre-existing dirty and untracked state.
  Evidence: `git status` shows many unrelated modified/untracked files.
  Consequence: only new isolated execution files will be added.
- Observation: Existing benchmark Slurm scripts use varied resources and are production benchmark entry points.
  Evidence: `benchmark/scripts/submit_comparison_batches.sbatch` and `benchmark/scripts/rosetta_output/submit_prism_analysis_all.sbatch`.
  Consequence: this design will not edit or wrap those scripts.

## Decision Log

- Decision: Use CSV with explicit `array_index` values 1–10 and task-local relative output paths.
  Rationale: explicit indexing avoids row-order ambiguity, and path constraints make isolation testable.
- Decision: Keep execution and submission concerns in one standalone Python runner, with `--dry-run` and explicit `--submit` modes.
  Rationale: local validation must not submit jobs accidentally; submission remains opt-in.
- Decision: Use Python `resource.getrusage(RUSAGE_CHILDREN)` for optional RSS measurement.
  Rationale: it is dependency-free and reports the child process maximum resident set size when available.
- Decision: Store scientific and scheduler retry IDs separately and set scientific pair success to an explicit unknown status.
  Rationale: retries and scheduler outcomes are operational metadata, not scientific evaluation results.

## Outcomes & Retrospective

The standalone runner, template, and focused tests are complete. Local shell syntax, Python compilation, three focused unit tests, and a CLI dry-run check pass. No Slurm job was submitted. The design remains intentionally separate from production manifests, evaluators, pipeline code, and project memory.

## Context and Orientation

The target repository is `/scratch/rshadi25/GitHub/PRISM-prescript`. General benchmark utilities live in `benchmark/scripts/`, and existing Slurm entry points live in `benchmark/scripts/` or `benchmark/jobs/`. This change is intentionally separate from manifests used by evaluators and production pipeline code.

The new runner will consume a CSV with columns `array_index`, `task_id`, `command`, `config_paths`, `input_paths`, `output_paths`, and optional `scientific_retry_id`. Paths in `config_paths` and `input_paths` are resolved relative to the manifest directory; relative output paths are resolved inside the task directory. Each row is executed with the task directory as its working directory.

## Plan of Work

Add `benchmark/scripts/isolated_kutem_runner.py` with manifest validation, task-directory derivation, deterministic file hashing, command execution, RSS collection, Slurm metadata capture, exit-record writing, dry-run planning, and opt-in `sbatch` submission. Add `benchmark/jobs/isolated_kutem_array.sbatch` containing the exact confirmed resource directives and invoking the runner for the current array index. Add a focused unittest module covering unique task paths, output containment, dry-run non-submission, and required exit-record fields.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, add the runner, template, tests, and this plan using `apply_patch`.
2. Run `bash -n benchmark/jobs/isolated_kutem_array.sbatch`.
3. Run the focused unittest module with the repository's available Python interpreter; the test must use only the standard library.
4. Run the runner's local `--dry-run` path against a temporary ten-row manifest and verify that no `sbatch` process is invoked and no task command runs.
5. Inspect the new-file diff only; do not stage, submit, or alter unrelated worktree state.

## Validation and Acceptance

Acceptance requires all of the following:

- The template declares `--partition=kutem`, `--account=kutem`, `--qos=kutem`, `--nodes=1`, `--ntasks=1`, `--cpus-per-task=2`, `--mem=2G`, `--time=00:05:00`, and `--array=1-10`.
- A valid manifest has exactly one row for each array index 1–10; missing or duplicate indices fail closed.
- Different indices resolve to different task directories under the run root, and relative output paths cannot escape them.
- `exit.json` contains timestamps, command/config/input/output hashes, return code, signal, RSS when available, Slurm IDs/resources, and separate retry IDs.
- The record states that scientific pair success is unknown and is not derived from batch completion.
- Focused tests pass; shell syntax passes; no Slurm job is submitted.

## Idempotence and Recovery

Each array index maps deterministically to `tasks/task_<index:04d>`. Re-running an index updates only that task directory's execution artifacts. Existing task directories are rejected unless `--allow-existing-task-dir` is supplied, preventing accidental mixing of retries. Scheduler retries should use a new run root or an explicit scheduler retry ID; scientific retries remain the manifest-level ID.

## Artifacts and Notes

Expected new files are:

- `benchmark/scripts/isolated_kutem_runner.py`
- `benchmark/jobs/isolated_kutem_array.sbatch`
- `benchmark/scripts/test_isolated_kutem_runner.py`

No Slurm output is expected because this task does not submit jobs. Local test temporary directories are created by the test framework and removed automatically.

## Interfaces and Dependencies

The runner uses Python 3 standard-library modules only: `argparse`, `csv`, `hashlib`, `json`, `os`, `pathlib`, `resource`, `shlex`, `signal`, `subprocess`, and `datetime`. The batch template expects `MANIFEST` and `RUN_ROOT` exported via `sbatch --export=ALL,...`; it invokes the runner with `SLURM_ARRAY_TASK_ID`. The submission helper uses the exact template path and explicit `sbatch` invocation but is never called in local dry-run mode.
