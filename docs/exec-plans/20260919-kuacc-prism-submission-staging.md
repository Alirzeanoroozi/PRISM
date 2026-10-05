# Stage KUACC-specific PRISM submission wrappers

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`,
`Decision Log`, and `Outcomes & Retrospective` current while implementing.

## Purpose / Big Picture

Provide an isolated, path-correct submission layer for the two existing KUACC
portable PRISM bundles. The result is visible as three wrapper scripts under a
new `kuacc_submission_20260919` directory in each bundle. Existing bundles,
canonical checkouts, raw inputs, and prior outputs remain unchanged.

## Progress

- [x] Verify KUACC identity, partitions, queue state, requested directories,
  bundle hashes, executables, and wrapper layout.
- [x] Identify hard-coded original-checkout paths and invalid partition/node
  directives in the supplied portable Slurm examples.
- [x] Add isolated KUACC CPU, GTalign, and scoring wrappers with fail-closed
  runtime checks.
- [x] Validate local syntax, hashes, remote staging, and non-submitting Slurm
  checks.
- [x] Prepare separate CPU and DockQ environments after resolving the
  NumPy/DockQ dependency conflict.
- [x] Submit and validate a bounded one-template BM55 smoke, preserving the
  initial failed roots and a separate fix-1 manifest.
- [x] Fix the C++ runtime path and isolated USalign executable path exposed by
  the first smoke.
- [x] Validate the existing `gtalign_env` binary and submit the bounded 946-
  prefix CPU comparison in parallel with a one-template GTalign environment
  probe.

## Surprises & Discoveries

- Observation: Both remote bundles already contain the same fourteen supplied
  Slurm scripts and four Python wrappers as the validated local bundle.
  Evidence: exact SHA-256 comparison on 2026-09-19.
- Observation: Those supplied scripts reference `/scratch/rshadi25/GitHub/PRISM*`,
  which is absent on KUACC, while the bundle source is under
  `/scratch/users/rshadi25/<bundle>/code/`.
- Observation: KUACC exposes `short`, `mid`, `long`, `longer`, and `ai` among
  other partitions; `cosbi`, `kutem`, and `kutem_gpu` are invalid for the
  supplied scripts, and pinned nodes `ag01`/`rk02` are absent.
- Observation: the existing `gtalign_env` is Python 3.10 with NumPy 2.2.6
  and lacks pandas, Biopython, FreeSASA, and DockQ. The system
  `/kuacc/users/rshadi25/bin/python` is not a valid PRISM runtime.
- Observation: The portable bundle omitted `optional/gtalign_gpu`, while
  `prism-new` contained the verified executable. The GTalign staging layer
  therefore requires a narrowly scoped optional-asset copy.
- Observation: KUACC accepts account `users` with the `users` QOS on `mid`;
  explicit QOS `mid` is rejected even though the partition advertises it.
- Observation: the supplied single environment failed while building DockQ
  2.1.3 with GCC 4.8.5 and also conflicts with DockQ's NumPy `<2` requirement.
  FreeSASA 2.2.1 built successfully.
- Observation: the first smoke could not load the bundled TMalign or GTalign
  binaries with KUACC's system `libstdc++`; USalign also resolved against the
  Python bin directory rather than the staged external-tools directory.
- Observation: after the fix, TMalign and USalign each completed 40/40 smoke
  alignments, while MultiProt completed 2/40. All three had zero transformed
  pairs under the current downstream filter.
- Observation: GTalign's embedded CUDA targets include `sm_75`, but live KUACC
  inventory exposes Tesla K20/K40/K80 devices; fix-1 reached CUDA and failed
  with `cudaErrorInvalidSymbol`.
- Observation: `/kuacc/users/rshadi25/.conda/envs/gtalign_env/bin/gtalign_gpu`
  is present and executable, but has the same `sm_75` target. Probe job
  `3121151` reproduced `cudaErrorInvalidSymbol` on KUACC, so the conda
  environment does not provide a compatible build.

## Decision Log

- Decision: Add a new wrapper layer rather than edit or overwrite the supplied
  bundle Slurm files. Rationale: preserve the portable evidence package and
  make rollback a directory removal/quarantine operation.
- Decision: Use account `users`, partition `mid`, QOS `users`, and no fixed
  node. Rationale: live `scontrol`/`sinfo` and `sbatch --test-only` acceptance on
  2026-09-19; the supplied alternatives are invalid or unauthorized.
- Decision: Default the wrappers to the validated versioned CPU and DockQ
  prefixes while retaining `RUN_PYTHON` and `DOCKQ_PYTHON` overrides.
  Rationale: the available system and legacy environments are incomplete.
- Decision: Keep the CPU pipeline at
  `/kuacc/users/rshadi25/.conda/envs/prism_portable_20260919` with NumPy 2.4.6
  and FreeSASA 2.2.1, and DockQ scoring at
  `/kuacc/users/rshadi25/.conda/envs/prism_dockq_20260919` with NumPy 1.26.4.
  Rationale: preserve the pipeline's pinned numerical stack while satisfying
  DockQ 2.1.3's NumPy constraint.
- Decision: Add the selected environment `lib/` to `LD_LIBRARY_PATH` in the
  isolated wrappers and patch only the two isolated bundle helper copies to
  resolve USalign from `pipeline_repo/external_tools`. Rationale: the first
  smoke provided direct loader/path evidence, and the original failed roots
  were preserved.
- Decision: Do not resubmit GTalign on KUACC. Rationale: the binary contains
  `sm_75` CUDA targets, while current eligible KUACC GPUs are Tesla K20/K40/K80;
  a different binary or GPU backend is required.
- Decision: Submit only the bounded CPU 946-prefix arms for the current
  comparison and retain GTalign as a one-template compatibility probe.
  Rationale: the CPU requests passed live scheduler checks, while GTalign's
  actual conda binary failed on the available hardware.

## Outcomes & Retrospective

The three wrappers were staged under the framework and under both bundle roots.
The missing portable-bundle GTalign executable and metadata were copied from
the verified `prism-new` optional directory. Local and remote syntax checks,
hash checks, bundle-reference checks, and six `sbatch --test-only` checks
passed; `squeue -u rshadi25` remained empty. The split environments now
validate with user-site packages disabled: CPU imports include NumPy 2.4.6,
pandas 2.3.3, Biopython 1.84, and FreeSASA 2.2.1; DockQ imports and its CLI
work in the NumPy 1.26.4 scoring prefix. The first one-template smoke exposed
and the fix-1 submission corrected the C++ runtime and USalign path issues.
Jobs `3121137`–`3121140` are persisted under
`/scratch/users/rshadi25/valar-remote-runs/20260919-kuacc-smoke-fix1/`:
TMalign and USalign completed 40/40 alignments, MultiProt completed 2/40,
and all three produced zero transformed pairs. GTalign failed with
`cudaErrorInvalidSymbol` on a KUACC GPU incompatible with its `sm_75` binary.
Scientific production execution and promotion remain blocked.
The CPU 946-prefix jobs `3121148`–`3121150` are currently running under
`/scratch/users/rshadi25/valar-remote-runs/20260919-kuacc-bm55-946/`; their
final scientific state is not yet known.

## Context and Orientation

The framework is `/home/rshadi25/valar-agent-framework`; its existing isolated
alignment-only submitters are under `scripts/prism-prescript/`. The portable
bundles are `/scratch/users/rshadi25/prism-new` and
`/scratch/users/rshadi25/portable_prism_bundle_20260916`. Each bundle has
`code/PRISM`, `code/PRISM-prescript`, `templates`, `datasets`, and `optional`.
The new wrappers use `run_cpu_pipeline.py`, `run_gtalign_pipeline.py`, and
`score_transformed_models.py` from the bundle and write only to caller-owned
run roots.

## Plan of Work

First validate the additive local wrappers. Then create new remote staging
directories and copy only the declared wrapper files with `rsync` after a dry
run. Finally verify remote hashes, shell syntax, required bundle paths, and
Slurm directives using `sbatch --test-only` on `mid` with account `users`.
Do not submit, monitor, cancel, or score a scientific production run without
an exact manifest and a fresh bounded authorization. The completed fix-1
smoke is evidence only, not a promotion result.

## Concrete Steps

1. From `/home/rshadi25/valar-agent-framework`, run `bash -n` on the three
   wrappers and `git diff --check`.
2. From the same directory, use `rsync -av --dry-run` for the declared
   wrapper directory to the framework and both new bundle staging paths.
3. Repeat the same `rsync` without `--dry-run` only after the dry-run file list
   is exactly the four new wrapper files.
4. On KUACC, create only
   `/scratch/users/rshadi25/<bundle>/kuacc_submission_20260919` and the
   designated Slurm log directory, then verify hashes and paths.
5. Run `sbatch --test-only --export=NONE` on each staged wrapper. Confirm
   `squeue -u rshadi25` remains empty. **Completed:** all six checks passed
   on `mid` with account/QOS `users`.

## Validation and Acceptance

Acceptance requires: all three local wrappers pass `bash -n`; remote hashes
match; both bundle roots resolve all referenced helper/tool files; directives
are accepted under `users`/`mid` with QOS `users`; no test job remains queued;
and both versioned environments pass import/CLI checks with user-site
packages disabled. These staging conditions passed. The bounded fix-1 smoke
also validated TMalign and USalign executable readiness, while exposing the
MultiProt filtering outcome and the KUACC GTalign hardware incompatibility.
Production execution remains a separate authorization and scientific-
validation step.

## Idempotence and Recovery

The staging directories are new and contain no raw or validated results.
Rerunning `rsync -av` is additive and idempotent. Do not use `--delete`. If
validation fails, leave the staged files for inspection and do not submit a
job; remove/quarantine only with explicit authorization.

## Artifacts and Notes

Local source: `scripts/prism-prescript/kuacc/`.
Remote staging targets:
`/scratch/users/rshadi25/prism-new/kuacc_submission_20260919/` and
`/scratch/users/rshadi25/portable_prism_bundle_20260916/kuacc_submission_20260919/`.
Remote Slurm logs: `/scratch/users/rshadi25/valar-remote-runs/prism-submissions/`.
Environment-preparation logs:
`/scratch/users/rshadi25/valar-remote-runs/env-prep-20260919/conda-create.log`,
`freesasa-repair.log`, `dockq-conda-create.log`, and `dockq-pip-install.log`.
Prepared environments:
`/kuacc/users/rshadi25/.conda/envs/prism_portable_20260919` and
`/kuacc/users/rshadi25/.conda/envs/prism_dockq_20260919`.

## Interfaces and Dependencies

Wrappers call the bundle Python interfaces with explicit `--pipeline-repo`,
`--pipeline-python`, `--timed-runner`, and run-root arguments. They require a
working Slurm `sbatch`, account/QOS authorization, the versioned CPU runtime
with `freesasa`, the bundle alignment tools, and the separate DockQ runtime
for scoring. No scheduler abstraction is introduced.

## Downstream refinement and scoring addendum — 2026-09-19

The active 946-prefix alignment roots are preserved and are not recomputed.
Three dependency-gated children were submitted after fresh `sbatch
--test-only` checks with KUACC's one-day `mid` limit:

- TMalign parent `3121148` -> successful child `3121164`
- MultiProt parent `3121149` -> corrected replacement child `3121167`
- USalign parent `3121150` -> corrected replacement child `3121168`

The additive wrappers are `refine_and_score.sbatch`,
`refine_and_score.py`, and `pyrosetta_refine_one.py`. They require the parent
`run_summary.json` to be completed and ranking-disabled, use external Rosetta
for the explicit refinement stage, probe PyRosetta through a supplied
interpreter, run the separate DockQ environment, and preserve atomic status,
timing, event, and per-candidate checkpoint files. No PRODIGY ranking is
included. A PyRosetta import failure is recorded as unavailable and does not
silently switch refiner identity.

The submission manifest is
`/scratch/users/rshadi25/valar-remote-runs/20260919-kuacc-bm55-946/downstream_submission_manifest.json`.
The final interpretation remains pending parent completion, child execution,
artifact reconciliation, and scientific validation. The scorer currently
scores transformed models and the additive scorer writes refined-stage DockQ
separately; both stages must remain explicitly separate in the final
comparison.

TMalign parent `3121148` completed with 37,840/37,840 successful alignment
records and 3 transformed pairs. Successful downstream child `3121164`
scored all 3 transformed rows, refined 2/3 through external Rosetta, and
scored 2/2 refined structures. The explicit PyRosetta probe returned
`ImportError: No module named pyrosetta`; no PyRosetta fallback was used. The
initial children `3121158` and `3121162` are retained as wrapper-failure
evidence, not scientific failures.

The first submitted MultiProt/USalign children (`3121159`, `3121160`) were
superseded before release because Slurm had already snapshotted the pre-fix
wrapper. Corrected replacements `3121167` and `3121168` were submitted after
fresh test-only checks; no prior job was canceled or modified.

USalign parent `3121150` completed with 37,840 successful alignments and zero
transformed pairs. Child `3121168` completed with valid zero-candidate
transformed/refined DockQ summaries and explicit PyRosetta unavailability.
MultiProt parent `3121149` completed with 2,621 successful and 35,219 failed
alignment records, yielding 132 transformed pairs; child `3121167` remains in
checkpointed external-Rosetta refinement.

## Batch-processing addendum — 2026-09-19

The batch path is implemented in `prepare_refinement_batches.py`,
`refine_one_external.py`, and `refine_batch.sbatch`. It uses one candidate per
array task, isolated Rosetta working directories, the same checkpoint keys as
the sequential aggregator, and a bounded array concurrency of 16. Manifests
were generated for all three roots (3, 132, and 0 candidates), and fresh
array `sbatch --test-only` checks passed without leaving test jobs queued.

The active MultiProt child `3121167` was not overlapped by an array because it
was already writing the same checkpoint/event namespace. Actual replacement of
that active job would have required explicit cancellation/replacement
authorization; the safe dependency-gated continuation preserved the existing
run instead.

A non-overlapping continuation was submitted without cancellation: array
`3121175` (`0-131%16`) waited for `3121167`, and aggregator `3121176` waited
for the array. All 132 array tasks completed with `resumed` status and the
aggregator exited `0:0`. Submission metadata is retained at
`/scratch/users/rshadi25/valar-remote-runs/20260919-kuacc-bm55-946/multiprot/downstream/batch_submission_manifest.json`.

Final validation found 132 refinement checkpoints (3 completed and 129
terminal-no-output), 132/132 transformed DockQ rows scored, and 3/3 refined
external-Rosetta DockQ rows scored. The authoritative parent timing was
3,233.03 s total: 2,741.59 s refinement, 472.82 s transformed DockQ, and
16.63 s refined DockQ. Ranking and PRODIGY were disabled. PyRosetta was
probed and recorded unavailable because the supplied interpreter raised
`ImportError: No module named pyrosetta`.

## Full BM55 all-template batch addendum — 2026-09-19

### Progress

- [x] Add and test the full-panel case-batch contract, including exact
  template/input hashes, query-source checks, attempt roots, resume states,
  per-case timing, and collision-safe aggregation.
- [x] Prepare 257-case manifests for TMalign, MultiProt, and USalign using
  all 19,855 entries from the frozen full template list.
- [x] Diagnose and preserve three pilot failure classes: Slurm spool-path
  resolution, attempt-root controller-log collision, and native-only case
  staging that omitted downloader query PDBs.
- [x] Correct case staging to include the native complex and the stripped
  four-character receptor/ligand source PDBs.
- [x] Submit full arrays and dependency-gated aggregators: `3121371`/`3121372`,
  `3121373`/`3121374`, and `3121375`/`3121376`.
- [ ] Complete execution, artifact validation, aggregation, refinement, and
  transformed/refined DockQ scientific validation.
- [ ] Resolve the PyRosetta transfer/environment authorization conflict.

### Contract evidence

The full list has 19,855 unique IDs and SHA-256
`4680d3eda8030861a40373cd193b0e8bef7c21a771a90553c4c49814e964b48d`. The
current PRISM alignment loop creates two query-by-template-chain records per
query case, yielding 79,420 raw records per case and 20,410,940 per aligner.
Ranking, PRODIGY, and refinement are disabled in the alignment arrays.

### Current risk

KUACC live QOS permits six 4-CPU tasks for this user at present, so `%16` is a
maximum rather than a promise; queued array work is pending with
`QOSMaxCpuPer...`. Do not cancel or duplicate the submitted arrays. Job
`3121364` is unrelated and must remain untouched. The remote
`scripts/kuacc-gpu-status.sh` required by the contract was not found, although
`squeue`, `sinfo`, and `scontrol` discovery succeeded.
