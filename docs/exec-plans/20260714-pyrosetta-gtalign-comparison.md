# Install PyRosetta and compare isolated refinement/alignment arms

This ExecPlan is living documentation for the optional PyRosetta arm and the
GTalign comparison.

## Purpose / Big Picture

Stage a separately versioned PyRosetta environment if an authorized Linux
wheel is available, expose it through an opt-in PRISM adapter, and compare it
with the existing external-Rosetta and GTalign/TM-align arms using identical
inputs, templates, evaluator rules, and task isolation. If licensing or
package access prevents installation, preserve that as a confirmed blocker
and do not alter the stable default pipeline.

## Progress

- [x] Inspect current environment and project constraints.
- [x] Query project NotebookLM sources and official installation references.
- [x] Stage an authorized PyRosetta wheel/environment.
- [x] Run PyRosetta import and one-pose refinement smoke tests.
- [x] Add a separate PyRosetta comparison adapter/manifest; execution is gated on environment availability.
- [x] Compare GTalign and TM-align on matched alignment inputs.
- [x] Freeze the approved 12-row stratified pilot, 946-template panel, arm/task manifests, and fail-closed collector.
- [ ] Run isolated benchmark shards and collect DockQ, iRMSD, interface size, runtime, and prediction counts (pilot arrays 1358990/1358991 are running; PyRosetta matched quality remains untested).

## Distributed validation update — 2026-07-16

- Stale Cosbi chain `1361684`/`1363282` was cancelled before restarting the comparison.
- The current 2x2 launcher and capacity-probe helpers are implemented in `benchmark/scripts/submit_variant_matrix.sh`, `benchmark/scripts/prepare_parallel_capacity_probe.py`, and `benchmark/scripts/prepare_variant_smoke_batch.py`.
- PyRosetta and GTAlign CPU/GPU preflights pass. The first GPU smoke failed only because a relative `gtalign_gpu` name was not on the compute-node PATH; the launcher now passes an absolute executable path.
- Corrected 25-template GTAlign GPU smoke jobs `1363478/1363479` completed on `ai08` with exit code 0 but zero transformations. Full 946-template three-case validation jobs `1363490/1363491` are queued. TM-align jobs and the eight-task KUTEM capacity probe remain queued behind the fully allocated `rk01` node.
- TM-align CPU smoke jobs `1363506/1363507` completed on AI nodes in 146 s and 125 s, respectively, and each produced one scoreable model. The shared difficult multichain case scored cross-DockQ means `0.004818` (external Rosetta) and `0.004877` (PyRosetta), with grouped iRMSD `25.566` and `25.969 Å`. Staging/scoring was tightened after detecting and rejecting an overlapping model chain-group interpretation.
- GTAlign+PyRosetta full-template validation `1363491` completed in 1113 s with 14/14 scoreable models after set-specific native routing. Cross-DockQ mean was `0.015760` and grouped iRMSD mean `15.662 Å`; GTAlign+external-Rosetta `1363490` remains in refinement.
- [x] Update report, findings ledger, and project memory.

## Surprises & Discoveries

- The current host is `ai11.kuvalar.ku.edu.tr` outside Slurm; heavy alignment/refinement must be submitted to KUTEM or a V100 GPU job.
- `gtalign_env` has Python 3.11.13 but no importable `pyrosetta` module.
- NotebookLM’s PRISM sources do not contain PyRosetta or GTalign deployment information.
- Official PyRosetta downloads provide licensed quarterly wheels/packages; installation requires an authorized distribution.
- The West PyRosetta mirror exposed a 1.659 GB cp311 wheel, but the transfer was cancelled at approximately 203 kB/s before completion. The East mirror failed certificate validation and a bounded trusted-host retry timed out. No wheel or environment was modified.
- On 2026-07-15, the authorized PyRosetta distribution became available in `gtalign_env` as `2026.3+releasequarterly.5e498f1409`.
- PyRosetta initialization succeeded, and the corrected one-pose refinement smoke produced a non-empty output PDB with total score `290.7085762721651`.
- KUTEM job 1356100 completed the real `1fgnH` versus `1kcaCH` smoke in a copied job-unique workspace: two paired chain records, TM-align 0.031 s, GTalign CPU 1.804 s, mean absolute TM-score difference 0.01811. This is an alignment-stage result only.
- Synthetic matched pilots completed: 10x10x50 residues (100 pairs), TM-align 2.349 s versus GTalign CPU 8.115 s; 50x50x50 residues (2500 pairs), TM-align 47.019 s versus GTalign CPU 46.543 s. These runs contain no DockQ/iRMSD evaluation.

## Decision Log

- Keep PyRosetta opt-in and isolated from `environment.yaml` until a licensed wheel, exact version, Python ABI, and import smoke are recorded.
- Use the existing external Rosetta refinement as the baseline; do not compare energy magnitudes across refinement backends.
- Compare GTalign and TM-align first at the alignment/pose boundary, then at downstream scores only for matched successful tasks.
- Treat all missing or malformed outputs as explicit failures and retain null structural metrics rather than imputing quality.

## Outcomes & Retrospective

PyRosetta is installed as a separate opt-in refinement arm. The import/init probe
and one-pose smoke pass after adapting the DockingProtocol setters for the
2026.3 API. The stable external Rosetta path remains the baseline; PyRosetta is
verified only at one-pose scale so far.

GTalign CPU is operational on KUTEM and produces paired alignment records for the
real PRISM smoke. It is not yet a verified replacement for TM-align in the full
pipeline: the real smoke has no downstream transformation/refinement/evaluation
stage, and the synthetic runtime pilots show parity only at the 50x50 size.

## Context and Orientation

The default runtime is `/home/rshadi25/.conda/envs/gtalign_env`; DockQ scoring
uses `/scratch/tmp/prism-dockq-env/bin/python`. Existing adapters are under
`src/`, batch jobs under `benchmark/jobs/`, and the retained comparison outputs
are under `tmp/agent/20260714-observational-score-replay-*`.

## Plan of Work

First obtain and hash an authorized PyRosetta package on a setup-capable host.
Create a separate environment under the user-managed environment root, probe
import/version/init, and run a one-pose adapter smoke. Only after that passes,
run the matched alignment/refinement pilot. GTalign and TM-align must receive
the same receptor/ligand/template files and use the same evaluator and output
budget. All batch tasks receive unique output roots and write exit/provenance
records.

## Concrete Steps

1. On a login/setup host, stage the official PyRosetta wheel into the project’s
   derived agent directory and record its SHA256, license/version, Python ABI,
   and download command.
2. Create `pyrosetta_prism` without modifying `gtalign_env`, install the wheel,
   and run `benchmark/scripts/probe_pyrosetta_environment.py`.
3. Run a one-pose PyRosetta refinement smoke using the same canonical pose used
   by the external-Rosetta smoke; record output hash, runtime, and score.
4. Add a separate arm manifest and isolated Slurm array for matched model rows.
5. Run GTalign/TM-align matched alignment shards and collect alignment runtime,
   normalized TM-scores, downstream DockQ/iRMSD, interface metrics, failures,
   and prediction counts.

## Validation and Acceptance

- PyRosetta import and `pyrosetta.init()` succeed in the pinned environment;
  `gtalign_env` contains `2026.3+releasequarterly.5e498f1409`.
- The adapter writes a distinct one-pose PDB and metadata sidecar.  Chain
  preservation and quality remain to be tested on matched benchmark shards.
- External Rosetta, PyRosetta, GTalign, and TM-align tasks use identical input
  hashes and evaluator mappings where comparison is claimed.
- Every task has an isolated directory, command, parameters, logs, output
  hashes, runtime, resources, and exit status.
- No score is promoted from a malformed or incomplete model.
- GTalign claims are based on paired runtime/quality data, not separate aggregate
  runs.

## Idempotence and Recovery

Never overwrite a staged wheel or validated result. Use a new run ID for each
installation or batch attempt. Failed environments remain as evidence; remove
only agent-created disposable caches after generating a cleanup manifest.

## Artifacts and Notes

Planned artifacts include `environment-pyrosetta.yaml`, a runtime manifest,
wheel hash record, PyRosetta probe JSON, one-pose smoke logs, arm/task manifests,
and a comparison report linked from `docs/pipeline-validation-report-20260714.md`.

## Interfaces and Dependencies

PyRosetta official quarterly wheel/package; Python 3.11 ABI if available;
Biopython 1.84; existing external Rosetta module; GTalign 0.19.0; TM-align;
DockQ; KUTEM CPU arrays; and, if a GPU implementation is actually used, a
fixed Tesla V100 Slurm allocation with recorded CUDA visibility and utilization.

## Validation checkpoint — 2026-07-16

The current 2x2 validation matrix has completed generation: GTAlign GPU with
external Rosetta (`1363490`, 1973 s) and PyRosetta (`1363491`, 1113 s) used the
946-template panel; TM-align smoke arms (`1363506`, 146 s; `1363507`, 125 s)
used 25 templates. All four exited 0. Canonical bijective scoring produced
14/14 valid GTAlign models in each refiner arm. External Rosetta means are
DockQ `0.054998` / best `0.055490` / grouped iRMSD `14.161 Å`; PyRosetta means
are `0.015760` / `0.016155` / `15.662 Å`. TM-align smoke scores are one model
per arm and remain plumbing observations. The missing rigid native `1gla.pdb`
was handled in a derived temporary native root assembled from Benchmark 5.5
bound receptor/ligand files; raw benchmark data were not changed. KUTEM probe
`1363475` is still pending behind an allocated `rk01`, so parallel capacity is
unresolved. Do not scale to full Benchmark 5.5 until equivalent TM-align
full-template validation and the native-data repair are complete.
