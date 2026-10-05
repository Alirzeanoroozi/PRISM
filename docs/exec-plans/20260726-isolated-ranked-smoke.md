# Validate Opt-in Ranking with an Isolated Paired Smoke

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`,
`Decision Log`, and `Outcomes & Retrospective` current as execution proceeds.

## Purpose / Big Picture

Demonstrate that the repaired opt-in ranking branch reduces refinement work in
a real current-pipeline run while leaving the stable unranked pipeline and its
artifacts untouched. Run one frozen case twice under stable thresholds: a
baseline arm that refines every transformed candidate and a ranked arm that
selects a deterministic top-K subset.

## Progress

- [x] Resolve the latest authoritative state from project memory and the scoped chronology graph.
- [x] Identify a stable-threshold case with multiple previously generated candidates.
- [x] Create and statically validate the isolated paired Slurm launcher.
- [x] Submit the initial top-2 job and wait for terminal Slurm and pipeline-stage status.
- [x] Submit the evidence-driven top-1 follow-up and wait for terminal status.
- [x] Compare candidate, selection, refinement, and provenance artifacts.
- [x] Update active memory and rebuild the scoped chronology graph.

## Surprises & Discoveries

- Observation: the repository-wide Graphify graph still uses the pre-#1504
  node-ID scheme and cannot provide authoritative chronology. The scoped graph
  under `docs/chronology/graphify-out/` correctly anchors the current state to
  2026-07-26, while detailed next actions remain authoritative in project
  memory.
- Observation: the one-template stable smoke (`1kcaCH`) does not produce a
  candidate for `1FGNH/1TFHA`, so it cannot exercise ranking selection.
- Observation: retained evidence for the same input pair shows four generated
  candidates from `1h5bAB`, `2f0xEH`, and `3lqmAB` under stable thresholds.
- Observation: job 1392708 produced two, rather than four, candidates with the
  current frozen source and inputs. Both baseline and top-2 ranked arms
  completed with two refined structures, so this run verifies ranking plumbing
  and paired reproducibility but cannot demonstrate resource reduction.
- Observation: external Rosetta produced different accepted-energy rows across
  paired arms despite identical pre-ranking candidates. Rosetta refinement is
  stochastic in this launcher, so the smoke supports an execution-load claim,
  not a paired model-quality claim.

## Decision Log

- Decision: use `1FGNH/1TFHA` and the three-template panel
  `1h5bAB,2f0xEH,3lqmAB` with unchanged stable thresholds.
  Rationale: this is an existing multi-candidate case and avoids treating the
  relaxed `1BPB/3K77` diagnostic as production evidence.
  Date/Author: 2026-07-26 / Codex
- Decision: compare independent baseline and ranked arms in a new Slurm run
  directory; make top-K an explicitly recorded `PRISM_SMOKE_TOP_K` parameter.
  Rationale: a paired run directly exposes resource reduction and prevents
  ranked outputs from overwriting baseline artifacts.
  Date/Author: 2026-07-26 / Codex
- Decision: copy the current `prism.py` and `src/` into each arm, stage local
  benchmark PDBs, and record hashes before execution.
  Rationale: the job must not download inputs at runtime or depend on later
  source edits while queued.
  Date/Author: 2026-07-26 / Codex

## Outcomes & Retrospective

Initial execution completed as Slurm job 1392708 (`COMPLETED`, exit `0:0`,
elapsed 00:02:30). Both arms returned zero, produced identical six-record
candidate audits with two generated candidates, completed all expected stages,
and generated two refined structures. Source and input manifests matched
byte-for-byte. The ranked arm completed its ranking stage and selected two of
two at top-2. This is valid end-to-end plumbing evidence, but not reduction
evidence; it motivated the completed top-1 follow-up below.

The top-1 follow-up completed as Slurm job 1392725 (`COMPLETED`, exit `0:0`,
elapsed 00:01:54). Both arms returned zero and had byte-identical source,
input, and six-record candidate-audit manifests. Each produced two transformed
candidates. Baseline refined two; ranking selected one and produced one refined
structure, choosing `1h5bAB/o1` (deterministic score 0.681002) over
`3lqmAB/o2` (0.568371). All expected stage events completed. Every refined PDB
contained 6,419 ATOM records across chains A and H. This satisfies the
resource-reduction acceptance contract while leaving stable launchers and
pipeline defaults unchanged.

The chronology builder was extended to inspect only bounded
`runs/<job-id>/results.tsv` aggregate files. Two focused tests pass, and the
rebuilt scoped graph records this run as `recorded:success` while avoiding the
copied per-arm source trees.

## Context and Orientation

Ranking is implemented in `prism.py` after transformation and before the
selected refiner. `src/candidate_selector.py` maps transformed pairs to the
fresh JSONL audit and applies deterministic baseline scoring. The repaired
branch is disabled by default and fail-open when audit coverage is incomplete.

The experiment launcher is
`tmp/agent/20260726-isolated-ranked-smoke/run_paired_ranked_smoke.sbatch`.
Each Slurm job writes only beneath
`tmp/agent/20260726-isolated-ranked-smoke/runs/<job-id>/`. Set
`PRISM_SMOKE_TOP_K` at submission time; its default is two and the resolved
value is recorded in `provenance/runtime.tsv`.

## Plan of Work

Stage identical source, input PDB, and template-interface snapshots for two
arms. Record source/input hashes, executable paths, Python version, module and
Slurm metadata, commands, candidate audits, and stage-status JSONL. Run the
baseline without ranking and the ranked arm with the requested top-K. Summarize
generated-candidate, transformed-pair, selected-pair, refined-structure, and
return-code counts without interpreting directory existence as success.

## Concrete Steps

1. Validate the launcher with `bash -n` and rerun the focused ranking tests.
2. Submit with `sbatch --export=ALL,PRISM_SMOKE_TOP_K=<K> tmp/agent/20260726-isolated-ranked-smoke/run_paired_ranked_smoke.sbatch`.
3. Monitor with `squeue -j <job-id>` and confirm terminal state with
   `sacct -j <job-id> --format=JobID,State,ExitCode,Elapsed,NodeList`.
4. Inspect `results.tsv`, both `status/stages.jsonl` files, both candidate
   audits, logs, structures, and provenance manifests.
5. Update this plan, active memory, and `docs/chronology/` only from terminal,
   reproducible evidence.

## Validation and Acceptance

Acceptance requires both pipeline commands to return zero; each stage-status
file to contain completed input, alignment, transformation, and refinement
events; the ranked arm additionally to contain a completed ranking event; and
the ranked selected/refinement-attempt count to equal the requested K when the
baseline has more than K generated transformed pairs. Input and source hashes
must be present for both arms. Any difference in pre-ranking candidate
generation is a failed paired comparison and must be investigated rather than
rationalized.

## Idempotence and Recovery

Every submission uses its Slurm job ID as a new output directory. The launcher
does not delete, modify, or reuse prior runs. A failed arm retains its log,
return code, audit, stage events, and partial outputs for diagnosis; rerun by
submitting a new job rather than editing or cleaning the failed run.

## Artifacts and Notes

- Stable operational guide: `docs/STABLE_PIPELINE.md`
- Ranking repair plan: `docs/exec-plans/20260725-ranking-integration-and-memory.md`
- Historical multi-candidate evidence:
  `tmp/agent/20260711-biological-ranking/current-1ahw-candidates.csv`
- Active memory: `.agents/skills/project-memory/references/`

## Interfaces and Dependencies

- Pipeline interpreter: `/home/rshadi25/.conda/envs/gtalign_env/bin/python`
- Rosetta module: `rosetta/2022.42`
- Aligner: stable default TMalign
- Refiner: stable default external Rosetta
- Ranking CLI: `--rank true --top-k "$PRISM_SMOKE_TOP_K" --rank-min-score 0.0`
