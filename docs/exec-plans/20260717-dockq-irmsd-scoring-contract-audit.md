# DockQ and iRMSD scoring-contract audit

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` up to date as work proceeds.

## Purpose / Big Picture

Establish whether current PRISM benchmark scores are scientifically valid by tracing model identities and partner chains from the batch manifest through refinement, staging, DockQ, and iRMSD. The result is an evidence-backed accept/reject decision for each scoring path; no production score will be changed or recomputed in place.

## Progress

- [x] Recover current project memory, benchmark decisions, and existing scoring-contract evidence.
- [x] Confirm the retained GTAlign/PyRosetta batch artifacts and its ad-hoc score output.
- [x] Trace each canonical and ad-hoc scoring path end-to-end, including chain and residue correspondence.
- [x] Compare the legacy iRMSD implementation with the guarded evaluator and assess the DockQ invocation mode.
- [x] Publish an accept/reject matrix and record durable conclusions.

## Surprises & Discoveries

- Observation: all 18 retained PyRosetta model filenames are explicit `unparseable_current_model_name` rows in the canonical stager, despite matching a batch input row.
  Evidence: `tmp/agent/20260717-benchmark55/scoring/gt_pyro/stage.csv`.
- Observation: the ad-hoc score used DockQ `--no_align` without recording the complete mapping or raw JSON; its reported per-interface iRMSD is not the grouped benchmark iRMSD.
  Evidence: `tmp/agent/20260717-benchmark55/scoring/score_final.py` and `scores_final.csv`.
- Observation: model `1p2cDF_1mlbAB_3lzt_o1` passes the strict raw chain and residue contract for `ABC:ABE`; legacy and paired iRMSD both equal `2.887 A` over 98 interface residues. In contrast, `2gk2AB_1fgnHL_1tfhA_o1` fails strict `--no_align` residue correspondence but legacy and paired sequence-aligned iRMSD agree at `17.584 A` over 129 interface residues.
  Evidence: direct read-only probes on 2026-07-17.

## Decision Log

- Decision: audit retained artifacts and source code before proposing any scoring repair.
  Rationale: benchmark decisions require a demonstrably valid mapping and correspondence contract.
  Date/Author: 2026-07-17 / Codex.

## Outcomes & Retrospective

The `score_final.py` path is rejected for confirmatory use: it invokes DockQ with automatic mapping and `--no_align`, preserves no raw JSON or complete selected mapping, and reports a best component interface rather than grouped benchmark metrics. The robust single-pair scorer is valid only when passed explicit complete groups; it enforces a strict correspondence check before allowing `--no_align` and retains DockQ JSON. The canonical bijective scorer is the acceptable benchmark path because it passes an explicit full mapping, keeps raw DockQ JSON, selects only requested receptor-ligand interfaces, and records grouped forward/reverse iRMSD. It cannot currently score PyRosetta output because the staging filename adapter rejects every current PyRosetta filename. Therefore no retained `gt_pyro` score is confirmatory until that adapter is repaired and the canonical scorer runs successfully on a compute node.

## Context and Orientation

The relevant scripts are `benchmark/scripts/stage_current_models_for_main_benchmark.py`, `benchmark/scripts/score_bijective_benchmark_models.py`, `benchmark/scripts/score_single_prism_pair.py`, `benchmark/scripts/standardized_evaluator.py`, and `benchmark/scripts/irmsd.py`. The retained smoke data live under `tmp/agent/20260717-benchmark55/`. The benchmark complex column defines native receptor and ligand groups; quality metrics are valid only when model partner groups map bijectively to those native groups and the metric's correspondence assumptions are met.

## Plan of Work

Trace the input manifest into transformation and refinement filenames, reconstruct the intended partner chain groups, and compare them with what each scorer passes to DockQ and iRMSD. Read DockQ's installed implementation for `--no_align` behavior. Test only deterministic, read-only microcases or disposable `/tmp` artifacts. Separate evidence for pipeline execution from evidence for metric validity.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, inspect `batch_0001/inputs.csv`, stage manifests, exit records, refinement metadata, and model chain IDs.
2. Read the full score and evaluator code; compare explicit mapping, auto-mapping, native filtering, and iRMSD correspondence behavior.
3. Run syntax/static checks and, where permitted, a minimal score reproduction writing only to `/tmp`.
4. Produce a score-path matrix with acceptance criteria and blockers.

## Validation and Acceptance

An accepted score path must retain the source model and native paths, benchmark row identity, explicit complete model:native chain mapping, raw DockQ JSON, requested cross-interface components, grouped iRMSD, and a successful correspondence validation. Any path without these records is exploratory only.

## Idempotence and Recovery

All inspection commands are read-only. Any temporary DockQ input/output belongs under `/tmp` and can be removed without affecting benchmark artifacts. No model, native, manifest, or validated result is overwritten.

## Artifacts and Notes

- `tmp/agent/20260717-benchmark55/scoring/gt_pyro/scores_final.csv`
- `tmp/agent/20260717-benchmark55/scoring/gt_pyro/stage.csv`
- `docs/exec-plans/20260716-score-contract-and-cleanup.md`

## Interfaces and Dependencies

The scorer uses `/scratch/tmp/prism-dockq-env/bin/python -m DockQ`; in the Codex sandbox DockQ's multiprocessing manager may fail due to socket permissions, so source and retained Slurm outputs are the primary evidence. The benchmark iRMSD script uses Biopython sequence alignment and index-coupled residue lists. PyRosetta outputs originate from `src/pyrosetta_refinement.py`, which uses `src.rosetta_refinement.combine_pdb` for partner assembly.
