# Set up and validate the FiberDock/MultiProt CLI pipeline

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` aligned with the actual run.

## Purpose / Big Picture

Set up the repository-local legacy CLI under `working_version/Multiprot-new/prism-fiberdock-cli`, run the smallest documented end-to-end case that the available runtime supports, and report whether each stage reaches completion.

## Progress

- [x] Inspect instructions, inputs, tools, runtime, and working-tree constraints.
- [x] Verify or build required local tools without changing global environments.
- [x] Run a minimal CLI smoke case and retain reproducible logs.
- [x] Run available tests or stage-level checks and classify failures.
- [x] Record durable setup/run findings in project memory if meaningful.

## Surprises & Discoveries

- The CLI tree already contains extracted FiberDock, MultiProt, NACCESS, and BeEM assets despite the README describing archive extraction.
- The login shell reports a read-only `.condarc`; commands should use non-login shells and avoid Conda writes.
- The system exposes Python 3.10, while the CLI source uses Python 2 syntax and requires Python 2.7.
- The existing `tmalignRosetta` environment provides Python 2.7.15, allowing the CLI to start without environment creation.
- The legacy NACCESS `accall` requires `libgfortran.so.3`, unavailable on this host; MultiProt and NMA terminate with `Bad system call` before useful output.
- Recompiling NACCESS `accall` with `conda run -n tmalignRosetta gfortran accall.f -o accall -O` resolves the NACCESS library failure; MultiProt remains blocked independently.

## Decision Log

- Decision: Preserve the existing dirty worktree and use a new job ID for all smoke outputs.
  Rationale: The repository contains extensive pre-existing benchmark work and generated artifacts.
  Date/Author: 2026-07-17 / Codex
- Decision: Prefer repository-local compilation/checks; do not install packages or modify global configuration without explicit approval.
  Rationale: The task is setup/run/test, and the project instructions prohibit unapproved environment mutation.
  Date/Author: 2026-07-17 / Codex

## Outcomes & Retrospective

The local setup and stage checks succeeded, and the CLI reached all stages for the one-pair smoke. NACCESS now produces valid surfaces. No refined model was produced because bundled MultiProt still fails with `Bad system call`, leaving zero usable transformations.

## Context and Orientation

The entry point is `working_version/Multiprot-new/prism-fiberdock-cli/prism.py`. Configuration is copied from `prism.ini` into `jobs/<jobId>/`, and tool paths are relative to that job work directory. Inputs include `pairlist_mmcif`, `pairlist_new`, and `template_list`; bundled templates and cached PDBs are expected by the configured relative paths.

## Plan of Work

Inspect the target tree and README, verify interpreter and executable compatibility, compile only missing local helper binaries if feasible, run a one-pair/small-template smoke case, then execute repository-provided tests or stage checks. Keep all outputs under a unique `jobs/` or `tmp/agent/` path and summarize blockers separately from scientific pipeline failures.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, inspect `README_CLI.md`, `README_setup.md`, `prism.ini`, inputs, binaries, and available Python runtimes.
2. From the CLI directory, extract `template.zip`, restore executable bits, patch the NACCESS wrapper to use its own directory, and probe BeEM/FiberDock/NACCESS/MultiProt/NMA. Completed; MultiProt/NMA and NACCESS remain blocked by host compatibility.
3. From the CLI directory, run `/home/rshadi25/.conda/envs/tmalignRosetta/bin/python prism.py smoke_pairs template_smoke codex_smoke_20260717` with a 120-second timeout and capture the log. Completed; exit 0 with no refined model.
4. Run Python 2 byte-compilation and the bundled BeEM conversion example. Completed successfully.

## Validation and Acceptance

Acceptance requires an observable setup result for each required tool, a CLI invocation that either completes or reaches a clearly identified runtime/pipeline stage, and targeted validation output. Achieved for setup and stage reachability; scientific end-to-end acceptance is blocked by missing `libgfortran.so.3` and host rejection of the bundled 32-bit binaries.

## Idempotence and Recovery

Use a fresh job ID for reruns. Existing jobs and cached PDBs are preserved. Do not delete or overwrite existing outputs. If a setup command changes a local generated binary, record it and rebuild from the local source/archive as needed.

## Artifacts and Notes

Artifacts: `working_version/Multiprot-new/prism-fiberdock-cli/jobs/codex_smoke_20260717/` and `tmp/agent/20260717-fiberdock-cli/smoke.log`. Setup added the extracted `working_version/Multiprot-new/prism-fiberdock-cli/template/` tree and changed the local NACCESS wrapper path resolution.

## Interfaces and Dependencies

The CLI depends on Python 2.7, MultiProt, FiberDock, NACCESS or POPS, Perl helpers, and optionally the BeEM C++ converter for mmCIF fallback. Network access may be needed only when uncached structures are requested.
