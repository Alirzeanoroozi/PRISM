# KUACC Rosetta, FiberDock, and MultiProt preflight

## Purpose

Prepare a bounded, isolated KUACC preflight for the PRISM-prescript MultiProt
alignment and Rosetta/FiberDock refinement dependencies. Confirm syntax,
paths, executable/runtime compatibility, input contracts, and one small smoke
execution before any production refinement or DockQ submission.

## Progress

- [x] Read the selected project contract and applicable execution/storage/
  verification procedures.
- [x] Confirmed the KUACC `rosetta/2022.42` module resolves the two Rosetta
  docking executables.
- [x] Located the KUACC MultiProt executable in the PRISM submission tree.
- [x] Inspect and stage the executable and FiberDock bundle; retain only the
  required legacy MultiProt runtime subset rather than the 297 MB workspace
  archive/database payload.
- [x] Complete the compact legacy runtime subset.
- [x] Verify FiberDock dependencies and execution directory contract.
- [x] Run isolated command and compute-node smoke tests on KUACC.
- [x] Validate evidence and report production submission status.

## Surprises

- The Rosetta module sets PATH entries but does not expose `ROSETTA3_DB`.
  The database path must therefore be supplied explicitly by the wrapper.
- The initial KUACC scan did not find a FiberDock module or installation in
  the searched application/user roots; this remains an open discovery item.
- The KUACC MultiProt executable is present under
  `/scratch/users/rshadi25/prism-new/code/PRISM-prescript/external_tools`.
- The full legacy MultiProt source workspace is 297 MB and contains large
  archive/database files. A partial transfer was stopped before use; the
  preflight will use a separately staged runtime-only subset and preserve the
  partial directory as untrusted evidence.
- The existing KUACC `prism-new` submission tree does not contain the
  `refine_one_external.py` wrapper expected by the Rosetta batch script. The
  wrapper must be staged explicitly before a future refinement submission.

## Decision Log

- Use a new run-scoped preflight root; do not modify the existing BM55 run,
  canonical PRISM submission tree, or source checkout.
- Do not submit production arrays and do not run DockQ in this phase.
- Treat command resolution, `--help`/version checks, and one tiny smoke run as
  separate evidence from Slurm submission or scientific validation.

## Outcomes

The bounded preflight passed. Production refinement remains unsubmitted until
the missing wrapper is staged together with the exact run manifest and valid
transformed PDB inputs.

## Context

- Source checkout: `/scratch/rshadi25/GitHub/PRISM-prescript`
- KUACC host: `rshadi25@login.kuacc.ku.edu.tr`
- Rosetta module: `rosetta/2022.42`
- Local MultiProt executable: `external_tools/multiprot.Linux`
- Local FiberDock bundle: `external_tools/fiberdock/`
- PRISM MultiProt adapter: `src/alignment_multiprot.py`
- PRISM FiberDock adapter: `src/fiberdock_refinement.py`
- PRISM Rosetta adapter: `src/rosetta_refinement.py`

## Plan

1. Inspect source and remote dependency trees and verify wrapper assumptions.
2. Build an explicit manifest with paths, hashes, permissions, and commands.
3. Stage only declared files under a new KUACC preflight namespace.
4. Run syntax/runtime checks and a bounded smoke case where inputs exist.
5. Retrieve evidence, validate artifacts, and document blockers before any
   production authorization is considered.

## Concrete Steps

- Validate Python syntax for the relevant adapters and staging utilities.
- Validate that the external-Rosetta batch wrapper is present in the staged
  bundle; absence in the existing submission tree is a submission blocker.
- Validate the Rosetta module and database path on KUACC.
- Validate MultiProt binary metadata and a minimal invocation.
- Validate all FiberDock companion scripts, libraries, and working-directory
  requirements; run a single-input probe only if suitable existing PDB inputs
  are available.
- Record return codes, elapsed time, stdout/stderr, hashes, and output paths.

## Validation

- Local: Python compile, shell syntax, executable/file checks, and manifest
  consistency.
- KUACC: module/tool discovery, `sbatch --test-only` only if a smoke job is
  needed, then one isolated smoke execution; use `squeue` and logs rather than
  `sacct`.
- Scientific validation is out of scope until the smoke contract passes and
  a separate production run is authorized.

## Results

- Login-node dependency preflight: PASS at 2026-09-27 12:56 TRT.
- Compute-node smoke job `3124094`: PASS on `be01.kuacc.ku.edu.tr`.
- Per-step smoke durations: MultiProt 0.0438 s; Rosetta prepack help
  0.2433 s; Rosetta docking help 0.1765 s; FiberDock 0.0312 s; NMA
  0.0126 s; Reduce 0.0145 s.
- Local adapter/wrapper Python compilation, shell syntax, focused contract
  tests, and `git diff --check`: PASS.
- No actual PDB alignment/refinement run: NOT TESTED because no valid PDB pair
  was present in the available KUACC/source roots.
- No production job or DockQ job submitted.

## Idempotence

The preflight destination must be new or explicitly empty only after a
read-only check. Existing production roots are never overwritten. Every
artifact is content-addressed in the manifest and all reruns use a new
timestamped preflight namespace.

## Artifacts

- Local plan: this file.
- Expected remote root:
  `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-kuacc-preflight-20260927`
- Expected evidence: manifest, module/tool reports, smoke logs, and output
  hashes under the remote root and a small retrieved evidence copy locally.

Retrieved evidence:
`evidence/prism-prescript/kuacc-preflight-20260927/`

## Interfaces

- Rosetta wrapper must pass an explicit database path and resolve both
  `docking_prepack_protocol.static.linuxgccrelease` and
  `docking_protocol.static.linuxgccrelease`.
- FiberDock must run from its own directory so its bundled `lib/` is found;
  all companion files are required.
- MultiProt must be invoked through the PRISM adapter with an explicit binary
  path and run-scoped working directory; no implicit current-directory or
  source-tree mutation is allowed.
