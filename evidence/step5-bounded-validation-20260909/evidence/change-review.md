# Step 5 evidence change review

## `current-pipeline-map.md`

- Purpose: record the maintained call path and separate direct source evidence
  from generic Graphify relationships.
- Critical regions: provider dispatch, disconnected `StructuralAligner`, score
  semantics, attrition locations, and refiner/evaluator boundaries.
- Correctness: direct imports in `prism.py` are treated as authoritative;
  Graphify is reported only as a temporary structural aid.
- Risk: the map is source/static evidence and cannot replace a fresh Slurm
  integration run.
- Verification: compared against `prism.py`, provider modules,
  `transformation.py`, `rosetta_refinement.py`, and `compare.py`.

## `tool-matrix.json` / `tool-matrix.tsv`

- Purpose: provide a reproducible comparison schema for bounded and retained
  combinations.
- Critical regions: exact command, executable version/provenance, scheduler
  identity, stage counts, failure reason, and output paths.
- Correctness: unknown scheduler fields and unrecorded hashes are `null`; rows
  that are not a frozen matched panel are explicitly `not_comparable`.
- Risk: historical rows cannot establish causal speed, ranking quality, or
  biological superiority.
- Verification: JSON parsing, JSONL parsing, 37-column TSV schema check, and
  row-by-row TSV/JSON consistency check passed.

## Stage ledgers

- `alignment-stage-ledger.json`: preserves 76,248 -> 2,997 -> 560 -> 71
  retained counts and records the 73,251 missing/unwritten remainder without
  inventing a provider failure reason.
- `no-drop-ledger.json`: distinguishes observed drops from unknown stages;
  clash rejection is 489, and the 1,877 residual is not relabeled as
  structural orphanage.
- `refiner-stage-ledger.json`: records retained Rosetta/FiberDock counts while
  exposing missing per-candidate return-code and score-gate provenance.
- `matched-candidate-ledger.json`: limits the retained 4-to-1 result to load
  reduction and records PRODIGY failure-preservation behavior.
- Risk: these are evidence ledgers, not implementation changes or validation
  of unrun arms.

## `error-inventory.json`, `final-validation.json`, and closeout files

- Purpose: append source-level risks, tool preflight, explicit status
  separation, and the current scheduler blocker to the durable run record.
- Critical regions: MultiProt Seccomp-2 restricted probe, disconnected
  StructuralAligner, nonuniform score contracts, Rosetta observability, DockQ
  provenance, and job `1657005` reconciliation.
- Correctness: no status is upgraded from `UNKNOWN`; isolated candidate code is
  not promoted; the canonical project fingerprint remains unchanged.
- Verification: full framework `pytest`, `compileall`, artifact parsing,
  `git diff --check` on the run, and canonical status fingerprint.
