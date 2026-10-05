# Reconstructed pre-consolidation project memory archive

**Status:** historical reference only; not the active project memory.

This document reconstructs durable context that was previously distributed
across the long project-memory notes. It was assembled on 2026-07-25 from
surviving repository documentation, execution plans, source code, and retained
artifacts. It cannot restore the former memory files verbatim. Where sources
conflict, the active memory under `.agents/skills/project-memory/references/`
and the newest dated evidence take precedence.

No results, benchmark outputs, raw inputs, or stable pipeline files were
modified to create this archive.

## Historical scope and evidence rules

- The repository BM5/5.5 tables are not an authoritative reproduction of the
  paper's 88-case BM3 cohort. `dataset_row_id` is the stable benchmark identity;
  normalized PDB pairs are audit fields only.
- Curated `r_u/l_u` archive members are source inputs and `r_b/l_b` members are
  native truth. Full PDB caches are audit comparators unless a separately
  validated chain-materialization step is recorded.
- Scheduler completion, a shell exit code, file counts, and scientific success
  are distinct. A run is scientifically complete only when its normal return,
  terminal stage records, refined-pose integrity, and evaluator results agree.
- Transformation halves are intermediate files, never scoreable docked models.

Primary surviving sources:

- `docs/pipeline-validation-report-20260714.md`
- `docs/validation/prism-pipeline-verification-20260718-blocked.md`
- `docs/exec-plans/20260718-prism-pipeline-verification-gates.md`
- `docs/ML_TRAINING.md`
- `docs/biological-ranking-pilot-20260711.md`
- `docs/legacy-multiprot-runtime.md`
- `docs/STABLE_PIPELINE.md`

## Historical pipeline record

### Current Python 3 pipeline

- Entry point: `prism.py`.
- Stages: input/download, surface extraction, alignment, transformation and
  filters, optional candidate ranking, selected refinement, optional scoring.
- Backends: NACCESS or FreeSASA; TM-align, GTalign, or MultiProt; external
  Rosetta, PyRosetta, or FiberDock.
- The current stable recipe is `docs/STABLE_PIPELINE.md`. Use the standard
  `gtalign_env` interpreter and isolated workspaces.
- The exact archived notes support nine confirmed combinations in the
  2026-07-21 controlled validation. Later 18-variant BM5.5 notes describe
  planned/submitted work with pending GPU and cancelled CPU jobs; they do not
  prove an all-18 completed benchmark. The normal operational default remains
  NACCESS + TM-align + external Rosetta.

### Legacy/reference pipeline

- `working_version/Multiprot-new/prism-fiberdock-cli/` is a Python 2
  compatibility/reference arm, not a target for current pipeline edits.
- Its direct historical execution is constrained by legacy 32-bit helpers,
  especially reduce.2 and compatible loader/library availability.
- An explicit reduce.3 substitution completed one controller-path diagnostic,
  but is not historical-equivalent and must not enter primary comparisons.

## Recovered experiments and conclusions

### Current versus historical comparisons

- Existing July current/legacy summaries are observational, not causal:
  source files, template panels, filters, output integrity, and scoring were
  not all matched.
- A current TM-align + Rosetta positive smoke exists for
  `1RGH_B + 1A19_B` using template `1b27AD`; the full observational run
  produced 143 Rosetta models. This is not a matched performance estimate.
- Root and historical working-version TM-align executables were recorded as
  byte-identical in `docs/pipeline-validation-report-20260714.md`; differences
  therefore arise downstream of the binary.
- No current-vs-legacy quality claim is valid until curated sources, template
  exposure, top-K budget, native-independent selection, and evaluator are
  frozen for both arms.

### FiberDock diagnostic and integrity finding

- A one-pair reduce.3 compatibility diagnostic generated final FiberDock
  artifacts after fixing a launcher parent-directory defect.
- The generated PDB collapsed both partners into chain `B` and reset residue
  numbering. It is invalid for the canonical two-partner scorer.
- A derived split-on-reset repair made that file scoreable (DockQ about 0.745),
  but it is diagnostic only: it does not prove intended chain identity,
  historical model recovery, or reduce.2 equivalence.
- The live issue is chain-preserving FiberDock output plus a permitted native
  reduce.2 runtime; do not use the repair or reduce.3 to open the primary arm.

### GTalign findings

- The previous GTalign pre-score failure was double filtering and was fixed by
  separating `PRISM_GTALIGN_PRE_SCORE` from the transformation TM threshold.
- Sparse broad-panel hits at TM-score ≥0.4 are expected for weakly similar
  targets/templates; they are not evidence of a parser failure.
- CPU and GPU GTalign output counts still differ materially under apparently
  matching settings. CPU results are not interchangeable with GPU results
  until an exact binary-version and argument audit is completed.
- GPU is the practical backend for 20K-template panels; CPU MultiProt/TM-align
  broad runs were too slow for routine full-benchmark execution.

### Filter, source-gate, and evaluator work

- Modern JSON and derived legacy protocol assets use incompatible encodings.
  Chain/orientation-aware normalization was implemented, but a large parity
  disagreement remains (recorded as 19,005 disagreements and 850 missing
  modern profiles in the 2026-07-18 verification evidence).
- `published_protocol` is intentionally fail-closed when required assets are
  absent; `geometry_only_experimental` is usable for wiring experiments but
  cannot substantiate a faithful historical-protocol claim.
- The source matrix was deterministic across AI/COSBI, but source authority
  remains blocked. The preserved policy retains 240 strict-eligible audit rows
  while 17 chain-contract rows require an authoritative decision.
- Canonical scoring records complete chain mappings, raw DockQ JSON hashes,
  and grouped iRMSD. Some difficult structures produced explicit DockQ runtime
  failures and remain excluded from quality denominators rather than dropped.

## Historical ranking and ML record

- Candidate auditing, deterministic baseline ranking, candidate-table building,
  DockQ/iRMSD label attachment, grouped splitting, tabular training, and a
  contact-model scaffold were implemented as optional components.
- The deterministic pilot for `1FGNH/1TFHA` against `1ahw` ranked `1h5bAB_o1`
  first (DockQ 0.031, iRMSD 16.937), but all four candidates were below the
  native-like DockQ threshold of 0.23. This is a one-complex observation.
- A five-complex grouped evaluation did not improve native-like top-1 success
  with the learned tabular reranker. Learned ranking stays disabled.
- 2026-07-25 corrected the runtime ranking integration: ranked runs create a
  new audit by default, the audit destination is resolved after CLI parsing,
  stale/partial audits fail open, and `--top-k` is validated. Focused tests:
  21 passed.

## Historical environment and execution notes

- Pipeline Python: `/home/rshadi25/.conda/envs/gtalign_env/bin/python`.
- DockQ Python: `/scratch/tmp/prism-dockq-env/bin/python`, DockQ 2.1.3 via
  `python -m DockQ`.
- External Rosetta: `module load rosetta/2022.42` and explicit
  `PRISM_ROSETTA_PREPACK`, `PRISM_ROSETTA_DOCK`, `PRISM_ROSETTA_DB` values.
- Run heavy work through Slurm; login nodes are for setup, inspection,
  submission, and collection. Use isolated task/run roots and retain command,
  environment, hashes, and job ID.
- Rescue settings (relaxed matches/TM/coverage/clash cutoffs) are diagnostic
  evidence only and must remain out of stable-production interpretation.

## Items intentionally not carried into active memory

- Superseded job-state snapshots, pending-job references, and transient Slurm
  controller errors.
- Repeated descriptions of the same GTalign pre-score, Rosetta-module, and
  Python-2 helper failures.
- Unpaired aggregate performance values presented without their provenance
  caveat.
- Proposed implementations that were superseded by later source-gate,
  evaluator, filter-normalization, or ranking decisions.

## Remaining recoverability limits

- The original long memory files were untracked local Markdown and were
  overwritten during consolidation; no verbatim pre-consolidation copy was
  available in Git.
- This archive preserves recoverable conclusions and evidence paths, not every
  historical shell command, log excerpt, or abandoned hypothesis.
- For a particular historical run, consult its retained `tmp/agent/` task
  directory and the execution plans before relying on this summary.
