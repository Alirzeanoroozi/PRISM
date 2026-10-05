# Isolated candidate changes

## USalign

Worktree: `/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-1e6873e93a4f4082a236f7348218ecd5`; HEAD `1026a7fc0609e19f19c48a68510e7ca39d67573d`;
feature-gated CLI branch and `src/alignment_usalign.py`; focused suite 20/20;
real probes jobs 1656964/1656966. The adapter preserves empty/failure records
and requires `-outfmt -1 -m -`. It is not on canonical `prism.py` and was not
promoted.

## PRODIGY

Worktree: `/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-4f0e2099ca254fde8dbd32b71a96e740`; HEAD `1026a7fc0609e19f19c48a68510e7ca39d67573d`;
feature-gated `src/prodigy_ranker.py`/candidate selector; focused suite 9/9;
real probes jobs 1656992/1656993. The candidate preserves no-contact groups
and uses explicit ranking states. It is not a quality-validated production
default and was not promoted.

## Shared stage contract

Worktree: `/scratch/tmp/prism-prescript-contract-20260913`; pinned baseline `1026a7fc0609e19f19c48a68510e7ca39d67573d`;
isolated Luna High candidate files `src/pipeline_contract.py`,
`src/alignment_result_adapter.py`, `benchmark/scripts/
build_alignment_event_ledger.py`, the additive `prism.py` stage-event
integration, and their tests. The focused contract/adapter/import/ledger suite
passed 34 tests and the complete isolated suite passed 37 tests. The candidate
also guards absent optional backends so default TMalign CLI import remains
usable while requested unavailable backends fail explicitly. Its adapter
preserves unknown, rejection, refinement-failure, and score-gate statuses as
explicit non-success records rather than treating a score as implicit success;
it also retains optional run/attempt identifiers and provider-specific numeric
metrics such as MultiProt RMSD. The ledger preserves one explicit event per
resolved manifest row, including missing/not-run records, and rejects missing
identity or duplicate candidate IDs before output. The candidate is not
integrated into the canonical checkout and does not claim provider,
transformation, ranking, refinement, or evaluator validation.

## Review risks before promotion

The canonical path still has nonuniform provider score semantics, incomplete
per-candidate refinement return codes, possible transformation audit gaps, and
an absent `.github/copilot-instructions.md` at the requested path. The isolated
contract remains partial: it lacks complete event/parent correlation fields,
does not yet wrap every pipeline stage, and the adapter/ledger are not wired
into provider writers. The ledger is tested against synthetic resolved
manifests only and does not repair the canonical raw-output hash semantics or
foreign-key lineage by itself. The candidate diffs must be rebased/reviewed against the dirty
canonical state in a new bounded worktree; this package deliberately makes no
source merge.
