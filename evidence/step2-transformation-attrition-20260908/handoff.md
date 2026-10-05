# Step 2 handoff

## Objective

Diagnose transformation-file attrition for the retained 1gte case and test the current transformation invariants without changing production code.

## Evidence inspected

- Project instructions and memory under `.agents/skills/project-memory/references/`.
- Current `prism.py`, `src/transformation.py`, `src/structural_aligner.py`, `src/alignment_gtalign.py`, `src/alignment_multiprot.py`, and focused transformation tests.
- Isolated worker diagnostics under `tmp/agent/worktrees/run-f9039c5bc245456e8c460b93e6f621af/`.
- Retained ledger and clash validation under `tmp/agent/20260908-prism-prescript-smoke-runtime-pin/project-workers/run-f9039c5bc245456e8c460b93e6f621af/`.
- Graphify code-only graph and explicit call-context query.

## Findings

- `OBSERVATION/EVIDENCE`: 2,997 alignment JSON sides are present from 76,248 potential sides; 73,251 are missing or unwritten at the alignment-prefilter boundary.
- `EVIDENCE`: 1,163 loaded sides lack their required partner; 1,834 are paired. The historical 1,877 figure is the residual `2,997 - 2*560`, not the same structural orphan definition.
- `EVIDENCE`: 2,997/2,997 transform fields and match/TM/coverage side gates pass in the retained run-era reconstruction.
- `EVIDENCE`: 560 orientation attempts materialize as 560 disk pairs; CA-clash filtering rejects 489 and retains 71.
- `UNKNOWN`: “574 transformable sides” cannot be reproduced without a frozen denominator/definition.
- `INFERENCE`: the retained attrition is dominated by upstream alignment-file availability and by an ambiguous count label, not by transformation thresholds, materialization, or the CA-clash filter.

## Decisions and rationale

- `REVIEWED`: do not change stable thresholds or production defaults.
- `REVIEWED`: do not merge the isolated diagnostic tools/tests; they are evidence-generating additions, not a demonstrated production fix.
- `REVIEWED`: use the `gtalign_env` USalign path explicitly; keep legacy MultiProt runtime conclusions separate from current pipeline behavior.

## Status distinctions

- `IMPLEMENTED`: isolated diagnostic ledger/clash tools and focused tests exist in the worker worktree.
- `TESTED`: worker diagnostic tests (12 in retained evidence), current focused transformation tests (15), and the bounded invariant probe passed.
- `VALIDATED`: retained 1gte ledger and CA-clash counts are validated within their historical geometry-only scope by Slurm evidence.
- `REVIEWED`: Graphify mapping, current-vs-run-era source difference, tool paths, and canonical Git preservation were reviewed.

## Risks and uncertainty

- Historical run-era source differs from the dirty canonical current source; no current full production rerun was performed.
- Geometry-only mode bypasses protocol hotspot/contact filtering.
- Live scheduler state could not be queried from this session.
- Copilot/Qwen attempt 5 reached local provider smoke but terminated at the 30-minute checkpoint; model output is not treated as validation.

## Recommended next action

Restore Slurm query access, freeze a current-source input manifest, and rerun only the ledger in an isolated worktree with an explicit definition for transformable sides. Do not alter thresholds until a current-source failure category is quantified.
