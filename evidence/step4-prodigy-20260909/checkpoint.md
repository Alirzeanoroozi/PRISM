# Step 4 checkpoint — optional PRODIGY

Timestamp: 2026-09-09 (Europe/Istanbul)

## Gate status

`TESTED` and `VALIDATED` within a bounded scope: the optional PRODIGY behavior and explicit failure states pass the isolated contract tests and a real local executable smoke on Slurm. The default-disabled path is unchanged. This is not a claim that the full PRISM scientific ranking quality has been validated.

## Evidence completed

- Local executable resolved at `/home/rshadi25/.conda/envs/gtalign_env/bin/prodigy`; local source reports version `2.4.0`.
- Direct Slurm probe `1656986` produced affinity `-65.827` from a retained paired fixture.
- Failure-preservation probe `1656992` reproduced `No contacts found for selection`; state sequence was `available`, `failed`, and both candidates were preserved.
- Corrected state smoke `1656993` produced `available`, `executed`, score state `executed`, one selected candidate, and affinity `-65.827`.
- TDD red-green evidence: unmodified isolated adapter failed all five state-contract tests; the minimal isolated state patch made all five pass.
- Isolated focused tests: `9 passed in 0.28s`.
- Canonical focused tests: `8 passed in 1.65s`; no canonical files were modified.
- Graphify deterministic extraction: 51 code files, 561 nodes, 1467 edges; graph and query limitation are recorded in `evidence/current-pipeline-map.md`.

## Boundaries and unknowns

- The explicit-state implementation is only in the isolated worktree; it has not been promoted to the dirty canonical checkout.
- The full current-tree PRISM run with `--rank true --rank-method prodigy`, native DockQ quality comparison, and promotion review were not run in this bounded step.
- Two broader candidate-selector tests are blocked by their canonical read-only working-directory fixture writes (`OSError 30`); this is recorded as harness evidence.
- PRODIGY does not establish biological validity by itself; affinity ordering remains an optional heuristic and needs downstream benchmark validation.

## Stop decision

Stop at the Step 4 boundary. Preserve the isolated candidate and all run evidence. Do not merge, reset, delete, or alter production defaults without a later authorized implementation/review step.
