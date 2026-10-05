# Step 4 handoff — optional PRODIGY

## Objective

Assess whether the locally installed PRODIGY can be an optional scoring stage without changing baseline behavior or silently dropping failed cases.

## Result

The executable and source are available locally and the existing ranking path is already opt-in. An isolated candidate adds an explicit state contract:

`available`, `executed`, `failed`, `skipped`, `not_configured`.

The candidate preserves all passed pairs when the executable is missing and preserves an affected candidate group when scoring fails. It keeps the existing PRODIGY thresholds and CLI argument order.

## Isolated candidate

Worktree: `/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-4f0e2099ca254fde8dbd32b71a96e740`

Files added there for review only:

- `src/prodigy_ranker.py`
- `src/candidate_selector.py` (minimal test shim)
- `tests/test_prodigy_ranker_existing.py`
- `tests/test_prodigy_state_contract.py`

The candidate is not merged or copied into the canonical checkout.

## Verification evidence

- Isolated focused tests: 9 passed.
- Real executable success: Slurm job `1656993`, affinity `-65.827`.
- Real failure preservation: Slurm job `1656992`, no-contact error, selected count remained 2.
- Canonical focused tests: 8 passed; canonical status count/hash remained unchanged.
- Graphify output and query limitations are recorded in `evidence/current-pipeline-map.md`.

## Next safe action

If implementation is authorized, review the isolated diff against the full canonical `src/candidate_selector.py`, add the state contract to the real selector integration, run the full focused suite from a writable run directory, then independently validate the complete PRISM path and benchmark quality before promotion. Do not claim full scientific validation from the adapter smoke alone.
