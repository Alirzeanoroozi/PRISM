# Step 3 handoff — USalign

## Objective

Determine whether the existing USalign installation can be integrated into
PRISM-prescript without changing default behavior.

## Result

USalign is `AVAILABLE` at
`/home/rshadi25/.conda/envs/gtalign_env/bin/USalign`, version `20241108`.
Compatibility is proven for a bounded real pair and a synthetic parser/failure
fixture. The default selector remains `tmalign` in the isolated candidate and
the canonical checkout was not modified.

## Candidate location

`/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-1e6873e93a4f4082a236f7348218ecd5`

Relevant candidate files:

- `src/alignment_usalign.py`: fail-closed parser and PRISM JSON writer.
- `prism.py`: opt-in `usalign` branch and `--usalign-path`; run-scoped output
  is passed directly to `transformer`.
- `tests/test_usalign_integration.py`: 17 focused tests.
- `tests/conftest.py`: isolated import shims for absent uncommitted modules.

## Verification

- Focused tests: `17 passed`.
- Full isolated tests: `20 passed` using `PYTHONPATH=.`.
- Compile check: passed with `PYTHONPYCACHEPREFIX=/tmp/prism-step3-pyc`.
- Live raw probe: Slurm job `1656964`, output hash
  `fa3c1bec91873b94c1ea89f08ed337ed100c87e83d29be7d0a27a4552d144e48`.
- Live adapter smoke: Slurm job `1656966`, output JSON hash
  `0601cf945843728b98986b5c9d88b271d378a4b11e99530d9005c235760f7a0f`.
- Transform check: 428 CA pairs, RMSD `0.7500643` after applying USalign's
  emitted transform.

## Integration contract

Invoke:

```text
USalign query.pdb interface.pdb -outfmt -1 -m -
```

Parse the real matrix rows as `row_index, translation, rotation[0:3]` and
retain the existing PRISM key direction: interface residue key to query
residue value. Missing tools, missing inputs, nonzero empty output, malformed
output, and missing transforms remain explicit empty records or a clear
executable error.

## Do not infer

This handoff does not prove full-panel scientific equivalence, full PRISM
pipeline success, or authorization to modify/merge the dirty canonical
checkout. Those require a separately bounded validation and explicit
promotion decision.
