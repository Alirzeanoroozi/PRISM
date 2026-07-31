# Changelog

## 2026-07-30

### Features

- Added an opt-in PRODIGY candidate-ranking adapter. It combines transformed
  receptor/ligand PDBs with collision-free chain IDs, invokes an externally
  installed `prodigy` executable, and retains score/input/command evidence.

- Added the canonical Phase 1 provenance core at `src/provenance`, covering
  contract/attempt identity, row-aware artifact observations, closeout, and a
  fail-closed consumer gate.
- Added end-to-end fixture coverage for mutation, missing artifacts, symlinks,
  duplicate keys, row mismatches, retries, secret redaction, and the CLI gate.

### Fixes

- Fixed PRODIGY command construction so the positional input PDB is passed
  before `--selection`; PRODIGY's argparse treats `--selection` as `nargs='+'`,
  so putting the input path after it made the path look like another chain
  group and caused return code 2.
- Existing provenance helpers now reuse canonical JSON and file hashing while
  preserving direct path-based CLI execution through a narrow repository-root
  import fallback.

### Learnings

- PRODIGY requires its own FreeSASA/NumPy environment; keeping it external
  preserves the verified DockQ environment and makes the scorer version
  explicit in ranking provenance.
- On the retained `5zngA,4eylA` / `1a0cCD` two-orientation case, opt-in
  PRODIGY top-1 ranking changes the forwarded candidate set from two
  candidates to one (`o1`, -65.827 kcal/mol versus `o2`, -65.274 kcal/mol).
  This is a selection/load-change observation, not a DockQ quality claim.

- Secret-bearing declaration values must be redacted before contract hashing;
  they are not allowed to become hidden identity inputs.
- A separate closeout view preserves the append-only ledger while making current
  artifact bytes and expected row identity explicit before consumption.
