# Changelog

## 2026-09-10

### Fixes

- Asset provenance now classifies consumed current/legacy interface and filter
  assets by resolved path and content hash, including copied or renamed files;
  interface-list coverage assets are included, and unresolved matches remain
  explicit instead of being treated as current.
- Ranking diagnostics now record the independent native-label source and
  SHA-256 and mark recovery provenance incomplete when either is absent.
- Alignment-dependent gate replay now fails closed to `unknown` when alignment
  return codes/statuses or raw-output hashes are incomplete, and ranking
  recovery is deferred when pre-ranking panels do not match exactly.

## 2026-09-09

### Features

- Added `--orientation {native,o1,o2}` to the current pipeline. `native`
  preserves the MultiProt-like dual-assignment behavior, while `o1` and `o2`
  provide fixed-orientation comparison runs.
- Added explicit template-panel selection (`--templates` or `--template-list`)
  and per-run transformation threshold controls. The selected input CSV now
  reaches both download and transformation stages.
- Added the `--scaffold-threshold` override with the established 5.0 Å default;
  the value is now carried in the typed threshold configuration and applied at
  surface extraction.
- Candidate-audit records now retain the resolved transformation thresholds,
  making native/o1/o2 yield differences traceable to their run configuration.
- Added `notebooks/pipeline_orientation_comparison.ipynb` with separate,
  default-disabled current-CLI cells for the native-default, `o1`, and `o2`
  arms.
- Added `src/stepwise_analysis.py` and extended the notebook with frozen-baseline,
  asset-parity, alignment-provenance, gate-ablation, clash-grid, refinement,
  ranking, and separate US-align preflight diagnostics.

### Fixes

- GTalign prefiltering now requires both reported TM-scores to meet the
  configured pre-score, avoiding acceptance of a one-sided low-score hit.
- The no-option orientation arm is now deterministically `native`, independent
  of inherited environment state. MultiProt records use their native
  match/coverage gate without first passing through generic TMalign gates.
- Candidate-audit statuses now distinguish missing alignments, protocol
  rejection, alignment-threshold rejection, transformation failure, and clash
  rejection.
- MultiProt interface coverage now uses the chain-specific template interface
  size as its denominator instead of the number of mapped residues.
- TMalign and MultiProt alignment JSONs now retain a raw-output hash and
  subprocess return code when available, so provenance gaps are visible before
  downstream interpretation.

### Learnings

- “Without orientation” is implemented as evaluating both implicit
  template-chain assignments and preserving their provenance, not as ignoring
  chain assignment. Contacts and transforms remain chain-specific.
- Explicit threshold configuration is applied at the transformation gates and
  remains separate from the MultiProt native match/coverage contract; relaxed
  diagnostic settings must not be treated as validated production defaults.

## 2026-08-23

### Features

- Added an executable transformation-gate ablation notebook that keeps every
  candidate visible across independent gates, cumulative diagnostic attrition,
  production-order replay, and leave-one-gate-out analysis.

### Learnings

- In the frozen `1FGNH`/`1TFHA` against `1kcaCH` demonstration, both
  orientations fail the TMalign score gate at the stable 0.5 threshold; one
  orientation would also fail the downstream C-alpha clash gate. Independent
  gate measurement is therefore necessary to distinguish early attrition from
  downstream geometry behavior.

- Added a selectable notebook case for `1gteA`/`1gteB` against `1h7xCD`.
  The case is prepared but cannot be interpreted until its four alignment JSON
  records are generated successfully.

### Features

- The selectable `1gteA`/`1gteB` against `1h7xCD` notebook case now has an
  opt-in runner that stages canonical assets, downloads the query PDB when
  needed, executes the current TMalign pipeline with no refinement, and feeds
  the generated alignment directory into the gate ledger.

## 2026-08-11

### Features

- Added an opt-in `legacy_compatible` MultiProt mode to the current adapter.
  It preserves legacy interface/query invocation order, optional `params.txt`,
  `Reference Molecule`, `Trans`, and up to three solver solutions.
- Promoted retained legacy solutions into the current transformation/filtering
  path, with collision-safe per-solution transformation filenames.

### Fixes

- Current MultiProt no longer discards valid legacy solutions because Kabsch
  cannot reconstruct a transform; the legacy `Trans`/reference-molecule
  transform is used only when the compatibility mode is explicitly selected.

### Learnings

- The 32-bit MultiProt binary is executable in the validated `kutem` Slurm
  runtime (`Seccomp: 0`) but is blocked by the direct shell runtime
  (`Seccomp: 2`); Python/Conda environment changes do not remove that kernel
  restriction.

## 2026-07-31

### Features

- Added `/PRISM`-compatible CLI aliases and runtime controls for input CSV,
  template limits, surface and GTalign spellings, MultiProt, PyRosetta,
  FiberDock, and explicit refinement gating.
- Both bare flags and legacy explicit boolean values remain accepted where the
  two entry points previously used different argparse conventions.

### Fixes

- Input CSV selection is now passed at call time instead of being limited to
  the downloader's import-time environment value.
- MultiProt executable and worker settings, PyRosetta output/init settings,
  and FiberDock root selection are now wired to their actual backend calls.

### Learnings

- CLI spelling parity is separable from default-policy parity. Prescript keeps
  its validated refinement-on default; `--no-refine` makes skipping explicit
  without silently changing existing commands.

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
