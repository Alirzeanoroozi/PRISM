# Bounded TMalign/MultiProt refinement comparison

## Purpose

Reuse retained BM55 transformed artifacts to compare high-confidence TMalign
and MultiProt candidates through both FiberDock and external Rosetta, then
score each refined model with the validated DockQ/iRMSD contract. Keep the
comparison bounded to 20 candidates per aligner (40 total) and preserve
candidate identity, timing, hashes, failures, and resumable checkpoints.

## Scope and assumptions

- Source run:
  `/scratch/users/rshadi25/valar-remote-runs/prism-prescript-bm55-full-20260919`
- Select the highest-confidence rows by the lower of the two side TM-scores,
  with mean TM-score and match count as deterministic tie breakers.
- Use existing transformed PDBs and native PDBs; do not recompute alignment.
- Run both refiners for every selected candidate where inputs validate.
- DockQ is run only on valid refined model outputs, with explicit native chain
  mapping and per-model raw JSON/hash provenance.

## Progress

- [x] Confirm retained TMalign/MultiProt transformed artifacts exist.
- [ ] Select and hash the bounded candidate manifest.
- [ ] Run one-candidate paired smoke for both aligners/refiners.
- [ ] Run remaining bounded refinement arrays.
- [ ] Run/validate DockQ for every refined output.
- [ ] Aggregate timing, energy, DockQ, failure, and no-drop tables.

## Safety and reproducibility

- Use a new run-scoped output root; never overwrite the parent alignment run.
- Use one candidate per resumable task and separate work directories for
  Rosetta and FiberDock because both legacy tools use shared output names.
- Keep transformed inputs, refined outputs, DockQ JSON, checkpoints, logs, and
  stage timing in the comparison root.
- No ranking, PRODIGY, PyRosetta, raw-data deletion, canonical promotion, or
  push is in scope.

## Validation gates

1. Candidate manifest contains exactly the selected count and both transformed
   files plus native PDB for every row.
2. One TMalign and one MultiProt candidate each pass through both refiners or
   fail with a recorded stage-specific reason before array expansion.
3. Refinement status and energy output are explicit for every candidate.
4. DockQ/iRMSD rows retain raw JSON hashes and explicit `scored`, `valid_unscored`,
   `score_failed`, or `not_run_no_model` states.
5. Aggregation proves no selected candidate was silently dropped.
