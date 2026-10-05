# Step 3 checkpoint — USalign

Timestamp: 2026-09-08 (Europe/Istanbul)

## State

- Step: `STEP 3 — USALIGN`
- Gate: `PASS_FOR_BOUNDED_COMPATIBILITY; CANONICAL_PROMOTION_DEFERRED`
- USalign classification: `AVAILABLE`
- Canonical project files modified by this step: `NO`
- Isolated candidate status: `IMPLEMENTED`, `TESTED`, `VALIDATED_BOUNDED`, `REVIEWED_ISOLATED_ONLY`

## Completed evidence

1. Recovered project instructions, memory, prior worker evidence, current
   canonical Git state, and the isolated candidate worktree.
2. Located the existing USalign wrapper in `gtalign_env`; verified version
   `20241108` and documented flags locally.
3. Ran bounded USalign probes in Slurm on cosbi. Job `1656964` produced the
   real transform output; job `1656966` ran the corrected adapter smoke.
4. Parsed a real alignment: 214 residues, TM-score `0.98359`, RMSD `0.71`,
   214 mapped residues, chain mapping `L -> A`.
5. Applied the emitted transform to 428 CA coordinates: RMSD
   `0.7500643081394642`; this supports the observed matrix convention.
6. Added focused TDD coverage in the isolated worktree. Full isolated test
   result: `20 passed` with `PYTHONPATH=.`.
7. Used Graphify and direct source inspection to establish that
   `structural_aligner.py` is not connected to the maintained `prism.py`
   entry point.

## Preserved failures and limits

- Job `1656965` is preserved as a wrapper error caused by omitted
  `PYTHONPATH`; corrected job `1656966` succeeded.
- `git diff --check` on the canonical checkout remains nonzero because of
  pre-existing dirty benchmark files. The Step 3 status count and hash are
  unchanged (`373`, `f63a1ccc...b1a8c2e`).
- The full real-project template panel and scientific equivalence against
  existing aligners were not run.
- No candidate files were copied, merged, or promoted into the canonical
  checkout.

## Next action

Stop at the Step 3 boundary. If promotion is later authorized, create a fresh
isolated worktree from the current canonical baseline, apply only the
feature-gated adapter, run a bounded PRISM fixture, and review the resulting
diff before any merge decision.
