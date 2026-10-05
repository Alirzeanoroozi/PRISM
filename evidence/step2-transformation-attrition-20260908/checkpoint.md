# Step 2 checkpoint — transformation attrition

Recorded: 2026-09-08 22:22 +03:00

State: `VALIDATED_FOR_RETAINED_1GTE_LEDGER`; `NOT_VALIDATED_FOR_CURRENT_FULL_PRODUCTION_RERUN`.

Completed:

- Recovered the retained 1gte run and its run-era source provenance.
- Reconstructed every requested attrition stage with counts; no stage was left unaccounted for in the retained case.
- Independently confirmed the retained CA-clash partition in Slurm evidence: 560 evaluated, 489 rejected, 71 passing, zero missing transformed pairs.
- Ran focused transformation tests in an isolated `/tmp` context: 15 passed.
- Ran a bounded current-source invariant probe: identity, translation, rotation, both orientation branches, missing reverse side, malformed transform, missing alignment file, and existing-output replacement all passed. The probe left the canonical Git status hash unchanged.
- Used Graphify code-only extraction and recorded the call map and executable paths.

Open/uncertain:

- `574 transformable sides` has no unique denominator in the retained evidence; it remains `UNKNOWN` rather than being promoted to a count.
- The ledger is a historical geometry-only diagnostic reconstruction, not a current full PRISM rerun and not a published-protocol quality result.
- Live Slurm status is unavailable from this login session because both `squeue` and `sacct` failed to connect; durable manifests record the Copilot attempts as failed at the walltime checkpoint (`138:0`) and the validation evidence is preserved.

Decision: no production code, default, threshold, or canonical output was changed. A minimal fix is not evidence-supported by this attrition case.

Next safe action: if a current-run causal fix is required, run a fresh bounded current-source ledger in a new isolated worktree after scheduler access is restored; first define the denominator for “transformable sides.”
