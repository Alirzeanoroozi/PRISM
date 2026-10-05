# PRISM-prescript Step 2 report — transformation attrition

Recorded: 2026-09-08 22:22 +03:00

## Outcome

The retained 1gte attrition funnel is reconstructed and validated within its stated historical scope. The dominant measured loss is before transformation: 73,251 of 76,248 potential sides have no alignment JSON, leaving 2,997. Of those loaded sides, all 2,997 pass the run-era transform-field and match/TM/coverage checks. Pairing leaves 560 complete orientation attempts; all 560 materialize; CA-clash filtering rejects 489 and leaves 71 final passes.

The earlier retained numbers are partly confirmed and partly mislabeled:

| Figure | Result | Classification |
|---|---:|---|
| alignment JSON files | 2,997 | `EVIDENCE/CONFIRMED` |
| orphan sides | 1,163 structural missing-partner sides; 1,877 residual unconsumed sides | `EVIDENCE/DEFINITION_MISMATCH` |
| transformable sides | 574 | `UNKNOWN/UNDEFINED_DENOMINATOR` |
| transformed pairs | 560 | `EVIDENCE/CONFIRMED` |
| clash rejections | 489 | `EVIDENCE/CONFIRMED` |
| final passes | 71 | `EVIDENCE/CONFIRMED` |

## Current-source invariant checks

In a temporary isolated execution context, the current `src/transformation.py` passed:

- identity, translation, and 90-degree rotation transforms;
- both `o1` and reverse `o2` orientation dispatch when all four side files exist;
- suppression of `o2` when its reverse partner is missing;
- rejection of malformed transform dimensions and missing alignment files;
- deterministic replacement of an existing output path.

The probe returned `all_pass: true`; the canonical Git status entry count and SHA-256 status hash were unchanged.

## Root-cause decision

No production fix is justified by this evidence. The data support an upstream alignment-file availability problem and an ambiguous historical count label, not a transformation threshold defect. The stable thresholds were not changed. The retained diagnostic tools and tests remain in the isolated worker worktree and were not merged or promoted.

## Tool/runtime status

Graphify 0.9.29 was used from `/home/rshadi25/.local/bin/graphify` to generate a 44-file code-only graph in `/tmp`; the explicit call query mapped the alignment-to-transformation flow. `USalign` is available in `gtalign_env` at `/home/rshadi25/.conda/envs/gtalign_env/bin/USalign` (Version 20241108). The project-local USalign path is absent. The legacy MultiProt binary is present at `external_tools/multiprot.Linux`, but its documented 32-bit/seccomp runtime blocker remains a separate execution-context issue.

## Qwen/Copilot status

The local provider smoke was recorded as passing, but the attempt-5 Copilot worker ended at its walltime checkpoint with exit `138:0`. Its retained artifacts are evidence inputs, not model-validation evidence. No production code was modified or merged.

## Acceptance gate

- Every requested attrition stage has a count and scope/status: **PASS for retained case**.
- Errors/blockers have reproducible evidence: **PASS for retained evidence; live scheduler status remains UNKNOWN**.
- Current alignment-to-transformation flow is mapped, including Graphify output: **PASS**.
- No project files were modified by this Step 2 work: **PASS; canonical status hash unchanged**.

This closes Step 2 at the retained-case validation boundary. A current full production rerun remains a separate, explicitly bounded step.
