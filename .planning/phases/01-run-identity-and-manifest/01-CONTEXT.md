# Phase 1: Run Identity and Manifest - Context

**Gathered:** 2026-07-29
**Mode:** standard
**Status:** Ready for planning

<domain>
## Phase Boundary

Define an immutable, inspectable identity and provenance contract for a PRISM
run before adding more execution evidence. The phase covers the run manifest,
input/configuration/runtime identity, materialized scientific-artifact ledger,
and pre-consumption validation fixture/gate. It does not redesign stage or
candidate lifecycle states, evaluation joins, Slurm orchestration, or the
researcher-facing evidence bundle; those are later phases.

</domain>

<decisions>
## Implementation Decisions

### Run identity and manifest shape
- Give every operational attempt a readable unique `run_id`, and compute a
  canonical manifest hash for the declared run contract. Do not use a
  content-addressed identity alone because repeated attempts must remain
  distinguishable.
- Store run/configuration identity in canonical JSON and artifact records in a
  row-oriented TSV ledger.
- Permit dirty repositories. Record the checked-out Git revision together
  with exact dirty-tree state fingerprints/content references rather than
  requiring a clean worktree.
- Preserve both raw user-declared selectors and normalized selectors consumed
  by execution.

### Artifact ledger
- Cover all materialized scientific artifacts: staged inputs, templates,
  intermediates, final models, natives, and evaluator artifacts. Missing or
  unavailable expected artifacts remain explicit records.
- For symlinked artifacts, hash target bytes and retain logical path, link
  target, resolved path, size, and target SHA256 metadata.
- Append artifact records as stages produce them, then perform a final
  closeout re-hash/validation for terminal completeness.
- Make an artifact record unique using durable dataset/selector row identity,
  scientific role, and run-relative path. Reject collisions and path-only
  joins.

### Validation and retry behavior
- Fail closed when a ledgered artifact hash changes before downstream
  consumption. Emit a precise mismatch record and retain the original ledger.
- Reject duplicate, missing, or mismatched dataset-row identities before the
  consuming stage. Preserve all observed records and diagnostics rather than
  silently normalizing or dropping them.
- Model corrections and retries as new immutable attempts linked to the
  original by parent/supersedes provenance; never repair the original ledger
  in place.
- Expose the same validator as a reusable library and a reproducible CLI gate.

### Runtime and configuration capture
- Capture a secret-safe allowlist of relevant environment/configuration keys,
  package versions, explicit seeds, resolved executables, and Slurm metadata;
  do not serialize the entire environment by default.
- Record structured command argv, working directory, effective relevant
  environment, and launcher/script hash rather than relying on a shell string
  alone.
- Use one execution-context schema for local and Slurm runs. Record explicit
  local context and nullable scheduler fields when Slurm is absent.
- Identify external tools by resolved path, reported version or probe output,
  and executable/script SHA256 when accessible.

### Agent's Discretion
- Exact JSON field names, canonicalization implementation, and TSV column
  ordering, provided they preserve the decisions above and remain stable.
- The concrete run-directory layout and whether the closeout validator is
  invoked by the main pipeline, a wrapper, or both.
- The precise representation of diff/content references, provided secrets
  are excluded and a dirty run can be distinguished from its Git HEAD.
- The smallest fixture data and CLI shape needed to prove changed-artifact
  and row-identity rejection before scoring.

</decisions>

<specifics>
## Specific Ideas

- Reuse the deterministic, standard-library provenance foundations in
  `benchmark/scripts/investigation_provenance.py` rather than introducing a
  competing capture format.
- Preserve the current operational distinction between a readable run
  directory and a stable scientific contract hash.
- Keep the stable NACCESS + TMalign + external-Rosetta defaults unchanged;
  Phase 1 adds identity/provenance around them rather than changing them.

</specifics>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

- `.planning/PROJECT.md`
- `.planning/REQUIREMENTS.md`
- `.planning/ROADMAP.md`
- `.planning/STATE.md`
- `benchmark/scripts/investigation_provenance.py`
- `src/candidate_audit.py`
- `prism.py`
- `docs/STABLE_PIPELINE.md`
- `.agents/skills/project-memory/references/summary.md`
- `.agents/skills/project-memory/references/decisions.md`
- `.agents/skills/project-memory/references/open_questions.md`

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `benchmark/scripts/investigation_provenance.py`: standard-library helpers
  already provide deterministic file hashing, secret-redacted environment and
  seed capture, package/executable resolution, Git provenance, and Slurm
  resource capture.
- `src/candidate_audit.py`: append-only JSONL candidate audit pattern with
  explicit terminal statuses and metadata can inform the artifact/event
  writer without changing candidate acceptance logic.
- `benchmark/scripts/collect_pipeline_verification_baseline.py`: rejects
  duplicate claim IDs and duplicate `(run_root, relative_path)` artifact keys
  when merging independently produced shards.

### Established Patterns
- `prism.py` currently generates timestamp/PID alignment and candidate-audit
  paths and has opt-in stage event recording through
  `PRISM_STAGE_STATUS_PATH`; the Phase 1 work should wrap or extend these
  hooks without changing stable defaults.
- Existing benchmark tooling uses explicit TSV manifests, JSON provenance,
  SHA256 values, isolated `tmp/agent/` run roots, and preserved failed/retry
  evidence.
- The project treats missing, failed, partial, unavailable, and non-scoreable
  records as explicit evidence and does not treat directory presence as
  scientific success.

### Integration Points
- The pipeline entry point `prism.py` is the natural source for run identity,
  command/config capture, and run-root initialization.
- Input/template staging, alignment/transformation/refinement output writers,
  and optional comparison output are the materialization boundaries where
  artifact records can be appended.
- Downstream scoring and collectors are the consuming gates that must reject
  changed hashes and row-identity mismatches before interpreting results.

</code_context>

<deferred>
## Deferred Ideas

- Stage and candidate lifecycle taxonomy, including timeout/partial/refiner
  failure semantics — Phase 2.
- Identity-safe score joins, chain mappings, interface scope, and evaluation
  denominator audits — Phase 3.
- Reproducible isolated launchers and bounded Slurm task recovery — Phase 4.
- Linked evidence bundles and the full regression gate — Phase 5.
- Content-addressed artifact reuse, queryable provenance indexes, workflow
  engine migration, and learned ranking quality — V2 or later.

</deferred>

---
*Phase: 01-run-identity-and-manifest*
*Context gathered: 2026-07-29*
