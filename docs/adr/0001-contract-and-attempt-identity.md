# ADR-0001: Separate declared contract identity from execution attempt

**Status:** Accepted
**Date:** 2026-07-29
**Scope:** Phase 1 — Run Identity and Manifest

## Context

PRISM needs to distinguish what a run declared from what happened during one operational execution. A single hash cannot safely represent selectors, source/tool inputs, Slurm placement, timestamps, evolving artifact observations, and retry history: those values change at different times and for different reasons.

The current pipeline also runs in a deliberately dirty worktree and has both stable production defaults and exploratory/legacy evidence. A clean-tree requirement or a full-worktree hash would either block normal work or absorb generated benchmark clutter into scientific identity.

## Decision

Model provenance as separate but linked objects:

- **Declared contract:** the canonical scientific declaration. Its hash covers normalized and raw selectors, explicit dataset/template/input inventories, declared source/configuration fingerprints, thresholds, backend choices, and tool executable/version fingerprints.
- **Execution attempt:** one operational run identified by a readable unique `run_id`. It records host/interpreter, timestamps, command execution, allowlisted environment, Slurm facts when present, status, and links to the declared contract hash.
- **Artifact observation:** a row/role/path-bound observation emitted as artifacts materialize. It is not part of the launch-time contract hash.
- **Run closure:** final expected-artifact validation and the whole-ledger digest, linked to the attempt and contract.

Use an explicit source inventory: hash Git HEAD plus tracked diff state and only declared untracked source/config/template files. Record other dirty paths as observations, not as implicit contract inputs.

An operational retry with the same declaration receives a new `run_id` and retains the same contract hash. A retry that changes declared inputs, source inventory, thresholds, backends, or tool identity receives a new contract hash and a parent/supersedes link.

## Consequences

### Positive

- Slurm job IDs, timestamps, host placement, and partial outputs no longer invalidate the declared scientific identity.
- Repeated attempts can be compared while preserving the original failure and retry history.
- Dirty worktrees remain usable without silently treating generated artifacts as source inputs.
- Downstream consumers can require a contract-backed closure without needing to interpret every runtime field.

### Negative

- Consumers must carry both `contract_hash` and `run_id` instead of treating one identifier as sufficient.
- The source/configuration inventory must be explicit and maintained by launchers or callers.
- The implementation needs a closure record in addition to the launch manifest.

## Rejected alternatives

- **Hash the entire attempt:** would change identity as Slurm/status/artifact facts evolve and would make retries indistinguishable from changed declarations.
- **Hash the full worktree:** impractical in this repository because generated benchmark outputs and scratch evidence are numerous and mutable.
- **Use Git HEAD only:** does not identify active tracked edits or declared untracked source/configuration files.

## Verification implications

- Two executions with identical declared inputs/configuration/source/tool inventory must produce the same contract hash even when run IDs and Slurm metadata differ.
- Changing a declared source/configuration file must change the contract hash.
- Changing only run ID, host, timestamps, or Slurm fields must not change the contract hash.
- A changed artifact must be diagnosed during closure/consumer validation rather than rewriting the contract.

## Related records

- `.planning/phases/01-run-identity-and-manifest/01-CONTEXT.md`
- `.planning/phases/01-run-identity-and-manifest/01-RESEARCH.md`
- `benchmark/scripts/investigation_provenance.py`
