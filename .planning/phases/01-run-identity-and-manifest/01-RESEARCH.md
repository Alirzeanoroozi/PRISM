# Phase 1: Run Identity and Manifest — Research

**Researched:** 2026-07-29
**Phase goal:** Define a durable, immutable identity for declared PRISM inputs,
configuration, runtime, tools, resources, and materialized artifacts.

## Don't Hand-Roll

| Problem | Recommended solution | Why | Provenance |
|---|---|---|---|
| File content hashing | Extend the existing chunked `hashlib.sha256()` helper in `benchmark/scripts/investigation_provenance.py`. | Python documents SHA-256 as a guaranteed constructor and hash objects as incremental byte digests; the current helper already reads bounded chunks. | [CITED: https://docs.python.org/3/library/hashlib.html] [VERIFIED: `investigation_provenance.py::sha256_file`] |
| Canonical manifest serialization | Reuse one canonical JSON function with sorted keys, compact separators, explicit UTF-8, and a documented non-finite-number policy; hash its exact UTF-8 bytes. | Python documents `sort_keys=True` for deterministic dictionary ordering and `separators=(',', ':')` for compact output. | [CITED: https://docs.python.org/3/library/json.html] [VERIFIED: `investigation_provenance.py::_canonical_json`] |
| Command capture/replay | Keep argv as a structured sequence and record cwd/effective allowlisted environment; retain a launcher/script digest separately. | Python recommends argument sequences over shell strings and fully qualified executables for reliability. | [CITED: https://docs.python.org/3/library/subprocess.html] [VERIFIED: `investigation_provenance.py::_redact_argv`] |
| Symlink identity | Use `Path.is_symlink()`/`readlink()` or `os.lstat(..., follow_symlinks=False)` for link metadata, and `Path.resolve()` plus normal file reads for target content. | The path entry and the bytes consumed are different facts. Python explicitly distinguishes following symlinks from inspecting the link itself. | [CITED: https://docs.python.org/3/library/pathlib.html] [CITED: https://docs.python.org/3/library/os.html#os.lstat] |
| Provenance vocabulary | Use a small local JSON/TSV schema with PROV-inspired entity/activity/derivation fields; do not add RDF/OWL as a Phase 1 dependency. | W3C PROV models entities, activities, agents, usage, generation, derivation, and invalidation, but allows domain-specific specializations. The project already uses JSON/TSV and needs a local CLI gate first. | [CITED: https://www.w3.org/TR/prov-o/] [CITED: https://www.w3.org/TR/prov-constraints/] [VERIFIED: existing JSON/TSV tooling] |
| Append-only retry history | Reuse `AppendOnlyLineage` and immutable-copy behavior as design patterns; create linked attempts instead of overwriting records. | The local lineage implementation rejects conflicting record IDs and terminal-event replacement, while immutable PDB copying refuses different bytes at an existing destination. | [VERIFIED: `benchmark/scripts/investigation_lineage.py::AppendOnlyLineage`] [VERIFIED: `benchmark/scripts/investigation_lineage.py::copy_pdb_immutable`] |

## Common Pitfalls

### Canonical JSON that is not actually canonical

**What goes wrong:** A manifest hash changes because dictionary ordering, whitespace, float handling, or path normalization changes even though the declared run is intended to be the same.

**Why:** `json.dumps` exposes choices such as `sort_keys`, `separators`, and `allow_nan`; a hash of an unspecified serialization is not a stable contract. [CITED: https://docs.python.org/3/library/json.html]

**How to avoid:** Centralize canonicalization, sort semantically unordered collections, reject or normalize non-finite numeric values, encode UTF-8 explicitly, and test byte-for-byte equality across input ordering variations.

### Treating a symlink path as either only a link or only a target

**What goes wrong:** A staged logical path appears unchanged while its target bytes changed, or the ledger loses that the pipeline consumed a symlinked staging view.

**Why:** `stat` follows symlinks by default, while `lstat` and `follow_symlinks=False` inspect the link entry. [CITED: https://docs.python.org/3/library/os.html#os.stat]

**How to avoid:** Record `path_kind`, logical relative path, link target text, resolved target path, target size, target SHA256, and explicit broken-link status. Hash target bytes for scientific identity while retaining link metadata for staging provenance.

### Assuming Git HEAD proves the worktree contents

**What goes wrong:** A run manifest points at a commit but omits tracked modifications, staged changes, or untracked input/configuration files.

**Why:** The existing `capture_git_provenance` records HEAD, porcelain status, and a `git diff HEAD` hash, but its diff command does not itself contain the bytes of untracked files. [VERIFIED: `investigation_provenance.py::capture_git_provenance`]

**How to avoid:** Preserve HEAD plus deterministic status/diff records and an explicit inventory/hash of declared untracked files that affect the run. Do not require a clean tree because this repository intentionally contains dirty user work.

### Using path-only artifact identity

**What goes wrong:** Two benchmark rows or scientific roles can produce the same relative filename, causing a downstream consumer to accept the wrong bytes.

**Why:** Existing baseline manifests key artifacts by `(run_root, relative_path)`, which is useful for a retained run but insufficient for the new row-level contract. [VERIFIED: `benchmark/scripts/build_pipeline_verification_baseline.py`]

**How to avoid:** Make the identity tuple include durable row/selector ID, scientific role, and run-relative path; reject duplicate keys before scoring. Keep the content hash as evidence, not as the only identity.

### Re-hashing only after successful completion

**What goes wrong:** A failed or timed-out stage leaves no auditable record of partial outputs, and a directory with stale files can look complete.

**Why:** Current pipeline output paths are created incrementally and current stage events are opt-in. [VERIFIED: `prism.py::record_stage_event`, existing `processed/` path construction]

**How to avoid:** Append artifact observations at materialization boundaries, record explicit missing/unavailable observations, and run a final closeout that validates expected identities and current bytes. Do not infer scientific success from directory presence.

### Leaking secrets through provenance

**What goes wrong:** An environment snapshot, config value, or command argument stores API tokens, passwords, or private values that should not be persisted.

**Why:** Runtime provenance is broad by nature, and secret names can occur in environment, nested config, and both `--key value` and `--key=value` argv forms.

**How to avoid:** Keep an allowlist, recursively redact secret-bearing names, redact both argv forms, test serialized output for sentinel secrets, and never fall back to full-environment capture by default. [VERIFIED: `investigation_provenance.py::_redact_value`, `::_redact_argv`, `capture_environment`]

### Confusing immutable contract identity with an operational attempt

**What goes wrong:** Re-running a corrected or failed job overwrites evidence, or identical declared inputs are treated as the same operational attempt.

**Why:** A content hash answers “what was declared,” while a run ID answers “which attempt produced these records.” W3C PROV also distinguishes entities, activities, and derivations rather than collapsing them into one identifier. [CITED: https://www.w3.org/TR/prov-o/] [VERIFIED: Phase 1 CONTEXT.md]

**How to avoid:** Keep readable per-attempt `run_id`, canonical contract hash, and explicit parent/supersedes links for retries. Never mutate the original ledger in place.

## Existing Patterns in This Codebase

- **Deterministic provenance capture:** `benchmark/scripts/investigation_provenance.py` already captures Git state, selected environment variables, seeds, package versions, executable paths, command redaction, config normalization, Slurm fields, file hashes, and stable JSON output. Extend it rather than creating a second runtime schema.
- **Append-only event identity:** `benchmark/scripts/investigation_lineage.py` defines fixed-schema lineage records, content-derived record IDs, terminal-event protection, and explicit failure reasons. Reuse its immutability and retry concepts; do not merge its candidate/ranking-specific fields into the Phase 1 artifact schema.
- **Existing artifact TSVs:** `benchmark/scripts/build_pipeline_verification_baseline.py` writes `artifact_manifest.tsv` rows with run root, relative path, SHA256, and byte size; `collect_pipeline_verification_baseline.py` rejects duplicate claim and artifact keys. Add row identity and role rather than replacing these conventions wholesale.
- **Template preflight:** `preflight_template_assets()` records duplicate counts, missing assets, format validity, size, and SHA256, including explicit missing rows. Its deterministic TSV writer is a useful model for artifact inventory output.
- **Pipeline integration hooks:** `prism.py` already derives a run ID for GTalign/candidate-audit paths and exposes opt-in stage event recording through `PRISM_STAGE_STATUS_PATH`. Integration should be opt-in or backward-compatible until the stable path is covered by fixtures.
- **Regression style:** `tests/test_investigation_provenance.py` uses temporary directories, controlled environments, deterministic output comparisons, and sentinel secret assertions. New manifest/ledger tests should follow this isolated pattern.

## Recommended Approach

Build the Phase 1 contract as a standard-library provenance module that consolidates `investigation_provenance.py` and the relevant immutable/TSV patterns, then expose a small CLI for manifest emission and pre-consumption validation. Keep the manifest JSON canonical and human-readable, keep the artifact ledger row-oriented, and make all records explicit about row identity, role, path kind, size, hash, and missing/error state. Integrate first through an isolated fixture and a wrapper/opt-in pipeline path; only then add the stable `prism.py` hook, preserving existing defaults and dirty-worktree behavior.

The minimum end-to-end proof should create a temporary run with raw and normalized selectors, a symlinked artifact, a missing expected artifact, and a retry link; it should emit deterministic manifest/ledger files, reject a modified artifact and a duplicated/mismatched row before a mock scoring consumer, and prove no secret sentinel appears in serialized provenance.

## Sources

- Python `hashlib`: https://docs.python.org/3/library/hashlib.html
- Python `json`: https://docs.python.org/3/library/json.html
- Python `subprocess`: https://docs.python.org/3/library/subprocess.html
- Python `pathlib`: https://docs.python.org/3/library/pathlib.html
- Python `os.stat`/`os.lstat`: https://docs.python.org/3/library/os.html#os.stat
- W3C PROV-O: https://www.w3.org/TR/prov-o/
- W3C PROV constraints: https://www.w3.org/TR/prov-constraints/

