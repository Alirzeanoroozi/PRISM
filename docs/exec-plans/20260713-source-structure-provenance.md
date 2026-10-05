# Build a source and structure provenance manifest

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`, `Decision Log`, and `Outcomes & Retrospective` aligned with the implementation.

## Purpose / Big Picture

Add an isolated benchmark provenance utility that reads the three checked-in benchmark CSVs, preserves each source row and its receptor/ligand/native chain context, inventories exact local selector files and exact benchmark archive members, and writes deterministic `source_manifest.tsv` and `structure_validation.tsv` artifacts. The utility must make source identity auditable when local files collide or archive members are unavailable.

## Progress

- [x] Route to PRISM-prescript and read project guidance/memory.
- [x] Inspect CSV schemas, local PDB layout, archive members, and Biopython availability.
- [x] Implement `benchmark/scripts/build_investigation_source_manifest.py`.
- [x] Add focused provenance tests under `tests/`.
- [x] Run focused tests and a bounded CLI smoke test.
- [x] Review changed regions and report findings/uncertainties.

## Surprises & Discoveries

- `T_Rigid.csv`, `T_medium.csv`, and `T_difficult.csv` use the same nine-column header and contain raw selectors such as `1FGN_LH`, `3LZT_`, and `1IK0_A(10)`.
- `benchmark/data/pdbs` contains 977 PDB files, including repeated basenames across difficulty directories and chainwise materializations; these are provenance collisions to retain, not deduplicate silently.
- `benchmark/originals/benchmark5.5.tgz` contains 1,084 non-metadata PDB members under `benchmark5.5/structures/`.
- The archive README reports “BENCHMARK VERSION 5.0” despite the archive path/name being `benchmark5.5`; the script will preserve paths/hashes and leave this version discrepancy visible.
- The system Python lacks Biopython, while `/home/rshadi25/.conda/envs/gtalign_env/bin/python` has Biopython 1.84. Focused validation will use the latter environment.

## Decision Log

- Decision: Keep selector-source and native-archive provenance as separate source roles and do not select one as a hidden fallback for the other.
  Rationale: The request requires source substitutions to be visible and native benchmark constituents have an exact archive naming contract.
  Date/Author: 2026-07-13 / Codex
- Decision: Emit one row per source candidate, including collisions and archive-member/extracted representations where applicable, with stable input-row order and role order.
  Rationale: A manifest that collapses candidates cannot explain duplicate basenames or source disagreements.
  Date/Author: 2026-07-13 / Codex
- Decision: Treat an absent archive prefix/member, missing local selector, or parser failure as an explicit unresolved/failed status; `--strict` returns nonzero after writing artifacts.
  Rationale: Provenance must fail closed without discarding the originating benchmark row.
  Date/Author: 2026-07-13 / Codex

## Outcomes & Retrospective

Implemented the standalone source/structure provenance track and focused tests. The test fixture proves exact local collision retention, archive-prefix recording, SHA-256 generation, Biopython metadata capture, deterministic reruns, and strict unresolved behavior. The real-repository `--limit 1` smoke emitted four source/validation rows for the first rigid case; all were resolved, parsed with Biopython 1.84, hashed, and byte-identical on rerun. No Slurm jobs or production pipeline code were touched.

The full three-dataset run was intentionally not launched in this focused turn; the bounded smoke and temporary-fixture tests cover the CLI and failure contracts without producing a broad derived artifact.

## Context and Orientation

The target repository is `/scratch/rshadi25/GitHub/PRISM-prescript`. The three inputs are `benchmark/data/T_Rigid.csv`, `benchmark/data/T_medium.csv`, and `benchmark/data/T_difficult.csv`. Selector structures are under `benchmark/data/pdbs`; native benchmark archive members are read directly from `benchmark/originals/benchmark5.5.tgz` rather than substituted with unrelated local files. New code belongs in `benchmark/scripts/`, and focused tests belong in `tests/`.

Each CSV row will receive a deterministic ID such as `rigid:000001`. The output will retain the raw CSV row as JSON in a TSV field, along with raw selector strings, parsed selector chains, native complex text, parsed native receptor/ligand chains, and difficulty. Selector roles are matched only to exact chainwise/local selector materializations. Native roles are matched only to exact archive members `<complex_id>_r_b.pdb` and `<complex_id>_l_b.pdb` using the native complex's PDB code.

## Plan of Work

Implement a standalone script using the standard library plus Biopython's PDB parser. Build deterministic indexes for local PDB candidates and archive members, expand each benchmark row into selector/native source records without deduplication, hash every available byte source, parse every available PDB, and write stable TSVs with explicit statuses. Add tests using temporary CSVs, colliding local chainwise files, and a tiny tar.gz archive. Validate the script directly and through `--limit 1`; do not launch Slurm or alter production pipeline code.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, add the standalone script and focused tests with `apply_patch`.
2. From the same directory, run `/home/rshadi25/.conda/envs/gtalign_env/bin/python -m pytest -q tests/test_build_investigation_source_manifest.py`.
3. From the same directory, run the script with `--repo-root . --output-dir <repo-local tmp/agent path> --limit 1`; inspect both TSV headers/rows and rerun to confirm byte-identical outputs.
4. Run a focused syntax/import check in the Biopython-capable environment and inspect `git diff --stat`/`git status --short` to ensure only the requested new paths plus the ExecPlan changed.

## Validation and Acceptance

- All three CSVs are read in deterministic order and every emitted row retains its dataset ID, raw selectors, chain assignments, native complex, difficulty, and source-row JSON.
- Duplicate local candidates remain separate rows and are labeled as collisions.
- Archive rows carry a nonempty exact `archive_prefix` whenever a matching member exists; missing archive data is explicit and never replaced by a local file.
- Both TSVs contain SHA-256 values for available sources and Biopython parser/version/structure metadata or explicit unresolved/parse-failed status.
- Repeated runs on unchanged inputs are byte-identical.
- `--strict` fails closed on unresolved/invalid structures after writing artifacts; `--limit` bounds rows for smoke tests.
- Focused tests pass without modifying existing production files, project memory, or unrelated tests.

## Idempotence and Recovery

The script creates/overwrites only the two named files in the requested output directory. Use a new `tmp/agent/<run-id>/` directory for smoke outputs. Rerunning with unchanged inputs is safe and deterministic. No source files, archive contents, production code, or memory files are modified.

## Artifacts and Notes

Implementation paths are `benchmark/scripts/build_investigation_source_manifest.py` and `tests/test_build_investigation_source_manifest.py`. Smoke artifacts are under `tmp/agent/20260713-source-structure-provenance-final-smoke/` and its deterministic rerun directory; each contains `source_manifest.tsv` and `structure_validation.tsv`.

## Interfaces and Dependencies

CLI: `--repo-root PATH`, `--output-dir PATH`, `--strict`, and `--limit INTEGER`. Runtime dependency: Biopython/PDBParser for structure validation. No Slurm, network access, project-module imports, or new package installation is required.
