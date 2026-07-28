# Conventions

## Python style

The current code uses readable module-level functions, lowercase
`snake_case` identifiers, uppercase constants, and docstrings for public or
contract-sensitive helpers. Newer benchmark modules use type annotations,
`pathlib.Path`, dataclasses, and `from __future__ import annotations`; older
pipeline modules are more permissive and use string paths. There is no
repository-wide formatter, linter, or static-type configuration.

## CLI and configuration

Command-line interfaces use `argparse` and are guarded by
`if __name__ == "__main__"`. Boolean CLI values are parsed explicitly in
`prism.py` via `parse_bool` to avoid Python’s `bool("false")` behavior. Runtime
variation is commonly exposed as `PRISM_*` environment variables, with source
defaults documented in `docs/STABLE_PIPELINE.md`.

When adding a setting, prefer an explicit CLI option or documented environment
variable, preserve the stable default, and record diagnostic overrides as
diagnostic rather than silently changing production behavior.

## Filesystem and data contracts

PDB and JSON filenames are treated as interfaces between stages. Preserve
chain-qualified IDs, orientation, template ID, `dataset_row_id`, and mapping
direction. Use isolated output roots and write derived data to new paths. Audit
and provenance records should retain status, reason, hashes, command metadata,
and explicit missing/failed states instead of dropping records.

Benchmark identity is row-based (`dataset_row_id`); do not collapse unrelated
rows by normalized PDB name. Ranking labels must join on row identity and model
SHA256, not on path alone. Report requested receptor–ligand cross-interface
DockQ separately from GlobalDockQ for multichain cases.

## Errors and subprocesses

The codebase mixes fail-fast exceptions with stage-level logging and explicit
failure rows. New validation and audit code should fail closed, preserve the
original exception/command where possible, and distinguish unavailable,
not-scoreable, failed, and successful outcomes. External subprocesses should
record return codes, stderr summaries, timeouts, and output paths; do not infer
scientific success merely from a surviving directory.

## Parallelism

Alignment and scoring use bounded thread/process pools or Slurm-internal worker
parallelism. NACCESS is serialized because its executable uses fixed filenames.
Avoid unbounded task submission, shared mutable output directories, and nested
parallelism that exceeds the allocated CPU count.

## Compatibility rules

The current default is NACCESS + TMalign + external Rosetta. GTalign, FreeSASA,
MultiProt, PyRosetta, FiberDock, ranking, and DockQ comparison are explicit
options. Keep the legacy Python 2 MultiProt/FiberDock tree separate and never
substitute a modern helper binary when the goal is historical equivalence.
