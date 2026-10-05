# KUACC PRISM-prescript dependency preflight

Date: 2026-09-27

Remote root:
`/scratch/users/rshadi25/valar-remote-runs/prism-prescript-kuacc-preflight-20260927`

Smoke job: `3124094`, executed on `be01.kuacc.ku.edu.tr` in KUACC `mid`.

## Passed

- `rosetta/2022.42` resolved both Rosetta docking executables and the explicit
  database directory.
- MultiProt emitted its expected usage contract.
- FiberDock, NMA, and Reduce emitted their expected usage contracts on a
  scheduled compute node.
- FiberDock had no missing shared libraries according to `ldd`.
- FiberDock Perl helpers passed `perl -c`.
- Source and staged hashes match for MultiProt, FiberDock, NMA, and Reduce;
  see `manifest/tool_hashes.sha256` and `manifest/compute_smoke_3124094.sha256`.

## Scope limits

This is a dependency/launcher smoke test, not a scientific refinement test:
no PDB pair was available in the staged/source KUACC roots, so no model was
aligned, refined, or scored. No production MultiProt, Rosetta, FiberDock, or
DockQ workload was submitted. The existing KUACC submission tree also lacked
`refine_one_external.py`; a checked copy was staged only under the isolated
preflight root.

The full legacy MultiProt workspace was not copied because it is 297 MB and
contains large archive/database payloads. The required runtime/configuration
subset is staged under `tools/multiprot_legacy_runtime`; a partial interrupted
copy under `tools/multiprot_legacy` is not trusted or used.
