# Concerns

## High-impact execution risks

- External tools are not fully pinned or installed by the repository. TMalign,
  GTalign, NACCESS, Rosetta, MultiProt, FiberDock, PyRosetta, DockQ, and
  FreeSASA each have separate runtime/version requirements.
- The current pipeline is strongly dependent on the repository root as the
  working directory and on relative paths such as `templates/`, `processed/`,
  and `external_tools/`. Running from another directory can fail or write to
  an unintended location.
- NACCESS and some legacy tools use fixed filenames or architecture-sensitive
  binaries. Parallel execution and restricted login-node/seccomp contexts can
  create false failures unless workspaces and execution context are controlled.

## Reproducibility and data risks

- Three environment recipes coexist and have different pinning. Historical
  installed environments and old launcher paths remain in the tree; use the
  verified interpreter and record versions rather than trusting directory
  names.
- PDB inputs are downloaded from a live archive without a repository lockfile
  or source checksum in the basic downloader. Benchmark workflows must freeze
  source inventories and hashes before comparison claims.
- The worktree contains extensive generated outputs, historical evidence, and
  user changes. Do not clean or overwrite them as part of routine maintenance.

## Pipeline design risks

- Several modules create directories and read environment variables at import
  time and maintain module-level mutable state (`passed_pairs`, template sizes,
  and configuration constants). Repeated in-process runs can leak state.
- External Rosetta refinement historically did not retain enough per-candidate
  subprocess return-code and score-gate detail to explain partial/missing
  canonical outputs. Add observability before interpreting attrition.
- FiberDock has a known output-contract issue: valid energy output may be
  written under `fiberdock_energies.ref` while existing parsing expects another
  filename. Fix only with an isolated fixture and regression test.
- MultiProt’s RMSD-derived compatibility score is not calibrated to the
  TMalign/GTalign TM-score threshold. Do not lower stable thresholds or make a
  quality claim from mixed score contracts without a paired calibration.

## Benchmark and scientific risks

- Chain order, raw selectors, curated native roles, and normalized PDB names
  can disagree. Preserve unresolved rows explicitly and use complete bijective
  mappings for confirmatory scores.
- GlobalDockQ can be dominated by a receptor-internal interface in multichain
  cases. Ranking labels must use the requested receptor–ligand cross-interface
  scope and report global metrics separately.
- Ranking is opt-in and mechanically validated, but current evidence is not
  sufficient to claim general quality improvement or enable it by default.
- Historical current-versus-legacy comparisons are observational unless inputs,
  template panel, evaluator, and refinement controls are frozen and matched.

## Security and maintainability

- The pipeline launches executable paths and builds subprocess environments from
  configuration. Validate executable paths and avoid logging credentials or
  untrusted command fragments.
- No authentication boundary exists because this is a local/HPC scientific
  tool, but downloaded PDB content and generated paths should still be treated
  as untrusted filesystem inputs.
- There is no repo-wide lint, type, coverage, or CI gate; regressions can land
  through untested external-tool paths. Focused contract tests and reproducible
  smoke/audit artifacts are the main safeguards.
