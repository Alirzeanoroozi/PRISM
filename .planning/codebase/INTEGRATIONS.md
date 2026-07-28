# Integrations

## External data services

- **PDB structure download:** `src/pdb_download.py` reads `inputs.csv` (or the
  `PRISM_INPUTS_CSV` override) and downloads four-character PDB entries over
  HTTPS from the PDBJ archive using `urllib.request`. It then materializes
  chain-qualified files under `processed/pdbs/`.
- **Native benchmark sources:** benchmark scripts stage curated/native PDBs and
  preserve row-level chain selectors. `dataset_row_id` is the durable identity
  for benchmark records; a normalized full PDB is not automatically a valid
  substitute for a chain-qualified source.
- There are no application databases, authentication providers, webhooks, or
  message queues in the maintained pipeline.

## Structural-analysis tools

- **NACCESS:** default surface extraction backend. It uses fixed working-file
  names, so `src/naccess_utils.py` serializes calls with a process-local lock
  and moves outputs into the processed workspace.
- **FreeSASA:** explicit alternative selected with `--surface_backend freesasa`
  and optionally `--freesasa_python`. It uses the project’s residue reference
  table and is not silently substituted for NACCESS.
- **TMalign:** default alignment backend, invoked from `src/alignment.py` and
  producing alignment JSON plus transformation data.
- **GTalign:** optional CPU/GPU backend, invoked from
  `src/alignment_gtalign.py`; output is isolated under
  `processed/alignment_gtalign/<run-id>/`.
- **MultiProt:** current experimental backend in `src/alignment_multiprot.py`;
  the separate 32-bit legacy binary and helper stack remain compatibility
  assets and require a suitable execution context.

## Refinement and scoring

- **External Rosetta:** `src/rosetta_refinement.py` invokes prepack and docking
  binaries. The stable setup requires `module load rosetta/2022.42` plus
  explicit `PRISM_ROSETTA_PREPACK`, `PRISM_ROSETTA_DOCK`, and
  `PRISM_ROSETTA_DB` settings.
- **PyRosetta:** explicit opt-in refiner. Its adapter records runtime and file
  metadata and does not act as a fallback for external Rosetta.
- **FiberDock:** explicit backend through `src/fiberdock_refinement.py` and
  `external_tools/fiberdock/`; it depends on helper binaries and declared
  energy/output-file contracts.
- **DockQ and iRMSD:** optional final evaluation in `src/compare.py` and the
  benchmark scripts. Current canonical scoring uses the repository’s verified
  DockQ 2.1.3 environment and preserves cross-interface and global scopes
  separately.

## HPC integration

Benchmark and reproducibility workflows use Slurm batch templates and runners
under `benchmark/jobs/` and `benchmark/scripts/`. Heavy alignment, refinement,
scoring, and large benchmark loops are intended for compute nodes. Login nodes
are used for inspection, setup, downloads, and submission. Internal worker
parallelism is used to reduce pressure on the cluster QoS job limit.

## Integration contracts

Most integrations communicate through files: CSV manifests, chain-qualified
PDBs, alignment JSON, transformed PDBs, refined PDBs, score CSV/JSON, JSONL
audits, and stage-status records. Filename conventions and relative paths are
therefore API contracts. A run should use an isolated working directory and
retain command, environment, input/output hashes, and scheduler metadata.
