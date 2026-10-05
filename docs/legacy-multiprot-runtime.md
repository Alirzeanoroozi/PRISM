# Legacy MultiProt runtime and smoke validation

The legacy PRISM arm is staged from the checked-out
`working_version/multiprot/external_tools` payload. The source tree and raw
benchmark files are not modified. Use
`benchmark/scripts/stage_legacy_tool_environment.py` to create a derived
environment and `benchmark/scripts/prepare_legacy_pipeline_smoke.py` to create
an isolated one-pair/one-template workspace.

The Python 2 compatibility runtime is:

- Python 2.7.15 from `/home/rshadi25/.conda/envs/tmalignRosetta/bin/python2.7`.
- NumPy 1.16.6 and PyMySQL 0.9.3 staged in a project-local site directory.
- `sitecustomize.py` exposes PyMySQL as the legacy `MySQLdb` import.
- MultiProt 1.6, POPS 1.5.3, and the checked-out FiberDock payload are hash
  recorded in `environment_manifest.json`.

The operational compatibility environment is the derived profile:

```sh
source tmp/agent/20260713-investigation-implementation/legacy-tool-environment-v5/activate.sh
python2 -c 'import numpy, MySQLdb; print(numpy.__version__)'
```

The historical profile is retained separately. Its derived wrapper is
relocated to the staged data directory, but the checked-out NACCESS binary
requires unavailable `libgfortran.so.3`; it is therefore blocked at surface
extraction. The operational profile explicitly uses
the repository's current NACCESS binary, which links against available
`libgfortran.so.5`. These profiles must not be pooled in scientific estimates.

Validation is performed with isolated KUTEM arrays:

```sh
sbatch benchmark/jobs/legacy_tool_probe_array.sbatch
sbatch benchmark/jobs/legacy_pipeline_smoke_array.sbatch
```

The first array independently probes NACCESS, POPS, MultiProt, and FiberDock.
The second compares the explicit NACCESS profiles using separate workspaces.
Each task writes `input_manifest.tsv`, `parameters.tsv`, `command.txt`, logs,
outputs, and `exit.json` under `tmp/agent/...`. A successful current-profile
controller smoke currently means `pipeline_plumbing_complete_no_candidates`:
all preprocessing, surface, alignment, transformation, and refinement setup
directories exist, but no candidate passed the historical filter, so no
FiberDock refinement was attempted. FiberDock's independent energy-only probe
does produce `resFile.ref`.
The environment manifest separately records `fiberdock_full_refinement=false`
because the bundled NMA/reduce helpers are 32-bit; the energy-only result does
not establish full refinement readiness.

The reference protocol lists Python/NumPy/MultiProt/NACCESS/FiberDock
requirements in `references/nprot.2011.367.md`. The exact legacy controller
also depends on deployment-specific database/template paths; the compatibility
snapshot replaces those boundaries with local event records and never performs
network, database, mail, or destructive cleanup actions.
