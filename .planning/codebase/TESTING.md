# Testing

## Framework and inventory

The repository primarily uses pytest with tests under `tests/`. A smaller set
of benchmark helpers uses `unittest.TestCase` directly. The suite covers current
pipeline parsing, coordinate/transformation behavior, surface backends,
candidate auditing/ranking, benchmark manifests, scoring contracts, lineage,
and Slurm runner planning. Specialized scripts under `tests/*_replacement` and
`tests/*_runtime` are exploratory or environment-dependent rather than a
single always-green unit suite.

## Test patterns

- Pure helpers are tested with small dictionaries, synthetic coordinates, and
  temporary directories.
- `pytest.raises`, parametrization, fixtures, and monkeypatch/mocks are common.
- External tools are generally mocked or tested through dry-run/contract
  adapters rather than invoked by every unit test.
- Tests that require unavailable PDB/template assets or licensed/optional
  runtimes use explicit skips; inspect skip reasons before interpreting a full
  pass as end-to-end validation.
- Benchmark contracts emphasize explicit failure rows, hashes, mapping
  bijections, output integrity, and provenance rather than only numeric values.

## Useful commands

From the repository root, the standard invocation is:

```bash
python -m pytest -q
```

The stable documentation also defines a focused current-pipeline command over
transformation, ranking, and reranker tests. Heavy pipeline smoke tests use
`benchmark/scripts/run_prism_pipeline_smoke.sh` and the documented
`gtalign_env` interpreter; benchmark scoring and alignment should run through
Slurm according to project HPC guidance.

## Validation boundaries

No repository-wide coverage threshold or CI configuration was found. There is
no configured formatter/linter gate. Passing unit tests does not validate
external binary availability, Rosetta output quality, PDB download behavior,
cluster scheduling, or biological correctness. Those require isolated smoke,
stage-status, output-integrity, and evaluator audits.

## Adding tests

Add the smallest focused regression test in `tests/` for current source changes.
For a new benchmark contract, test both success and explicit failure or
not-scoreable paths, including row identity and hash behavior. Keep production
inputs and validated result artifacts read-only; use `tmp_path` or a new
`tmp/agent/<run-id>/` workspace for derived fixtures.
