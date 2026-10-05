# Stable PRISM pipelines

## Stable operational recipe

The only maintained recipe documented as stable is TMalign + NACCESS + external
Rosetta, using the current `gtalign_env` interpreter and defaults from
`docs/STABLE_PIPELINE.md`:

```bash
PRISM_PIPELINE_PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python \
  bash benchmark/scripts/run_prism_pipeline_smoke.sh
```

The smoke is an execution/wiring check and may legitimately produce zero
candidates. It is not scientific validation.

## Bounded candidates

* USalign: isolated candidate only; use the transform-producing equivalent of
  `USalign A.pdb B.pdb -outfmt -1 -m -`. Expected output must include parsed
  match pairs, rotation matrix, translation, score fields, return code, and an
  explicit failure/empty record. It is not promoted or stable.
* PRODIGY: isolated opt-in ranking candidate only, version 2.4.0. Affinity is a
  ranking feature, not a DockQ metric. No stable command is claimed until a
  same-set top-k experiment measures overhead, savings, regret, and quality.
* GTalign CPU/GPU: current opt-in implementations pass retained software/smoke
  checks, but CPU/GPU score and candidate-set parity is unresolved. No exact
  current-panel stable benchmark command is claimed here.
* MultiProt, FiberDock, FreeSASA, and PyRosetta remain supported/experimental
  or legacy choices with explicit evidence gaps; MultiProt's proxy score must
  remain separate from TM-score.

Promotion threshold for every provider: focused regression tests, deterministic
small smoke, real Slurm run, stage/no-drop manifests, reproducibility metadata,
quantitative output checks, matched cross-arm comparison, and independent review.
