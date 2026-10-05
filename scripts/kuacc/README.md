# KUACC portable PRISM submitters

These wrappers are an isolated submission layer for the two existing portable
bundles. They do not modify the bundle source tree or canonical PRISM outputs.
The default bundle root is the parent directory of this directory, so a copy
under `.../prism-new/kuacc_submission_20260919/` uses that bundle. Set
`PRISM_BUNDLE_ROOT` explicitly when the wrappers are staged elsewhere.

The live KUACC association verified on 2026-09-19 is account `users` with the
`mid` partition and `users` QOS. The wrappers deliberately avoid the unavailable
`cosbi`, `kutem`, and `kutem_gpu` partitions and do not pin a node. GTalign
requests one generic GPU; the exact eligible node is selected by Slurm.

Runtime defaults and overrides:

- `RUN_PYTHON` defaults to
  `/kuacc/users/rshadi25/.conda/envs/prism_portable_20260919/bin/python`,
  which contains the validated PRISM CPU dependencies, including `freesasa`.
  Override it only with another validated environment.
- The wrappers prepend the selected environment's `lib/` to
  `LD_LIBRARY_PATH`; the bundled TMalign, USalign, and GTalign executables
  require a newer C++ runtime than the KUACC system default.
- CPU wrapper: `DATASET_DIR`, `TEMPLATE_LIMIT`, `RUN_ROOT`, and `SCENARIO`.
  `ALIGNER` defaults to `tmalign` and accepts `tmalign`, `multiprot`, or
  `usalign`.
- GTalign wrapper: `DATASET_DIR`, `TEMPLATE_LIMIT`, `RUN_ROOT`, and
  `SCENARIO`; `GTALIGN_BIN` defaults to the bundle's optional executable.
  `GTALIGN_BIN` can point to a conda-provided executable, but its embedded CUDA
  target must match the allocated GPU generation; `gtalign_env` on KUACC is
  currently `sm_75` and is not compatible with the available Tesla K20/K40/K80
  nodes.
- Scoring wrapper: `RUN_ROOT`; `DOCKQ_PYTHON` defaults to
  `/kuacc/users/rshadi25/.conda/envs/prism_dockq_20260919/bin/python`.

The default template list is `templates/template_panel_full.txt`. Use a new
run root for every arm and preserve the completed run root before scoring.
Slurm stdout/stderr are placed under
`/scratch/users/rshadi25/valar-remote-runs/prism-submissions/`.

The 2026-09-19 fix-1 smoke also patched the isolated bundle helper copies so
USalign resolves from `code/PRISM-prescript/external_tools/USalign`, where the
executable is staged. The pre-fix helper copies are retained in the KUACC
smoke evidence directory.

Before any real submission, validate without creating a job:

```bash
sbatch --test-only --export=NONE run_cpu.sbatch
sbatch --test-only --export=NONE run_gtalign.sbatch
sbatch --test-only --export=NONE score_models.sbatch
```

`--test-only` validates scheduler directives only. It does not validate the
runtime environment or scientific outputs. The older
`/kuacc/users/rshadi25/bin/python` and `gtalign_env` remain untouched; the
versioned prefixes above are the validated runtime pair.
