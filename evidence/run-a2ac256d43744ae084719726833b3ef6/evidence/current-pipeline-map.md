# PRISM-prescript current alignment-to-transformation map

Status: `REVIEWED` for the current canonical dirty checkout; this is a
read-only map and does not assert that the pipeline is production-validated.

## Entry point and dispatch

`prism.py:64-162` runs the stages in this order:

```text
inputs.csv / PDB download
  -> surface extraction
  -> structural alignment
  -> transformation/filtering
  -> optional ranking (baseline or PRODIGY)
  -> selected refinement
  -> optional comparison
```

At `prism.py:109-141`, current dispatch is:

```text
--aligner tmalign  -> src.alignment.align(..., processed/alignment_tmalign/<run_id>)
--aligner multiprot -> src.alignment_multiprot.align_multiprot(...,
                         processed/alignment_multiprot/<run_id>,
                         multiprot_path, mode, params, solutions)
anything else      -> src.alignment_gtalign.align_gtalign(...,
                         processed/alignment_gtalign/<run_id>)
```

The `anything else` branch is a catch-all, not a USalign integration point.
The parser at `prism.py:235` rejects `usalign` before dispatch, with a
reproducible exit code `2`.

## Alignment output contract into transformation

`prism.py:152-158` passes the selected `alignment_output_dir` to
`src.transformation.transformer()`.

`src/transformation.py:43-80` then:

1. reads `inputs.csv`;
2. loads each template interface manifest;
3. normalizes query IDs;
4. opens `<alignment_dir>/<query>_<template>_<chain>.json`;
5. constructs both interface orientations in `process_pair_for_template()`;
6. applies aligner-specific score/match gates and optional published-protocol
   hotspot/contact checks;
7. calls `create_transformed_pair()` for accepted alignment pairs.

The transformation stage therefore requires every future aligner to emit the
existing JSON fields (`match_count`, `match_dict`, `translation`,
`rotation_mat`, `tm_score`) plus a clear aligner/score contract. It does not
discover an adapter dynamically.

## Existing aligner paths

- TMalign: `src/alignment.py`; default executable
  `external_tools/TMalign`; current dirty tree writes run-scoped outputs.
- GTalign: `src/alignment_gtalign.py`; caller supplies CPU/GPU executable
  path; output is parsed into the same downstream transform fields.
- MultiProt: `src/alignment_multiprot.py`; default
  `external_tools/multiprot.Linux`, an ELF32 executable. Retained logs show
  `Bad system call` under a restricted runtime. The current login shell is
  `Seccomp: 0`, so it is not a reproduction of that restricted context.
- USalign: no canonical adapter, no `structural_aligner.py`, no parser choice,
  and no dispatch branch found. This is an unresolved disconnected path, not
  a validated absence of any external USalign binary elsewhere.

## Optional PRODIGY path

When `--rank true --rank-method prodigy` is selected after transformation,
`prism.py:170-184` calls `select_top_candidates()` and the opt-in
`src/prodigy_ranker.py`. The repository contains the vendored `prodigy/`
source, but `prodigy` is not on `PATH`; no live scorer validation was run.
Baseline ranking remains the default.

## Graphify evidence

Bounded in-memory Graphify AST extraction over `prism.py`, `src/alignment.py`,
`src/alignment_multiprot.py`, `src/alignment_gtalign.py`, and
`src/transformation.py` produced 70 nodes. Extracted edges include:

- `prism.main -> src_alignment_align` (call, `prism.py:L111`)
- `prism.main -> src_alignment_multiprot_align_multiprot` (call,
  `prism.py:L114`)
- `prism.main -> src_alignment_gtalign_align_gtalign` (call,
  `prism.py:L131`)
- `prism.main -> src_transformation_transformer` (call,
  `prism.py:L153`)
- `src_transformation_transformer ->
  src_transformation_process_pair_for_template` (call,
  `src/transformation.py:L58`)

The retained `graphify-out.backup-20260729-154649/graph.json` is available as
read-only historical evidence. The current `graphify-out/graph.json` is
missing, and no rebuild or overwrite was performed.

## Evidence boundary

This map establishes source-level flow and reproducible blockers. It does not
validate a live USalign parser, a successful MultiProt execution under the
restricted runtime, or PRODIGY execution.
