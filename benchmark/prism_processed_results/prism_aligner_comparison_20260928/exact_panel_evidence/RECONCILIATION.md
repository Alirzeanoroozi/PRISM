# VALAR reconciliation checkpoint — 2026-09-13

This package preserves the two evidence namespaces while treating them as one
selected underlying project.

* PRISM-prescript job `1657005`: live accounting shows `CANCELLED by 1365446`,
  zero runtime, and no durable output. The old Step 5 `blocked/unknown`
  description is superseded for this job.
* prism-refactoring run `run-2688ab679e0f47409f2a19e182483a75`: job `1657891`
  failed with exit `127:0` and retry `1657894` failed with exit `1:0` on
  `ai22`. The direct manifest still points at retry state `PENDING`/`WAITING_FOR_JOB`,
  but checkpoints, logs, handoff, and accounting show no worker completion
  report. The discrepancy is retained, not rewritten.
* Slurm controllers were up at this checkpoint, but current AI resources were
  fragmented by unrelated user allocations. No replacement was submitted.
* New isolated runner job `1659355` failed before staging because of a shell
  initialization bug (`/etc/bashrc` with `set -u`); fixed GTalign reruns are
  separately manifested as jobs `1659356` and `1659358`. Checked-panel
  staging job `1659357` failed closed on a missing interface and remains
  preserved. The exact USalign panels completed as jobs `1659360` (946-prefix)
  and `1659361` (19,062-entry calculated panel), using the transform-producing
  `-outfmt -1 -m -` invocation and eight CPU workers.
  The matched TMalign references completed as jobs `1659448` (946-prefix) and
  `1659449` (19,062-entry calculated panel), using eight CPU workers and
  per-side JSONL/no-drop records.

Authoritative inputs: `evidence/prism-prescript/`,
`evidence/prism-refactoring/run-2688ab679e0f47409f2a19e182483a75/`, and the
live accounting output recorded during this goal.
