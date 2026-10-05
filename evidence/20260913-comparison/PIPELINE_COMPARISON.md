# PRISM pipeline comparison

## Stage-by-stage interpretation

The matrix uses the common stage order `alignment → transformation/filtering →
ranking → refinement → evaluation`. Empty fields are unknown, not zero. A zero
candidate count is only a result when the stage manifest explicitly records it.

### Alignment

The retained logical 946 panel is the only same-panel alignment comparison with
complete pair counts: both TMalign and GTalign produced 1,892/1,892 pair records.
The runtime artifact does not record whether GTalign was CPU or GPU, so it cannot
answer the GPU-vs-CPU question. The historical large run suggests GTalign GPU is
operationally attractive (906 s) but it used a 19,855-entry historical list and
is not a clean comparison to current 19,948/19,062 panels. The bounded USalign
probe established the transform-producing invocation; the new exact-panel
USalign runs extend that evidence to alignment-only timing.

The new isolated alignment-only runs provide a cleaner current-panel search
comparison. On the checked-list 946-template prefix, GTalign GPU took 7.89 s
total (4.81 s search), TMalign took 12.31 s total (11.48 s search), and USalign
took 22.23 s total (21.50 s search). On the calculated panel, 19,062 source
entries yielded 19,058 materialized templates: GTalign GPU took 144.62 s total
(73.96 s search), TMalign took 157.18 s total (144.95 s search), and USalign
took 452.11 s total (439.97 s search), each CPU reference using eight workers.
These are alignment-stage throughput observations only: GTalign emits a
batched search output whereas the pairwise runners record one result per
interface, so candidate-set agreement and downstream quality are still
unmeasured.

### Transformation and filtering

The validated retained 1gte ledger contains 2,997 alignment JSONs, 560 complete
orientation attempts, 489 clash rejections, and 71 retained geometry passes.
It also records 1,163 structural missing-partner sides and 1,877 residual
unconsumed sides. The cause of 73,251 missing/unwritten potential sides is not
known. This is a historical geometry-only ledger, not a provider comparison.

### Ranking

Bounded evidence shows candidate-load reduction but not speedup. In the retained
current-tree example, baseline forwarded 4 candidates and refined 2; a top-1
ranked arm forwarded 1 and refined 1, but wall time was 62.7 s versus 56.3 s.
PRODIGY selection and no-contact preservation are tested; DockQ retention,
regret, and end-to-end savings remain unknown. MultiProt's RMSD-derived proxy is
kept out of TM-score comparisons.

### Refinement and evaluation

External Rosetta and PyRosetta have historical successful outputs, but their
per-candidate failure observability and evaluator parity are incomplete. DockQ
2.1.3 is wired in the repository-local scoring environment. Historical GTalign
GPU quality is high in the scored subset, but only 2/10 BM5.5 pairs were directly
scorable in one report because of chain mapping; this limits causal conclusions.

## Scalability and reliability

Historical TMalign generated 714,780 alignment JSONs for the large run, whereas
GTalign GPU generated 19,561 alignment hits and 60/64 transforms in the cited
arms. This indicates a likely alignment/search and filesystem/output bottleneck,
but stage-separated timing and candidate-set agreement for an exact common panel
are missing. The highest-value optimization is therefore a clean, batched,
stage-instrumented GTalign GPU versus USalign versus TMalign run at both exact
panels, with reusable preprocessing and no-drop manifests.

## Comparison rule

Only rows with matching query set, template manifest/hash, surface, thresholds,
ranking, refiner, evaluator, source revision, and complete stage manifests may
support a causal claim. All other rows remain historical, observational, or
not comparable in the CSV.

## Final decision question

The current evidence supports only a provisional answer: GTalign GPU + external
Rosetta is the fastest and best observed large-library arm, TMalign + NACCESS +
external Rosetta is the recommended reproducibility-first default, and PRODIGY
should remain opt-in pending same-set top-k DockQ/regret and end-to-end timing.
Before promotion, exact current-panel counts/hashes, matched inputs, stage
timings, USalign compatibility, GPU/CPU parity, complete evaluator mapping,
and independent review are still required.
