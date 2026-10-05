# Step 4 report — optional PRODIGY assessment

## Outcome

`PRODIGY` is locally available and technically suitable for an optional ranking stage. The existing PRISM CLI gates it behind `--rank true --rank-method prodigy`; defaults remain baseline ranking disabled. The bounded assessment supports retaining the optional feature, subject to later full-pipeline and scientific-quality validation.

## Tool and input findings

The executable is `/home/rshadi25/.conda/envs/gtalign_env/bin/prodigy`. Local source is under `prodigy/src/prodigy_prot`, version `2.4.0`, commit `6c459597af54329b65b039ad2458c0976e09c43c`. Local CLI help and source establish a positional PDB/mmCIF input, `--selection` chain groups, and defaults of 5.5 Angstrom, 0.05 accessibility threshold, and 25.0 Celsius. The positional input must precede `--selection` because the option consumes multiple values.

The adapter constructs a combined complex with disjoint left/right chain namespaces, invokes PRODIGY, parses quiet affinity output, and ranks more-negative values first. These are local source/runtime findings, not external documentation claims.

## Error inventory

1. `PRODIGY-001`: executable and local API available; validated by direct Slurm invocation.
2. `PRODIGY-002`: command ordering, chain grouping, contact semantics, output ordering, and defaults established locally.
3. `PRODIGY-003`: ranking is feature-gated; baseline defaults remain unchanged.
4. `PRODIGY-004`: retained paired evidence shows top-k selection changes when opt-in ranking is enabled.
5. `PRODIGY-005`: isolated explicit state contract passes five state tests and a real success smoke.
6. `PRODIGY-006`: no-contact failure is reproduced and preserves the full candidate group.
7. `PRODIGY-007`: disabled behavior and canonical focused tests remain unchanged.
8. `PRODIGY-008`: two broader tests are blocked by relative fixture writes from the canonical read-only cwd.
9. `PRODIGY-009`: full-tree quality, DockQ comparison, and promotion remain unknown.

Full evidence, classifications, hashes, job IDs, and status distinctions are in `evidence/error-inventory.json`.

## State contract

The isolated candidate records `available`, `executed`, `failed`, `skipped`, and `not_configured` in `prodigy-state.jsonl`. Missing configuration returns the original candidates. A score failure returns the complete affected group rather than silently dropping it. No PRODIGY threshold was changed.

## Verification and limitations

The isolated focused suite passed 9 tests. Slurm job `1656993` exercised the local executable successfully with affinity `-65.827`. Slurm job `1656992` exercised a real no-contact failure and preserved both candidates. Canonical focused tests passed 8 tests, and the canonical Git status remained at 373 entries with the recorded unchanged status hash.

Graphify was used to locate and query the tool/path graph. It identified the PRODIGY adapter’s internal calls and the `prism.py` ranking call site, but did not produce a direct AST edge between the two modules; direct source inspection therefore remains authoritative for the complete route.

This step does not establish that PRODIGY affinity is a scientifically superior ranking criterion, nor does it validate a full production run. The isolated candidate must be reviewed against the real selector before any authorized promotion.
