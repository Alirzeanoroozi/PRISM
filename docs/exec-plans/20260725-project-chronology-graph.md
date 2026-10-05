# Build an Auditable Project Chronology Graph

This ExecPlan is a living document. Keep `Progress`, `Surprises & Discoveries`,
`Decision Log`, and `Outcomes & Retrospective` up to date as work proceeds.

## Purpose / Big Picture

Produce a small, path-qualified Graphify graph that represents the sequence of
PRISM-prescript work without overwriting the repository-wide graph. The graph
must preserve `tmp/agent` as evidence while excluding copied environments,
vendored dependencies, caches, and prior Graphify corpora.

## Progress

- [x] Audit the existing graph and identify the lack of temporal metadata.
- [x] Define the include/exclude boundary with project memory.
- [x] Add a deterministic chronology-corpus builder.
- [x] Generate and validate the separate directed Graphify graph.
- [x] Record durable output paths and limitations in project memory.

## Surprises & Discoveries

- Observation: the existing `graphify-out/graph.json` has 208,846 nodes and
  337,510 links, but no temporal links or temporal node metadata.
  Evidence: direct JSON audit on 2026-07-25.
- Observation: `tmp/agent` is the dominant source of graph noise because a
  minority of run directories contain copied environments and dependencies.
  Evidence: prior node audit found 128,904 vendor-like paths among 172,831
  `tmp` nodes.
- Observation: a dated run-directory prefix proves the calendar day but not
  the order of events within that day.
  Evidence: most run IDs have `YYYYMMDD-...`, not a complete timestamp or
  terminal event sequence.

## Decision Log

- Decision: create `docs/chronology/` as a new derived artifact rather than
  modifying `graphify-out/`.
  Rationale: preserve the old graph and ensure an independently auditable
  chronology source.
  Date/Author: 2026-07-25 / Codex + user direction.
- Decision: model dated run directories and dated project documents as events;
  derive event state only from recorded status artifacts.
  Rationale: chronology needs explicit, reproducible temporal relations rather
  than inference from graph traversal.
  Date/Author: 2026-07-25 / Codex + user direction.
- Decision: date active-memory snapshots from their file modification time,
  while retaining filename/heading dates for historical documents.
  Rationale: active memory is a living snapshot and should appear on its actual
  revision date rather than remain fixed at its first consolidation date.
  Date/Author: 2026-07-26 / Codex + user direction.

## Outcomes & Retrospective

Completed 2026-07-25: `tools/build_project_chronology.py` produced 156 events
over 22 calendar dates, with 183 nodes and 334 directed edges. Graphify's
multigraph diagnostic found no missing endpoints, self-loops, duplicate edges,
or edge collapse. The original repository-wide `graphify-out/graph.json` hash
was unchanged.

2026-07-26 maintenance: active-memory dates now derive from file modification
time, so chronology advances when the operational memory is revised. The active
memory gained a compact verification/debugging path and the legacy seccomp
execution boundary; these are durable instructions, not run logs.

2026-07-26 environment correction: the chronology builder imports Graphify and
must use the interpreter recorded at `graphify-out/.graphify_python`; it is not
part of the pipeline `gtalign_env` environment.

2026-07-28 maintenance: after the matched MultiProt/TMalign calibration,
native-gate replay, matched refiner replay, and active-memory revision, the
chronology builder produced 168 events, 199 path-qualified nodes, 360 directed
edges, and 23 explicit `before` relations.
Graphify integrity checks found no dangling endpoints, self-loops, duplicate
edges, or endpoint collapse. Same-day ordering remains intentionally
unasserted.

## Context and Orientation

- Source evidence: `tmp/agent/`, `docs/exec-plans/`,
  `docs/memory-archive/`, and `.agents/skills/project-memory/references/`.
- Builder: `tools/build_project_chronology.py`.
- Derived output: `docs/chronology/`, including `manifest.tsv`, generated
  event notes, Graphify extraction JSON, and `graphify-out/`.
- Stable pipeline files, benchmark outputs, and the existing `graphify-out/`
  are read-only inputs to this task.

## Plan of Work

Create a deterministic inventory of dated project evidence. Each dated run
directory or dated project document becomes an event with path-qualified ID,
recorded date, evidence count, and conservative status. Generate explicit
`before` edges between verified calendar dates, plus `occurred_on`,
`has_status`, and a `current_as_of` anchor. Do not infer a same-day event
order. Use Graphify's
build/report/export components to create a directed graph from this
deterministic extraction.

## Concrete Steps

1. From `/scratch/rshadi25/GitHub/PRISM-prescript`, run
   `$(cat graphify-out/.graphify_python) tools/build_project_chronology.py --root . --out docs/chronology`.
2. Inspect `docs/chronology/manifest.tsv` for source count, exclusion reasons,
   source paths, date origin, and observed status.
3. Run `graphify query` against
   `docs/chronology/graphify-out/graph.json` for the most recent event and its
   preceding event.
4. Run `graphify diagnose multigraph` against the same graph.

## Validation and Acceptance

- The original `graphify-out/graph.json` hash is unchanged.
- `docs/chronology/manifest.tsv` includes dated `tmp/agent` runs and excludes
  vendor/cache/Graphify-derived paths by an explicit reason.
- The directed graph has at least one `before` relation and every event has a
  path-qualified source location.
- A query using the chronology graph returns dated event nodes rather than
  vendored/test-environment nodes.

## Idempotence and Recovery

- The builder rewrites only `docs/chronology/`, a derived output directory.
- Rerunning from unchanged evidence produces the same event order and IDs.
- If a generated result is unsuitable, remove only `docs/chronology/`; source
  evidence and the existing Graphify graph remain unaffected.

## Artifacts and Notes

- `docs/chronology/manifest.tsv` is the audit trail.
- `docs/chronology/events/*.md` are source-linked event summaries.
- `docs/chronology/graphify-out/graph.json` is the directed graph.

## Interfaces and Dependencies

- Python standard library for inventory/metadata parsing.
- Installed Graphify Python package for graph construction, clustering, JSON
  export, diagnostics, and report generation.
- No network, package installation, Slurm, or pipeline execution is required.
