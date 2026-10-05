# Architecture Research

**Domain:** Reproducible HPC protein–protein docking pipelines
**Researched:** 2026-07-28
**Confidence:** HIGH for the existing PRISM boundaries; MEDIUM for future workflow-engine adoption

## Component Boundaries

### System Overview

```text
┌──────────────────────────────────────────────────────────────┐
│ Research run definition                                      │
│ inputs • selectors • templates • thresholds • tool versions  │
└──────────────────────────────┬───────────────────────────────┘
                               │ frozen manifest
┌──────────────────────────────▼───────────────────────────────┐
│ Current PRISM stage adapters                                 │
│ input → surfaces → alignment → transform/filter → refine     │
└──────────────┬──────────────────────┬────────────────────────┘
               │ artifacts             │ status/provenance
┌──────────────▼──────────────┐  ┌────▼─────────────────────────┐
│ PDB/PDBx, JSON, PDB outputs │  │ append-only ledger + hashes   │
└──────────────┬──────────────┘  └────┬─────────────────────────┘
               │                       │
┌──────────────▼──────────────────────▼────────────────────────┐
│ Evaluation and evidence layer                               │
│ DockQ/iRMSD • mapping audit • completeness • report bundle  │
└──────────────────────────────────────────────────────────────┘
```

### Component Responsibilities

| Component | Responsibility | Typical Implementation |
|-----------|----------------|------------------------|
| Run definition | Freeze input rows, chain roles, templates, thresholds, backend, environment, and output root | Versioned YAML/JSON/TSV manifest plus Git revision |
| Stage adapters | Invoke current Python functions and external binaries with explicit contracts | Existing `prism.py` and `src/*` modules; add a shared execution/status wrapper |
| Artifact store | Keep source, intermediate, final, and evaluator files isolated and addressable | Run-scoped directories, stable filenames, SHA256 inventory |
| Status ledger | Explain every stage/candidate outcome | JSONL/TSV append-only records with status, reason, return code, and timestamps |
| Evaluation adapter | Compute and preserve scoped DockQ/iRMSD results | DockQ JSON, mapping table, interface scope, evaluator version, hashes |
| Evidence bundle | Give researchers a completion verdict and navigable artifacts | Markdown summary plus machine-readable tables and raw-output links |
| Scheduler runner | Map independent tasks to Slurm without losing identity | Array task IDs or bounded internal workers; per-task `exit.json` |

## Data Flow

### Primary Data Flow

```text
inputs.csv + frozen template manifest
      → source materialization and chain validation
      → surface artifacts
      → aligner JSON/transforms
      → transformed candidate PDBs + filter audit
      → refined model PDBs + refiner status
      → DockQ/iRMSD JSON/TSV + mapping audit
      → evidence bundle + completion verdict
```

### Data Flow Description

| Flow | Source | Destination | Format | Notes |
|------|--------|-------------|--------|-------|
| Source capture | RCSB/PDB archive or curated native files | Isolated `processed/pdbs/` and source inventory | PDB/mmCIF + TSV/JSON | RCSB exposes entry, entity, assembly, and chain-level API objects; retain retrieval metadata ([RCSB APIs](https://www.rcsb.org/docs/programmatic-access/web-apis-overview)). |
| Alignment | Query surfaces + template interfaces | Alignment stage | JSON + transform files | Include backend, command, version, match contract, and raw stdout/stderr. |
| Transformation | Alignment records + template assets | Candidate workspace | PDB + JSONL audit | Preserve filter mode, thresholds, orientation, match coverage, and clash result. |
| Refinement | Candidate partner PDBs | Final model workspace | PDB + status/score metadata | Require return code, expected output, score gate, and output hash. |
| Evaluation | Final models + curated native complexes | Scored evidence | DockQ JSON, score tables, mapping records | DockQ can emit per-interface data and JSON; do not collapse scope prematurely ([DockQ](https://github.com/wallnerlab/DockQ)). |
| Batch execution | Task manifest | Slurm task directories | `exit.json`, logs, artifacts | Use `SLURM_ARRAY_TASK_ID` or explicit task IDs; cap array concurrency and retain job metadata ([Slurm](https://slurm.schedmd.com/job_array.html)). |

### Edge-case flow

Failures remain rows in the ledger: missing source, unavailable alignment,
filter rejection, timeout, partial refiner output, failed score, no native
interface, and explicit `not_scoreable`. A scheduler cancellation overrides
surviving artifacts in the completion verdict.

## Build Order

Suggested implementation sequence based on dependencies:

| Order | Component | Dependencies | Rationale |
|-------|-----------|--------------|-----------|
| 1 | Run manifest and identity schema | None | Define the row, source, command, backend, and output identities before emitting more evidence. |
| 2 | Shared stage execution/status wrapper | #1 | Every adapter must write consistent started/completed/failed/timeout records. |
| 3 | Artifact inventory and hash ledger | #1, #2 | Make joins and rerun comparisons content-aware rather than path-only. |
| 4 | Refiner and evaluator output contracts | #2, #3 | Refinement and scoring are the highest-risk external boundaries; record return codes and explicit gates. |
| 5 | Completeness/mapping audit | #3, #4 | Validate denominators, row identity, native mappings, and scope before interpretation. |
| 6 | Researcher evidence bundle | #4, #5 | Present the trustworthy, human-readable completion verdict after machine checks exist. |
| 7 | Bounded Slurm runner improvements | #1–#5 | Scale the already-defined contracts; do not scale ambiguous output semantics. |
| 8 | Workflow-engine prototype | #1–#7 | Only then compare Snakemake/Nextflow migration against equivalent current behavior. |

## Integration Points

### External Integrations

| Integration | Type | Protocol | Auth | Notes |
|------------|------|----------|------|-------|
| RCSB PDB | Data API/download | HTTPS REST/GraphQL/files | None for public data | Use chain/entity semantics and record retrieved bytes/metadata; RCSB offers JSON APIs and BinaryCIF model subsets. |
| NACCESS/FreeSASA | Local executable/library | Filesystem/subprocess | None | NACCESS fixed filenames require serialized or isolated execution. |
| TMalign/GTalign/MultiProt | Local executable adapters | Subprocess + JSON/text/files | None | Backend-specific score contracts must not be conflated. |
| Rosetta/PyRosetta/FiberDock | Local executable/library | Subprocess/Python API | License/runtime dependent | Capture versions, paths, return codes, and output gates. Rosetta’s official docking docs identify protocol source and integration tests ([RosettaDock](https://docs.rosettacommons.org/docs/latest/application_documentation/docking/docking-protocol)). |
| Slurm | Scheduler | `sbatch`/`squeue`/`sacct` | Cluster account/QOS | Job completion is operational state, not scientific validity. |

### Internal Boundaries

| Boundary | Left Side | Right Side | Contract |
|----------|-----------|------------|----------|
| Input identity | CSV/native manifest | Downloader/stager | `dataset_row_id`, raw selectors, normalized chain IDs, source hash |
| Alignment | External aligner | Transformer | Alignment JSON with backend, match records, score contract, transform, raw-output identity |
| Candidate | Transformer | Refiner | Candidate PDB pair, orientation/template metadata, filter/audit status |
| Refiner | External backend | Evaluator | Canonical model path, return code, score-gate decision, model hash |
| Evaluation | DockQ/iRMSD | Evidence bundle | Scoped metrics, model/native hashes, chain mapping, evaluator version, status |

## Recommended Project Structure

```text
src/
├── pipeline/          # stage orchestration and shared contracts (future)
├── provenance/        # manifests, hashes, ledgers, completion verdicts
├── adapters/          # external aligner, refiner, and evaluator adapters
├── eval/              # DockQ/iRMSD and mapping logic
└── ranking/           # opt-in candidate selection and offline evaluation
benchmark/
├── scripts/           # manifest builders, runners, collectors, audits
├── jobs/              # Slurm templates and resource profiles
└── schemas/           # durable TSV/JSON contracts (future)
tests/                 # unit, contract, and isolated integration fixtures
docs/                  # stable operations, decisions, and evidence reports
results/<case>/<run>/  # derived artifacts; never overwrite raw inputs
```

This structure is a recommendation for incremental extraction, not a mandate
to reorganize the existing repository in one migration.

---
*Architecture research for: reproducible HPC docking pipelines*
*Researched: 2026-07-28*
