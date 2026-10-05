# Stack Research

**Domain:** Reproducible protein–protein docking and structural-bioinformatics pipelines
**Researched:** 2026-07-28
**Confidence:** HIGH for the existing tool contracts and official workflow capabilities; MEDIUM for adopting a new workflow engine in this repository

## Recommended Stack

### Core Technologies

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| Python | 3.11.x | Current pipeline and benchmark adapters | Matches the validated `gtalign_env` runtime and the pinned `environment.yaml`; avoids a disruptive language migration. |
| Biopython | 1.84–1.85 | PDB/mmCIF parsing and residue/chain operations | Already embedded across `src/` and benchmark code; preserves the project’s structural-data abstractions. |
| Slurm | Site-managed | HPC scheduling and resource accounting | The project already targets Slurm; official job arrays expose task IDs, bounded concurrency, and dependency semantics needed for auditable batches ([Slurm job arrays](https://slurm.schedmd.com/job_array.html)). |
| Git + immutable run manifests | Repository version | Code, configuration, input/output hashes, and provenance | NIH and DFG guidance emphasize version control, metadata, documentation, and traceability for research software ([NIH guidance](https://datascience.nih.gov/tools-and-analytics/best-practices-for-sharing-research-software-faq), [DFG principles](https://www.dfg.de/en/basics-topics/digital-topics/research-software/principles)). |

### Supporting Libraries

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| DockQ | 2.1.3 | Native-vs-model docking quality and interface metrics | Use for final scored artifacts; preserve JSON, chain mappings, interface scope, and version. The official tool supports per-interface output and JSON export ([DockQ repository](https://github.com/wallnerlab/DockQ)). |
| NumPy | 1.26.4 | Coordinates, transforms, and numerical metrics | Keep the validated project pin for deterministic array behavior. |
| pandas | 2.3.3 | Input and benchmark tables | Keep row-level `dataset_row_id` and explicit status columns through every transformation. |
| FreeSASA | 2.2.1 | Optional surface-area backend | Use only when explicitly selected; retain NACCESS for the stable legacy-compatible path. |
| Snakemake | 9.23.1 documentation line | Candidate orchestration/provenance layer | Evaluate for future manifest-driven benchmark DAGs and reports; do not replace the stable `prism.py` path before a paired migration proof. Snakemake documents scalable execution, software environments, reports, and provenance ([documentation](https://snakemake.readthedocs.io/en/stable/)). |

### Development Tools

| Tool | Purpose | Notes |
|------|---------|-------|
| `pytest` | Focused unit and contract tests | Test pure adapters, status contracts, hashes, mappings, and explicit failures without requiring every external binary. |
| Conda environment recipes | Python and scientific dependency setup | Prefer the maintained Python 3.11 recipe; record the exact environment and binary paths in run manifests. |
| `sbatch`, `squeue`, `sacct` | Reproducible HPC submission and monitoring | Record job and array IDs, partition, resources, node, command, and work directory. |
| `sha256sum`/Python hashing | Artifact identity | Hash source inputs, staged inputs, raw outputs, canonical outputs, and evaluator JSON before joining records. |

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| Existing Python stage adapters plus manifests | Nextflow 26.x | Use if portability across clusters/clouds and container-first dataflow become primary; Nextflow explicitly targets portable HPC workflows and implicit parallelism ([Nextflow](https://nextflow.io/)). |
| Snakemake for a future DAG layer | Nextflow | Prefer Nextflow for multi-platform deployment and channel-based dataflow; prefer Snakemake for Python-native rules, current repo integration, and report/provenance features. |
| PDB/PDBx source capture through RCSB services | Ad-hoc live URL construction only | Use RCSB REST/GraphQL/Search APIs when metadata, chain identity, or batch selection is needed; RCSB exposes entry, entity, assembly, and chain-level objects ([RCSB APIs](https://www.rcsb.org/docs/programmatic-access/web-apis-overview)). |
| DockQ 2.1.3 with project mapping adapters | Unscoped scalar docking scores | Use alternatives only for sensitivity analyses; canonical ranking labels must retain requested cross-interface scope and mapping provenance. |

## What NOT to Use

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| A new workflow-engine rewrite as the first reliability change | It expands migration risk while the current pipeline already has working stage boundaries and evidence contracts. | Harden current Python adapters and add a thin manifest/runner layer first. |
| Live, unversioned downloads as benchmark truth | RCSB data and identifiers evolve; unrecorded retrieval dates, selectors, and hashes make reruns ambiguous. | Freeze source inventories, chain selectors, retrieval metadata, and file hashes; use RCSB APIs where metadata is required. |
| Global DockQ alone for multichain ranking | Official DockQ reports multiple interfaces and mappings; a global aggregate can hide the requested receptor–ligand behavior. | Store per-interface/cross-interface results plus GlobalDockQ separately. |
| Unbounded local multiprocessing or unconstrained Slurm arrays | It can exceed memory, file-system, or QoS limits and obscure which task failed. | Use bounded workers, isolated workspaces, array concurrency limits, and explicit task status. |
| Replacing NACCESS/TMalign/Rosetta defaults with experimental backends | It would invalidate existing comparisons before quality and provenance are independently established. | Keep experimental aligners/refiners opt-in and record their contracts distinctly. |

## Versions

### Version Compatibility

| Package A | Compatible With | Notes |
|-----------|-----------------|-------|
| Python 3.11.13 | NumPy 1.26.4, pandas 2.3.3, Biopython 1.84 | Matches `environment.yaml`; confirm exact host interpreter before a run. |
| Python 3.11 | DockQ 2.1.3, FreeSASA 2.2.1 | Existing maintained recipe; keep DockQ’s virtual-environment entry path intact. |
| Rosetta 2022.42 | Current external-Rosetta adapter | Requires site module loading and explicit database/prepack/dock paths. |
| Snakemake 9.x | Conda/Slurm-backed rules | Treat as a future orchestration dependency; validate network filesystem persistence before adoption. Snakemake warns about metadata-file scaling and SQLite lock contention on shared filesystems ([provenance guidance](https://snakemake.readthedocs.io/en/v9.19.0/executing/provenance.html)). |

### Installation

```bash
# Existing project environment recipe
conda env create -f environment.yaml

# Python dependencies are intentionally kept in the repository recipes.
# Rosetta and licensed PyRosetta require separately authorized installation.

# Optional future workflow evaluation
conda install -c bioconda -c conda-forge snakemake
```

**Confidence notes:** version values for this repository come from local
recipes and validated project memory; workflow-engine feature claims come from
the official documentation linked above. No workflow-engine migration is
recommended until an isolated prototype proves provenance, failure recovery,
and Slurm behavior on this cluster.

---
*Stack research for: reproducible protein–protein docking pipelines*
*Researched: 2026-07-28*
