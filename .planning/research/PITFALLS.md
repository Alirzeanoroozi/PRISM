# Pitfalls Research

**Domain:** Reproducible HPC protein–protein docking and structural evaluation
**Researched:** 2026-07-28
**Confidence:** HIGH for the observed repository risks and official workflow guidance; MEDIUM for frequency estimates across the wider ecosystem

## Common Mistakes

| # | Mistake | Severity | Frequency | Impact |
|---|---------|----------|-----------|--------|
| 1 | Treating scheduler completion or surviving output directories as scientific success | CRITICAL | Frequent in batch research | Partial or canceled runs enter score tables and inflate completion or quality claims. |
| 2 | Joining models, natives, and scores by path or normalized PDB name instead of row and chain identity | CRITICAL | Frequent in multichain benchmarks | Metrics are assigned to the wrong complex, interface, or partner orientation. |
| 3 | Comparing methods with different inputs, template panels, thresholds, mappings, or evaluators | CRITICAL | Frequent in retrospective comparisons | Observational differences are misreported as aligner/refiner quality effects. |
| 4 | Recording only scalar DockQ/global aggregates | HIGH | Common when reports are optimized for one table | Requested cross-interface behavior, mapping ambiguity, and no-native-interface cases disappear. |
| 5 | Running external tools without capturing versions, commands, stderr, return codes, and output hashes | HIGH | Frequent in HPC prototypes | Failures cannot be reproduced or diagnosed after the job workspace changes. |
| 6 | Scaling an ambiguous pipeline before hardening its contracts | HIGH | Common when large benchmark work is urgent | Parallel failures multiply, shared files collide, and missing rows become hard to reconstruct. |

### Mistake Details

**1. Scheduler success is mistaken for scientific success**
- **What happens:** A Slurm job exits zero after a subprocess warning, partial
  output, or a canceled dependency, and the collector treats the directory as
  a valid prediction.
- **Why it happens:** Operational status and scientific completion are modeled
  as one boolean.
- **Example:** A Rosetta process leaves a raw PDB but fails the canonical score
  gate; a batch summary counts it as refined.
- **Fix cost:** High after publication or benchmark aggregation; low if a stage
  ledger is designed first.

**2. Identity is reduced to filenames**
- **What happens:** Chain-qualified selectors, model/native mapping, or
  duplicate benchmark rows collapse into one record.
- **Why it happens:** PDB IDs are convenient keys, while biological roles are
  row-level and interface-specific.
- **Example:** The same normalized entry appears with different receptor and
  ligand selectors but receives one shared score.
- **Fix cost:** Very high once labels and derived tables are distributed.

**3. Method comparisons are not paired**
- **What happens:** One backend uses a different panel, fallback, filter, or
  evaluator, yet score means are compared directly.
- **Why it happens:** Existing artifacts are easier to aggregate than a frozen
  paired replay.
- **Example:** Current TMalign + Rosetta and legacy MultiProt + FiberDock are
  compared despite different runtime and output contracts.
- **Fix cost:** High; requires new controlled runs and evidence ledgers.

**4. Interface scope is lost**
- **What happens:** A multimer score is reported as if it represented the
  requested receptor–ligand interface.
- **Why it happens:** Global values are simpler to rank and visualize.
- **Example:** DockQ’s multi-interface output is reduced to one aggregate even
  when a receptor-internal interface dominates it.
- **Fix cost:** Medium if raw JSON remains; high if only the scalar survives.

**5. External execution is opaque**
- **What happens:** A missing or malformed model is attributed to biology when
  it was caused by a path, module, architecture, timeout, or score parser.
- **Why it happens:** Logs are scattered and return codes are not retained in
  the result table.
- **Example:** A FiberDock energy file exists under a filename the parser does
  not inspect.
- **Fix cost:** Medium for new runs; high for historical runs without logs.

**6. Parallelism comes before contracts**
- **What happens:** Multiple workers share fixed tool files or append to
  ambiguous outputs, creating nondeterministic corruption.
- **Why it happens:** Throughput is optimized before workspace and identity
  isolation.
- **Example:** NACCESS fixed filenames collide, or array tasks write the same
  model path.
- **Fix cost:** High after large runs; low when per-task roots are designed up
  front.

## Warning Signs

Early indicators that a project is heading toward common pitfalls:

| Warning Sign | Indicates | Action |
|-------------|-----------|--------|
| A result table has no explicit failed/not-scoreable rows | Mistake 1 or 2 | Stop aggregation; compare task manifest to result ledger and restore missing states. |
| Joins use only PDB ID, basename, or path | Mistake 2 | Require `dataset_row_id`, chain selectors, model hash, and native hash in the join key. |
| Two methods have different source/template/evaluator manifests | Mistake 3 | Label the comparison observational or rerun under a frozen paired contract. |
| Only a global DockQ column is retained | Mistake 4 | Preserve raw JSON and requested interface rows before publishing summaries. |
| Logs lack executable version or return code | Mistake 5 | Add command/runtime capture and rerun a small diagnostic case. |
| Multiple tasks share a work directory or fixed filenames | Mistake 6 | Isolate task roots, serialize unsafe tools, and cap worker/array concurrency. |

## Prevention Strategies

Proactive measures to avoid the mistakes above:

| Strategy | Prevents | When to Apply | How |
|----------|----------|---------------|-----|
| Completion contract | #1 | Before any collector or scorer | Require terminal stage status, expected output integrity, evaluator result, and audit pass; scheduler state is only one input. |
| Durable identity schema | #2 | At manifest creation | Store row ID, raw selectors, normalized chains, template/orientation, model/native hashes, and mapping direction. |
| Frozen comparison manifest | #3 | Before a method comparison | Freeze inputs, templates, exclusions, candidate budget, thresholds, backend versions, random seeds where possible, and evaluator. |
| Scope-preserving score schema | #4 | During evaluator integration | Store GlobalDockQ, each requested cross-interface result, mapping, status, and original evaluator JSON. DockQ documents per-interface and JSON output ([DockQ](https://github.com/wallnerlab/DockQ)). |
| Execution metadata envelope | #5 | In every subprocess wrapper | Capture command, executable hash/version, environment whitelist, stdout/stderr paths, return code, timeout, output inventory, and hashes. |
| Isolated bounded execution | #6 | Before scale-out | Give every task a unique run root; use Slurm array limits or bounded internal workers; record array/task IDs ([Slurm](https://slurm.schedmd.com/job_array.html)). |
| Versioned research artifacts | All | Before publication or handoff | Keep code, workflow, environment, metadata, and documentation accessible and versioned, consistent with NIH/DFG research-software guidance ([NIH](https://datascience.nih.gov/tools-and-analytics/best-practices-for-sharing-research-software-faq), [DFG](https://www.dfg.de/en/basics-topics/digital-topics/research-software/principles)). |

## Domain-Specific Patterns

### Patterns That Look Right But Aren't

| Pattern | Why It Seems Good | Actual Problem | Better Approach |
|---------|-------------------|----------------|-----------------|
| “One row per PDB entry” | Simple benchmark table | Chain role and biological pair identity can differ within an entry | Use durable row and chain-role identity; retain normalized PDB as a secondary field. |
| “One aggregate score per model” | Easy ranking and plotting | Multichain interfaces and no-native-interface states are hidden | Store scoped interface rows and derive aggregates only as labeled summaries. |
| “Use a workflow engine and reproducibility is solved” | Snakemake/Nextflow provide useful provenance/scaling primitives | Engine metadata cannot repair ambiguous scientific contracts or wrong input joins; shared filesystems can also introduce persistence contention ([Snakemake provenance](https://snakemake.readthedocs.io/en/v9.19.0/executing/provenance.html)) | Define scientific identity/status schemas first, then evaluate an engine against them. |

### Patterns That Look Wrong But Work

| Pattern | Why It Seems Bad | Why It Actually Works | When to Use |
|---------|------------------|----------------------|-------------|
| Append-only JSONL alongside normalized tables | Duplicates information | Preserves event history and failed attempts while tables provide queryable summaries | Stage/candidate ledgers where retries and partial work matter. |
| One large internally parallelized Slurm job | Appears less modular than many jobs | Reduces scheduler/QOS pressure while retaining bounded workers and one resource envelope | Many independent short tasks on a cluster with tight job limits. |
| Keep legacy and current pipelines separate | Appears to duplicate tooling | Prevents incompatible contracts from contaminating quality claims | Historical compatibility or causal comparison investigations. |

---
*Pitfalls research for: reproducible HPC docking pipelines*
*Researched: 2026-07-28*
