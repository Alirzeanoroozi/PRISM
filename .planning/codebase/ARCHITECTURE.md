# Architecture

## System shape

The maintained system is a synchronous, file-oriented Python pipeline rather
than a web service. `prism.py` is the orchestration entry point and invokes
stage functions in a fixed order:

```text
input CSV/PDB download
  -> optional template analysis/generation
  -> surface extraction
  -> structural alignment
  -> transformation, protocol/geometry filtering, clash checks
  -> optional deterministic candidate ranking
  -> selected refinement backend
  -> optional DockQ/iRMSD comparison
```

Each stage consumes files created by the previous stage and returns a small
Python collection, commonly a list of pair paths or alignment records. Stage
start/completion/failure events can be appended as JSONL through
`PRISM_STAGE_STATUS_PATH`.

## Main layers

- **Orchestration:** `prism.py` parses CLI arguments, sets selected environment
  values, chooses the aligner/refiner, and records stage events.
- **Input/materialization:** `src/pdb_download.py` normalizes PDB-plus-chain
  identifiers, downloads source structures, and creates chain-filtered PDBs.
- **Template preparation:** `src/analyse_pdbs.py`, `template_filtering.py`,
  `template_generate.py`, `hotspot.py`, `interface.py`, and `contact.py`
  derive template interfaces, contacts, hotspots, and manifests.
- **Alignment adapters:** `alignment.py`, `alignment_gtalign.py`, and
  `alignment_multiprot.py` normalize different external tools into alignment
  JSON consumed by transformation.
- **Geometric filtering:** `transformation.py` applies aligner-specific match
  contracts, optional published-protocol assets, coordinate transforms, and
  clash thresholds; `candidate_audit.py` records provenance.
- **Refinement adapters:** Rosetta, PyRosetta, and FiberDock modules turn
  transformed partner structures into final candidate models.
- **Evaluation/ranking:** `compare.py`, `src/eval/`, `candidate_ranker.py`,
  `candidate_selector.py`, and benchmark scripts score, label, audit, and
  select candidates.

## Data and state flow

The primary filesystem state is under `processed/` and `templates/` during a
run. Inputs and templates are staged into an isolated workspace; alignments,
transformation intermediates, refinement outputs, and evaluation records are
named using PDB IDs, chains, templates, and orientations. Benchmark workflows
add manifest-driven tables, hashes, native complexes, score artifacts, and
lineage/audit records under `benchmark/` or `tmp/agent/`.

## Abstraction boundaries

External tools are wrapped by small adapter functions rather than a shared
plugin framework. The common boundary is file format plus tuple/dict metadata:
alignment JSON exposes matches and transforms, transformed PDB pairs are passed
to refiners, and output filename parsing reconstructs comparison metadata.
Candidate audits and stage-status files provide the main observability layer.

## Legacy boundary

`working_version/Multiprot-new/prism-fiberdock-cli/` is a Python 2
MultiProt/FiberDock reference implementation. It is not part of the maintained
Python 3 execution path and should not be modified when changing current
pipeline behavior. Current FiberDock support is an explicit Python 3 adapter,
not proof of historical equivalence.

## Architectural risks

The orchestration and many modules use module-level mutable state, import-time
directory creation, relative paths, and environment variables read at import
time. These choices make isolated work directories, process boundaries, and
explicit provenance important when adding features or running parallel jobs.
