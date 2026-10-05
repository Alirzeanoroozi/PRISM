# PRISM-prescript current alignment-to-PRODIGY map

Status labels are deliberate: `IMPLEMENTED` means code exists in the isolated candidate; `TESTED` means a deterministic check ran; `VALIDATED` means the behavior was exercised with the local executable or Slurm fixture; `REVIEWED` means source, memory, Git state, and evidence were inspected. These labels are not interchangeable.

## Flow

```text
prism.py
  -> selected alignment backend
       tmalign | multiprot | gtalign
  -> processed/alignment* 
  -> src.transformation.transformer(...)
  -> passed_pairs
  -> optional ranking gate (--rank true)
       rank_method=baseline  -> baseline candidate selector
       rank_method=prodigy   -> PRODIGY adapter
                                  combine left/right PDBs
                                  rename chains into disjoint groups
                                  run local PRODIGY CLI
                                  parse quiet affinity
                                  rank lower/more-negative values first
                                  forward top_k per receptor/ligand group
  -> selected pairs
  -> refinement stage
```

The alignment stages are upstream producers of transformed candidate pairs; PRODIGY is downstream of transformation and upstream of refinement. It does not replace alignment or transformation.

## Source and tool locations

| Component | Location / command | Evidence classification |
|---|---|---|
| CLI and stage orchestration | `/scratch/rshadi25/GitHub/PRISM-prescript/prism.py` | `OBSERVATION` from source inspection |
| Candidate dispatch | `src/candidate_selector.py`, called from `prism.py` | `OBSERVATION` from source inspection |
| PRODIGY adapter | `src/prodigy_ranker.py` | `OBSERVATION`; isolated candidate contains explicit state contract |
| Local executable | `/home/rshadi25/.conda/envs/gtalign_env/bin/prodigy` | `EVIDENCE`; wrapper hash recorded in `error-inventory.json` |
| Local source | `/scratch/rshadi25/GitHub/PRISM-prescript/prodigy/src/prodigy_prot` | `EVIDENCE`; version `2.4.0`, commit recorded in inventory |
| Graphify graph | `/scratch/tmp/prism-step4-graphify-20260909/graphify-out/graph.json` | `EVIDENCE`; isolated code-only extraction |

## CLI contract

The local CLI accepts one positional PDB/mmCIF input followed by `--selection`. Because `--selection` consumes one or more values, the adapter must emit the input path before the selection arguments:

```text
prodigy -q --distance-cutoff 5.5 --acc-threshold 0.05 \
  --temperature 25.0 combined.pdb --selection LEFT RIGHT
```

The adapter combines the two input structures and gives the two sides disjoint chain namespaces. PRODIGY calculates contacts between the selected groups. Quiet output supplies the affinity value used for ordering; lower, more-negative values are preferred.

## Explicit optional states in the isolated candidate

| State | Trigger | Candidate disposition |
|---|---|---|
| `available` | executable resolves before scoring | continue |
| `executed` | all candidate groups score and selection completes | return selected top-k |
| `failed` | timeout, nonzero output, malformed result, or no contacts | preserve the complete affected group |
| `skipped` | no candidates are supplied | return empty input |
| `not_configured` | executable is absent | preserve all passed pairs |

The state log is `processed/ranking/prodigy/prodigy-state.jsonl` in a run. The default CLI remains disabled (`rank=false`, `rank_method=baseline`), so the normal path does not invoke PRODIGY.

## Graphify qualification

Graphify found the PRODIGY adapter node and its local calls (`_executable_is_available`, `_append_state`, `_group_key`, and `score_candidate`) and found the `prism.py` ranking-stage import/call site. Its shortest-path query did not establish a direct `prism.py -> prodigy_ranker.py` edge; this is recorded as a graph extraction limitation, not treated as evidence of a disconnected runtime. Direct source inspection is authoritative for the complete route through `candidate_selector.py`.
