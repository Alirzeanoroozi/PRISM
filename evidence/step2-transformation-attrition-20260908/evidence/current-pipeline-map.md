# Current alignment-to-transformation map

Classification: `EVIDENCE` for observed source relationships; `INFERENCE` only where the relationship is interpreted as pipeline flow.

## Graphify command and output

Graphify was found at `/home/rshadi25/.local/bin/graphify`, version `0.9.29`. Because the repository's `graphify-out/graph.json` was absent, a code-only graph was built outside the project:

```bash
GRAPHIFY_TMP=$(mktemp -d /tmp/prism-prescript-graphify-step2.XXXXXX)
graphify extract /scratch/rshadi25/GitHub/PRISM-prescript/src \
  --code-only --no-cluster --out "$GRAPHIFY_TMP"
```

Observed output: 44 code files, 501 nodes, 1,148 edges. Graph hash: `7c2909034fa0921ae38549bbdfaae89d36f5c8012d716e1510613e8e0b70fdec`.

The most useful query was:

```bash
graphify query "Which functions connect structural alignment outputs to transformation.py, especially load_alignment, process_pair_for_template, create_transformed_pair, and pair_has_acceptable_clashes?" \
  --context call --budget 2200 --graph "$GRAPHIFY_TMP/graphify-out/graph.json"
```

Graphify reported these extracted call edges:

```text
transformer() -> process_pair_for_template()
process_pair_for_template() -> load_filter_assets(), _protocol_hotspots(), _protocol_contacts(), _alignment_variants(), alignment_passes_thresholds(), write_audit_record()
alignment_passes_thresholds() -> hotspot_analysis(), alignment_score_passes()
create_transformed_pair() -> apply_tm_transform()
pair_has_acceptable_clashes() -> read_ca_coordinates(), distance_calculator()
```

The shortest Graphify paths between the module nodes were only shared-import paths through `json`; the explicit call-context query is therefore the stronger architectural evidence:

```text
structural_aligner.py --imports--> json <--imports-- transformation.py
alignment_gtalign.py  --imports--> json <--imports-- transformation.py
alignment_multiprot.py --imports--> json <--imports-- transformation.py
```

## Source flow

1. `prism.py` selects the alignment backend and creates a run-scoped alignment directory. The GTalign path is `processed/alignment_gtalign/<run_id>/`; the MultiProt path is `processed/alignment_multiprot/<run_id>/`.
2. `prism.py` calls `src.transformation.transformer(templates, alignment_dir=alignment_output_dir)` after alignment.
3. `transformer()` reads `inputs.csv`, loads template interface JSON, stores chain sizes, and dispatches every receptor/ligand pair to `process_pair_for_template()`.
4. `process_pair_for_template()` loads the two sides for orientation `o1` and the swapped two sides for orientation `o2`. Missing alignment files are converted to zero-score fallback records, which fail the match/TM/coverage gate without raising out of the pair loop.
5. `alignment_passes_thresholds()` applies the stable match count, score-contract, hotspot, and coverage checks. The retained ledger's run-era case used `geometry_only_experimental`, so the contact/hotspot gate was diagnostic no-op evidence only.
6. `create_transformed_pair()` resolves canonical PDB IDs, writes `_L.pdb` and `_R.pdb`, checks transform success and file existence, then applies the CA-clash gate.
7. `pair_has_acceptable_clashes()` counts cross-pair CA distances strictly below 3.0 Å and rejects at the fifth clash. The retained case produced 71 final passes from 560 materialized attempts.

## Tool/path routing

- Python: `/home/rshadi25/.conda/envs/gtalign_env/bin/python`.
- USalign: `/home/rshadi25/.conda/envs/gtalign_env/bin/USalign` and lowercase alias `/home/rshadi25/.conda/envs/gtalign_env/bin/usalign`; bounded `-h` probe reports Version 20241108.
- Project-local `external_tools/USalign`: absent; do not infer that the environment executable is automatically routed by the current CLI.
- Retained MultiProt binary: `/scratch/rshadi25/GitHub/PRISM-prescript/external_tools/multiprot.Linux`; project memory classifies it as legacy and records the 32-bit/seccomp runtime blocker.
- Structural-alignment abstraction: `src/structural_aligner.py`; current retained-run alignment source: `src/alignment_gtalign.py`; transformation: `src/transformation.py`.

## Scope limitation

The Graphify graph was intentionally limited to `src/` and was not treated as proof of runtime behavior. Runtime counts come from the retained artifacts and the bounded Slurm reconstruction, not from Graphify or model confidence.
