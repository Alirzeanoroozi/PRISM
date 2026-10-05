# PRODIGY candidate ranking

PRODIGY is an opt-in ranking scorer for transformed PRISM candidates. It is
not a replacement for the stable NACCESS + TMalign + external-Rosetta path and
it is not enabled by `--rank` alone.

## Login-node setup

The repository is cloned at `prodigy/`. Install it into a separate environment
because the current verified DockQ environment uses a different NumPy/FreeSASA
contract:

```bash
cd /scratch/rshadi25/GitHub/PRISM-prescript
python -m venv /scratch/rshadi25/tmp/prism-prodigy-env
/scratch/rshadi25/tmp/prism-prodigy-env/bin/python -m pip install -e ./prodigy
```

Alternatively, use a user-managed conda environment. Do not modify the
verified DockQ environment in `benchmark/prism_processed/env/`.

## Pipeline use

```bash
PRISM_PRODIGY_EXECUTABLE=/scratch/rshadi25/tmp/prism-prodigy-env/bin/prodigy \
PRISM_RANK=true \
PRISM_RANK_METHOD=prodigy \
PRISM_TOP_K=1 \
python prism.py
```

Equivalent CLI options are `--rank --rank-method prodigy
--prodigy-executable <path> --top-k 1`.

For each transformed receptor/ligand pair, PRISM writes a chain-renamed
combined PDB, PRODIGY stdout/stderr, and a JSON score record under
`processed/ranking/prodigy/`. The score is the predicted affinity in
kcal·mol⁻¹; lower (more negative) values rank first. The exact command,
input hash, return code, and failure status are retained.

If any candidate in a receptor/ligand group cannot be scored, PRISM preserves
the entire group instead of silently promoting a partial ranking. PRODIGY
scores are exploratory ranking evidence and must not be treated as native
DockQ labels or as evidence of improved biological quality without an
independent, frozen native-complex evaluation.

## Paired smoke result

Evidence retained under
`tmp/agent/20260730-prodigy-ranking-paired-test/summary-corrected.json`:

- Case: `5zngA,4eylA`
- Template/orientations: `1a0cCD`, `o1` and `o2`
- Without ranking: both transformed candidates are forwarded.
- With `rank_method=prodigy` and `top_k=1`: only `o1` is forwarded.
- PRODIGY affinities: `o1 = -65.827 kcal/mol`, `o2 = -65.274 kcal/mol`.

The adapter originally failed this smoke because the command placed the input
PDB after `--selection`; PRODIGY's `--selection` parser consumes all following
tokens as selection values. The fixed order is:

```bash
prodigy -q --distance-cutoff 5.5 --acc-threshold 0.05 --temperature 25.0 \
  combined.pdb --selection A,B C,D
```
