# PRISM-prescript Step 3 pipeline map

## Canonical maintained path

The current canonical `prism.py` path is:

```text
input/download
  -> template load or generation
  -> surface extraction
  -> selected alignment provider
       tmalign      -> processed/alignment_tmalign/<run_id>/
       multiprot    -> processed/alignment_multiprot/<run_id>/
       gtalign      -> processed/alignment_gtalign/<run_id>/
  -> src.transformation.transformer(templates, alignment_dir=...)
  -> optional ranking
  -> optional refinement
  -> optional comparison
```

The canonical checkout is dirty by design and was not edited in Step 3. Its
recorded selector remains `tmalign` by default. The canonical source imports
`src.alignment`, `src.alignment_gtalign`, and `src.alignment_multiprot`; it has
no `StructuralAligner` import or call.

## Isolated USalign candidate path

```text
--aligner usalign --usalign-path <existing executable>
  -> align_usalign(query surface PDB, template interface PDB)
  -> USalign query.pdb interface.pdb -outfmt -1 -m -
  -> parse alignment triple, score, 5-column matrix rows
  -> build match_dict: interface residue -> query residue
  -> processed/alignment_usalign/<run_id>/*.json
  -> existing transformer(alignment_dir=<run-scoped directory>)
```

The candidate is in the isolated worktree
`/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/worktrees/run-1e6873e93a4f4082a236f7348218ecd5`.
It adds no shared `processed/alignment` alias for USalign. The parser keeps
the synthetic four-column fixture compatibility while accepting the observed
USalign 20241108 five-column matrix rows.

## StructuralAligner connectivity result

Graphify 0.9.29 was run on the project `src` tree. Graph hash:
`7c2909034fa0921ae38549bbdfaae89d36f5c8012d716e1510613e8e0b70fdec`.

The recorded query was:

```text
graphify path 'StructuralAligner' 'transformer()' --graph /tmp/prism-prescript-graphify-step2.shM8dh/graphify-out/graph.json
```

Graphify returned a generic four-hop path through the shared `json` import:

```text
StructuralAligner <--contains-- structural_aligner.py
  --imports--> json <--imports-- transformation.py --contains--> transformer()
```

This is not a call edge from `prism.py` to `StructuralAligner`. The direct
source inspection is authoritative for entry-point connectivity: the
maintained CLI calls `transformer` directly after the provider branch.

## Tool and transform evidence

- Existing executable: `/home/rshadi25/.conda/envs/gtalign_env/bin/USalign`.
- Reported version: `20241108`.
- Default `-outfmt -1` output had alignment and TM-score but no matrix.
- `-m -` emitted rows `m, t[m], u[m][0], u[m][1], u[m][2]` and the equation
  `X=t[0]+u[0][0]*x+u[0][1]*y+u[0][2]*z` (and analogous Y/Z rows).
- A real Slurm probe parsed 214 aligned residues and mapped second-file chain
  `L` residues to first-file chain `A` residues. Applying the parsed matrix to
  428 CA coordinates gave RMSD `0.7500643`.

## Boundaries

This map proves the provider contract and one real end-to-end adapter case. It
does not claim full-panel scientific equivalence or canonical promotion.
