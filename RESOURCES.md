# Trustworthy PRISM provenance Resources

## Knowledge

- [Python `hashlib` documentation](https://docs.python.org/3/library/hashlib.html)
  Use for the guarantees and incremental behavior of SHA-256 file hashing.
- [Python `json` documentation](https://docs.python.org/3/library/json.html)
  Use for deterministic serialization choices such as `sort_keys` and `separators`.
- [Python `subprocess` documentation](https://docs.python.org/3/library/subprocess.html)
  Use for structured argv, executable resolution, and replayable command capture.
- [Python `pathlib` documentation](https://docs.python.org/3/library/pathlib.html)
  Use for resolving paths and distinguishing symlink metadata from target files.
- [W3C PROV-O](https://www.w3.org/TR/prov-o/)
  Use as the conceptual vocabulary for entities, activities, agents, usage, generation, and derivation.
- [W3C PROV constraints](https://www.w3.org/TR/prov-constraints/)
  Use for understanding identity and lifetime semantics of changing entities.
- [Phase 1 context](.planning/phases/01-run-identity-and-manifest/01-CONTEXT.md)
  The authoritative local record of decisions for this project phase.
- [Phase 1 ADRs](docs/adr/0001-contract-and-attempt-identity.md)
  Local decisions separating contract identity, attempts, artifacts, and closure.
- [Slurm `sbatch` documentation](https://slurm.schedmd.com/sbatch.html)
  Primary reference for batch scripts, `#SBATCH` directives, resource requests, and submission behavior.
- [Slurm `sinfo` documentation](https://slurm.schedmd.com/sinfo.html)
  Primary reference for inspecting partitions, node states, resources, and formatted summaries.
- [Slurm `scontrol` documentation](https://slurm.schedmd.com/scontrol.html)
  Primary reference for detailed partition, node, and job inspection.
- [Slurm GRES documentation](https://slurm.schedmd.com/gres.html)
  Primary reference for typed GPU/resource requests such as `--gres=gpu:tesla_v100:1`.

## Wisdom

- Not yet collected. A structural-biology researcher reviewing whether the evidence vocabulary matches real interpretation decisions would be valuable.

## Gaps

- No single external source defines the exact PRISM artifact schema; local contracts and fixtures must supply that domain-specific detail.
- The cluster's account/QOS policy is site-specific; always confirm it with `sacctmgr` and `sbatch --test-only` before a long run.
