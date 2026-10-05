# Bounded goal

In the assigned isolated editing worktree, diagnose and improve exactly two PRISM-prescript runtime pipeline failures: (1) the current prism.py rejects --aligner usalign and the structural_aligner.py USalign path is disconnected; (2) MultiProt is blocked on the observed Seccomp/32-bit runtime combination. Recover project instructions and memory first. Preserve the canonical dirty checkout and all existing data. Do not download or build USalign or PRODIGY, do not use hosted/remote models, do not force-bypass Seccomp, and do not claim live tool validation when the executable/runtime is unavailable. Implement only bounded, test-backed code changes in the isolated worktree; if a blocker cannot be safely removed, make it explicit and fail closed with reproducible diagnostics. Keep stable TMalign/NACCESS/external-Rosetta defaults unchanged.

## Success criteria

- Canonical PRISM-prescript checkout is not modified; all code changes are confined to the assigned isolated worktree.
- USalign CLI/adapter behavior is either boundedly integrated with an explicit executable/configuration contract or has a clear fail-closed runtime boundary; synthetic tests cover the reachable behavior.
- MultiProt Seccomp/32-bit incompatibility is represented by explicit, reproducible diagnostics and tests without setting PRISM_MULTIPROT_FORCE or claiming compatibility.
- Focused deterministic tests and CLI/static checks run; failures, unavailable binaries, and validation limits are preserved in the handoff.
- The worker records IMPLEMENTED, TESTED, VALIDATED, and REVIEWED separately and leaves a durable checkpoint, report, handoff, and decision evidence.
