# PyRosetta installation status

## Confirmed

- Current pipeline environment: `/home/rshadi25/.conda/envs/gtalign_env`, Python 3.11.13; `importlib.util.find_spec("pyrosetta")` returns `None`.
- Official quarterly candidate: `pyrosetta-2026.3+releasequarterly.5e498f1409-cp311-cp311-linux_x86_64.whl`, approximately 1.659 GB.
- West official mirror download started but reached approximately 203 kB/s, with an estimated 1:58:45 remaining; it was cancelled before producing a completed wheel.
- East official mirror initially failed TLS verification (`SSLCertVerificationError`); the explicit trusted-host retry reached the wheel but timed out after 120 seconds without a completed artifact.
- No existing PyRosetta wheel or installation was found under the searched user/project/dataset paths.
- No environment was modified and no unverified package was installed.

## Installation path implemented

Once an authorized wheel is staged, run:

```bash
benchmark/scripts/install_pyrosetta_prism.sh \
  /path/to/pyrosetta-2026.3+releasequarterly.5e498f1409-cp311-cp311-linux_x86_64.whl \
  <sha256>
```

The script creates `/home/rshadi25/.conda/envs/pyrosetta_prism`, refuses to
overwrite an existing environment, validates the Python ABI and optional
SHA256, installs with `--no-index --no-deps`, and imports PyRosetta.

The separate pipeline entrypoint is
`benchmark/scripts/run_pyrosetta_refinement.py`; the isolated array job is
`benchmark/jobs/pyrosetta_refinement_array.sbatch`.

## Sources

- Official downloads and quarterly wheel instructions: <https://www.pyrosetta.org/downloads>
- Official licensing information: <https://www.pyrosetta.org/home/licensing-pyrosetta>
- GTalign source/release information: <https://github.com/minmarg/gtalign_alpha>
