#!/home/rshadi25/.conda/envs/boltz_env/bin/python
"""Run the existing DockQ extension with a NumPy-1.x-compatible interpreter."""

import sys

import numpy  # noqa: F401  # load NumPy 1.26 before DockQ's compiled extension

sys.path.append("/home/rshadi25/.conda/envs/gtalign_env/lib/python3.11/site-packages")

from DockQ.DockQ import main


if __name__ == "__main__":
    main()
