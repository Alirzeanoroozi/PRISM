#!/usr/bin/env bash

# Fill these in with the local installation paths for the external tools that
# the upstream MaSIF preprocessing scripts expect.

export APBS_BIN=/path/to/apbs/bin/apbs
export MULTIVALUE_BIN=/path/to/apbs/share/apbs/tools/bin/multivalue
export PDB2PQR_BIN=/path/to/pdb2pqr/pdb2pqr
export REDUCE_HET_DICT=/path/to/reduce/reduce_wwPDB_het_dict.txt
export PYMESH_PATH=/path/to/PyMesh
export MSMS_BIN=/path/to/msms/msms

# Optional convenience additions:
export PATH=$PATH:/path/to/reduce
