#!/bin/bash
#SBATCH --partition=ai
#SBATCH --qos=ai
#SBATCH --account=ai
#SBATCH --job-name=mp-native-replay
#SBATCH --output=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-native-gate-replay/replay_%j.out
#SBATCH --error=/scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-native-gate-replay/replay_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=02:00:00

set -euo pipefail
source /opt/ohpc/admin/lmod/lmod/init/bash
module load rosetta/2022.42
cd "$SLURM_SUBMIT_DIR"
PYTHON=/home/rshadi25/.conda/envs/gtalign_env/bin/python
export PRISM_ROSETTA_PREPACK="${PRISM_ROSETTA_PREPACK:-/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static/docking_prepack_protocol.static.linuxgccrelease}"
export PRISM_ROSETTA_DOCK="${PRISM_ROSETTA_DOCK:-/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static/docking_protocol.static.linuxgccrelease}"
export PRISM_ROSETTA_DB="${PRISM_ROSETTA_DB:-/opt/ohpc/pub/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database/}"
"$PYTHON" benchmark/scripts/run_multiprot_native_gate_replay.py \
  /scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260722-pipeline-validation/multiprot_external \
  /scratch/rshadi25/GitHub/PRISM-prescript/tmp/agent/20260728-multiprot-native-gate-replay/results
