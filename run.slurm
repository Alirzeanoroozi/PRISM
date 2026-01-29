#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=kutem_gpu
#SBATCH --account=kutem
#SBATCH --gres=gpu:tesla_a100:1
#SBATCH --time=72:00:00
#SBATCH --mem=40G
#SBATCH --output=logs/comp.out
#SBATCH --error=logs/comp.out

# -p kutem_gpu -A kutem --qos=kutem --nodelist=rk02

# Initialise environment and modules
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/bin/activate bg
export LD_LIBRARY_PATH=${CONDA_BASE}/lib

nvidia-smi

module load rosetta/2022.42
python prism.py