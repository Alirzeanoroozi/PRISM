#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=cosbi
#SBATCH --account=cosbi
#SBATCH --time=72:00:00
#SBATCH --mem=40G
#SBATCH --output=logs/comp.out
#SBATCH --error=logs/comp.out

# Initialise environment and modules
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/bin/activate bg
export LD_LIBRARY_PATH=${CONDA_BASE}/lib

# zip -r templates.zip templates/

# module load rosetta/2022.42
echo "Running PRISM"
# python -u prism.py --generate_templates True
python prism.py
