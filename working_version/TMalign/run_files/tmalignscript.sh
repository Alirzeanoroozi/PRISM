#!/bin/bash
#SBATCH --job-name=joblist_001
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --qos=users
#SBATCH --account=users
##SBATCH --exclude=ai[01-14]
#SBATCH --time=7-00:00:00
#SBATCH --output=out-%x.%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rshadi25@ku.edu.tr
eval "$(conda shell.bash hook)"
conda activate tmalignRosetta
module load rosetta/2022.42
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$script_dir/.." && pwd)
cd "$project_root"
python2 prism.py "$SLURM_JOB_NAME" "$script_dir/template_default" "$SLURM_JOB_NAME"
