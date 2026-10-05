#!/bin/bash
#SBATCH --job-name=prism-fiberdock
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --qos=users
#SBATCH --account=users
##SBATCH --exclude=ai[01-14]
#SBATCH --time=7-00:00:00
#SBATCH --output=out-%x.%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fcankara20@ku.edu.tr

# --- environment ---------------------------------------------------------
# FiberDock/MultiProt build: no Rosetta / no TMalign needed, just Python 2.
# Reusing the conda env that already provides python2 (swap if you have a
# dedicated py2 env for this build).
eval "$(conda shell.bash hook)"
conda init bash
conda activate tmalignRosetta

# --- run -----------------------------------------------------------------
# EDIT this to wherever you uploaded prism-fiberdock-cli on the cluster:
cd /scratch/users/fcankara20/hpc_run/prism-fiberdock-cli

# Same pattern as tmalignscript.sh: the job name doubles as the pair_list
# filename and the jobId; template_default is the template list.
python2 prism.py $SLURM_JOB_NAME template_default $SLURM_JOB_NAME
