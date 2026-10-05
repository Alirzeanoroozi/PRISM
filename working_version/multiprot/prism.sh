#!/bin/sh
#
#SBATCH -p cosbi
#SBATCH --ntasks=1

cd /cosbi/web/apps/prism/

python2 prism.py pair_list template_list ${1} > jobs/${1}/debug.log;
