**VALAR runs**

To do 

==Move the scripts and files per user to Cluster==

// Fix paths //

1. rsync -av --progress PATH_TO_prism-fiberdock-cli
  hdemirel22@172.20.240.205:/scratch/hdemirel22/hpc_run/prism-fiberdock-cli/




2. Transfer input files



 ==Set up environment once for each account==


3. ./setup_cluster.sh 


4. Run 

module load anaconda3/2025.06
conda create -n prism-fiberdock -c conda-forge python=2.7  

4. Fix sh script path

for file in fatma_*; do sbatch -J "$file" prism.sh; done

[For a smoke run]
sbatch -J "pairlist_deneme" prism.sh

[smoke run on login node]
conda activate prism-fiberdock
head -20 template_default > template_smoke # will not produce result probably
echo "1cew 2ghuD" > smoke_pairs
python2 prism.py smoke_pairs template_smoke smoketest
find jobs/smoketest -name '*.HB' -size 0     # must print nothing
