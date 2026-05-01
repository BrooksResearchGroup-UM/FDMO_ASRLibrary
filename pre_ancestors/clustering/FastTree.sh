#!/bin/bash 
#SBATCH --job-name=FastTree
#SBATCH --time=10:00:00
#SBATCH -p gpuA5500
##SBATCH --output=/dev/null
#SBATCH --chdir=/home/cdchiang/for_azamh/vae
# Remember to change the directory above!!!

# Load modules
module load fasttree/fasttree

FastTree -lg -cat 20 -gamma /home/cdchiang/for_azamh/vae/data/PF01494_proc_output_5_20230920_2/PF01494_MSA_truncated.fasta > /home/cdchiang/for_azamh/vae/data/PF01494_proc_output_5_20230920_2/PF01494_MSA_truncated_tree_LG_G20.newick

exit
















