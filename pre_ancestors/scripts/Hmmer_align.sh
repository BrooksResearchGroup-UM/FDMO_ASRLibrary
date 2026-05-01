#!/bin/bash 
#SBATCH --job-name=Hmmer_align
#SBATCH --time=02:00:00
#SBATCH -p gpu
#SBATCH --chdir=/home/cdchiang/for_azamh/vae/data
# Remember to change the directory above!!!

# Load modules

module load hmmer/3.3.2

# Run
# seed file: PF01494_seed.hmm
# input example: PF01494_input_3.fasta
# output example: hmmer_test.fasta
hmmalign --outformat afa PF01494_seed.hmm PF01494_input_3_20230920_2.fasta > PF01494_MSA_4_20230920_2.fasta 

exit
















