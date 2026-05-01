#!/usr/bin/bash
#SBATCH --job-name=esmfold
#SBATCH --partition=gpuA5500
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=50GB
#SBATCH --array=170
#SBATCH --output=/home/azamh/esm_screen/slurm/esmfold/esmfold-%a.out

# Activate conda enviornment
source ~/.bashrc
conda activate esmfold

# Path to fold.py script
fold_py=/home/azamh/esm_screen/esm/scripts/fold.py

# Path to model directory
model_dir=/home/azamh/esm_screen/esm_pdb

# Array id
id=$SLURM_ARRAY_TASK_ID
fasta_path=/home/azamh/esm_screen/fasta/$id.fasta


echo "Array id: $id"
echo "Running $fold_py"
echo "Model directory: $model_dir"
echo "Fasta path: $fasta_path"

# Run fold_py
python $fold_py -i $fasta_path -o $model_dir
