#!/bin/bash
#SBATCH --job-name=02_Prepare_GRN         # Job name
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --partition=batch_cpu                 # Partition name
#SBATCH --qos=3h
#SBATCH --mem=100G                       # Total memory
#SBATCH --cpus-per-task=1                # 1 CPU
#SBATCH --output=02_Prepare_GRN.out
#SBATCH --error=02_Prepare_GRN.err

source ~/.bashrc
conda activate pando

cd /projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando
Rscript ./02_Prepare_GRN.R