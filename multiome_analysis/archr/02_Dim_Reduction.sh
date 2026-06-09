#!/bin/bash
#SBATCH --job-name=Dim_Reduction_02         # Job name
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --partition=batch_cpu                 # Partition name
#SBATCH --qos=1d
#SBATCH --mem=250G                       # Total memory
#SBATCH --cpus-per-task=1                # 1 CPU
#SBATCH --output=Dim_Reduction_02.out
#SBATCH --error=Dim_Reduction_02.err

source ~/.bashrc
conda activate archr_to_signac

cd /projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/archr
Rscript ./02_Dim_Reduction.R

