#!/bin/bash
#SBATCH --job-name=06_chromVAR_Deviations         # Job name
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --partition=batch_cpu                 # Partition name
#SBATCH --qos=3h
#SBATCH --mem=90G                       # Total memory
#SBATCH --cpus-per-task=1                # 1 CPU
#SBATCH --output=06_chromVAR_Deviations.out
#SBATCH --error=06_chromVAR_Deviations.err

source ~/.bashrc
conda activate archr_to_signac

cd /projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/archr
Rscript ./06_chromVAR_Deviations.R