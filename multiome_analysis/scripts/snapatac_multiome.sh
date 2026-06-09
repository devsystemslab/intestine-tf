#!/bin/bash
#SBATCH --job-name=SnapATACMultiome         # Job name
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --partition=batch_cpu                 # Partition name
#SBATCH --qos=1d
#SBATCH --mem=250G                       # Total memory
#SBATCH --cpus-per-task=1                # 1 CPU
#SBATCH --output=/home/adaml9/logs/SnapATACMultiome.out
#SBATCH --error=/home/adaml9/logs/SnapATACMultiome.err

# Load environment
source /home/adaml9/.bashrc
conda activate snapatac2

python snapatac_multiome.py