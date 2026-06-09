#!/bin/bash
#SBATCH --job-name=pyscenic         # Job name
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --partition=batch_cpu                 # Partition name
#SBATCH --mem=50G                       # Total memory
#SBATCH --qos=1d
#SBATCH --cpus-per-task=5                # 1 CPU
#SBATCH --output=pyscenic.out
#SBATCH --error=pyscenic.err

# Load modules and activate conda environment
source ~/.bashrc
conda activate scenic_protocol

# Change to the data directory
cd /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/anndata_objects/tf_ko_panel_control_pyscenic

# Loop through all processed_input_*.loom files
for loom in processed_input_bin*.loom; do
    bin_name=$(basename "$loom" .loom)
    echo "Running pyscenic grn on $loom"
    
    pyscenic grn "$loom" allTFs_hg38.txt \
        -o "${bin_name}_adj.csv" \
        --num_workers ${SLURM_CPUS_PER_TASK}
done