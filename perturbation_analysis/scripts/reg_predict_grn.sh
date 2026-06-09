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

# Run pyscenic command
pyscenic ctx adj.tsv \
    /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/external_data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/external_data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/external_data/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/external_data/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/external_data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname processed_input.loom \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 20

    