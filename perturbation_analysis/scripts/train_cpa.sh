#!/bin/bash
#SBATCH --job-name=train_cpa       # Job name
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --partition=batch_gpu                 # Partition name
#SBATCH --gres=gpu:5                     # Request 5 GPU
#SBATCH --qos=1d                     
#SBATCH --mem=200G                       # Total memory
#SBATCH --cpus-per-task=5                # 5 CPU
#SBATCH --output=train_cpa.out     # Standard output
#SBATCH --error=train_cpa.err      # Standard error

source ~/.bashrc
conda activate cpa

python train_cpa.py \
  --adata /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/anndata_objects/tf_ko_panel_atlas.h5ad \
  --outdir /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/models/cpa \
  --max_epochs 100 \
  --condition_column condition \
  --control_label Control