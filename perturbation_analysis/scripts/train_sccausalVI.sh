#!/bin/bash
#SBATCH --job-name=train_sccausalVI       # Job name
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --partition=batch_gpu                 # Partition name
#SBATCH --gres=gpu:5                     # Request GPU
#SBATCH --qos=1d                     
#SBATCH --mem=200G                       # Total memory
#SBATCH --cpus-per-task=5                # 5 CPU
#SBATCH --output=train_sccausalVI.out     # Standard output
#SBATCH --error=train_sccausalVI.err      # Standard error

source ~/.bashrc
conda activate sccausalvi

python train_sccausalVI.py \
  --adata /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/anndata_objects/tf_ko_panel_atlas.h5ad \
  --outdir /projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/models/scCausalVI \
  --n_salient_latent 10 \
  --n_background_latent 10 \
  --max_epochs 50 \
  --condition_column condition \
  --control_label Control