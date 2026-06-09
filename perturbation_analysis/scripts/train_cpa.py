import os
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import cpa
import torch
import argparse

model_params = {
    "n_latent": 64,
    "recon_loss": "nb",
    "doser_type": "linear",
    "n_hidden_encoder": 128,
    "n_layers_encoder": 2,
    "n_hidden_decoder": 512,
    "n_layers_decoder": 2,
    "use_batch_norm_encoder": True,
    "use_layer_norm_encoder": False,
    "use_batch_norm_decoder": False,
    "use_layer_norm_decoder": True,
    "dropout_rate_encoder": 0.0,
    "dropout_rate_decoder": 0.1,
    "variational": False,
    "seed": 6977,
}

trainer_params = {
    "n_epochs_kl_warmup": None,
    "n_epochs_pretrain_ae": 30,
    "n_epochs_adv_warmup": 50,
    "n_epochs_mixup_warmup": 0,
    "mixup_alpha": 0.0,
    "adv_steps": None,
    "n_hidden_adv": 64,
    "n_layers_adv": 3,
    "use_batch_norm_adv": True,
    "use_layer_norm_adv": False,
    "dropout_rate_adv": 0.3,
    "reg_adv": 20.0,
    "pen_adv": 5.0,
    "lr": 0.0003,
    "wd": 4e-07,
    "adv_lr": 0.0003,
    "adv_wd": 4e-07,
    "adv_loss": "cce",
    "doser_lr": 0.0003,
    "doser_wd": 4e-07,
    "do_clip_grad": True,
    "gradient_clip_value": 1.0,
    "step_size_lr": 10,
}

def main():
    parser = argparse.ArgumentParser(description="Train a scCausalVI model")
    parser.add_argument(
        "--adata",
        type=str,
        required=True,
        help="Path to the input AnnData (.h5ad) file",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Directory to save the trained model",
    )
    parser.add_argument(
        "--max_epochs",
        type=int,
        default=500,
        help="Maximum number of training epochs (default: 500)",
    )
    parser.add_argument(
        "--condition_column",
        type=str,
        default="condition",
        help="obs column that defines control vs. target cells (default: 'condition')",
    )
    parser.add_argument(
        "--control_label",
        type=str,
        default="Control",
        help="Label in obs[condition_column] that defines background/control cells (default: 'Control')",
    )

    args = parser.parse_args()
    print("CUDA available:", torch.cuda.is_available())

    # Load data
    adata = sc.read(args.adata)
    
    # Add control column
    adata.obs["control"] = (adata.obs[args.condition_column] == args.control_label) * 1

    # Setup CPA
    cpa.CPA.setup_anndata(adata,
                      perturbation_key=args.condition_column,
                      control_group=args.control_label,
                      is_count_data=True,
                      max_comb_len=1,
                     )

    # Create model
    model = cpa.CPA(adata=adata,
                **model_params,
               )
    
    # Train model
    os.makedirs(args.outdir, exist_ok=True)
    
    use_gpu = torch.cuda.is_available()
    model.train(max_epochs=args.max_epochs,
            use_gpu=use_gpu,
            batch_size=512,
            plan_kwargs=trainer_params,
            early_stopping_patience=5,
            check_val_every_n_epoch=5,
            save_path=args.outdir,
           )
    print(f"Model saved to {args.outdir}")
    
    # Save the AnnData with latent representations
    latent_outputs = model.get_latent_representation(adata, batch_size=2048)

    print("Write AnnData with latent representations to disk")
    latent_outputs.write(os.path.join(args.outdir, "adata_with_latent.h5ad"))


if __name__ == "__main__":
    main()
