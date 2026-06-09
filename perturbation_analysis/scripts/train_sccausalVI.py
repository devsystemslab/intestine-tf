import os
import scanpy as sc
import scvi
import anndata as ad
import numpy as np
import pandas as pd
from scCausalVI import scCausalVIModel
import torch
import argparse
import jax

torch.set_float32_matmul_precision('medium')

print(jax.__version__)
print(jax.devices())


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
        "--n_salient_latent",
        type=int,
        default=10,
        help="Number of salient latent dimensions (default: 10)",
    )
    parser.add_argument(
        "--n_background_latent",
        type=int,
        default=10,
        help="Number of background latent dimensions (default: 10)",
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

    scvi.settings.seed = 0
    print("Last run with scvi-tools version:", scvi.__version__)
    print("CUDA available:", torch.cuda.is_available())

    # Load data
    adata = sc.read(args.adata)

    # Setup scCausalVI
    scCausalVIModel.setup_anndata(adata, condition_key=args.condition_column, layer="counts")

    # Map condition labels to integers
    condition2int = (
        adata.obs.groupby(args.condition_column, observed=False)["_scvi_condition"]
        .first()
        .to_dict()
    )

    model = scCausalVIModel(
        adata,
        condition2int=condition2int,
        control=args.control_label,
        n_background_latent=args.n_background_latent,
        n_te_latent=args.n_salient_latent,
        n_layers=2,
        n_hidden=128,
        use_mmd=True,
        mmd_weight=10,
        norm_weight=0.2,
    )

    # Get group indices list
    conditions = adata.obs[args.condition_column].unique().tolist()
    group_indices_list = [
        np.where(adata.obs[args.condition_column] == group)[0] for group in conditions
    ]

    # Train model
    model.train(group_indices_list, max_epochs=args.max_epochs, use_gpu=True)

    # Save model
    os.makedirs(args.outdir, exist_ok=True)
    model.save(args.outdir, overwrite=True)
    print(f"Model saved to {args.outdir}")

    # Save the AnnData with latent representations
    adata.obsm["latent_bg"], adata.obsm["latent_t"] = model.get_latent_representation()

    print("Write AnnData with latent representations to disk")
    adata.write(os.path.join(args.outdir, "adata_with_latent.h5ad"))


if __name__ == "__main__":
    main()
