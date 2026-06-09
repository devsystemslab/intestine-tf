import argparse
import os
import numpy as np
import scanpy as sc
import scvi
import torch


def main():
    parser = argparse.ArgumentParser(description="Train a ContrastiveVI model with scvi-tools")
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

    # Setup ContrastiveVI
    scvi.external.ContrastiveVI.setup_anndata(adata, layer="counts")

    contrastive_vi_model = scvi.external.ContrastiveVI(
        adata,
        n_salient_latent=args.n_salient_latent,
        n_background_latent=args.n_background_latent,
        use_observed_lib_size=False,
    )

    # Define background and target indices
    background_indices = np.where(adata.obs[args.condition_column] == args.control_label)[0]
    target_indices = np.where(adata.obs[args.condition_column] != args.control_label)[0]

    print(f"Background cells: {len(background_indices)}, Target cells: {len(target_indices)}")

    # Train
    contrastive_vi_model.train(
        background_indices=background_indices,
        target_indices=target_indices,
        early_stopping=True,
        max_epochs=args.max_epochs,
    )

    # Save model
    os.makedirs(args.outdir, exist_ok=True)
    contrastive_vi_model.save(args.outdir, overwrite=True)
    print(f"Model saved to {args.outdir}")

    # Save the AnnData with latent representations
    adata.obsm["salient_rep"] = contrastive_vi_model.get_latent_representation(
        adata, representation_kind="salient"
    )
    adata.obsm["shared_rep"] = contrastive_vi_model.get_latent_representation(
        adata, representation_kind="background"
    )

    print("Write AnnData with latent representations to disk")
    adata.write(os.path.join(args.outdir, "adata_with_latent.h5ad"))


if __name__ == "__main__":
    main()
