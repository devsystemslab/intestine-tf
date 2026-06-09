# Intestine TF

This repository contains the analysis code accompanying our publication on "Arrayed cell state reporters and perturbation screens in human organoids map the regulatory landscape of intestinal epithelial fate". It brings together the main computational workflows used for the study, including comparative analyses, FACS-based measurements, VASA, multiome and ArchR processing, perturbation analyses, and reporter analyses.

## Repository layout

- `comparative_analysis/`: comparative analyses and plotting notebooks.
- `facs_analysis/`: FACS preprocessing and fluorescence ratio analyses.
- `multiome_analysis/`: multiome processing, ArchR workflows, assembly notebooks, and Pando-based regulatory network analyses.
- `perturbation_analysis/`: perturbation screens, annotation, mapping, and downstream statistical analyses.
- `reporter_analysis/`: reporter expression analyses.
- `vasa_analysis/`: VASA loading, trajectory analysis, and transcription factor selection.
- `global.R`: shared project paths and analysis parameters used by several R workflows.

## Notes on reproducibility

The repository is organized as a working analysis codebase rather than a packaged software project. Several scripts and notebooks rely on local file paths, cluster execution wrappers, and environment-specific software installations, so paths and runtime settings will need to be adapted before rerunning the workflows in a new environment.

## Data availability

Raw and processed data can be made available upon reasonable request.
