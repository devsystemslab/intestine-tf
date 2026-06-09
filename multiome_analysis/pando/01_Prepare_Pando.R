# Title: "Prepare for Pando"
# Author: "Lukas Adam"

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

## Load libraries
library(ArchR)
library(Seurat)
library(Signac)
library(dplyr)
library(purrr)
library(tidyr)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
`%notin%` <- Negate(`%in%`)

# Source global definitions
source("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global.R")
set.seed(0)

#Set/Create Working Directory to Folder
output_dir <- paste0(io$archr.project.dir, "/Multiome_v3")
if (!dir.exists(output_dir)) {dir.create(output_dir)}
setwd(output_dir)

#Set Threads to be used
addArchRThreads()

# Things to set for subseting projects. Also find and replace proj name for subsample with desired project name (or similar)
new_project_save_name = "EEC_only_multiome"
subscript = "human"

# Add default hg38 genome
addArchRGenome("hg38")
geneAnnotation <- getGeneAnnotation()
genomeAnnotation <- getGenomeAnnotation()

library(ArchR)
# Load project
proj <- loadArchRProject(new_project_save_name)

rna_mat <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
rna_counts <- assay(rna_mat)
rna_features <- rowData(rna_mat)
rownames(rna_counts) <- rna_features$name

peak_mat <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
atac_counts <- assay(peak_mat)
atac_features <- gr %>% as.data.frame() %>% unite(seqnames, start, end, sep="-") %>% pull()
rownames(atac_counts) <- atac_features

# Create RNA assay
rna_assay <- CreateAssayObject(counts = rna_counts)

# Create ATAC assay (use peaks as features)
atac_assay <- CreateAssayObject(counts = atac_counts)

# Initialize combined object
seurat_obj <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  meta.data = as.data.frame(proj@cellColData)
)

# Add ATAC assay
seurat_obj[["ATAC"]] <- atac_assay
DefaultAssay(seurat_obj) <- "RNA"

# Add embeddings
seurat_obj[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(proj@embeddings$UMAP$df),
  key = "UMAP_",
  assay = "RNA"
)

seurat_obj[["pca_harmony"]] <- CreateDimReducObject(
  embeddings = as.matrix(proj@embeddings$UMAP$df),
  key = "PCA_Harmony_",
  assay = "RNA"
)

seurat_obj[["umap_harmony"]] <- CreateDimReducObject(
  embeddings = as.matrix(proj@embeddings$UMAP$df),
  key = "UMAP_Harmony_",
  assay = "RNA"
)

saveRDS(seurat_obj, file = paste0(output_dir, "/seurat_object_for_pando.rds"))
