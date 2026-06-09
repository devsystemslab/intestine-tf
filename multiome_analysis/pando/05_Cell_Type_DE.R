# Title: "Identify cell type-specific DE genes"
# Author: "Lukas Adam"

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Memory and compatibility options
options(
  future.globals.maxSize = 3e+09,
  Seurat.object.assay.version = "v3"
)
## Load libraries
library(Seurat)
library(Signac)
library(dplyr)
library(purrr)
library(ggplot2)
library(parallel)
library(SCpubr)
`%notin%` <- Negate(`%in%`)

source("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global.R")

# Read in the dataset
data.seurat <- readRDS(paste0(io$outdir.processed, "/parse.multiome.merged.cont3.rds"))

# Subset to only eecs
cells.keep <- data.seurat@meta.data %>% 
                dplyr::filter(final_annotation %in% c("EEC Progenitors", "EC Cells", "X Cells", "K Cells", "D Cells", "I/N Cells")) %>%
                rownames()

data.seurat <- subset(data.seurat, cells = cells.keep)
data.seurat$final_annotation %>% unique()

# First we have to set the RNA assay as DefaultAssay
DefaultAssay(data.seurat) <- "RNA"

# Get the counts matrix
counts <- GetAssayData(data.seurat, assay = "RNA", slot = "counts")

# Hot fix for PERCC1 gene name, since we used 2020 genome version for mapping
rownames(counts)[rownames(counts) == "AL032819.3"] <- "PERCC1"

# Create the RNA assay again with fixed gene names
data.seurat[["RNA"]] <- CreateAssayObject(counts = counts)

# Normalize the RNA assay
data.seurat <- NormalizeData(data.seurat, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

# Set the identities to the final annotation
Idents(data.seurat) <- data.seurat$final_annotation

# Identify DE genes for all cell types
all.markers <- FindAllMarkers(object = data.seurat, assay = "RNA", slot = "data")

# Plot DE genes as heatmap
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

p <- VlnPlot(data.seurat, features = unique(top10$gene))

ggsave("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/test.pdf", p,
       width = 20, height = 20)