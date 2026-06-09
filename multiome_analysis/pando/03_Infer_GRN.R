# Title: "Infer GRN using Pando"
# Author: "Lukas Adam"

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Memory and compatibility options
options(
  future.globals.maxSize = 3e+09,
  Seurat.object.assay.version = "v3"
)
library(reticulate)
reticulate::use_python("/pmount/projects/site/pred/ihb-g-deco/USERS/adaml9/miniforge3/envs/pando/bin/python", required = TRUE)

## Load libraries
library(Seurat)
library(Signac)
library(dplyr)
library(purrr)
library(tidyverse)
library(Pando)
library(presto)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
library(ggraph)
library(ggrepel)
library(foreach)
library(doParallel)
`%notin%` <- Negate(`%in%`)

source("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global.R")

# Set color palette
ct_colors <- c(
  # Proliferative axis
  "Stem Cells"       = "#1f77b4",  # deep blue
  "TA Cells"         = "#4c91c6",  # medium blue
  "Cycling Cells"    = "#7aaed6",  # light blue
  
  # Absorptive lineage
  "Enterocytes"      = "#2ca02c",  # green
  "BEST4+ Enterocytes" = "#1abc9c",# teal
  
  # Secretory lineage
  "Goblet Cells"     = "#ff7f0e",  # orange
  
  # EEC progenitors
  "EEC Progenitors"  = "#9467bd",  # purple
  
  # EEC subtypes
  "EC Cells"         = "#d62728",  # red
  "D Cells"          = "#17becf",  # cyan
  "X Cells"          = "#e377c2",  # magenta
  "I/N Cells"        = "#7f7f7f",  # grey
  "K Cells"          = "#bcbd22"   # olive
)

# Set figdir 
figdir <- paste0(io$archr.project.dir, "/Multiome_v3/figures/pando")
# Create figure directory recursively
if (!dir.exists(figdir)) {dir.create(figdir, recursive = TRUE)}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

for(n_features in c(2000, 3000, 4000, 5000)) {
  for(peak_to_gene_method in c("GREAT")) {
    for(method in c("glm")) {
      
      input_file <- paste0(io$outdir.processed, "/pando_grn_object_eecs_custom", n_features, ".rds")
      output_file <- paste0(io$outdir.processed, "/pando_grn_object_eecs_fit_custom_", tolower(peak_to_gene_method), "_", method, "_", n_features, ".rds")

      grn_object <- readRDS(file = input_file)

      # Infer gene regulatory network
      grn_object <- infer_grn(grn_object, peak_to_gene_method = peak_to_gene_method, parallel = T, method = method)
      
      saveRDS(grn_object, file = output_file)
      
    }
  }
}