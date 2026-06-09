# Title: "Compute Gene Scores using ArchR"
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
library(tibble)
library(GenomicFeatures)
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
new_project_save_name = "EEC_multiome"
subscript = "human"

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 0) Load previously defined archr project & seurat object containing EEC multiome dataset

# Add default hg38 genome
addArchRGenome("hg38")
geneAnnotation <- getGeneAnnotation()
genomeAnnotation <- getGenomeAnnotation()

proj <- loadArchRProject(new_project_save_name)

markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "CellType",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

saveRDS(markersGS, paste0(output_dir, "markersGenesAll.summarizedExperiment.rds"))

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

mat <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")

# Convert the matrix to a sparse matrix
library(Matrix)

# Extract components
matSparse <- assay(mat)  # Sparse matrix of gene scores
genes <- mat@elementMetadata$name  # Gene names
cellMetaData <- colData(mat)  # Cell metadata

# Save objects as RDS files
saveRDS(matSparse, paste0(output_dir, "GeneScoreMatrixSparse.rds"))  # Sparse matrix
saveRDS(genes, paste0(output_dir, "GeneNames.rds"))                 # Gene names
saveRDS(cellMetaData, paste0(output_dir, "CellMetaData.rds"))       # Cell metadata

# Save matSparse as mtx
Matrix::writeMM(matSparse, file = paste0(output_dir, "GeneScoreMatrixSparse.mtx"))

# Save genes as tsv
write.table(genes, file = paste0(output_dir, "GeneNames.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save cell metadata as tsv
write.table(as.data.frame(cellMetaData), file = paste0(output_dir, "CellMetaData.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


