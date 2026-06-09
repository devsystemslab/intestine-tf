# Title: "Initialize object using ArchR"
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
new_project_save_name = "EEC_multiome"
subscript = "human"

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# 0) Load data

data.seurat.multiome <- readRDS(paste0(io$outdir.processed, "/parse.multiome.merged.cont3.rds"))
seurat_obj <- data.seurat.multiome
DefaultAssay(seurat_obj) <- "ATAC"

# Extract the fragment file path from the Signac object
fragment_files <- lapply(seurat_obj@assays$ATAC@fragments, function(x) x@path)

# Extract the sample names
sample.names <- map_chr(fragment_files, function(x) gsub(".tsv.gz", "", basename(x)))

# Add default hg38 genome
addArchRGenome("hg38")
geneAnnotation <- getGeneAnnotation()
genomeAnnotation <- getGenomeAnnotation()

# Define Arrow files based on fragment data
ArrowFiles <- createArrowFiles(
  inputFiles = as.vector(unlist(fragment_files)),
  sampleNames = sample.names,  # Adjust name as needed
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  excludeChr = c("chrY")
)

ArrowFiles <- list.files(".", pattern=".arrow")

# Create ArchR project with Arrow files and metadata
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = new_project_save_name,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  copyArrows = FALSE
)
saveArchRProject(ArchRProj = proj, outputDirectory = new_project_save_name, load = FALSE, overwrite = FALSE)