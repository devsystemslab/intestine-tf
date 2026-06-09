# Title: "Call Peaks in Human using ArchR"
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
library(purrr)
library(lubridate)
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
# 0) Load previously defined archr project & seurat object containing EEC multiome dataset

# Add default hg38 genome
addArchRGenome("hg38")
geneAnnotation <- getGeneAnnotation()
genomeAnnotation <- getGenomeAnnotation()

library(ArchR)
set.seed(1)
library(parallel)

pathToMacs2 <- tools$macs2_path

# Load project again
proj <- loadArchRProject(new_project_save_name)

# Get genome size
genomeSize <- sum(width(proj@genomeAnnotation$chromSizes))
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "CellType", force=TRUE)
proj <- addReproduciblePeakSet(
     ArchRProj = proj, groupBy = "CellType", force = TRUE, 
        pathToMacs2 = pathToMacs2, maxPeaks = 250000, genomeSize = genomeSize
    )

#Add Peak Matrix
proj <- addPeakMatrix(ArchRProj = proj, force = TRUE)
saveArchRProject(ArchRProj = proj, outputDirectory = new_project_save_name, load = FALSE, overwrite = TRUE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Load project
proj <- loadArchRProject(new_project_save_name)

# Export consensus peaks
peaks <- getPeakSet(proj)

# Export peaks as BED file
bed_file <- paste0(output_dir, "/consensus_peaks.bed")
rtracklayer::export.bed(peaks, con = bed_file)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Load project
proj <- loadArchRProject(new_project_save_name)

# Identify marker peaks using ArchR
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "CellType",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "binomial"
)

saveRDS(markersPeaks, paste0(output_dir, "markersPeaksAll.summarizedExperiment.binomial.rds"))
