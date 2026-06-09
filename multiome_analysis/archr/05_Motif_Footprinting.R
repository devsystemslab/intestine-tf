# Title: "Perform Dim Reduction using ArchR"
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
library(AnnotationForge)
library(GenomicFeatures)
library(org.Hsapiens.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
`%notin%` <- Negate(`%in%`)

# Source global definitions
source("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global.R")
# Source colors
source("/projects/site/pred/ihb-g-deco/USERS/adaml9/phd_frameworks/articulate/plotting/color_palettes.R")

set.seed(0)

#Set/Create Working Directory to Folder
output_dir <- paste0(io$archr.project.dir, "/Multiome_v1")
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

# Load project
proj <- loadArchRProject(new_project_save_name)

motifPositions <- getPositions(proj)

motifs <- c("NEUROG3", "RFX6", "PERCC1", "HHEX", "ONECUT3", "ATOH1", "HES6", 
            "ASCL1", "PROX1", "^GFI1", "PAX6", "SOX4", "ST18", "HES4", 
            "^ARX")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "CellType"
)

outfile.dir <- paste(io$figdir, "Multiome", "DA", "ArchR", sep="/", collapse="")
if(!dir.exists(outfile.dir)){
  dir.create(outfile.dir, recursive=TRUE)
}
# Get the current date and time using lubridate
current_datetime <- format(now(), "%Y%m%d_%H%M%S")

# Construct the output file path with date and time
outfile.path <- file.path(outfile.dir, paste0("ArchR_Motif_Footprints_", current_datetime, ".pdf"))

pdf(outfile.path, width=10, height=10)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
dev.off()