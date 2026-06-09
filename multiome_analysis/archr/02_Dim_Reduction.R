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

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

proj <- loadArchRProject(new_project_save_name)

# Load cell multiome dataset
data.seurat.multiome <- readRDS(paste0(io$outdir.processed, "/parse.multiome.merged.cont3.rds"))
DefaultAssay(data.seurat.multiome) <- "RNA"


data.seurat.multiome@meta.data <- data.seurat.multiome@meta.data %>%
  rownames_to_column(var = "cellID") %>%
  separate(
    cellID,
    into = c("cell_id", "sample_id"),
    sep = "(?<=-\\d(?:_\\d)?)-"
  ) %>%
  mutate(cell_id = str_replace(cell_id, "-\\d(_\\d)?$", "-1")) %>%
  unite(newcellID, c("sample_id", "cell_id"), sep = "#", remove = FALSE) %>%
  column_to_rownames(var = "newcellID")


cellNames <- intersect(rownames(data.seurat.multiome@meta.data), rownames(proj@cellColData))

proj <- subsetCells(ArchRProj = proj, cellNames = cellNames)

# Define cell annotations
cell_annotations <- data.seurat.multiome@meta.data %>% 
                        rownames_to_column(var="cellID") %>% 
                        dplyr::select(cellID, final_annotation) %>% 
                        dplyr::rename(c("cellType" = "final_annotation")) %>% 
                        dplyr::filter(cellID %in% cellNames)

# Add annotations
proj <- addCellColData(
  ArchRProj = proj,
  data = as.character(cell_annotations$cellType),
  cells = cell_annotations$cellID,
  name = "CellType",
  force = TRUE
)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 4, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.1, 0.2, 0.4), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 15000, 
    dimsToUse = 1:30,
    force = TRUE
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)

p <- plotEmbedding(proj, name = "CellType")

ggsave(
    filename = paste0("ArchR_UMAP_CellType_", subscript, ".pdf"), 
    plot = p, 
    device = "pdf", 
    path = output_dir, 
    width = 6, 
    height = 5
)

saveArchRProject(ArchRProj = proj, outputDirectory = new_project_save_name, load = FALSE, overwrite = TRUE)


