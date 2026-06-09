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

# In addition to these features we also add all of ours TFs 
tfs <- c(
  "Prox1", "Etv4", "Hmgb3", "FOXN3", "NR2E3", "BHLHE40", "Atoh1", "Foxa1", "Foxa2",
  "Pax4", "Sox4", "Znf800", "Klf4", "Neurod2", "Pdx1", "Hhex", "Lmx1a",
  "Spdef", "Hnf4g", "Neurog3", "Tox3", "hEts1", "RFX3",  "ARID5B", "CXXC4",
  "ZNF326", "spib", "Bambi", "DACH2", "THRB", "ZBTB18", "ZBTB7C", "tcf7",
  "ID1", "ID3", "Pou2af3", "ZNF608", "ZNF704", "MXD1", "TEAD1", "Isl1", "ZBTB46",
  "ZNF445", "HMX2",  "Arnt2", "GATA4", "ETV1", "PRDM16",
  "Pou2af2", "Mnx1", "HIF1A", "LCORL", "MEIS3", "ZNF711", "Hes1", "Klf4", "HES4",
  "PERCC1", "NPAS2", "FOS", "GLIS3", "GTF2IRD1", "JAZF1", "KLF12", "SETBP1", 
  "Insm1", "Lmx1b", "Rfx6", "HES2", "Neurod1", "Mnx1", "Nkx2-2", "CSRNP3",
  "L3MBTL3", "RORA", "Arx", "Isl1", "Pax6", "Runx1t1", "Hdac9", "Myt1l", "Gfi1",
  "Zcchc12", "EGR4", "ZKSCAN1", "Nkx6-1", "Pou2af3", "Insm1", "Lmx1b", "Percc1"
)
tfs <- toupper(tfs) %>% unique()

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

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

# Rename POU2F2 in the dataset to POU2AF2
rownames(counts)[rownames(counts) == "POU2F2"] <- "POU2AF2"

# Rename POU2F3 in the dataset to POU2AF3
rownames(counts)[rownames(counts) == "POU2F3"] <- "POU2AF3"

# Rename ETS1 in the dataset to HETS1
rownames(counts)[rownames(counts) == "ETS1"] <- "HETS1"

# Subset to protein coding genes
protein_coding_genes <- read_tsv("/projects/site/pred/ihb-g-deco/PUBLIC_DATA/DB/protein_coding_genes.tsv") %>% 
                            dplyr::filter(!is.na(external_gene_name)) %>% 
                            pull(external_gene_name) 
print(paste0("Number of protein coding genes: ", length(protein_coding_genes)))

# Subset counts to protein coding genes
counts <- counts[rownames(counts) %in% protein_coding_genes, ]

# Remove genes with less than 3 cells with non-zero expression
genes.keep <- rowSums(counts > 0) >= 3
counts <- counts[genes.keep, ]
print(paste0("Number of genes after filtering: ", nrow(counts)))

# Remove genes with less than 10 total counts across all cells
genes.keep <- rowSums(counts) >= 10
counts <- counts[genes.keep, ]
print(paste0("Number of genes after filtering: ", nrow(counts)))

# Create a new RNA assay with the filtered counts matrix
data.seurat[["RNA"]] <- CreateAssayObject(counts = counts)

# We have to split the RNA assay by sample
data.seurat[["RNA"]] <- split(data.seurat[["RNA"]], f = data.seurat$sample)

# Normalize the RNA assay
data.seurat <- NormalizeData(data.seurat, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

# We now join the layers back together
data.seurat <- JoinLayers(data.seurat, assays = "RNA")

saveRDS(data.seurat, file = paste0(io$outdir.processed, "/parse.multiome.merged.cont3.eecs_for_pando.rds"))

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

for (nfeatures in c(2000, 3000, 4000, 5000)) {

  # Read in the dataset
  data.seurat <- readRDS(file = paste0(io$outdir.processed, "/parse.multiome.merged.cont3.eecs_for_pando.rds"))

  # Set the variable features for the RNA assay, which will be used for the GRN inference
  data.seurat <- FindVariableFeatures(data.seurat, assay = "RNA", selection.method = "vst", nfeatures = nfeatures)

  # Get the TFs that are not present in the dataset
  tfs_not_present <- tfs[!tfs %in% rownames(data.seurat)]
  print(paste0("Number of TFs not present in the dataset: ", length(tfs_not_present)))
  print("TFs not present in the dataset:")
  print(tfs_not_present)

  # Filter to only those TFs that are present in the dataset
  tfs_used <- tfs[tfs %in% rownames(data.seurat)]
  print(paste0("Number of TFs in the dataset: ", length(tfs_used)))

  # Now create the set of features
  var.features <- rownames(data.seurat)[!is.na(data.seurat@assays$RNA@meta.data$var.features)]
  var.features <- unique(c(var.features, tfs_used))

  print(paste0("Number of variable features: ", length(var.features)))

  # We also have to select the regions to which
  # we want to constrain the network
  regions <- phastConsElements20Mammals.UCSC.hg38

  main_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
  keep_peaks <- as.logical(seqnames(granges(data.seurat[["ATAC"]])) %in% main_chroms)

  annotation <- Annotation(data.seurat[["ATAC"]])
  annotation <- annotation[annotation@seqnames %in% main_chroms]

  counts <- data.seurat[["ATAC"]]@counts[keep_peaks, ]

  ranges <- granges(data.seurat[["ATAC"]])[keep_peaks]
  fragments <- data.seurat[["ATAC"]]@fragments

  data.seurat[["ATAC"]] <- CreateChromatinAssay(
    counts = counts,
    sep = c("-", "-"),
    ranges = ranges,
    annotation = annotation,
    fragments = fragments,
  )

  motif_annotations <- readRDS(file = "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr/motif_PFMatrixList.rds")

  motif_summary <- readRDS(file = "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr/motif_summary.rds")

  motif2tf_use <- motif_summary %>%
                      as.data.frame() %>%
                      dplyr::filter(name %in% var.features) %>% 
                      arrange(name) 
  motifs_use <- motif_annotations[rownames(motif2tf_use)]                 
  motif2tf_use <- motif2tf_use %>% rownames_to_column("motif") %>%
                      dplyr::select(motif, name) %>%
                      dplyr::rename(tf = name)

  grn_object <- initiate_grn(data.seurat, peak_assay = "ATAC", rna_assay = "RNA")

  # Scan candidate regions for TF binding motifs
  grn_object <- find_motifs(
      grn_object,
      pfm = motifs_use, 
      motif_tfs = motif2tf_use,
      genome = BSgenome.Hsapiens.UCSC.hg38
  )

  saveRDS(grn_object, file = paste0(io$outdir.processed, "/pando_grn_object_eecs_custom", nfeatures, ".rds"))
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

#grn_object <- readRDS(file = paste0(io$outdir.processed, "/pando_grn_object_eecs_custom.rds"))

#var.features <- grn_object@data@assays$RNA@meta.data$var.features
#var.features <- var.features[!is.na(var.features)]

#c("SST", "GHRL", "PYY", "PPY", "GCG", "NTS", "TAC1", "SCT", "MLN", "ARX") %in% var.features