# Title: "Export ArchR to h5ad"
# Author: "Lukas Adam"

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

## Load libraries
library(SeuratObject)
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

# Add default hg38 genome
addArchRGenome("hg38")
geneAnnotation <- getGeneAnnotation()
genomeAnnotation <- getGenomeAnnotation()

library(ArchR)
# Load project
proj <- loadArchRProject(new_project_save_name)

# Get peak matrix
peakMat <- getMatrixFromProject(proj, "PeakMatrix")

# Get project output directory
output.dir <- proj@projectMetadata$outputDirectory

# Create new directory for h5ad export
h5ad.dir <- file.path(output.dir, "h5ad_export")
if (!dir.exists(h5ad.dir)) {dir.create(h5ad.dir)}

# Write matrix in mtx format
Matrix::writeMM(assays(peakMat)$PeakMatrix, file = file.path(h5ad.dir, "peak_matrix.mtx"))

# Save cell metadata
cell_metadata <- as.data.frame(proj@cellColData)
cell_metadata$cell_id <- rownames(cell_metadata)
write.table(cell_metadata, file = file.path(h5ad.dir, "cell_metadata.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

# Write features
peak_ranges <- rowRanges(peakMat)
peak_names <- paste0(seqnames(peak_ranges), ":", start(peak_ranges), "-", end(peak_ranges))
features_df <- data.frame(peak_names)
write.table(features_df, file = file.path(h5ad.dir, "features.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Load motif matches
motif_mat <- readRDS(proj@peakAnnotation[[1]]$Matches)
motif_mat <- assays(motif_mat)[[1]]
rownames(motif_mat) <- peak_names
colnames(motif_mat) <- proj@peakAnnotation[[1]]$motifSummary[colnames(motif_mat),]$name
# write sparse matrix
Matrix::writeMM(motif_mat, file = file.path(h5ad.dir, "motif_matrix.mtx"))

# Save motif names
motif_names <- data.frame(motif_name = colnames(motif_mat))
write.table(motif_names, file = file.path(h5ad.dir, "motif_names.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

motifMat <- getMatrixFromProject(proj, "MotifMatrix")
# Write matrix in mtx format
Matrix::writeMM(assays(motifMat)$z, file = file.path(h5ad.dir, "chromvar_matrix.mtx"))

# Save motif names
motif_names <- data.frame(motif_name = rownames(motifMat))
write.table(motif_names, file = file.path(h5ad.dir, "chromvar_motif_names.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Export motif position weight matrices
motif_pwms <- proj@peakAnnotation[[1]]$motifs
motif_pwms <- motif_pwms[names(motif_pwms) %in% colnames(motif_mat)]

library(TFBSTools)

# compute motif similarity matrix
sim <- TFBSTools::PWMsimilarity(motif_pwms, type="Pearson")




# Get TSS for each gene
library("biomaRt")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
tss <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", values=proj@geneAnnotation$genes$symbol %>% unname(),
             mart = ensembl)

tss %>% dplyr::select(external_gene_name, chromosome_name, transcription_start_site) %>%
  dplyr::rename(GeneSymbol = external_gene_name,
                GeneChr = chromosome_name,
                GenePos = transcription_start_site) %>%
  write.table(file = file.path(h5ad.dir, "gene_tss.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")


