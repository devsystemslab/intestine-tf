# Title: "Run chromVAR using ArchR"
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

# Load project
proj <- loadArchRProject(new_project_save_name)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)
proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(ArchRProj = proj, outputDirectory = new_project_save_name, load = FALSE, overwrite = TRUE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

proj <- loadArchRProject(new_project_save_name)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

source("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/archr/utils.R")

markerMotifs <- getFeatures(proj, useMatrix = "MotifMatrix")
markerMotifs

seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "CellType")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

# Extract the chromVAR matrix
chromVAR_matrix <- as.matrix(seZ@assays@data$MotifMatrix)
rownames(chromVAR_matrix) <- proj@peakAnnotation[[1]]$motifSummary[rowData(seZ)$name,]$name

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

tf_list_df <- read.csv("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/external/data/cisbp/TF_Information_all_motifs.txt", sep="\t", header=TRUE)
tf_list_df <- tf_list_df %>% dplyr::filter(TF_Species == "Homo_sapiens")


# Vector of target TFs
tfs <- c(
  "PROX1", "ETV4", "HMGB3", "FOXN3", "NR2E3", "BHLHE40", "ATOH1", "FOXA1", "FOXA2",
  "PAX4", "SOX4", "ZNF800", "KLF4", "NEUROD2", "PDX1", "HHEX", "LMX1A", "SPDEF",
  "HNF4G", "NEUROG3", "TOX3", "ETS1", "RFX3", "ARID5B", "CXXC4", "ZNF326", "SPIB",
  "BAMBI", "DACH2", "THRB", "ZBTB18", "ZBTB7C", "TCF7", "ID1", "ID3", "POU2AF3",
  "ZNF608", "ZNF704", "MXD1", "TEAD1", "ISL1", "ZBTB46", "ZNF445", "HMX2", "ARNT2",
  "GATA4", "ETV1", "PRDM16", "POU2AF2", "MNX1", "HIF1A", "LCORL", "MEIS3", "ZNF711",
  "HES1", "HES4", "PERCC1", "NPAS2", "FOS", "GLIS3", "GTF2IRD1", "JAZF1",
  "KLF12", "SETBP1", "INSM1", "LMX1B", "RFX6", "HES2", "NEUROD1", "NKX22",
  "CSRNP3", "L3MBTL3", "RORA", "ARX", "PAX6", "RUNX1T1", "HDAC9", "MYT1L",
  "GFI1", "ZCCHC12", "EGR4", "ZKSCAN1", "NKX61"
)

print("Genes missing from ArchR:")
print(setdiff(tfs, rownames(chromVAR_matrix) %>% unique()))

print("Genes missing from CISBP:")
print(setdiff(tfs, tf_list_df$TF_Name %>% unique()))

# Subset chromVAR matrix to those motifs
matched_motifs <- intersect(tfs, rownames(chromVAR_matrix) %>% unique())
chromVAR_subset <- chromVAR_matrix[matched_motifs, ]
chromVAR_subset <- chromVAR_subset[, c("Stem Cells", "TA Cells", "Cycling Cells",
                                       "Enterocytes", "BEST4+ Enterocytes", 
                                       "Goblet Cells", "EEC Progenitors", 
                                       "EC Cells", "X Cells", "D Cells", "I/N Cells", "K Cells")]

group_by_max <- function(expr) {
  # expr: data.frame or matrix with genes as rows, groups as columns
  max_group <- apply(expr, 1, function(x) colnames(expr)[which.max(x)])
  ordered <- c()
  for (group in colnames(expr)) {
    ordered <- c(ordered, rownames(expr)[max_group == group])
  }
  return(ordered)
}

chromVAR_subset <- chromVAR_subset[group_by_max(chromVAR_subset), ]
# zscore normalize rows to make deviations comparable
chromVAR_subset <- t(scale(t(chromVAR_subset)))

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Get the 11 colors from RdBu
colors <- rev(brewer.pal(11, "RdBu"))
breakpoints <- seq(0, 1, length.out = 11)
col_fun <- colorRamp2(breakpoints, colors)

ht <- .ArchRHeatmap(t(chromVAR_subset), 
                    limits = c(-2, 2),
                    showRowDendrogram=FALSE, 
                    showColDendrogram=FALSE, 
                    clusterCols = FALSE,
                    clusterRows = FALSE,
                    labelCols=TRUE, 
                    labelRows=TRUE,
                    padding=5, 
                    color=col_fun)

outfile.dir <- paste(io$figdir, "Multiome", "DA", "ArchR", sep="/", collapse="")
if(!dir.exists(outfile.dir)){
  dir.create(outfile.dir, recursive=TRUE)
}
# Get the current date and time using lubridate
current_datetime <- format(now(), "%Y%m%d_%H%M%S")

# Construct the output file path with date and time
outfile.path <- file.path(outfile.dir, paste0("chromVAR_TF_KO_Panel_Heatmap_", current_datetime, ".pdf"))

pdf(outfile.path, width=20, height=5)
draw(ht)
dev.off()

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Convert chromVAR_matrix to long form
chromVAR_long <- chromVAR_matrix %>%
  as.data.frame() %>%
  rownames_to_column("motif") %>%
  pivot_longer(-motif, names_to = "cell_type", values_to = "score")

# Function to compute exclusivity for one cell type
compute_exclusivity <- function(target, df_long) {
  df_long %>%
    group_by(motif) %>%
    mutate(
      max_other = max(score[cell_type != target], na.rm = TRUE),
      exclusive_score = score[cell_type == target] - max_other
    ) %>%
    ungroup() %>%
    filter(cell_type == target) %>%
    arrange(desc(exclusive_score)) %>%
    slice_head(n = 10) %>%
    mutate(target_cell = target)
}

# Run for all cell types and combine
all_top_exclusive <- map_dfr(unique(chromVAR_long$cell_type),
                             compute_exclusivity,
                             df_long = chromVAR_long)

# Subset original chromVAR_matrix to top motifs per cell type
top_motifs <- unique(all_top_exclusive$motif)
chromVAR_matrix_top_n <- chromVAR_matrix[top_motifs, ]

# Optional: preview results
print(all_top_exclusive)
print(chromVAR_matrix_top_n)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Get the 11 colors from RdBu
colors <- rev(brewer.pal(11, "RdBu"))
breakpoints <- seq(0, 1, length.out = 11)
col_fun <- colorRamp2(breakpoints, colors)

ht <- .ArchRHeatmap(t(chromVAR_matrix_top_n), 
                    limits = c(-5, 5),
                    showRowDendrogram=TRUE, 
                    showColDendrogram=TRUE, 
                    clusterCols = TRUE,
                    clusterRows = TRUE,
                    labelCols=TRUE, 
                    labelRows=TRUE,
                    padding=5, 
                    color=col_fun)

outfile.dir <- paste(io$figdir, "Multiome", "DA", "ArchR", sep="/", collapse="")
if(!dir.exists(outfile.dir)){
  dir.create(outfile.dir, recursive=TRUE)
}
# Get the current date and time using lubridate
current_datetime <- format(now(), "%Y%m%d_%H%M%S")

# Construct the output file path with date and time
outfile.path <- file.path(outfile.dir, paste0("chromVAR_Heatmap_", current_datetime, ".pdf"))

pdf(outfile.path, width=20, height=5)
draw(ht)
dev.off()

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

markerMotifs <- getFeatures(proj, useMatrix = "MotifMatrix")

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- grep("NANOG", markerMotifs, value = TRUE)

name <- "NANOG_433"
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "MotifMatrix", 
    name = paste("z:", name), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

outfile.dir <- paste(io$figdir, "Multiome", "DA", "ArchR", sep="/", collapse="")
if(!dir.exists(outfile.dir)){
  dir.create(outfile.dir, recursive=TRUE)
}
# Get the current date and time using lubridate
current_datetime <- format(now(), "%Y%m%d_%H%M%S")

# Construct the output file path with date and time
outfile.path <- file.path(outfile.dir, paste0("ArchR_chromVAR_UMAP_", name, "_", current_datetime, ".pdf"))

ggsave(outfile.path, p, width=6, height=6)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################


corGSM_MM <- correlateMatrices(
    ArchRProj = proj,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

library(ggrepel)

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta)) +
  geom_point(aes(color = TFRegulator)) + 
  geom_text_repel(data=corGSM_MM[corGSM_MM$TFRegulator=="YES",] %>% head(10), aes(label=MotifMatrix_matchName)) + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

outfile.dir <- paste(io$figdir, "Multiome", "DA", "ArchR", sep="/", collapse="")
if(!dir.exists(outfile.dir)){
  dir.create(outfile.dir, recursive=TRUE)
}
# Get the current date and time using lubridate
current_datetime <- format(now(), "%Y%m%d_%H%M%S")

# Construct the output file path with date and time
outfile.path <- file.path(outfile.dir, paste0("ArchR_pos_tf_regulators_", name, "_", current_datetime, ".pdf"))

ggsave(outfile.path, p, width=6, height=6)
