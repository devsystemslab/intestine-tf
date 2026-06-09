# Title: "Trajectory analysis using ArchR"
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
library(patchwork)
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
new_project_save_name = "EEC_only_multiome"
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
# Load project
proj <- loadArchRProject("EEC_multiome")

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Subset to EECs
cellNames <- proj@cellColData %>% 
    as.data.frame() %>%
    dplyr::filter(CellType %in% c("EEC Progenitors", "EC Cells", "D Cells", "I/N Cells", "X Cells", "K Cells")) %>%
    rownames()

proj <- subsetCells(ArchRProj = proj, cellNames = cellNames)

saveArchRProject(ArchRProj = proj, outputDirectory = new_project_save_name, load = FALSE, overwrite = TRUE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

proj <- loadArchRProject(new_project_save_name)

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

proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"   # or another batch variable
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)

saveArchRProject(ArchRProj = proj, outputDirectory = new_project_save_name, load = FALSE, overwrite = TRUE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

library(readr)

proj <- loadArchRProject(new_project_save_name)

# Load embeddings
pca_harmony <- read_tsv("/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr/X_pca_harmony.tsv") %>% as.data.frame() 
umap_harmony <- read_tsv("/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr/X_umap_harmony.tsv") %>% as.data.frame()

# Ensure the rownames are correct
rownames(pca_harmony) <- pca_harmony[,1]
pca_harmony <- pca_harmony[,-1]
colnames(pca_harmony) <- paste0("pca_harmony#Dim", colnames(pca_harmony))

rownames(umap_harmony) <- umap_harmony[,1]
umap_harmony <- umap_harmony[,-1]
colnames(umap_harmony) <- paste0("umap_harmony#Dim", colnames(umap_harmony))

pca_harmony <- pca_harmony %>% 
    rownames_to_column(var = "cellID") %>%
    separate(
        cellID,
        into = c("cell_id", "sample_id"),
        sep = "(?<=-\\d(?:_\\d)?)-"
    ) %>%
    mutate(cell_id = str_replace(cell_id, "-\\d(_\\d)?$", "-1")) %>%
    unite(newcellID, c("sample_id", "cell_id"), sep = "#", remove = FALSE) %>%
    column_to_rownames(var = "newcellID") %>% 
    dplyr::select(-c("sample_id", "cell_id"))

umap_harmony <- umap_harmony %>% 
    rownames_to_column(var = "cellID") %>%
    separate(
        cellID,
        into = c("cell_id", "sample_id"),
        sep = "(?<=-\\d(?:_\\d)?)-"
    ) %>%
    mutate(cell_id = str_replace(cell_id, "-\\d(_\\d)?$", "-1")) %>%
    unite(newcellID, c("sample_id", "cell_id"), sep = "#", remove = FALSE) %>%
    column_to_rownames(var = "newcellID") %>% 
    dplyr::select(-c("sample_id", "cell_id"))

# Match order of cells to ArchR project
pca_harmony <- pca_harmony[getCellNames(proj), , drop = FALSE]
umap_harmony <- umap_harmony[getCellNames(proj), , drop = FALSE]

proj@embeddings$pca_harmony <- SimpleList(df = pca_harmony, params = list())
proj@embeddings$umap_harmony <- SimpleList(df = umap_harmony, params = list())

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

library(Matrix)
library(SummarizedExperiment)

base_dir <- "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr/"

# Read in matrix
counts <- readMM(file.path(base_dir, "gene_exp_matrix.mtx"))

# Convert from dgTMatrix → dgCMatrix
counts <- as(counts, "dgCMatrix")

# Read barcodes and features
barcodes <- read.delim(file.path(base_dir, "gene_exp_barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)
features <- read.delim(file.path(base_dir, "gene_exp_features.tsv"), header = FALSE, stringsAsFactors = FALSE)

# Read metadata 
metadata <- read.delim(file.path(base_dir, "cell_metadata.tsv"), header = TRUE, stringsAsFactors = FALSE)
rownames(metadata) <- rownames(umap_harmony)
metadata$CellType <- as.factor(metadata$final_annotation)

# Assign row and column names
counts <- t(counts)
rownames(counts) <- features[[1]]
colnames(counts) <- rownames(umap_harmony)

# Drop duplicated gene symbols first
gene_uniq <- geneAnnotation$genes[!duplicated(geneAnnotation$genes$symbol)]

# Match by unique symbol
idx <- match(features$V1, gene_uniq$symbol)

# Keep only matched
matched_features <- features$V1[!is.na(idx)]
matched_annotation <- gene_uniq[idx[!is.na(idx)]]

counts <- counts[matched_features, ]

seRNA <- SummarizedExperiment(
  assays = list(counts = counts),
  colData = metadata,
  rowRanges = matched_annotation
)

proj <- addGeneExpressionMatrix(
    input = proj,
    seRNA = seRNA,
)

saveArchRProject(ArchRProj = proj, outputDirectory = new_project_save_name, load = FALSE, overwrite = TRUE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

proj <- loadArchRProject(new_project_save_name)

proj <- addTrajectory(
    ArchRProj = proj, 
    name = "XTraj", 
    groupBy = "CellType",
    trajectory = c("EEC Progenitors", "X Cells"), 
    embedding = "umap_harmony", 
    reducedDims = "pca_harmony",
    force = TRUE
)

proj <- addTrajectory(
    ArchRProj = proj, 
    name = "ECTraj", 
    groupBy = "CellType",
    trajectory = c("EEC Progenitors", "EC Cells"), 
    embedding = "umap_harmony", 
    reducedDims = "pca_harmony",
    force = TRUE
)

proj <- addTrajectory(
    ArchRProj = proj, 
    name = "KTraj", 
    groupBy = "CellType",
    trajectory = c("EEC Progenitors", "K Cells"), 
    embedding = "umap_harmony", 
    reducedDims = "pca_harmony",
    force = TRUE
)

proj <- addTrajectory(
    ArchRProj = proj, 
    name = "DTraj", 
    groupBy = "CellType",
    trajectory = c("EEC Progenitors", "D Cells"), 
    embedding = "umap_harmony", 
    reducedDims = "pca_harmony",
    force = TRUE
)

proj <- addTrajectory(
    ArchRProj = proj, 
    name = "ITraj", 
    groupBy = "CellType",
    trajectory = c("EEC Progenitors", "I/N Cells"), 
    embedding = "umap_harmony", 
    reducedDims = "pca_harmony",
    force = TRUE
)

proj <- addImputeWeights(proj)

saveArchRProject(ArchRProj = proj, outputDirectory = new_project_save_name, load = FALSE, overwrite = TRUE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

library(RColorBrewer)

RdBu_r_250 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(250)
PuOr_r_250 <- colorRampPalette(rev(brewer.pal(11, "PuOr")))(250)
RdPu_250 <- colorRampPalette(brewer.pal(11, "RdPu"))(250)
BuPu_250 <- colorRampPalette(brewer.pal(11, "BuPu"))(250)

proj <- loadArchRProject(new_project_save_name)

# Export motif accessibility values as matrix
motif_vals <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")
motif_mat <- assays(motif_vals)[["z"]]
write.csv(as.data.frame(as.matrix(motif_mat)), file = "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr/motif_accessibility_matrix.csv")

# Export motif PFMatrixList
motif_annotations <- proj@peakAnnotation[[1]]$motifs
saveRDS(motif_annotations, file = "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr/motif_PFMatrixList.rds")

motifSummary <- proj@peakAnnotation[[1]]$motifSummary
saveRDS(motifSummary, file = "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr/motif_summary.rds")



p2 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "MotifMatrix", 
    name = "z:RFX3_723", 
    embedding = "umap_harmony"
)

pdf("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/test.pdf")
print(p2)
dev.off()

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

library(TFBSTools)
library(ggseqlogo)
library(ggrastr)

motif_mat <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")
mat_motif <- assays(motif_mat)[["z"]]
imputed_mat_motif <- imputeMatrix(
    mat = mat_motif,    # or "GeneExpressionMatrix" depending on your project
    imputeWeights = getImputeWeights(proj)
)
rownames(imputed_mat_motif) <- rownames(motif_mat)
colnames(imputed_mat_motif) <- colnames(motif_mat)

# Export imputed motif accessibility values as csv
write.csv(as.data.frame(as.matrix(imputed_mat_motif)), file = "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr/Multiome_v3/EEC_multiome/h5ad_export/imputed_motif_accessibility_matrix.csv")

expr_mat <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
mat_expr <- assays(expr_mat)[["GeneExpressionMatrix"]]
imputed_mat_expr <- imputeMatrix(
    mat = mat_expr,    # or "GeneExpressionMatrix" depending on your project
    imputeWeights = getImputeWeights(proj)
)
rownames(imputed_mat_expr) <- rownames(expr_mat)
colnames(imputed_mat_expr) <- colnames(expr_mat)

tf_queries <- c("LMX1B", "RFX3", "PAX6", "ARX", "PROX1", "TEAD1", "SOX4", "PAX4", "ISL1", "PLAGL1", "FEV", "NPAS2")

for(tf_query in tf_queries){

    # Extract UMAP coordinates
    umap_df <- getEmbedding(ArchRProj = proj, embedding = "umap_harmony", returnDF = TRUE)
    umap_df$cellNames <- rownames(umap_df)

    # Extract motif accessibility values (chromVAR z-scores)
    motif_feature <- grep(paste0("z:", tf_query), getFeatures(proj, useMatrix = "MotifMatrix"), value = TRUE)[1]
    motif_vec <- imputed_mat_motif[gsub("z:", "", motif_feature), ]
    motif_logo <- proj@peakAnnotation[[1]]$motifs[[gsub("z:", "", motif_feature)]]
    # convert PWM (log-odds) to probability matrix (PFM)
    motif_pfm <- apply(motif_logo@profileMatrix, 2, function(col) motif_logo@bg * 2^col)
    
    gene_idx <- which(rowData(expr_mat)$name == tf_query)
    imputed_expr_vec <- imputed_mat_expr[gene_idx, ]
    imputed_expr_vec <- log2(imputed_expr_vec + 1)
    expr_vec <- mat_expr[gene_idx, ]
    expr_vec <- log2(expr_vec + 1)

    # Combine into one dataframe
    plot_df <- umap_df %>%
    mutate(
        imputed_motif_accessibility = motif_vec[cellNames],
        gene_expression = expr_vec[cellNames],
        imputed_gene_expression = imputed_expr_vec[cellNames],
    )

    p1 <- ggplot(plot_df %>% arrange(imputed_motif_accessibility), aes(`umap_harmony#Dim0`, `umap_harmony#Dim1`, color = imputed_motif_accessibility)) +
    geom_point_rast(size = 0.3) +
    scale_color_gradientn(
        colors = BuPu_250,
        oob = scales::squish,
        limits=c(min(plot_df$imputed_motif_accessibility),max(plot_df$imputed_motif_accessibility))
    ) + 
    ggtitle(paste(tf_query, "imputed motif accessibility")) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5))

    p2 <- ggplot(plot_df %>% arrange(gene_expression), aes(`umap_harmony#Dim0`, `umap_harmony#Dim1`, color = gene_expression)) +
    geom_point_rast(size = 0.3) +
    scale_color_gradientn(
        colors = RdPu_250,
        limits=c(min(plot_df$gene_expression),max(plot_df$gene_expression))
    ) +
    ggtitle(paste(tf_query, "expression")) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5))

    p3 <- ggplot(plot_df %>% arrange(imputed_gene_expression), aes(`umap_harmony#Dim0`, `umap_harmony#Dim1`, color = imputed_gene_expression)) +
    geom_point_rast(size = 0.3) +
    scale_color_gradientn(
        colors = RdPu_250,
        limits=c(min(plot_df$imputed_gene_expression),max(plot_df$imputed_gene_expression))
    ) +
    ggtitle(paste(tf_query, "imputed expression")) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5))

    p4 <- ggseqlogo(motif_pfm) + 
        ggtitle(paste(tf_query, "motif logo")) + 
        theme(plot.title = element_text(hjust = 0.5))

    p <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4, widths = c(1,1,1,1))
    outfile.dir <- paste(io$figdir, "Multiome", "DA", "ArchR", sep="/", collapse="")
    if(!dir.exists(outfile.dir)){
    dir.create(outfile.dir, recursive=TRUE)
    }
    # Get the current date and time using lubridate
    current_datetime <- format(now(), "%Y%m%d_%H%M%S")  
    # Construct the output file path with date and time
    outfile.path <- file.path(outfile.dir, paste0("ArchR_umap_motif_accessibility_", tf_query, "_", current_datetime, ".pdf"))
    ggsave(outfile.path, plot = p, width = 20, height = 5)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

labelMarkers <- tfs <- c(
  "PROX1", "ETV4", "HMGB3", "FOXN3", "NR2E3", "BHLHE40", "ATOH1", "FOXA1", "FOXA2",
  "PAX4", "SOX4", "ZNF800", "KLF4", "NEUROD2", "PDX1", "HHEX", "LMX1A", "SPDEF",
  "HNF4G", "NEUROG3", "TOX3", "ETS1", "RFX3", "ARID5B", "CXXC4", "ZNF326", "SPIB",
  "BAMBI", "DACH2", "THRB", "ZBTB18", "ZBTB7C", "TCF7", "ID1", "ID3", "POU2AF3",
  "ZNF608", "ZNF704", "MXD1", "TEAD1", "ISL1", "ZBTB46", "ZNF445", "HMX2", "ARNT2",
  "GATA4", "ETV1", "PRDM16", "POU2AF2", "MNX1", "HIF1A", "LCORL", "MEIS3", "ZNF711",
  "HES1", "HES4", "PERCC1", "NPAS2", "FOS", "GLIS3", "GTF2IRD1", "JAZF1",
  "KLF12", "SETBP1", "INSM1", "LMX1B", "RFX6", "HES2", "NEUROD1", "NKX22",
  "CSRNP3", "L3MBTL3", "RORA", "ARX", "PAX6", "RUNX1T1", "HDAC9", "MYT1L",
  "GFI1", "ZCCHC12", "EGR4", "ZKSCAN1", "NKX6-1"
)
useMatrices <- c("GeneExpressionMatrix", "PeakMatrix", "MotifMatrix")
trajNames <- c("XTraj", "ECTraj", "KTraj", "DTraj", "ITraj")

for(trajName in trajNames){
    trajGEM <- getTrajectory(ArchRProj = proj, name = trajName, useMatrix = "GeneExpressionMatrix", log2Norm = TRUE)
    trajMM <- getTrajectory(ArchRProj = proj, name = trajName, useMatrix = "MotifMatrix", log2Norm = FALSE)

    corGSM_MM <- correlateTrajectories(trajGEM, trajMM, corCutOff = 0, varCutOff1 = 0, varCutOff2 = 0)

    df <- as.data.frame(corGSM_MM[["correlatedMappings"]])
    df <- df %>% dplyr::filter(grepl("z:",name2)) %>% dplyr::filter(!grepl("-AS", name1))
    df <- df %>% 
        dplyr::filter(abs(Correlation) > 0.05) %>%
        dplyr::filter(VarAssay1 > 0.1 & VarAssay2 > 0.1) 

    trajGEM2 <- trajGEM[df$name1, ]
    trajMM2 <- trajMM[df$name2, ]

    trajCombined <- trajGEM2
    assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGEM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

    combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
    rowOrder <- match(rownames(combinedMat), rownames(trajGEM2))

    ht1 <- plotTrajectoryHeatmap(trajCombined, pal = PuOr_r_250, varCutOff = 0, rowOrder = rowOrder)
    ht2 <- plotTrajectoryHeatmap(trajGEM2, pal = RdPu_250, varCutOff = 0, rowOrder = rowOrder, labelMarkers=labelMarkers)
    ht3 <- plotTrajectoryHeatmap(trajMM2, pal = RdBu_r_250, varCutOff = 0, rowOrder = rowOrder)

    outfile.dir <- paste(io$figdir, "Multiome", "DA", "ArchR", sep="/", collapse="")
    if(!dir.exists(outfile.dir)){
    dir.create(outfile.dir, recursive=TRUE)
    }
    # Get the current date and time using lubridate
    current_datetime <- format(now(), "%Y%m%d_%H%M%S")

    # Construct the output file path with date and time
    outfile.path <- file.path(outfile.dir, paste0("ArchR_", trajName, "_integrated_traj_", current_datetime, ".pdf"))

    pdf(outfile.path, width=14, height=8)
    ComplexHeatmap::draw(ht1 + ht2 + ht3)
    dev.off()

    pdf(paste0("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/", trajName, "_traj.pdf"),
            width = 14, height = 8)
    ComplexHeatmap::draw(ht1 + ht2 + ht3)
    dev.off()
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################


for(trajName in trajNames){
    for(useMatrix in useMatrices){
        if(useMatrix=="PeakMatrix"){
            traj <- getTrajectory(ArchRProj = proj, name = trajName, useMatrix = useMatrix, log2Norm = TRUE)

            p <- plotTrajectoryHeatmap(traj, pal = PuOr_r_250)
        } else {
            traj <- getTrajectory(ArchRProj = proj, name = trajName, useMatrix = useMatrix, log2Norm = FALSE)

            p <- plotTrajectoryHeatmap(traj, pal = RdBu_r_250)
        }

        outfile.dir <- paste(io$figdir, "Multiome", "DA", "ArchR", sep="/", collapse="")
        if(!dir.exists(outfile.dir)){
        dir.create(outfile.dir, recursive=TRUE)
        }
        # Get the current date and time using lubridate
        current_datetime <- format(now(), "%Y%m%d_%H%M%S")

        # Construct the output file path with date and time
        outfile.path <- file.path(outfile.dir, paste0("ArchR_", trajName, "_", useMatrix, "_", current_datetime, ".pdf"))

        pdf(outfile.path, width=5, height=8)
        print(p)
        dev.off()
        
    }
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

traj  <- getTrajectory(ArchRProj = proj, name = "KTraj", useMatrix = "MotifMatrix", log2Norm = FALSE)

p1 <- plotTrajectoryHeatmap(traj, pal = paletteContinuous(set = "solarExtra"))


outfile.dir <- paste(io$figdir, "Multiome", "DA", "ArchR", sep="/", collapse="")
if(!dir.exists(outfile.dir)){
  dir.create(outfile.dir, recursive=TRUE)
}
# Get the current date and time using lubridate
current_datetime <- format(now(), "%Y%m%d_%H%M%S")

# Construct the output file path with date and time
outfile.path <- file.path(outfile.dir, paste0("ArchR_Traj_X_cells_", current_datetime, ".pdf"))

pdf(outfile.path, width=5, height=8)
print(p1)
dev.off()

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Load proj
proj <- loadArchRProject(new_project_save_name)

# Export proj as bigwg
bw <- getGroupBW(proj, groupBy = "CellType", tileSize=10)

# Get fragments from object
getFragmentsFromProject(proj)