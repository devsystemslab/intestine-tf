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
library(Matrix)
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

####################################################################################################
## Helper: convert sparse matrix to long format (only non-zero entries)
####################################################################################################

sparse_to_long <- function(mat, row_var = "row", col_var = "col", value_var = "value") {
  nz <- summary(mat)

  tibble(
    !!row_var := rownames(mat)[nz$i],
    !!col_var := colnames(mat)[nz$j],
    !!value_var := nz$x
  ) %>%
    filter(.data[[value_var]] == 1)
}

####################################################################################################
## Load GRN object and extract peaks
####################################################################################################

grn_object <- readRDS(
  file = paste0(io$outdir.processed, "/pando_grn_object_eecs_custom.rds")
)
params <- Params(grn_object)
gene_annot <- Signac::Annotation(GetAssay(grn_object, params$peak_assay))
regions <- NetworkRegions(grn_object)
peak_to_gene_method <- "GREAT"

distances <- c(1000, 10000, 100000)
only_tss <- TRUE

for (d in distances) {

  upstream   <- d
  downstream <- d

  peaks_near_gene <- find_peaks_near_genes(
    peaks      = regions@ranges,
    method     = peak_to_gene_method,
    genes      = gene_annot,
    upstream   = upstream,
    downstream = downstream,
    only_tss   = only_tss
  )

  # Write sparse matrix for peaks near genes
  prefix <- paste0(
    io$tabdir,
    "/pando_peak_gene_links_",
    upstream, "bp"
  )

  writeMM(peaks_near_gene, paste0(prefix, ".mtx"))
  writeLines(rownames(peaks_near_gene), paste0(prefix, "_rows.txt"))
  writeLines(colnames(peaks_near_gene), paste0(prefix, "_cols.txt"))

  df_genes <- sparse_to_long(
    peaks_near_gene,
    row_var = "peak",
    col_var = "gene"
  )

  df_motifs <- sparse_to_long(
    grn_object@grn@regions@motifs@data,
    row_var = "peak",
    col_var = "motif"
  ) %>%
    left_join(
      grn_object@grn@regions@motifs@motif.names %>%
        enframe() %>%
        unnest() %>%
        dplyr::rename(motif = "name", motif.name = "value"),
      by = "motif"
    ) %>%
    dplyr::select(-motif, -value)

  merged <- df_genes %>%
    full_join(df_motifs, by = "peak") %>%
    mutate(
      upstream   = upstream,
      downstream = downstream,
      only_tss   = only_tss
    )

  outfile <- paste0(
    io$tabdir,
    "/pando_peak_gene_motif_links_",
    upstream, "bp.csv"
  )

  write.csv(merged, outfile, row.names = FALSE)
}

####################################################################################################
## Inspect links for LMX1B -> PAX4
####################################################################################################

# Sanity check. These regions should be bound by LMX1B motif
df_motifs %>% 
  dplyr::filter(peak %in% c("chr7-127809012-127809923",
                            "chr7-127864990-127866121",
                            "chr7-128017605-128018549",
                            "chr7-128123281-128124222")) %>% 
  dplyr::filter(motif.name=="LMX1B")

grn_object@grn@networks$glm_network@coefs %>%
  dplyr::filter(tf=="LMX1B") %>% 
  dplyr::filter(target=="PAX4")

####################################################################################################
## Merge peak–gene and peak–motif tables
####################################################################################################

merged <- read.csv(file = paste0(io$tabdir, "/pando_peak_gene_motif_links.csv"))

merged %>% 
  dplyr::filter(peak %in% c("chr7-127809012-127809923",
                            "chr7-127864990-127866121",
                            "chr7-128017605-128018549",
                            "chr7-128123281-128124222"))

