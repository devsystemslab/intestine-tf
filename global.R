## Define I/O ##
io <- list()
io$indir.facs <- "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/facs_data/20240312"
io$outdir.processed <- "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/processed"
io$outdir.atac <- "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac"
io$figdir <- "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/figures"
io$tabdir <- "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tables"
io$archr.project.dir <- "/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/atac/archr"

## Define Options ##
opts <- list()
opts$nCount_ATAC_max <- 100000 
opts$nCount_RNA_max <- 25000 
opts$nCount_ATAC_min <- 1000 
opts$nCount_RNA_min <- 1000 
opts$nFeature_RNA_min <- 1000
opts$nFeature_ATAC_min <- 1000
opts$nucleosome_signal_max <- 2 
opts$TSS_enrichment_min <- 1
opts$percent_mito_max <- 30

## Define Tools ##
tools <- list()
tools$macs2_path <- "/home/adaml9/scratch/miniforge3/envs/pando/bin/macs2"

## Define ProjectParams ##
project.params <- list()

project.params$ct.colors <- c(
  "cycling cells" = "#aeaf18",
  "GAU1+ cells" = "#dddddd",
  "Early EECs" = "#d3d19a",
  "K cells" = "#58B6D7",
  "I cells" = "#629BD2",
  "I/K cells" = "#73cdd1",
  "Late I/K cells" = "#73cdd1",
  "Early I/K cells" = "#E6F4F1",
  "Early X cells" = "#c6c6e0",
  "X cells" = "#85a7cf",
  "Late X cells" = "#85a7cf",
  "Early ECs" = "#f8bbaa",
  "ECs" = "#ee8865",
  "Late ECs" = "#ee8865",
  "PCSK5+ cells" = "#b6c425",
  "D cells" = "#6ac077",
  "SMOC2+ cells" = "#ABBE84",
  "L cells" = "#dbc8e0",
  "N cells" = "#cac8e0",
  "TPH1+ ECs" = "#e48865",
  "TAC1+ ECs" = "#ff8865",
  "Goblet cells" = "#FFA040",
  "Paneth cells I" = "#24c2df",
  "Paneth cells II" = "#1c97ad",
  "Stem cells" = "#ABBE84",
  "TA cells" = "#9AADA6",
  "Enterocytes" = "#15AD73"
)

project.params$species.colors <- c(
  "Human" = "#3498db",
  "Mouse" = "#e74c3c"
)


project.params$condition.bmp.colors <- c(
  "Control" = "#24c2df",
  "BMP" = "#15AD73"
)