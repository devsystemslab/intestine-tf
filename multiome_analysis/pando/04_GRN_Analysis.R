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
library(ggrepel)
library(igraph)
library(ggraph)
library(tidygraph)

set.seed(10)

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
  "Goblet Cells"     = "#f57f20",  # orange
  
  # EEC progenitors
  "EEC Progenitors"  = "#d3d19a",  # purple
  
  # EEC subtypes
  "EC Cells"         = "#ee8764",  
  "D Cells"          = "#6bc077",  # cyan
  "X Cells"          = "#85a7ce",  # magenta
  "I/N Cells"        = "#639cd2",  # grey
  "K Cells"          = "#59b7d8"   # olive
)

# Set figdir 
figdir <- paste0(io$archr.project.dir, "/Multiome_v3/figures/pando")
# Create figure directory recursively
if (!dir.exists(figdir)) {dir.create(figdir, recursive = TRUE)}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

source("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/utils.R")

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Loop through all GRN objects, find modules, and save results

for(n_features in c(2000, 3000, 4000, 5000)) {
  for(peak_to_gene_method in c("GREAT")) {
    for(method in c("glm")) {
      print(paste0("Processing ", peak_to_gene_method, " ", method))
      
      # Load input file
      input_file <- paste0(io$outdir.processed, "/pando_grn_object_eecs_fit_custom_", tolower(peak_to_gene_method), "_", method, "_", n_features, ".rds")
    
      # Load GRN
      grn_object <- readRDS(file = input_file)

      # Load modules
      grn_object@grn@networks[[paste0(method, "_network")]] <- find_modules(
          grn_object@grn@networks[[paste0(method, "_network")]], 
          p_thresh = 0.1,
          nvar_thresh = 2, 
          min_genes_per_module = 1, 
          rsq_thresh = 0.01
      )
      # Save TF-target-peak interactions
      df <- grn_object@grn@networks[[paste0(method, "_network")]]@coefs 
      df %>% 
        write.table(
          file = paste0(io$tabdir, "/pando_grn_eecs_tf_target_peak_interactions_", tolower(peak_to_gene_method), "_", method, "_", n_features, ".tsv"), 
          sep = "\t", quote = FALSE, row.names = FALSE
        )
      

      # Save tf modules
      df <- grn_object@grn@networks[[paste0(method, "_network")]]@modules@meta
      print(paste0("Number of TF modules: ", length(unique(df$tf))))
      print(paste0("Number of TF-target pairs: ", nrow(df)))
      
      df %>% write.table(
                  paste0(io$tabdir, "/pando_grn_eecs_tf_modules_", tolower(peak_to_gene_method), "_", method, "_", n_features, ".tsv"), 
                  sep = "\t", quote = FALSE, row.names = FALSE)

      # Save GRN object with modules
      saveRDS(grn_object, file = paste0(io$outdir.processed, "/pando_grn_object_eecs_fit_custom_", tolower(peak_to_gene_method), "_", method, "_", n_features, "_with_modules.rds"))
    } 
  }
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Load GRN object with highest number of TF modules and TF-target pairs
n_features <- 3000
peak_to_gene_method <- "GREAT"
method_base <- "glm"
method <- paste0(method_base, "_network")

grn_object <- readRDS(file = paste0(io$outdir.processed, "/pando_grn_object_eecs_fit_custom_", tolower(peak_to_gene_method), "_", method_base, "_", n_features, "_with_modules.rds"))
# Get TFs inside grn_object
grn_tfs <- unique(grn_object@grn@networks[[method]]@coefs$tf) 
# Get features inside grn_object
grn_features <- names(grn_object@grn@networks[[method]]@features)

# Plot distribution of estimate 
p <- grn_object@grn@networks[[method]]@coefs %>% 
    ggplot(aes(x=estimate)) + 
    geom_histogram(bins=50) + 
    theme_classic() + 
    labs(title = "Distribution of effect sizes for TF-target-peak interactions")

ggsave(filename = paste0(figdir, "/pando_grn_eecs_estimate_distribution_", tolower(peak_to_gene_method), "_", method_base, ".pdf"), plot = p, width = 5, height = 4)


# Plot distribution of p-values
p <- grn_object@grn@networks[[method]]@coefs %>% 
    ggplot(aes(x=pval)) + 
    geom_histogram(bins=50) + 
    theme_classic() + 
    labs(title = "Distribution of p-values for TF-target-peak interactions")

ggsave(filename = paste0(figdir, "/pando_grn_eecs_pval_distribution_", tolower(peak_to_gene_method), "_", method_base, ".pdf"), plot = p, width = 5, height = 4)


# Plot scatter plot of p-value vs. estimate
p <- grn_object@grn@networks[[method]]@coefs %>% 
    ggplot(aes(x=estimate, y=-log10(pval))) + 
    geom_point(alpha=0.5) + 
    theme_classic() + 
    labs(title = "TF-target-peak interactions: Effect size vs. significance")
  
ggsave(filename = paste0(figdir, "/pando_grn_eecs_pval_vs_estimate_", tolower(peak_to_gene_method), "_", method_base, ".pdf"), plot = p, width = 5, height = 4)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Save TF-target-peak interactions
grn_object@grn@networks[[method]]@coefs %>% 
  write.table(
    file = paste0(io$tabdir, "/pando_grn_eecs_tf_target_peak_interactions.tsv"), 
    sep = "\t", quote = FALSE, row.names = FALSE
  )

# Save tf modules
write.table(grn_object@grn@networks[[method]]@modules@meta, 
            paste0(io$tabdir, "/pando_grn_eecs_tf_modules.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

p <- plot_gof(grn_object, point_size=3)

ggsave(filename = paste0(figdir, "/pando_grn_eecs_gof_plot.pdf"), plot = p, width = 5, height = 4)

p <- plot_module_metrics(grn_object)

ggsave(filename = paste0(figdir, "/pando_grn_eecs_module_metrics_plot.pdf"), plot = p, width = 7, height = 5)

# Check individual modules
grn_object@grn@networks[[method]]@coefs %>% 
    dplyr::filter(tf == "LMX1B")

grn_object@grn@networks[[method]]@modules@meta %>% 
    dplyr::filter(tf == "LMX1B") 

"LMX1B" %in% c(grn_object@grn@networks[[method]]@features %>% names())

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

source("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/utils.R")

# Define which TFs to label
label_tfs <- c("NEUROG3", "LMX1A", "FEV", "LMX1B", "PROX1", "ARX", 
               "RFX3", "HHEX", "PAX4", "ISL1", "PAX6", "GATA4", "TEAD1")

# Instead of plotting gene x gene network, we plot the TFxTF network
gene_graph <- get_tf_network_graph(grn_object, 
    graph_name='umap_graph',
    features = c(grn_tfs, grn_features)
)

gene_graph <- gene_graph %>%
    tidygraph::activate(nodes) %>%
    mutate(tf = name %in% grn_tfs)

gene_graph <- gene_graph %>%
        tidygraph::activate(edges) %>%
        mutate(from_node=tidygraph::.N()$name[from], to_node=tidygraph::.N()$name[to]) %>%
        tidygraph::activate(nodes) %>%
        mutate(centrality=tidygraph::centrality_pagerank())

graph.de <- presto::wilcoxauc(grn_object@data,
          group_by = "final_annotation", 
          seurat_assay = "RNA",
          assay = "data") %>% 
          group_by(feature) %>%
          slice_max(logFC, n=1, with_ties=FALSE) %>%
          ungroup()

gene_graph <- gene_graph %>%
                tidygraph::activate(nodes) %>%
                left_join(graph.de, by = c("name" = "feature")) %>%
                mutate(log10_padj=-log10(padj))

p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
p <- p + geom_edge_diagonal(aes(color=corr), width=0.4, alpha=0.8)
p <- p + geom_node_point(aes(fill=group, size=centrality), shape=21, stroke=0.1) 
p <- p + theme_void() 
p <- p + geom_label_repel(
            data=gene_graph %>% 
                    tidygraph::activate("nodes") %>% 
                    as.data.frame %>% 
                    dplyr::filter(name %in% label_tfs),
            aes(x=UMAP_1, y=UMAP_2, label=name),
            size=2, box.padding = 2, max.overlaps = Inf
        )
p <- p + scale_fill_manual(values = ct_colors)
p <- p + scale_size(range = c(1, 8))
p <- p + scale_edge_colour_distiller(
    palette = "RdBu",
    limits = c(-1, 1),
    direction = -1   # reverses the palette (so it's RdBu_r)
  ) +
      guides(fill = guide_legend(override.aes = list(size=5)))

ggsave(filename = paste0(figdir, "/pando_grn_eecs_tf_network_graph.pdf"), plot = p, width = 5, height = 5)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Define which TFs to label
label_tfs <- c("NEUROG3", "LMX1A", "FEV", "LMX1B", "PROX1", "ARX", 
               "RFX3", "HHEX", "PAX4", "ISL1", "PAX6", "GATA4", "TEAD1")

# Instead of plotting gene x gene network, we plot the TFxTF network
grn_object <- Pando::get_network_graph(grn_object, 
    graph_name='umap_graph',
    features = c(grn_tfs, grn_features)
)

gene_graph <- grn_object@grn@networks[[method]]@graphs$umap_graph %>%
    tidygraph::activate(nodes) %>%
    mutate(tf = name %in% grn_tfs)

gene_graph <- gene_graph %>%
        tidygraph::activate(edges) %>%
        mutate(from_node=tidygraph::.N()$name[from], to_node=tidygraph::.N()$name[to]) %>%
        tidygraph::activate(nodes) %>%
        mutate(centrality=tidygraph::centrality_pagerank())

graph.de <- presto::wilcoxauc(grn_object@data,
          group_by = "final_annotation", 
          seurat_assay = "RNA",
          assay = "data") %>% 
          group_by(feature) %>%
          slice_max(logFC, n=1, with_ties=FALSE) %>%
          ungroup()

gene_graph <- gene_graph %>%
                tidygraph::activate(nodes) %>%
                left_join(graph.de, by = c("name" = "feature")) %>%
                mutate(log10_padj=-log10(padj))

p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
p <- p + geom_edge_diagonal(aes(color=corr), width=0.4, alpha=0.8)
p <- p + geom_node_point(aes(fill=group, size=centrality), shape=21, stroke=0.1) 
p <- p + theme_void() 
p <- p + geom_label_repel(
            data=gene_graph %>% 
                    tidygraph::activate("nodes") %>% 
                    as.data.frame %>% 
                    dplyr::filter(name %in% label_tfs),
            aes(x=UMAP_1, y=UMAP_2, label=name),
            size=2, box.padding = 2, max.overlaps = Inf
        )
p <- p + scale_fill_manual(values = ct_colors)
p <- p + scale_size(range = c(1, 8))
p <- p + scale_edge_colour_distiller(
    palette = "RdBu",
    limits = c(-1, 1),
    direction = -1   # reverses the palette (so it's RdBu_r)
  ) +
      guides(fill = guide_legend(override.aes = list(size=5)))

ggsave(filename = paste0(figdir, "/pando_grn_eecs_network_graph.pdf"), plot = p, width = 5, height = 5)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

grn_object <- get_network_graph(
    grn_object, 
    graph_name = 'full_graph', 
    umap_method = 'none',
    features = c(grn_tfs, grn_features)
)

tf_queries <- c("INSM1") #c("LMX1A", "PAX4", "INSM1", "RFX3", "RFX6", "ISL1", "PAX6", "TEAD1", "GATA4", "NEUROG3", "ARX", "SOX4")
for (tf_query in tf_queries) {

  # Try catch to handle cases where the TF is not present in the graph
  tryCatch({
    print(paste0("Processing TF: ", tf_query))
    # Get subgraph of TF and its target genes
    grn_object <- get_tf_network(grn_object, tf=tf_query, graph='full_graph', keep_all_edges = TRUE)

  }, error = function(e) {
    message(paste0("Error processing TF: ", tf_query, " - ", e$message))
    next
  })

  # Extract subgraph for the queried TF
  subgraph <- NetworkGraph(grn_object, graph='tf_graphs')[[tf_query]]

  # Get distances
  dists <- distances(subgraph, v = tf_query)
  keep_nodes <- V(subgraph)[dists[1, ] < 3]
  
  # Subset to the induced graph
  sub_g <- induced_subgraph(subgraph, vids = keep_nodes)

  # Convert igraph to tidygraph for ggraph plotting
  sub_g <- as_tbl_graph(sub_g)

  # Label nodes with gene names, but only for TFs (not target genes) to avoid overcrowding the plot
  sub_g <- sub_g %>%
    activate(nodes) %>%
    mutate(label = ifelse(name %in% grn_object@grn@networks[[method]]@coefs$tf, name, ""))

  p <- ggraph(sub_g, layout = "tree", root = tf_query, circular = TRUE) +
      geom_edge_diagonal(aes(colour=dir), alpha = 0.1) +
      scale_edge_colour_distiller(
        palette = "RdBu",
        direction = -1   # reverses the palette (so it's RdBu_r)
      ) +
      geom_node_point(size = 1, alpha = 0.1) +
      geom_node_text(aes(label = label), repel = TRUE, size = 3) +
      theme_graph()

  ggsave(filename = paste0(figdir, "/pando_grn_eecs_tf_subgraph_", tf_query, ".pdf"), 
         plot = p, width = 10, height = 10, device = cairo_pdf)
}

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################


# Inspect interactions between RFX6 and RFX3
grn_object@grn@networks[[method]]@coefs %>% 
    dplyr::filter((tf == "RFX6" & target == "RFX3") | (tf == "RFX3" & target == "RFX6")) %>%
    print()

grn_object@grn@networks[[method]]@modules@meta %>% 
    dplyr::filter((tf == "RFX6" & target == "RFX3") | (tf == "RFX3" & target == "RFX6")) %>%
    print()






############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Modify base GRN object to include only edges with
# TFs and target genes that were perturbed in our screen

# Load perturbation data
library(dplyr)
library(readr)
library(tidygraph)

# Load + keep only the annotations of interest
tf_ko_de_results <- read_tsv(
  file.path("/projects/site/pred/ihb-intestine-evo/lukas_area/eec_fate_project/tf_ko_screen/panel/tables",
            "tf_ko_panel_contrastiveVI_de_results_filtered.tsv")
) %>%
  rename(annotation = group) %>%
  filter(annotation %in% c("EC Cells", "D Cells", "X Cells", "K Cells", "I/N Cells",
                           "Early EEC Progenitors", "Late EEC Progenitors"))

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Build perturbation TF-TF graph
perturb_edges <- tf_ko_de_results %>%
  filter(significant != "NS") %>%
  transmute(
    perturbed_tf = toupper(test_condition),
    target_gene  = toupper(feature),
    effect       = coef
  ) %>%
  group_by(perturbed_tf, target_gene) %>%
  summarise(effect_median = median(effect, na.rm = TRUE), .groups = "drop") %>%
  mutate(effect_clipped = pmax(pmin(effect_median, 2), -2))

# TF list from perturbation screen
perturbed_tfs <- unique(perturb_edges$perturbed_tf)

# keep TF→TF edges only
perturb_tf_edges <- perturb_edges %>%
  filter(perturbed_tf %in% perturbed_tfs,
         target_gene %in% perturbed_tfs)

# build tidygraph (TF perturbation network)
perturb_tf_graph <- as_tbl_graph(
  perturb_tf_edges %>% dplyr::select(from = perturbed_tf, to = target_gene, weight = effect_clipped),
  directed = TRUE
)

perturb_node_df <- perturb_tf_graph %>%
  activate(nodes) %>%
  mutate(node_id = row_number()) %>%
  as_tibble() %>%
  dplyr::select(node_id, gene = name)

tf_perturb_score <- perturb_tf_graph %>%
  activate(edges) %>%
  as_tibble() %>%
  group_by(to) %>%
  summarise(tf_effect = median(weight, na.rm = TRUE)) %>%
  rename(node_id = to) %>%
  left_join(perturb_node_df, by = "node_id") %>%
  dplyr::select(tf = gene, tf_effect)

# Build multiome TF-TF graph
multiome_grn <- grn_object@grn@networks$glm_network@graphs$full_graph %>%
  activate(nodes) %>%
  mutate(is_tf = name %in% unique(grn_object@grn@networks$glm_network@coefs$tf))

multiome_tf_graph <- multiome_grn %>%
  # keep only TF → TF edges
  activate(edges) %>%
  filter(.N()$is_tf[from] & .N()$is_tf[to]) %>%
  # compute degree in this filtered graph
  activate(nodes) %>%
  mutate(tf_degree = centrality_degree()) %>%
  # keep TF nodes with at least 1 TF neighbor
  filter(is_tf & tf_degree > 1)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# TF master regulator rankings
library(tidygraph)
library(dplyr)
library(graph)

# Work on node table
tf_rank_table <- multiome_tf_graph %>%
  activate(nodes) %>%
  mutate(
    out_degree        = centrality_degree(mode = "out"),
    in_degree         = centrality_degree(mode = "in"),
    out_strength      = centrality_degree(mode = "out", weights = estimate),
    in_strength       = centrality_degree(mode = "in",  weights = estimate),
    pagerank_directed = centrality_pagerank(directed = TRUE, weights = abs(estimate)),
    betweenness       = centrality_betweenness(directed = TRUE, weights = abs(estimate)),
    closeness         = centrality_closeness()
  ) %>%
  as_tibble() %>%
  rename(tf = name)

tf_rank_table <- tf_rank_table %>%
  inner_join(tf_perturb_score, by = "tf")

rank_norm <- function(x) {
  r <- rank(-x, ties.method = "average")
  r / max(r)
}

tf_rank_table <- tf_rank_table %>%
  mutate(
    r_out_strength   = rank_norm(out_strength),
    r_pagerank       = rank_norm(pagerank_directed),
    r_tf_effect      = rank_norm(tf_effect),

    # weighted combo; bias slightly to perturbation evidence
    master_regulator_score = tf_effect
  ) %>%
  arrange(desc(master_regulator_score))

tf_rank_table %>%
  select(tf, tf_effect, out_strength, pagerank_directed, master_regulator_score) %>%
  arrange(desc(master_regulator_score)) %>%
  print(n = 50)


############################################################################################################################
# Perform PageRank mass sweep analysis

library(Matrix)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidygraph)
library(igraph)
library(ggraph)


### helpers -------------------------------------------------------------

ppr_out_from_S <- function(G, S_names, alpha = 0.15) {
  S_names <- intersect(S_names, V(G)$name)
  pers <- rep(0, vcount(G))
  pers[match(S_names, V(G)$name)] <- 1 / length(S_names)
  pr <- page_rank(G, damping = 1 - alpha, personalized = pers)$vector
  pr # returns full vector, we'll sum over all nodes later if needed
}

ppr_out_mass <- function(G, S_names, alpha = 0.15) {
  S_names <- intersect(S_names, V(G)$name)
  pers <- rep(0, vcount(G))
  pers[match(S_names, V(G)$name)] <- 1 / length(S_names)
  pr <- page_rank(G, damping = 1 - alpha, personalized = pers)$vector
  sum(pr)
}

coverage_k <- function(G, S_names, k = 3) {
  S_idx <- match(S_names, V(G)$name)
  dist_to_S <- distances(G, to = S_idx, mode = "all")
  covered <- apply(dist_to_S, 1, function(x) any(x <= k))
  mean(covered)
}

avg_dist_to_S <- function(G, S_names) {
  S_idx <- match(S_names, V(G)$name)
  dist_to_S <- distances(G, to = S_idx, mode = "all")
  dmin <- apply(dist_to_S, 1, min)
  mean(dmin[is.finite(dmin)])
}

max_dist_to_S <- function(G, S_names) {
  S_idx <- match(S_names, V(G)$name)
  dist_to_S <- distances(G, to = S_idx, mode = "all")
  dmin <- apply(dist_to_S, 1, min)
  max(dmin[is.finite(dmin)])
}

conductance_global <- function(G, S_names) {
  S_names <- intersect(S_names, V(G)$name)
  S_idx <- match(S_names, V(G)$name)
  cut_edges <- sum(sapply(E(G), function(e) {
    vs <- ends(G, e)
    (vs[1] %in% S_names) != (vs[2] %in% S_names)
  }))
  volS <- sum(degree(G)[S_idx])
  cut_edges / volS
}

### main panel ----------------------------------------------------------

global_extrapolation_panel <- function(G, S_names, alpha = 0.15) {
  S_names <- intersect(S_names, V(G)$name)

  tibble(
    # seed coverage of graph
    coverage_k2 = coverage_k(G, S_names, k = 2),
    coverage_k3 = coverage_k(G, S_names, k = 3),

    # how close whole graph is to S
    mean_dist = avg_dist_to_S(G, S_names),
    max_dist = max_dist_to_S(G, S_names),

    # outward signal spreading
    ppr_out = sum(ppr_out_from_S(G, S_names, alpha)),

    # S "opens" to graph (higher = better for extrapolation)
    conductance = conductance_global(G, S_names)
  )
}


g <- as.igraph(multiome_tf_graph)
S_names <- intersect(tf_perturb_score$tf, V(g)$name)

panel_global <- global_extrapolation_panel(g, S_names)

panel_global


# convert to tidygraph for plotting
tg <- as_tbl_graph(multiome_tf_graph) %>%
  mutate(
    in_S = name %in% S_names
  )

g <- ggraph(tg, layout = "fr") +
  geom_edge_link(alpha = 0.05) +
  geom_node_point(aes(color = in_S), size = 2) +
  geom_node_text(
    aes(label = name),
    size = 3, box.padding = 1, max.overlaps = Inf
  ) +
  scale_color_manual(values = c("grey70", "red")) +
  theme_void() +
  labs(title = "Perturbed TFs (red) inside Multiome TF Network")

ggsave("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/test6.pdf")


############################################################################################################################
# Perform message passing / diffusion of perturbation signals through the GRN

library(Matrix)
library(tidygraph)
library(dplyr)
library(igraph)
library(ggraph)

# Propagate effect score through GRN edges
A <- multiome_tf_graph %>%
  activate(edges) %>%
  mutate(estimate = abs(estimate)) %>%
  as_adjacency_matrix(attr = "estimate", sparse = TRUE)

node_names <- multiome_tf_graph %>%
  activate(nodes) %>%
  pull(name)

# initialize zero for all TFs in GRN
x <- rep(0, length(node_names))
names(x) <- node_names

# fill values for TFs that have perturb signals
overlap <- intersect(tf_perturb_score$tf, node_names)
x[overlap] <- tf_perturb_score$tf_effect[match(overlap, tf_perturb_score$tf)]

deg_out <- Matrix::rowSums(abs(A))                  # abs handles signed edges
L <- Diagonal(x = deg_out) - A                      # directed Laplacian
I <- Diagonal(n = length(x))

alpha <- 0.2  # diffusion strength (smoother if 0.05–0.3 range)

y <- solve(I + alpha * L, x)                        # diffused vector
names(y) <- names(x)

multiome_tf_graph <- multiome_tf_graph %>%
  activate(nodes) %>%
  mutate(
    perturb_seed = x[name],
    perturb_diffused = y[name]
  )

signal.df <- multiome_tf_graph %>%
  activate(nodes) %>%
  as.data.frame() %>% 
  dplyr::filter(perturb_seed == 0) %>% 
  dplyr::filter(perturb_diffused != 0) %>%
  mutate(diff = perturb_diffused - perturb_seed) %>%
  arrange(perturb_diffused)

p <- ggplot(signal.df, aes(x=diff)) +
  geom_histogram(bins=50) +
  theme_classic() +
  labs(x = "Diffused TF Perturbation", y = "Number of TFs")
ggsave("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/test5.pdf")







p <- ggraph(multiome_tf_graph, layout = "fr") +
  geom_edge_link(alpha = 0.1) +
  geom_node_point(aes(color = perturb_diffused), size = 3) +
  scale_color_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-2, 2)
  ) +
  theme_void() +
  labs(color = "Diffused TF\nPerturbation")

ggsave("/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/test5.pdf")

############################################################################################################################










# TFs to keep (present in perturbed edges)
keep_tfs <- unique(edge_coef$from)

# Filter your gene_graph by TFs present in g and edges supported by g,
#    then JOIN the coef onto the edges by *names*, not by integer indices.
gene_graph_filtered <- gene_graph %>%
  # remove TF nodes not present in perturbed set
  activate(nodes) %>%
  filter(!(tf & !toupper(name) %in% keep_tfs)) %>%
  # move to edges and create name-based endpoints
  activate(edges) %>%
  mutate(from_name = .N()$name[from],
         to_name   = .N()$name[to]) %>%
  # enforce edges being in the perturbed set
  semi_join(edge_coef, by = c("from_name" = "from", "to_name" = "to")) %>%
  # attach coef; use it as weight (fallback to 0 if a gap remains)
  left_join(edge_coef, by = c("from_name" = "from", "to_name" = "to")) %>%
  # clean up helper columns
  select(-from_name, -to_name) %>%
  # drop isolates
  activate(nodes) %>%
  convert(to_subgraph, !node_is_isolated()) %>%
  mutate(centrality=tidygraph::centrality_pagerank())

# Make it explicitly directed
gene_graph_directed <- gene_graph_filtered %>% to_directed()

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

p <- ggraph(gene_graph_directed, x = UMAP_1, y = UMAP_2) +
  geom_edge_diagonal(
    aes(color = coef), 
    width = 0.1,
    arrow = arrow(length = unit(2, 'mm'), type = "closed"),
    end_cap = circle(1, 'mm'),
    arrow.fill = NA  # will make arrowheads use the same color as the edge
  ) +
  theme_graph() + 
  geom_node_point(
    aes(fill = group, size = centrality),
    shape = 21, stroke = 0.1
  ) +
  geom_label_repel(
    data = gene_graph_directed %>%
      tidygraph::activate("nodes") %>%
      as.data.frame() %>%
      dplyr::filter(name %in% keep_tfs),
    aes(x = UMAP_1, y = UMAP_2, label = name),
    size = 5, box.padding = 2, max.overlaps = Inf
  ) +
  scale_fill_manual(values = ct_colors) +
  scale_size(range = c(1, 8)) +
  scale_edge_colour_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-2, 2)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_void()

ggsave(filename = "/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/test2.pdf", plot = p, width = 10, height = 5)

ggsave(filename = paste0(figdir, "/pando_grn_eecs_network_graph_perturbation_informed.pdf"), plot = p, width = 10, height = 10)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)

target_gene <- "ARX"

# Convert to igraph
g_ig <- as.igraph(gene_graph_directed)

# Safety check
if (!(target_gene %in% V(g_ig)$name)) {
  stop(paste("Gene", target_gene, "not found in graph vertices."))
}

# 1-hop neighborhood
neighbors_vertices <- igraph::ego(g_ig, order = 2, nodes = target_gene, mode = "all")[[1]]

# Create induced subgraph
subgraph_ig <- igraph::induced_subgraph(g_ig, vids = neighbors_vertices)

# Convert back to tidygraph
subgraph_gene <- as_tbl_graph(subgraph_ig)

# 🔹 Remove isolated/unconnected nodes
subgraph_gene <- subgraph_gene %>%
  activate(nodes) %>%
  filter(!node_is_isolated())

p <- ggraph(subgraph_gene, layout = 'fr') +
  geom_edge_diagonal(
    aes(color = coef),
    width = 0.5,
    arrow = arrow(length = unit(2, 'mm'), type = "closed"),
    end_cap = circle(1, 'mm'),
    arrow.fill = NA
  ) +
  #geom_node_point(aes(fill = group), shape = 21, size=5, stroke = 0.1) +
  geom_node_text(
    aes(label = name),  # sugiyama supplies x/y
    size = 5, box.padding = 2, max.overlaps = Inf
  ) +
  scale_fill_manual(values = ct_colors) +
  scale_size(range = c(1, 8)) +
  scale_edge_colour_distiller(palette = "RdBu", direction = -1, limits = c(-2, 2)) +
  theme_void()

#ggsave(
#  filename = "/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/test.pdf",
#  plot = p, width = 10, height = 10
#)

ggsave(filename = "/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/test3.pdf", plot = p, width = 15, height = 10)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)

target_gene <- "LMX1B"

# Convert to igraph
g_ig <- as.igraph(gene_graph_directed)

# Safety check
if (!(target_gene %in% V(g_ig)$name)) {
  stop(paste("Gene", target_gene, "not found in graph vertices."))
}

# 1-hop neighborhood
neighbors_vertices <- igraph::ego(g_ig, order = 2, nodes = target_gene, mode = "all")[[1]]

# Create induced subgraph
subgraph_ig <- igraph::induced_subgraph(g_ig, vids = neighbors_vertices)

# Convert back to tidygraph
subgraph_gene <- as_tbl_graph(subgraph_ig)

# 🔹 Remove isolated/unconnected nodes
subgraph_gene <- subgraph_gene %>%
  activate(nodes) %>%
  filter(!node_is_isolated())

p <- ggraph(subgraph_gene, layout = 'fr') +
  geom_edge_diagonal(
    aes(color = coef),
    width = 0.5,
    arrow = arrow(length = unit(2, 'mm'), type = "closed"),
    end_cap = circle(1, 'mm'),
    arrow.fill = NA
  ) +
  #geom_node_point(aes(fill = group), shape = 21, size=5, stroke = 0.1) +
  geom_node_text(
    aes(label = name),  # sugiyama supplies x/y
    size = 5, box.padding = 2, max.overlaps = Inf
  ) +
  scale_fill_manual(values = ct_colors) +
  scale_size(range = c(1, 8)) +
  scale_edge_colour_distiller(palette = "RdBu", direction = -1, limits = c(-2, 2)) +
  theme_void()

#ggsave(
#  filename = "/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/test.pdf",
#  plot = p, width = 10, height = 10
#)

ggsave(filename = "/projects/site/pred/ihb-g-deco/USERS/adaml9/intestine_fate/notebooks/global_analysis/multiome/standalone_v2/pando/test4.pdf", plot = p, width = 5, height = 5)

