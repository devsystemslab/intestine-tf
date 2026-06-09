find_modules <- function(
    object,
    p_thresh = 0.05,
    rsq_thresh = 0.1,
    nvar_thresh = 10,
    min_genes_per_module = 5,
    xgb_method = c('tf', 'target'),
    xgb_top = 50,
    verbose = TRUE
){
    fit_method <- NetworkParams(object)$method
    xgb_method <- match.arg(xgb_method)

    models_use <- gof(object) %>%
        filter(rsq>rsq_thresh & nvariables>nvar_thresh) %>%
        pull(target) %>%
        unique()

    modules <- coef(object) %>%
        filter(target %in% models_use)

    if (fit_method %in% c('cv.glmnet', 'glmnet')){
        modules <- modules %>%
            filter(estimate != 0)
    } else if (fit_method == 'xgb'){
        modules <- modules %>%
            group_by_at(xgb_method) %>%
            top_n(xgb_top, gain) %>%
            mutate(estimate=sign(corr)*gain)
    } else {
        modules <- modules %>%
            filter(ifelse(is.na(padj), T, padj<p_thresh))
    }

    modules <- modules %>%
        group_by(target) %>%
        mutate(nvars=n()) %>%
        group_by(target, tf) %>%
        mutate(tf_sites_per_gene=n()) %>%
        group_by(target) %>%
        mutate(
            tf_per_gene=length(unique(tf)),
            peak_per_gene=length(unique(region))
        ) %>%
        group_by(tf) %>%
        mutate(gene_per_tf=length(unique(target))) %>%
        group_by(target, tf)

    if (fit_method %in% c('cv.glmnet', 'glmnet', 'xgb')){
        modules <- modules %>%
            reframe(
                estimate=sum(estimate),
                n_regions=peak_per_gene,
                n_genes=gene_per_tf,
                n_tfs=tf_per_gene,
                regions=paste(region, collapse=';')
            )
    } else {
        modules <- modules %>%
            reframe(
                estimate=sum(estimate),
                n_regions=peak_per_gene,
                n_genes=gene_per_tf,
                n_tfs=tf_per_gene,
                regions=paste(region, collapse=';'),
                pval=min(pval),
                padj=min(padj)
            )
    }

    modules <- modules %>%
        distinct() %>%
        arrange(tf)

    module_pos <- modules %>%
        filter(estimate>0) %>%
        group_by(tf) %>% filter(n()>min_genes_per_module) %>%
        group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>%
        map(function(x) x$target)

    module_neg <- modules %>%
        filter(estimate<0) %>%
        group_by(tf) %>% filter(n()>min_genes_per_module) %>%
        group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>%
        map(function(x) x$target)

    regions_pos <- modules %>%
        filter(estimate>0) %>%
        group_by(tf) %>% filter(n()>min_genes_per_module) %>%
        group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>%
        map(function(x) unlist(str_split(x$regions, ';')))

    regions_neg <- modules %>%
        filter(estimate<0) %>%
        group_by(tf) %>% filter(n()>min_genes_per_module) %>%
        group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>%
        map(function(x) unlist(str_split(x$regions, ';')))

    module_feats <- list(
        'genes_pos' = module_pos,
        'genes_neg' = module_neg,
        'regions_pos' = regions_pos,
        'regions_neg' = regions_neg
    )

 
    module_meta <- dplyr::select(modules, tf, target, everything())
    object@modules@meta <- module_meta
    object@modules@features <- module_feats
    object@modules@params <- list(
        p_thresh = p_thresh,
        rsq_thresh = rsq_thresh,
        nvar_thresh = nvar_thresh,
        min_genes_per_module = min_genes_per_module
    )
    return(object)
}


#' Compute UMAP embedding
get_umap <- function(
    x,
    n_pcs = 30,
    ...
){
    if (ncol(x)>100){
        pca_mat <- irlba::prcomp_irlba(x, n=n_pcs)$x
        rownames(pca_mat) <- rownames(x)
        x <- as.matrix(pca_mat)
    }
    umap_tbl <- uwot::umap(x, ...) %>%
        {colnames(.) <- c('UMAP_1', 'UMAP_2'); .} %>%
        as_tibble(rownames='gene')
    return(umap_tbl)
}

get_tf_network_graph <- function(
    object,
    network = DefaultNetwork(object),
    graph_name = 'module_graph',
    rna_assay = 'RNA',
    rna_layer = 'data',
    features = NULL,
    random_seed = 111,
    verbose = TRUE,
    ...
){
    modules <- NetworkModules(object, network=network)

    if (length(modules@params)==0){
        stop('No modules found, please run `find_modules()` first.')
    }

    if (is.null(features)){
        features <- NetworkFeatures(object, network=network)
    }

    rna_expr <- t(LayerData(object, assay = rna_assay, layer = rna_layer))
    features <- intersect(features, colnames(rna_expr))

    message("Computing gene-gene correlation")
    rna_expr <- rna_expr[, features]
    gene_cor <- sparse_cor(rna_expr)

    gene_cor_df <- gene_cor %>%
        as_tibble(rownames = "source") %>%
        pivot_longer(!source, names_to = "target", values_to = "corr")

    modules_use <- modules@meta %>%
        filter(
            target %in% colnames(rna_expr),
            tf %in% colnames(rna_expr)
        )

    gene_net <- modules_use %>%
        select(tf, target, everything()) %>%
        group_by(target) %>%
        left_join(gene_cor_df, by = c("tf" = "source", "target")) %>%
        mutate(corr = replace_na(corr, 0))

    # target x TF regulatory matrix
    reg_mat <- gene_net %>%
        select(target, tf, estimate) %>%
        pivot_wider(
            names_from = tf,
            values_from = estimate,
            values_fill = 0
        ) %>%
        column_to_rownames("target") %>%
        as.matrix()

    message("Computing weighted regulatory factor")

    # target x TF weight matrix
    reg_factor_mat <- abs(reg_mat) + 1

    # target x TF coexpression-weighted regulatory matrix
    coex_target_tf <- gene_cor[
        rownames(reg_factor_mat),
        colnames(reg_factor_mat),
        drop = FALSE
    ] * sqrt(reg_factor_mat)

    # transpose so UMAP is computed on TFs instead of target genes
    coex_mat <- t(coex_target_tf)

    print('Computing UMAP embedding')
    set.seed(random_seed)
    coex_umap <- get_umap(as.matrix(coex_mat), ...)

    print('Getting network graph')
    gene_graph <- as_tbl_graph(gene_net) %>%
        activate(edges) %>%
        mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
        mutate(dir=sign(estimate)) %>%
        activate(nodes) %>%
        mutate(centrality=centrality_pagerank())

    print(gene_graph)
    print(coex_umap)
    gene_graph <- gene_graph %>%
        inner_join(coex_umap, by=c('name'='gene'))

    return(gene_graph)
}
