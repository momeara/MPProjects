

library(plyr)
library(tidyverse)
library(MPStats)
library(monocle3)
library(Rmagic)

source("parameters.R")
source("scripts/get_dataset.R")

datasets <- readr::read_tsv("raw_data/datasets_20210208.tsv", col_types = readr::cols())
marker_genes <- readr::read_tsv(parameters$marker_genes_fname)

# on chip control
load("intermediate_data/cds_CZ_2.Rdata")

# remove genes with zero expression (36,601 -> 27,412)
cds_CZ_2 <- cds_CZ_2[Matrix::rowSums(exprs(cds_CZ_2)) != 0, ]

gene_set <- marker_genes %>%
    dplyr::filter(gene_set == "kupffer")

cds_CZ_2_magic <- cds_CZ_2 %>%
    SingleCellExperiment::counts() %>%
    Matrix::t()

cds_CZ_2_magic@Dimnames <- list(
    cds_CZ_2_magic@Dimnames[[1]],
    rowData(cds_CZ_2)$gene_short_name)

cds_CZ_2_magic <- cds_CZ_2_magic %>%
    Rmagic::magic(
        genes = gene_set$gene,
        knn_dist.method = "cosine")


expression_data <- as.matrix(cds_CZ_2_magic$result) %>% Matrix::t()
cell_metadata <- colData(cds_CZ_2) %>% as.data.frame()
gene_metadata <- rowData(cds_CZ_2) %>%
    as.data.frame() %>%
    dplyr::semi_join(gene_set, by = c("gene_short_name" = "gene"))
rownames(expression_data) <- rownames(gene_metadata)

cds_CZ_2_magic_2 <- monocle3::new_cell_data_set(
    expression_data = expression_data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata)

reducedDim(cds_CZ_2_magic_2, "UMAP") <- reducedDim(cds_CZ_2, "UMAP")


plot <- cds_CZ_2_magic_2 %>%
    monocle3::plot_cells(genes = unique(gene_set$gene)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
    filename = paste0(
        "product/figures/CZ_2/umap_magic_cupffer_", "cupffer", "_", MPStats::date_code(), ".png"),
    width = 15,
    height = 14,
    plot = plot)




expand.grid(
    a = seq(1, 5, length.out = 3),
    b = seq(1, 2, length.out = 3)) %>%
    dplyr::rowwise() %>%
    dplyr::do({
    params <- .
    cat("Generating UMAP for params a=", params$a, ", b=", params$b, "\n", sep ="")
    cds_CZ_2 <- cds_CZ_2 %>%
        monocle3::reduce_dimension(
            preprocess_method = "PCA",
            a = params$a,
            b = params$b,
            verbose = TRUE)

    # cds_CZ_2 <- cds_CZ_2 %>%
    #     monocle3::cluster_cells(
    #         resolution = 1e-4,
    #         verbose = TRUE)
    
    plot <- cds_CZ_2 %>%
        monocle3::plot_cells()
    
    ggplot2::ggsave(
        filename = paste0(
            "product/figures/CZ/umap_magic_a=", params$a, ",b=", params$b, "_20210218.pdf"))
    })

##
# look for specific bio-markers in each cluster
marker_test_res <- cds_CZ %>%
    monocle3::top_markers(
        group_cells_by = "cluster",
        reference_cells = 1000,
        cores = 8)

top_specific_markers <- marker_test_res %>%
    dplyr::filter(fraction_expressing >= 0.10) %>%
    dplyr::group_by(cell_group) %>%
    dplyr::top_n(5, pseudo_R2)

top_specific_markers %>%
    readr::write_tsv("product/figures/CZ/top_specific_markers_5e-4_20210218.tsv")


top_specific_marker_ids <- top_specific_markers %>% pull(gene_id) %>% unique()

plot <- cds_CZ %>%
    monocle3::plot_genes_by_group(
        top_specific_marker_ids,
        group_cells_by = "cluster",
        ordering_type = "maximal_on_diag",
        max.size = 3)

ggplot2::ggsave("product/figures/CZ/top_speciic_markers_heatmap_5e-4_20210218.pdf")


