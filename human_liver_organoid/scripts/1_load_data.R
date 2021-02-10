library(plyr)
library(tidyverse)
library(MPStats)
library(googlesheets4)
library(hdf5r)
library(monocle3)
library(Matrix)

source("parameters.R")


# gather the dataset summary and ids from google drive
# HLOs/datasets
datasets <- googlesheets4::read_sheet(
    ss = parameters$datasets_googlesheet_id)

datasets %>%
    readr::write_tsv("raw_data/datasets_20210208.tsv")

datasets <- readr::read_tsv("raw_data/datasets_20210208.tsv")

# Load the output of the 10x genomics CellRanger pipeline stored in turbo
# https://cole-trapnell-lab.github.io/monocle3/docs/starting/#10x-output

# to begin we'll load the off-chip and on-chip datasets and
# and do a quick umap embedding as a sanity check


# off chip control
cds_CZ <- monocle3::load_cellranger_data(
    paste(parameters$data_base_dir, datasets$data_path[1], sep = "/"))

# we'll use the monocle3 defaults,
# expect set the Louvian clustering resolution to 1e-4 to get ~8 clusters
cds_CZ <- cds_CZ %>%
    monocle3::preprocess_cds(num_dim = 100) %>%
    monocle3::reduce_dimension(preprocess_method = "PCA") %>%
    monocle3::cluster_cells(
        resolution = 1e-4,
        verbose = TRUE)

plot <- cds_CZ %>%
    monocle3::plot_cells()

ggplot2::ggsave(
    filename = "product/figures/CZ/umap_20210208.pdf")


# on chip control
cds_CZ_1 <- monocle3::load_cellranger_data(
    paste(parameters$data_base_dir, datasets$data_path[2], sep = "/"))

cds_CZ_1 <- cds_CZ_1 %>%
    monocle3::preprocess_cds(num_dim = 100) %>%
    monocle3::reduce_dimension(preprocess_method = "PCA") %>%
    monocle3::cluster_cells(
        resolution = 1e-4,
        verbose = TRUE)

plot <- cds_CZ_1 %>%
    monocle3::plot_cells()

ggplot2::ggsave(
    filename = "product/figures/CZ_1/umap_20210208.pdf")


#######


# this adds a 'sample' column to colData(cds_CZ_CZ_1)
# indicating which dataset the cell came from ['CZ', 'CZ_1']
cds_CZ_CZ_1 <- monocle3::combine_cds(
    cds_list = list(
        CZ = cds_CZ,
        CZ_1 = cds_CZ_1))

cds_CZ_CZ_1 <- cds_CZ_CZ_1 %>%
    monocle3::preprocess_cds(num_dim = 100)

# this creates a new reducedDims object 'Aligned'
cds_CZ_CZ_1 <- cds_CZ_CZ_1 %>%
    monocle3::align_cds(
        num_dim = 100,
        alignment_group = "sample",
        alignment_k = 200,
        verbose = TRUE)

cds_CZ_CZ_1 <- cds_CZ_CZ_1 %>%
    monocle3::reduce_dimension(preprocess_method = "Aligned")

cds_CZ_CZ_1 <- cds_CZ_CZ_1 %>%
    monocle3::cluster_cells(
        resolution = 5e-4,
        verbose = TRUE)

plot <- cds_CZ_CZ_1 %>%
    monocle3::plot_cells(
        color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/CZ_CZ_1/umap_aligned_k=200_20210208.pdf")


plot <- cds_CZ_CZ_1 %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/CZ_CZ_1/umap_aligned_k=200_clusters_5e-4_20210208.pdf")

plot <- cds_CZ_CZ_1 %>%
    monocle3::plot_cells(color_cells_by = "partition")
ggplot2::ggsave(
    filename = "product/figures/CZ_CZ_1/umap_aligned_k=200_partitions_20210208.pdf")


# look for specific bio-markers in each cluster
marker_test_res <- cds_CZ_CZ_1 %>%
    monocle3::top_markers(
        group_cells_by = "cluster",
        reference_cells = 1000,
        cores = 8)

top_specific_markers <- marker_test_res %>%
    dplyr::filter(fraction_expressing >= 0.10) %>%
    dplyr::group_by(cell_group) %>%
    dplyr::top_n(5, pseudo_R2)

top_specific_markers %>%
    readr::write_tsv("product/figures/CZ_CZ_1/top_specific_markers_5e-4_20210210.tsv")


top_specific_marker_ids <- top_specific_markers %>% pull(gene_id) %>% unique()

plot <- cds_CZ_CZ_1 %>%
    monocle3::plot_genes_by_group(
        top_specific_marker_ids,
        group_cells_by = "cluster",
        ordering_type = "maximal_on_diag",
        max.size = 3)

ggplot2::ggsave("product/figures/CZ_CZ_1/top_speciic_markers_heatmap_5e-4_20210210.pdf")

######

# Combine all on-chip datasets

# annoyingly, Monocle3 wants the matrix data to be gzipped
# but the cellrange just has it uncompressed files for a few of the datasets,
# so just make a compressed version in the same directory where needed
for (i in c(2:7)) {
    matrix_dir <- paste(
        parameters$data_base_dir,
        datasets$data_path[i],
        "outs",
        "filtered_feature_bc_matrix",
        sep = "/")
    cat("# compressing matrix files for path '", matrix_dir, "'", sep = "", sep = "\n")
    cat(paste0("gzip -c ", matrix_dir, "/barcodes.tsv > ", matrix_dir, "/barcodes.tsv.gz"), sep = "", "\n")
    cat(paste0("gzip -c ", matrix_dir, "/features.tsv > ", matrix_dir, "/features.tsv.gz"), sep = "", "\n")
    cat(paste0("gzip -c ", matrix_dir, "/matrix.mtx > ", matrix_dir, "/matrix.mtx.gz"), sep = "", "\n")
}

cds_CZ_x <- monocle3::combine_cds(
    cds_list = list(
        Acetaminophen = monocle3::load_cellranger_data(paste(parameters$data_base_dir, datasets$data_path[2], sep = "/")),
        Control = monocle3::load_cellranger_data(paste(parameters$data_base_dir, datasets$data_path[3], sep = "/")),
        Fialuridine = monocle3::load_cellranger_data(paste(parameters$data_base_dir, datasets$data_path[4], sep = "/")),
        Tenofovir = monocle3::load_cellranger_data(paste(parameters$data_base_dir, datasets$data_path[5], sep = "/")),
        Inarigivir = monocle3::load_cellranger_data(paste(parameters$data_base_dir, datasets$data_path[6], sep = "/")),
        Tenofovir_Inarigivir  = monocle3::load_cellranger_data(paste(parameters$data_base_dir, datasets$data_path[7], sep = "/"))))


cds_CZ_x <- cds_CZ_x %>%
    monocle3::preprocess_cds(num_dim = 100)


cds_CZ_x <- cds_CZ_x %>%
    monocle3::align_cds(
        num_dim = 100,
        alignment_group = "sample",
        alignment_k = 200,
        verbose = TRUE)


cds_CZ_x <- cds_CZ_x %>%
    monocle3::reduce_dimension(preprocess_method = "Aligned")

cds_CZ_x <- cds_CZ_x %>%
    monocle3::cluster_cells(
        resolution = 1e-5,
        verbose = TRUE)

plot <- cds_CZ_x %>%
    monocle3::plot_cells(
        color_cells_by = "sample")
ggplot2::ggsave(
    filename = "product/figures/CZ_x/umap_aligned_k=200_1e-5_20210208.pdf")

plot <- cds_CZ_x %>%
    monocle3::plot_cells()
ggplot2::ggsave(
    filename = "product/figures/CZ_x/umap_aligned_k=200_clusters_1e-5_20210208.pdf")



###
# look for specific bio-markers in each cluster
marker_test_res <- cds_CZ_x %>%
    monocle3::top_markers(
        group_cells_by = "cluster",
        reference_cells = 1000,
        cores = 8)

top_specific_markers <- marker_test_res %>%
    dplyr::filter(fraction_expressing >= 0.10) %>%
    dplyr::group_by(cell_group) %>%
    dplyr::top_n(5, pseudo_R2)


plot <- cds_CZ_x %>%
    monocle3::plot_genes_by_group(
        top_specific_marker_ids,
        group_cells_by = "cluster",
        ordering_type = "maximal_on_diag",
        max.size = 3)

ggplot2::ggsave(
    "product/figures/CZ_x/top_specific_markers_heatmap_1e-5_20210210.pdf",
    width = 9,
    height = 13)


cluster_by_sample <- tibble::tibble(
    cluster_label = cds_CZ_x@clusters[['UMAP']]$clusters,
    sample = colData(cds_CZ_x)$sample) %>%
    dplyr::mutate(
        sample_name = sample %>%
            dplyr::recode_factor(
                "CZ_1" = "Acetaminophen",
                "CZ_2" = "Control",
                "CZ_3" = "Fialuridine",
                "CZ_4" = "Tenofovir",
                "CZ_5" = "Inarigivir",
                "CZ_6" = "Tenofovir_Inarigivir"))


cluster_by_sample %>%
    dplyr::count(cluster_label, sample_name) %>%
    tidyr::pivot_wider(
        names_from = sample_name,
        values_from = n)
