library(plyr)
library(tidyverse)
library(MPStats)
library(monocle3)

source("parameters.R")
source("scripts/get_dataset.R")

datasets <- readr::read_tsv("raw_data/datasets_20210208.tsv", col_types = readr::cols())


# to begin we'll load the off-chip and on-chip datasets and
# and do a quick umap embedding as a sanity check

# off chip control
cds_CZ <- get_dataset("2603-CZ")



# expect set the Louvian clustering resolution to 1e-4 to get ~8 clusters
cds_CZ <- cds_CZ %>%
    monocle3::preprocess_cds(num_dim = 100) %>%
    monocle3::reduce_dimension(preprocess_method = "PCA") %>%
    monocle3::cluster_cells(
        resolution = 1e-4,
        verbose = TRUE)
save(cds_CZ, file = "intermediate_data/cds_CZ.Rdata")
    



# on chip control
cds_CZ_2 <- get_dataset("2602-CZ-2")

cds_CZ_2 <- cds_CZ_2 %>%
    monocle3::preprocess_cds(num_dim = 100) %>%
    monocle3::reduce_dimension(preprocess_method = "PCA") %>%
    monocle3::cluster_cells(
        resolution = 1e-4,
        verbose = TRUE)
save(cds_CZ_2, file = "intermediate_data/cds_CZ_2.Rdata")



#######

# compare on chip with off chip
# this adds a 'sample' column to colData(cds_CZ_CZ_2)
# indicating which dataset the cell came from ['CZ', 'CZ_2']
cds_CZ_CZ_2 <- monocle3::combine_cds(
    cds_list = list(
        CZ = cds_CZ,
        CZ_2 = cds_CZ_2))

cds_CZ_CZ_2 <- cds_CZ_CZ_2 %>%
    monocle3::preprocess_cds(num_dim = 100)

# this creates a new reducedDims object 'Aligned'
cds_CZ_CZ_2 <- cds_CZ_CZ_2 %>%
    monocle3::align_cds(
        num_dim = 100,
        alignment_group = "sample",
        alignment_k = 200,
        verbose = TRUE)

cds_CZ_CZ_2 <- cds_CZ_CZ_2 %>%
    monocle3::reduce_dimension(preprocess_method = "Aligned")

cds_CZ_CZ_2 <- cds_CZ_CZ_2 %>%
    monocle3::cluster_cells(
        resolution = 5e-4,
        verbose = TRUE)
save(cds_CZ_CZ_2, file = "intermediate_data/cds_CZ_CZ_2.Rdata")



# look for specific bio-markers in each cluster
marker_test_res <- cds_CZ_CZ_2 %>%
    monocle3::top_markers(
        group_cells_by = "cluster",
        reference_cells = 1000,
        cores = 8)

top_specific_markers <- marker_test_res %>%
    dplyr::filter(fraction_expressing >= 0.10) %>%
    dplyr::group_by(cell_group) %>%
    dplyr::top_n(5, pseudo_R2)

top_specific_markers %>%
    readr::write_tsv("product/figures/CZ_CZ_2/top_specific_markers_5e-4_20210210.tsv")


top_specific_marker_ids <- top_specific_markers %>% pull(gene_id) %>% unique()

plot <- cds_CZ_CZ_2 %>%
    monocle3::plot_genes_by_group(
        top_specific_marker_ids,
        group_cells_by = "cluster",
        ordering_type = "maximal_on_diag",
        max.size = 3)

ggplot2::ggsave("product/figures/CZ_CZ_2/top_speciic_markers_heatmap_5e-4_20210210.pdf")


#####################
# combine some of the on-chip treatments
cds_CZ_treat <- monocle3::combine_cds(
    cds_list = list(
        Control = get_dataset("2602-CZ-2"),
        Fialuridine = get_dataset("2602-CZ-3"),
        Tenofovir = get_dataset("2602-CZ-4"),
        Tenofovir_Inarigivir = get_dataset("2602-CZ-6")))


cds_CZ_treat <- cds_CZ_treat %>%
    monocle3::preprocess_cds(num_dim = 100)


cds_CZ_treat <- cds_CZ_treat %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        verbose = TRUE)

# resolution: 1e-5
#   Modularity is 0.734636987448508
#   Quality is 2116543.19454
#   Significance is 868460.182810843
#   Number of clusters is 12
# resolution: 1e-06;
#   Modularity is 0.446991090433721;
#   Quality is 2125276.611134;
#   Significance is 318007.466408565;
#   Number of clusters is 8
cds_CZ_treat <- cds_CZ_treat %>%
    monocle3::cluster_cells(
        resolution = 1e-5,
        verbose = TRUE)

save(cds_CZ_treat, file = "intermediate_data/cds_CZ_treat.Rdata")



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
        Acetaminophen = get_dataset("2602-CZ-1"),
        Control = get_dataset("2602-CZ-2"),
        Fialuridine = get_dataset("2602-CZ-3"),
        Tenofovir = get_dataset("2602-CZ-4"),
        Inarigivir = get_dataset("2602-CZ-5"),
        Tenofovir_Inarigivir = get_dataset("2602-CZ-6")))


cds_CZ_x <- cds_CZ_x %>%
    monocle3::preprocess_cds(num_dim = 100)


# Don't align because they were all run in the same lane
#cds_CZ_x <- cds_CZ_x %>%
#    monocle3::align_cds(
#        num_dim = 100,
#        alignment_group = "sample",
#        alignment_k = 200,
#        verbose = TRUE)


cds_CZ_x <- cds_CZ_x %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        verbose = TRUE)

cds_CZ_x <- cds_CZ_x %>%
    monocle3::cluster_cells(
        resolution = 1e-5,
        verbose = TRUE)
save(cds_CZ_x, file = "intermediate_data/cds_CZ_x.Rdata")

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

top_specific_marker_ids <- top_specific_markers %>% pull(gene_id) %>% unique()

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


cluster_by_sample %>%
    dplyr::count(cluster_label, sample_name) %>%
    tidyr::pivot_wider(
        names_from = sample_name,
        values_from = n)


##################
# Load Ouchi dataset

cds_Ouchi2019 <- get_dataset("GSM3731527")
cds_Ouchi2019 <- cds_Ouchi2019 %>%
    monocle3::preprocess_cds(num_dim = 100)

cds_Ouchi2019 <- cds_Ouchi2019 %>%
    monocle3::reduce_dimension(
        n_epochs = 2000,
        preprocess_method = "PCA",
        umap.n_neighbors = 500,
        cores = 2,
        a = 10,
        b = 1,
        verbose = TRUE)


cds_Ouchi2019 <- cds_Ouchi2019 %>%
    monocle3::cluster_cells(
        resolution = 1e-4,
        verbose = TRUE)
save(cds_Ouchi2019, file = "intermediate_data/cds_Ouchi2019.Rdata")


##################
# Ouchi2019 + on_chip control

cds_Ouchi2019_CZ_2 <- monocle3::combine_cds(
    cds_list = list(
        Ouchi2019 = get_dataset("GSM3731527"),
        Control = get_dataset("2602-CZ-2")))

cds_Ouchi2019_CZ_2 <- cds_Ouchi2019_CZ_2 %>%
    monocle3::preprocess_cds(num_dim = 100)


cds_Ouchi2019_CZ_2 <- cds_Ouchi2019_CZ_2 %>%
    monocle3::align_cds(
        num_dim = 100,
        alignment_group = "sample",
        alignment_k = 200,
        verbose = TRUE)


cds_Ouchi2019_CZ_2 <- cds_Ouchi2019_CZ_2 %>%
    monocle3::reduce_dimension(preprocess_method = "Aligned")

cds_Ouchi2019_CZ_2 <- cds_Ouchi2019_CZ_2 %>%
    monocle3::cluster_cells(
        resolution = 1e-5,
        verbose = TRUE)
save(cds_Ouchi2019_CZ_2, file = "intermediate_data/cds_Ouchi2019_CZ_2.Rdata")



##################
# Ouchi2019 and CZ_x

cds_Ouchi2019_CZ_x <- monocle3::combine_cds(
    cds_list = list(
        Ouchi2019 = get_dataset("GSM3731527"),
        Acetaminophen = get_dataset("2602-CZ-1"),
        Control = get_dataset("2602-CZ-2"),
        Fialuridine = get_dataset("2602-CZ-3"),
        Tenofovir = get_dataset("2602-CZ-4"),
        Inarigivir = get_dataset("2602-CZ-5"),
        Tenofovir_Inarigivir = get_dataset("2602-CZ-6")))


cds_Ouchi2019_CZ_x <- cds_Ouchi2019_CZ_x %>%
    monocle3::preprocess_cds(num_dim = 100)


cds_Ouchi2019_CZ_x <- cds_Ouchi2019_CZ_x %>%
    monocle3::align_cds(
        num_dim = 100,
        alignment_group = "sample",
        alignment_k = 200,
        verbose = TRUE)


cds_Ouchi2019_CZ_x <- cds_Ouchi2019_CZ_x %>%
    monocle3::reduce_dimension(preprocess_method = "Aligned")

cds_Ouchi2019_CZ_x <- cds_Ouchi2019_CZ_x %>%
    monocle3::cluster_cells(
        resolution = 1e-5,
        verbose = TRUE)
save(cds_Ouchi2019_CZ_x, file = "intermediate_data/cds_Ouchi2019_CZ_x.Rdata")



###############################
# On Chip Hepatocyte clusters #
###############################
cds_CZ_x_hepatocyte <- cds_CZ_x
cds_CZ_x_hepatocyte <- cds_CZ_x_hepatocyte %>%
    MPStats::compute_qc_covariates()

cds_CZ_x_hepatocyte <- cds_CZ_x_hepatocyte[
    # genes that are observed in at least 5 cells
    SingleCellExperiment::counts(cds_CZ_x_hepatocyte) %>% Matrix::rowSums() >= 5,
    # hepatocyte cluster
    (monocle3::clusters(cds_CZ_x_hepatocyte, reduction_method = "UMAP") %in% c(1, 4, 7, 9, 11)) &
    # cells that have at least 1000 reads which are at most 30% mitochondrial
    (SummarizedExperiment::colData(cds_CZ_x_hepatocyte)[["count_depth"]] >= 10000) &
    (SummarizedExperiment::colData(cds_CZ_x_hepatocyte)[["mt_fraction"]] < .3)]

cds_CZ_x_hepatocyte <- cds_CZ_x_hepatocyte %>%
    monocle3::preprocess_cds(num_dim = 100)

cds_CZ_x_hepatocyte <- cds_CZ_x_hepatocyte %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        verbose = TRUE)

save(cds_CZ_x_hepatocyte, file = "intermediate_data/cds_CZ_x_hepatocyte.Rdata")


###
# look for specific bio-markers in each cluster
marker_test_res <- cds_Ouchi2019_CZ_x %>%
    monocle3::top_markers(
        group_cells_by = "cluster",
        reference_cells = 1000,
        cores = 8)

top_specific_markers <- marker_test_res %>%
    dplyr::filter(fraction_expressing >= 0.10) %>%
    dplyr::group_by(cell_group) %>%
    dplyr::top_n(5, pseudo_R2)

top_specific_marker_ids <- top_specific_markers %>% pull(gene_id) %>% unique()

plot <- cds_Ouchi2019_CZ_x %>%
    monocle3::plot_genes_by_group(
        top_specific_marker_ids,
        group_cells_by = "cluster",
        ordering_type = "maximal_on_diag",
        max.size = 3)

ggplot2::ggsave(
    "product/figures/Ouchi2019_CZ_x/top_specific_markers_heatmap_1e-5_20210210.pdf",
    width = 9,
    height = 13)


cluster_by_sample <- tibble::tibble(
    cluster_label = cds_Ouchi2019_CZ_x@clusters[["UMAP"]]$clusters,
    sample = colData(cds_Ouchi2019_CZ_x)$sample) %>%


cluster_by_sample %>%
    dplyr::count(cluster_label, sample_name) %>%
    tidyr::pivot_wider(
        names_from = sample_name,
        values_from = n)
