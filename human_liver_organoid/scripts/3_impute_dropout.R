


library(plyr)
library(tidyverse)
library(MPStats)
library(monocle3)
library(Rmagic)

source("parameters.R")
source("scripts/get_dataset.R")

datasets <- readr::read_tsv("raw_data/datasets_20210208.tsv", col_types = readr::cols())


# off chip control
cds_CZ <- get_dataset("2603-CZ")

cds_CZ_magic <- cds_CZ %>%
    SingleCellExperiment::counts() %>%
    Rmagic::magic()

cds_CZ <- monocle3::new_cell_data_set(
    expression_data = as.matrix(cds_CZ_magic$result) %>% Matrix::t(),
    cell_metadata = rowData(cds_CZ),
    gene_metadata = colData(cds_CZ) %>% as.data.frame())

cds_CZ <- cds_CZ %>%
    monocle3::preprocess_cds(num_dim = 300) %>%
    monocle3::reduce_dimension(preprocess_method = "PCA") %>%
    monocle3::cluster_cells(
        resolution = 1e-4,
        verbose = TRUE)

cds_CZ <- cds_CZ %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        a = 10,
        b = 1,
        verbose = TRUE)

plot <- cds_CZ %>%
    monocle3::plot_cells()

ggplot2::ggsave(
    filename = "product/figures/CZ/umap_pca300_a=10,b=1_magic_20210218.pdf")

