library(plyr)
library(tidyverse)
library(fuzzyjoin)
library(ggplot2)
library(readxl)
library(MPStats)
library(arrow)
library(rdist)

cell_feature_columns <- readr::read_tsv(
    "product/cell_feature_columns_TS_202008.tsv")

cell_metadata_columns <- readr::read_tsv(
    "product/cell_metadata_columns_TS_202008.tsv")


###################
# Plot distortion #
###################

cell_features <- dplyr::bind_cols(
    arrow::read_parquet(
        file = "product/TS2_2M_Cell_MasterDataTable.parquet"),
    arrow::read_parquet(
        file = "intermediate_data/UMAP_embedding_TS2_2M_epochs=2000_20200901/umap_embedding.parquet"))

# sample disjoint sets of source and target cells
n_cells <- 2000
cell_indices <- cell_features %>% nrow() %>% sample(size = 2*n_cells)
source_cell_indices <- cell_indices[1:n_cells]
target_cell_indices <- cell_indices[(n_cells+1):(2*n_cells)]

data <- data.frame(
    native = rdist::cdist(
        X = cell_features %>%
            dplyr::slice(source_cell_indices) %>%
            dplyr::select(tidyselect::one_of(cell_feature_columns$feature)),
        Y = cell_features %>%
            dplyr::slice(target_cell_indices) %>%
            dplyr::select(tidyselect::one_of(cell_feature_columns$feature))) %>%
        as.vector() %>%
        rank(na.last = TRUE, ties.method = "average"),
    embedded = rdist::cdist(
        X = cell_features %>%
            dplyr::slice(source_cell_indices) %>%
            dplyr::select(UMAP_1, UMAP_2),
        Y = cell_features %>%
            dplyr::slice(target_cell_indices) %>%
            dplyr::select(UMAP_1, UMAP_2)) %>%
        as.vector() %>%
        rank(na.last = TRUE, ties.method = "average")) %>%
    dplyr::mutate(
        native = native / (n_cells * n_cells),
        embedded = embedded / (n_cells * n_cells))

cor(data$native, data$embedded, method="spearman")

p <- ggplot2::ggplot(data = data) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
             x = embedded,
             y = native),
        size = .005,
        alpha = .01) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "blue") +
    ggplot2::scale_x_continuous(
         "Embedded cell-cell distance rank",
         expand = c(0,0)) +
    ggplot2::scale_y_continuous(
        "Native cell-cell distance rank",
        expand = c(0,0)) +
    ggplot2::ggtitle(
        label = "Native vs Embedding cell-cell distances: TS2 Dataset",
        subtitle = "4M sampled distances between cells")

ggplot2::ggsave(
    filename = "intermediate_data/UMAP_embedding_TS2_2M_epochs=2000_20200901/figures/sp_distance_rank_correlation.png",
    width = 8.5,
    height = 8.5)
