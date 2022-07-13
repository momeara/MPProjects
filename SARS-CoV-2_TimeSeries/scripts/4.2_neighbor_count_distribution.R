
library(tidyverse)
library(arrow)

source("scripts/compute_spatial_features.R")

cell_features_TS2PL1 <- arrow::read_parquet(
    file = "product/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet") %>%
  add_well_coordinates()

cell_feature_columns <- readr::read_tsv(
  "product/cell_feature_columns_TS_202008.tsv")

cell_metadata_columns <- readr::read_tsv(
  "product/cell_metadata_columns_TS_202008.tsv")

# generate neighbor graph
neighbor_graph <- cell_features_TS2PL1 %>%
  dplyr::filter(time_point == "Uninfected") %>%
  compute_neighbor_graph_by_well(verbose = TRUE)


# plot distribution neighbor count distribution per field
neighbor_distribution <- neighbor_graph %>%
  dplyr::count(
    plate_id, row, column,
    cell1_well_index,
    name = "n_neighbors") %>%
  dplyr::group_by(
    plate_id, row, column,
    n_neighbors) %>%
  dplyr::summarize(
    neighbor_count_density = dplyr::n()) %>%
  dplyr::ungroup()


plot <- ggplot2::ggplot(data = neighbor_distribution) +
  ggplot2::theme_bw() +
  ggplot2::geom_line(
    mapping = ggplot2::aes(
      x = n_neighbors,
      y = neighbor_count_density,
      group = paste0("Row:", row, " Col:", column)),
    alpha = 0.7) +
  ggplot2::scale_x_continuous(
    "Number of Neighbors < 100 pixels",
    breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    expand = c(0, 0),
    limits = c(1, 12)) +
  ggplot2::scale_y_continuous("Count") +
  ggplot2::theme(panel.grid.minor.x = element_blank())

ggplot2::ggsave(
  filename = "product/figures/TS2/spatial/neighbor_distribution_uinfected_by_well_20220419.pdf",
  plot = plot,
  width = 5,
  height = 3)

ggplot2::ggsave(
  filename = "product/figures/TS2/spatial/neighbor_distribution_uinfected_by_well_20220419.png",
  plot = plot,
  width = 5,
  height = 3)
