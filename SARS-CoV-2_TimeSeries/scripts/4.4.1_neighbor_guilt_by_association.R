
library(tidyverse)
library(arrow)
library(EGAD)

source("scripts/compute_spatial_features.R")

cell_features_TS2PL1 <- arrow::read_parquet(
  file = "~/Documents/maom_lab/MPProjects/SARS-CoV-2_TimeSeries/product/TS2PL1_Cell_MasterDataTable.parquet") %>%
  compute_well_coordinates()

cell_feature_columns <- readr::read_tsv(
  "product/cell_feature_columns_TS_202008.tsv")

cell_metadata_columns <- readr::read_tsv(
  "product/cell_metadata_columns_TS_202008.tsv")

cell_features_rel <- cell_features_TS2PL1 %>%
  dplyr::filter(time_point == "Uninfected")


# generate neighbor graph
neighbor_graph <- cell_features_rel %>%
  compute_neighbor_graph_by_well(verbose = TRUE)


feature_gba <- cell_feature_columns %>%
  dplyr::rowwise() %>%
  tidyr::map_drc(function(row){
    feature <- row$feature[1]
    cat("Cell feature: ", feature, "\n", sep = "")
    cell_features_rel %>%
      dplyr::select(
        plate_id,
        row,
        column,
        cell1_well_index = well_index,
        feature_value1 = tidyselect::any_of(feature)) %>%
      dplyr::left_join(
        neighbor_graph,
        by = c("plate_id", "row", "column", "cell1_well_index")) %>%
      dplyr::left_join(
        cell_features_rel %>%
          dplyr::select(
            plate_id,
            row,
            column,
            cell2_well_index = well_index,
            feature_value2 = tidyselect::any_of(feature)),
        by = c("plate_id", "row", "column", "cell2_well_index")) %>%
      dplyr::group_by(cell1_well_index) %>%
      dplyr::summarize(
        {{feature}} := mean(feature_value2, na.rm = TRUE),
        .groups = "drop")
  }) %>%
  dplyr::ungroup()




