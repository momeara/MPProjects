
library(tidyverse)
library(arrow)

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
rm(cell_features_TS2PL1)
gc()

# generate neighbor graph
neighbor_graph <- cell_features_rel %>%
  compute_neighbor_graph_by_well(verbose = TRUE)


# plot distribution neighbor count distribution per field
neighbor_distribution <- neighbor_graph %>%
  dplyr::count(
    plate_id, row, column,
    cell1_number,
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
      color = Image_Metadata_WellID)) +
  ggplot2::scale_x_continuous(
    "Number of Neighbors < 100 pixels",
    breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
  ggplot2::scale_y_continuous("Count") +
  ggplot2::scale_color_discrete("Well") +
  ggplot2::theme(panel.grid.minor.x = element_blank())

ggplot2::ggsave(
  filename = "product/figures/TS2/spatial/neighbor_distribution_uinfected_by_well_20220418.pdf",
  plot = plot,
  width = 5,
  height = 3)


##################



# compute per-feature neighbor correlation
cell_features_rel <- cell_features_TS2PL1 %>%
  dplyr::filter(time_point == "Uninfected")

spatial_correlation <- cell_feature_columns %>%
  dplyr::rowwise() %>%
  dplyr::do({
    cell_feature_column <- .
    cat("Cell feature: ", cell_feature_column$feature, "\n", sep = "")
    neighbor_graph %>%
      dplyr::left_join(
        cell_features_rel %>%
          dplyr::select(
            plate_id,
            row,
            column,
            Image_Metadata_FieldID,
            cell1_number = Cells_Number_Object_Number,
            feature_value1 = tidyselect::any_of(cell_feature_column$feature)),
        by = c("plate_id", "row", "column", "Image_Metadata_FieldID", "cell1_number")) %>%
      dplyr::left_join(
        cell_features_rel %>%
          dplyr::select(
            plate_id,
            row,
            column,
            Image_Metadata_FieldID,
            cell2_number = Cells_Number_Object_Number,
            feature_value2 = tidyselect::any_of(cell_feature_column$feature)),
        by = c("plate_id", "row", "column", "Image_Metadata_FieldID", "cell2_number")) %>%
      dplyr::group_by(Image_Metadata_FieldID) %>%
      dplyr::summarize(
        feature_name = cell_feature_column$feature,
        rank_correlation = cor(feature_value1, feature_value2, method="spearman"))
  })

spatial_correlation_summary <- spatial_correlation %>%
  dplyr::group_by(feature_name) %>%
  dplyr::summarize(
    rank_correlation_mean = mean(rank_correlation),
    rank_correlation_sd = sd(rank_correlation)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(rank_correlation_mean))


# plot an example of a high neighbor correlation feature
feature <- "Nuclei_Intensity_MeanIntensityEdge_Hoe"
feature <- "Cells_Intensity_MeanIntensity_NP"
feature <- "Cells_RadialDistribution_ZernikeMagnitude_NP_0_0"
feature <- "Cytoplasm_Intensity_MeanIntensityEdge_NP"
feature <- "Nuclei_Intensity_MeanIntensityEdge_Hoe"
plot_data <- neighbor_graph %>%
  dplyr::left_join(
    cell_features_rel %>%
      dplyr::select(
        plate_id,
        row,
        column,
        Image_Metadata_FieldID,
        cell1_number = Cells_Number_Object_Number,
        feature_value1 = tidyselect::any_of(feature)),
    by = c("plate_id", "row", "column", "Image_Metadata_FieldID", "cell1_number")) %>%
  dplyr::left_join(
    cell_features_rel %>%
      dplyr::select(
        plate_id,
        row,
        column,
        Image_Metadata_FieldID,
        cell2_number = Cells_Number_Object_Number,
        feature_value2 = tidyselect::any_of(feature)),
    by = c("plate_id", "row", "column", "Image_Metadata_FieldID", "cell2_number"))

plot <- ggplot2::ggplot(
  data = plot_data %>%
    dplyr::ungroup() %>%
    dplyr::filter(row == 1, column == 18)) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = feature_value1,
      y = feature_value2),
    size = .5,
    alpha = .1,
    shape = 16) +
  ggplot2::facet_wrap(facets = dplyr::vars(Image_Metadata_FieldID)) +
  ggplot2::coord_fixed() +
  ggplot2::ggtitle(paste0("Neighbor correlation: ", feature)) +
  ggplot2::scale_x_continuous(paste0("Query cell")) +
  ggplot2::scale_y_continuous(paste0("Neighbor cell"))

ggplot2::ggsave(
  filename = paste0("product/figures/TS2/spatial/correlation_uninfected_row1_column18_by_field_", feature, "_202204018.pdf"),
  plot = plot,
  width = 15,
  height = 15,
  useDingbats = FALSE)

ggplot2::ggsave(
  filename = paste0("product/figures/TS2/spatial/correlation_uninfected_row1_column18_by_field_", feature, "_202204018.png"),
  plot = plot,
  width = 9,
  height = 9)

############
feature <- "Cells_AreaShape_Zernike_8_8"
plot_data <- neighbor_graph %>%
  dplyr::left_join(
    cell_features_rel %>%
      dplyr::select(
        plate_id,
        row,
        column,
        Image_Metadata_FieldID,
        cell1_number = Cells_Number_Object_Number,
        feature_value1 = tidyselect::any_of(feature)),
    by = c("plate_id", "row", "column", "Image_Metadata_FieldID", "cell1_number")) %>%
  dplyr::left_join(
    cell_features_rel %>%
      dplyr::select(
        plate_id,
        row,
        column,
        Image_Metadata_FieldID,
        cell2_number = Cells_Number_Object_Number,
        feature_value2 = tidyselect::any_of(feature)),
    by = c("plate_id", "row", "column", "Image_Metadata_FieldID", "cell2_number"))

plot <- ggplot2::ggplot(
  data = plot_data %>%
    dplyr::ungroup()) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = feature_value1,
      y = feature_value2),
    size = .5,
    alpha = .1,
    shape = 16) +
  #ggplot2::facet_wrap(facets = dplyr::vars(row, column)) +
  ggplot2::coord_fixed() +
  ggplot2::ggtitle(paste0("Neighbor correlation: ", feature)) +
  ggplot2::scale_x_continuous(paste0("Query cell")) +
  ggplot2::scale_y_continuous(paste0("Neighbor cell"))

ggplot2::ggsave(
  filename = paste0("product/figures/TS2/spatial/correlation_uninfected_", feature, "_202204018.pdf"),
  plot = plot,
  width = 7,
  height = 7,
  useDingbats = FALSE)

ggplot2::ggsave(
  filename = paste0("product/figures/TS2/spatial/correlation_uninfected_", feature, "_202204018.png"),
  plot = plot,
  width = 7,
  height = 7)



