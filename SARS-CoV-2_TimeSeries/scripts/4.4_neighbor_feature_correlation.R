
library(tidyverse)
library(arrow)

source("scripts/compute_spatial_features.R")

cell_features_TS2PL1 <- arrow::read_parquet(
  file = "product/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet") %>%
  compute_well_coordinates()

cell_feature_columns <- readr::read_tsv(
  "product/cell_feature_columns_TS_202008.tsv")

cell_metadata_columns <- readr::read_tsv(
  "product/cell_metadata_columns_TS_202008.tsv")


cell_features_rel <- cell_features_TS2PL1 %>%
  dplyr::filter(time_point == "Uninfected")

# generate neighbor graph
neighbor_graph <- cell_features_rel %>%
  compute_neighbor_graph_by_well()


spatial_correlation <- cell_feature_columns %>%
  dplyr::rowwise() %>%
  dplyr::do({
    feature <- .$feature
    cat("Cell feature: ", feature, "\n", sep = "")

    smoothed_feature <- cell_features_rel %>%
      dplyr::group_by(plate_id, row, column) %>%
      dplyr::do({
        data <- .
        cat(
          "plate id: ", data$plate_id[1], " ",
          "row: ", data$row[1], " ",
          "column: ", data$column[1], " ",
          "n cells: ", nrow(data), "\n",
          sep = "")
        model <- mgcv::gam(
          formula = as.formula(paste0(feature, "~ s(well_x, well_y)")),
          data = data,
          method = "REML")
        data.frame(
          well_index = data$well_index,
          feature_value = data[[feature]],
          feature_value_smoothed = residuals(model),
          deviance_explained = summary(model)$dev.expl)
      }) %>%
      dplyr::ungroup()

    feature_graph <- neighbor_graph %>%
      dplyr::left_join(
        smoothed_feature %>%
          dplyr::rename(
            cell1_well_index = well_index,
            cell1_feature_value = feature_value,
            cell1_feature_value_smoothed = feature_value_smoothed),
        by = c("plate_id", "row", "column", "cell1_well_index")) %>%
      dplyr::left_join(
        smoothed_feature %>%
          dplyr::rename(
            cell2_well_index = well_index,
            cell2_feature_value = feature_value,
            cell2_feature_value_smoothed = feature_value_smoothed),
        by = c("plate_id", "row", "column", "cell2_well_index")) %>%
      dplyr::group_by(plate_id, row, column) %>%
        dplyr::summarize(
          feature = feature,
          rank_correlation = cor(
            x = cell1_feature_value,
            y = cell2_feature_value,
            method = "spearman"),
          rank_correlation_smoothed = cor(
            x = cell1_feature_value_smoothed,
            y = cell2_feature_value_smoothed,
            method = "spearman"),
          .groups = "drop") %>%
      dplyr::left_join(
        smoothed_feature %>%
          dplyr::distinct(plate_id, row, column, deviance_explained),
        by = c("plate_id", "row", "column"))
  })

spatial_correlation_summary <- spatial_correlation %>%
  dplyr::group_by(feature) %>%
  dplyr::summarize(
    rank_correlation_mean = mean(rank_correlation),
    rank_correlation_sd = sd(rank_correlation),
    rank_correlation_smoothed_mean = mean(rank_correlation_smoothed),
    rank_correlation_smoothed_sd = sd(rank_correlation_smoothed),
    smooth_deviance_explained_mean = mean(deviance_explained),
    smooth_deviance_explained_sd = sd(deviance_explained)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(rank_correlation_smoothed_mean))


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


