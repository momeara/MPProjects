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
  dplyr::filter(time_point == "48 hours") %>%
  dplyr::filter(row == 2, column == 3)
rm(cell_features_TS2PL1)
gc()

base_predictor <- cell_feature_columns$feature %>%
  purrr::map_dfc(function(feature){
    cat("Cell feature: ", feature, "\n", sep = "")
    model <- lm(
      formula = as.formula(paste0(feature, " ~ .")),
      data = cell_features_rel %>%
        dplyr::select(tidyselect::one_of(cell_feature_columns$feature)))

    tibble::tibble(
      {{feature}} := predict(model))
  }) %>%
  dplyr::ungroup()

well_smooths <- cell_features_rel %>%
  dplyr::left_join(
    cell_features_rel %>%
      dplyr::group_by(plate_id, row, column) %>%
      dplyr::do({
        well_data <- .
        cat(
          "plate_id:", well_data$plate_id[1], " ",
          "row:", well_data$row[1], " ",
          "column:",  well_data$column[1], " ",
          "n cells:", nrow(data), "\n",
          sep = "")
        well_data %>%
          dplyr::select(
            plate_id, row, column, well_index) %>%
          dplyr::bind_cols(
            cell_feature_columns$feature %>%
              purrr::map_dfc(function(feature){
              cat("Cell feature: ", feature, "\n", sep = "")
                tryCatch({
                  model <- mgcv::gam(
                    formula = as.formula(paste0(feature, "~ s(well_x, well_y)")),
                    data = well_data,
                    method = "REML")
                  tibble::tibble(
                    {{feature}} := as.vector(predict(model)))
                }, error = function(m){
                  cat("Error:", m$message, "\n", sep = "")
                  tibble::tibble(
                  {{feature}} := NA)
                }) %>%
              })) %>%
      }) %>%
      dplyr::ungroup(),
    by = c("plate_id", "row", "column", "well_index"))


# generate neighbor graph
neighbor_graph <- cell_features_rel %>%
  compute_neighbor_graph_by_well(verbose = TRUE)


predictors <- cell_feature_columns$feature %>%
  purrr::map_dfr(function(feature){
    cat("Cell feature: ", feature, "\n", sep = "")
    cell_features_rel %>%
      dplyr::select(
        plate_id,
        row,
        column,
        well_index = well_index,
        feature_value = tidyselect::any_of(feature)) %>%
      dplyr::bind_cols(
        base_predictor %>%
          dplyr::select(
            feature_value_base = tidyselect::any_of(feature))) %>%
      dplyr::bind_cols(
        well_smooths %>%
          dplyr::select(
            feature_value_gam_pred = tidyselect::any_of(feature))) %>%
      dplyr::left_join(
        neighbor_graph %>%
          dplyr::rename(
            well_index = cell1_well_index,
            neighbor_well_index = cell2_well_index),
        by = c("plate_id", "row", "column", "well_index")) %>%
      dplyr::left_join(
        cell_features_rel %>%
          dplyr::select(
            plate_id,
            row,
            column,
            neighbor_well_index = well_index,
            feature_value_neighbor = tidyselect::any_of(feature)),
        by = c("plate_id", "row", "column", "neighbor_well_index")) %>%
      dplyr::group_by(plate_id, row, column, well_index) %>%
      dplyr::summarize(
        feature = {{feature}},
        feature_value = feature_value[1],
        feature_value_base = feature_value_base[1],
        feature_value_gam_pred = feature_value_gam_pred[1],
        feature_value_neighbor_mean = mean(feature_value_neighbor),
        .groups = "drop")
  })

predictors_summary <- predictors %>%
  dplyr::filter(!is.na(feature_value_neighbor_mean)) %>%
  dplyr::group_by(feature) %>%
  dplyr::summarize(
    rank_correlation_value_base = cor(feature_value, feature_value_base, method = "spearman"),
    rank_correlation_value_neighbor = cor(feature_value, feature_value_neighbor_mean, method = "spearman"),
    rank_correlation_value_gam_pred = cor(feature_value, feature_value_gam_pred, method = "spearman"),
    rank_correlation_smoothed_neighbor = cor(feature_value - feature_value_gam_pred, feature_value_neighbor_mean, method = "spearman"),
    rank_correlation_base_neighbor = cor(feature_value_base, feature_value_neighbor_mean, method = "spearman"),
    rank_correlation_residual_neighbor = cor(feature_value_base - feature_value, feature_value_neighbor_mean, method = "spearman"),
    .groups = "drop") %>%
  dplyr::arrange(desc(rank_correlation_smoothed_neighbor))


feature <- "Cytoplasm_Intensity_MedianIntensity_Hoe"
plot <- GGally::ggpairs(
  data = predictors %>%
    dplyr::filter(feature == {{feature}}),
  columns = c(
    "feature_value",
    "feature_value_gam_pred",
    "feature_value_neighbor_mean"),
  lower = list(
    continuous = GGally::wrap(
      "points",
      alpha = 0.3,
      size=0.1),
    combo = GGally::wrap(
      "dot",
      alpha = 0.4,
      size=0.2)))
plot

feature <- "Cytoplasm_Intensity_MedianIntensity_Hoe"
plot_data <- predictors %>%
  dplyr::filter(feature == {{feature}}) %>%
  dplyr::filter(!is.na(feature_value_neighbor_mean))


plot1 <- ggplot2::ggplot(data = plot_data) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = feature_value_gam_pred,
      y = feature_value),
    alpha = .4,
    size = .8,
    shape = 16) +
  ggplot2::geom_smooth(
    mapping = ggplot2::aes(
      x = feature_value_gam_pred,
      y = feature_value),
    formula = y ~ x,
    method = "gam") +
  ggplot2::scale_x_continuous("Smoothed prediction") +
  ggplot2::scale_y_continuous("Feature value")

plot2 <- ggplot2::ggplot(data = plot_data) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = feature_value_neighbor_mean,
      y = feature_value),
    alpha = .4,
    size = .8,
    shape = 16) +
  ggplot2::geom_smooth(
    mapping = ggplot2::aes(
      x = feature_value_neighbor_mean,
      y = feature_value),
    formula = y ~ x,
    method = "gam") +
  ggplot2::scale_x_continuous("Mean neighbor prediction") +
  ggplot2::scale_y_continuous("Feature value")

plot3 <- ggplot2::ggplot(data = plot_data) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = feature_value_neighbor_mean,
      y = feature_value - feature_value_gam_pred),
    alpha = .4,
    size = 0.8,
    shape = 16) +
  ggplot2::geom_smooth(
    mapping = ggplot2::aes(
      x = feature_value_neighbor_mean,
      y = feature_value - feature_value_gam_pred),
    formula = y ~ x,
    method = "gam") +
  ggplot2::scale_x_continuous("Mean neighbor prediction") +
  ggplot2::scale_y_continuous("Smoothed feature value")

plot <- plot1 + plot2 + plot3 +
  patchwork::plot_annotation(
    title = "Spatial regression",
    subtitle = paste0("SARS-CoV-2 (48hpi): ", feature))

ggplot2::ggsave(
  filename = paste0("product/figures/gam_vs_neighbor_row2_col3_48hpi_feature_", feature, "_20220425.pdf"),
  plot = plot,
  width = 10,
  height = 3.5)
ggplot2::ggsave(
  filename = paste0("product/figures/gam_vs_neighbor_row2_col3_48hpi_feature_", feature, "_20220425.png"),
  plot = plot,
  width = 10,
  height = 3.5)


# y=feature | y=feature  | y=feature - gam
# y=gam     | x=neighbor | x=neighbor


ggplot2::ggplot2(data = plot_data) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    x = feature_value_prediction)


feature <- "Cytoplasm_Intensity_MedianIntensity_Hoe"
ggplot2::ggplot(
  base_vs_neighbor_predictors %>%
  dplyr::filter(feature == {{feature}})) +
  ggplot2::theme_bw() +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      y = feature_value,
      x = feature_value_neighbor_mean),
    size = .7,
    alpha = .6) +
  ggplot2::geom_smooth(
    mapping = ggplot2::aes(
      y = feature_value,
      x = feature_value_neighbor_mean),
    formula = y ~ x,
    method = "lm") +
  #ggplot2::coord_fixed() +
  ggplot2::ggtitle(
    label = "Prediction by neighbor voting",
    subtitle = paste0("SARS-CoV-2 (48hpi): ", feature)) +
  ggplot2::scale_x_log10("Mean neighbor feature value") +
  ggplot2::scale_y_log10("Feature value")

