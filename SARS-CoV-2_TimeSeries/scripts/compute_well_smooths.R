
library(tidyverse)


# generate neighbor graph
neighbor_graph <- cell_features_rel %>%
  compute_neighbor_graph_by_well(verbose = TRUE)


spatial_correlation <- cell_feature_columns %>%
  dplyr::rowwise() %>%
  dplyr::do({
    feature <- .$feature
    cat("Cell feature: ", feature, "\n", sep = "")

    smoothed_feature <- cell_features_rel %>%
      dplyr::group_by(Plate_ID, row, column) %>%
      dplyr::do({
        data <- .
        cat(
          "plate id: ", data$Plate_ID[1], " ",
          "row: ", data$row[1], " ",
          "column: ", data$column[1], " ",
          "n cells: ", nrow(data), "\n",
          sep = "")
        tryCatch({
          model <- mgcv::gam(
            formula = as.formula(paste0(feature, "~ s(well_x, well_y)")),
            data = data,
            method = "REML")
          data.frame(
            well_index = data$well_index,
            feature_value = data[[feature]],
            feature_value_smoothed = residuals(model),
            deviance_explained = summary(model)$dev.expl)
        }, error = function(m){
          cat("Error:", m$message, "\n", sep = "")
          data.frame(
            well_index = data$well_index,
            feature_value = data[[feature]],
            feature_value_smoothed = NA,
            deviance_explained = NA)
        })
      }) %>%
      dplyr::ungroup()

    feature_graph <- neighbor_graph %>%
      dplyr::left_join(
        smoothed_feature %>%
          dplyr::rename(
            cell1_well_index = well_index,
            cell1_feature_value = feature_value,
            cell1_feature_value_smoothed = feature_value_smoothed),
        by = c("Plate_ID", "row", "column", "cell1_well_index")) %>%
      dplyr::left_join(
        smoothed_feature %>%
          dplyr::rename(
            cell2_well_index = well_index,
            cell2_feature_value = feature_value,
            cell2_feature_value_smoothed = feature_value_smoothed),
        by = c("Plate_ID", "row", "column", "cell2_well_index")) %>%
      dplyr::group_by(Plate_ID, row, column) %>%
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
          dplyr::distinct(Plate_ID, row, column, deviance_explained),
        by = c("Plate_ID", "row", "column"))
  })

spatial_correlation_summary <- spatial_correlation %>%
  dplyr::left_join(
    cell_features %>%
      dplyr::distinct(Plate_ID, row, column, concentration_label),
    by = c("Plate_ID", "row", "column")) %>%
  dplyr::group_by(feature, concentration_label) %>%
  dplyr::summarize(
    rank_correlation_mean = mean(rank_correlation),
    rank_correlation_sd = sd(rank_correlation),
    rank_correlation_smoothed_mean = mean(rank_correlation_smoothed),
    rank_correlation_smoothed_sd = sd(rank_correlation_smoothed),
    smooth_deviance_explained_mean = mean(deviance_explained),
    smooth_deviance_explained_sd = sd(deviance_explained),
    .groups = "drop") %>%
  dplyr::arrange(desc(rank_correlation_mean))



########

cell_features_rel <- cell_features %>%
  dplyr::filter(Concentration > 0)

# plot an example of a high neighbor correlation feature
feature <- "Intensity_MeanIntensityEdge_Virus"
plot_data <- neighbor_graph %>%
  dplyr::group_by(Plate_ID, row, column) %>%
  dplyr::filter(cell1_well_index < cell2_well_index) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(
    cell_features_rel %>%
      dplyr::select(
        Plate_ID,
        row,
        column,
        cell1_well_index = well_index,
        feature_value1 = tidyselect::any_of(feature)),
    by = c("Plate_ID", "row", "column", "cell1_well_index")) %>%
  dplyr::inner_join(
    cell_features_rel %>%
      dplyr::select(
        Plate_ID,
        row,
        column,
        cell2_well_index = well_index,
        feature_value2 = tidyselect::any_of(feature)),
    by = c("Plate_ID", "row", "column", "cell2_well_index")) %>%
  dplyr::mutate(
    concentration_label = paste0("Nic: ", concentration_label, " uM"))

plot <- ggplot2::ggplot(data = plot_data) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position ="bottom") +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = feature_value1,
      y = feature_value2),
    size = .9,
    alpha = .5,
    shape = 16) +
  ggplot2::geom_smooth(
    mapping = ggplot2::aes(
      x = feature_value1,
      y = feature_value2,
      group = paste(row, column)),
    method = "lm",
    alpha = .2,
    size = .5,
    formula = y ~ x) +
  ggplot2::coord_fixed() +
  ggplot2::facet_wrap(
    facets = dplyr::vars(concentration_label),
    nrow = 1) +
  ggplot2::ggtitle(paste0("Neighbor correlation: ", feature)) +
  ggplot2::scale_x_continuous(paste0("Query cell")) +
  ggplot2::scale_y_continuous(paste0("Neighbor cell"))
plot

ggplot2::ggsave(
  filename = paste0("product/figures/correlation_by_nic_", feature, "_202204018.pdf"),
  plot = plot,
  width = 8,
  height = 4,
  useDingbats = FALSE)

ggplot2::ggsave(
  filename = paste0("product/figures/correlation_by_nic_", feature, "_202204018.png"),
  plot = plot,
  width = 8,
  height = 4)
