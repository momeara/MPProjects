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



cell_features_rel <- cell_features_TS2PL1 %>%
  dplyr::filter(row == 1, column == 18)


field_location <- cell_features_rel %>%
  dplyr::group_by(
    plate_id,
    row,
    column,
    Image_Metadata_FieldID,
    is_control,
    time_point) %>%
  dplyr::do({
    data <- .
    cat(
      "plate_id: ", data$plate_id[1],
      "row: ", data$row[1],
      "column: ", data$column[1],
      "field: ", data$Image_Metadata_FieldID[1], "\n")
    cell_feature_columns %>%
      dplyr::rowwise() %>%
      dplyr::do({
        feature <- .$feature
        model <- lm(
          formula = paste0(feature, " ~ 0 + Cells_AreaShape_Center_X + Cells_AreaShape_Center_Y"),
          data = data)
        data.frame(
          plate_id = data$plate_id[1],
          row = data$row[1],
          column = data$column[1],
          Image_Metadata_FieldID = data$Image_Metadata_FieldID[1],
          is_control = data$is_control[1],
          time_point = data$time_point[1],
          feature = feature,
          coef_Cells_AreaShape_Center_X = model$coefficients["Cells_AreaShape_Center_X"],
          coef_Cells_AreaShape_Center_Y = model$coefficients["Cells_AreaShape_Center_Y"],
          r.squared = summary(model)$r.squared)
      })
  }) %>%
  dplyr::ungroup()



#################3
#feature <- "Nuclei_Intensity_MeanIntensityEdge_ConA"
feature <- "Nuclei_Intensity_MeanIntensityEdge_Hoe"
plot <- ggplot2::ggplot(cell_features_rel) +
  ggplot2::theme_bw() +
  ggplot2::facet_grid(
    rows = dplyr::vars(field_row),
    cols = dplyr::vars(field_column)) +
  ggplot2::coord_fixed() +
  ggplot2::scale_x_continuous(
    "Field Pixels X",
    limits = c(0, 2560),
    expand = c(0, 0),
    breaks = c(0, 500, 1000, 1500, 2000)) +
  ggplot2::scale_y_continuous("Field Pixels Y", limits = c(0, 2160), expand = c(0, 0)) +
  ggplot2::scale_color_viridis_c("Feature rank in well") +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::theme(panel.spacing = ggplot2::unit(0, "cm")) +
  ggplot2::guides (color = ggplot2::guide_colourbar(barwidth = 15)) +
  ggplot2::ggtitle(paste0("Well distribution for feature: ", feature)) +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = Cells_AreaShape_Center_X,
      y = Cells_AreaShape_Center_Y,
      color = rank(Nuclei_Intensity_MeanIntensityEdge_ConA)),
    size = 3,
    shape = 16,
    alpha = .8)

ggplot2::ggsave(
  filename = paste0("product/figures/TS2/spatial/well_distribution_uninfected_row1_column18_by_field_", feature, "_202204018.pdf"),
  plot = plot,
  width = 9,
  height = 9,
  useDingbats = FALSE)

ggplot2::ggsave(
  filename = paste0("product/figures/TS2/spatial/well_distribution_uninfected_row1_column18_by_field_", feature, "_202204018.png"),
  plot = plot,
  width = 9,
  height = 9)



#####################

model <- mgcv::gam(
  formula = Nuclei_Intensity_LowerQuartileIntensity_NP ~ s(well_x, well_y),
  data = cell_features_rel,
  method = "REML")

smooth_fit <- mgcViz::getViz(model) %>%
  mgcViz::sm(1) %>%
  plot()

plot_data <- cell_features_rel %>%
  dplyr::mutate(
    residuals = model %>% residuals()) %>%
  dplyr::mutate(
    field_row = dplyr::case_when(
      Image_Metadata_FieldID %in% c("0007", "0008", "0009") ~ 1,
      Image_Metadata_FieldID %in% c("0004", "0005", "0006") ~ 2,
      Image_Metadata_FieldID %in% c("0001", "0002", "0003") ~ 3),
    field_column = dplyr::case_when(
      Image_Metadata_FieldID %in% c("0007", "0004", "0001") ~ 1,
      Image_Metadata_FieldID %in% c("0008", "0005", "0002") ~ 2,
      Image_Metadata_FieldID %in% c("0009", "0006", "0003") ~ 3),
    well_x = (field_column-1)*2560 + Cells_AreaShape_Center_X,
    well_y = (3-field_row)*2160 + Cells_AreaShape_Center_Y)


feature <- "Nuclei_Intensity_LowerQuartileIntensity_NP"
plot <- ggplot2::ggplot(plot_data) +
  ggplot2::theme_bw() +
  ggplot2::coord_fixed() +
  ggplot2::scale_x_continuous(
    "Well Pixels X",
    limits = c(0, 2560*3),
    expand = c(0, 0)) +
  ggplot2::scale_y_continuous("Well Pixels Y", limits = c(0, 2160*3), expand = c(0, 0)) +
  ggplot2::scale_color_viridis_c("Feature rank in well") +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::theme(panel.spacing = ggplot2::unit(0, "cm")) +
  ggplot2::guides (color = ggplot2::guide_colourbar(barwidth = 15)) +
  ggplot2::ggtitle(paste0("GAM residauls for feature: ", feature)) +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = well_x,
      y = well_y,
      color = rank(residuals)),
      #color = rank(Nuclei_Intensity_LowerQuartileIntensity_NP)),
    size = 3,
    shape = 16,
    alpha = .8) +
  ggplot2::geom_contour(
    data = smooth_fit$ggObj$data,
    mapping = ggplot2::aes(
      x = x,
      y = y,
      z = tz),
    color = "orange",
    size = 1)
plot

ggplot2::ggsave(
  filename = paste0("product/figures/TS2/spatial/well_distribution_48hpi_row2_column3_gam_residual_", feature, "_202204018.pdf"),
  plot = plot,
  width = 9,
  height = 9,
  useDingbats = FALSE)

ggplot2::ggsave(
  filename = paste0("product/figures/TS2/spatial/well_distribution_48hpi_row2_column3_gam_residual_", feature, "_202204018.png"),
  plot = plot,
  width = 9,
  height = 9)



