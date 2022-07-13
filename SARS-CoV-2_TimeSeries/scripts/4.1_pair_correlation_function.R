

library(tidyverse)
library(spatstat)

source("scripts/compute_spatial_features.R")

cell_features_TS2PL1 <- arrow::read_parquet(
  file = "product/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet") %>%
  compute_well_coordinates()

ripleys_k <- cell_features_TS2PL1 %>%
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
    tryCatch({
      pcf <- spatstat.core::Kest(
        data %>%
          dplyr::select(
            field_x,
            field_y) %>%
          spatstat.geom::as.ppp(
            W = spatstat.geom::owin(
              xrange = c(0, field_width_pixels),
              yrange = c(0, field_height_pixels))),
        correction = "isotropic",
        rmax = min(field_width_pixels, field_height_pixels)/2) %>%
        spatstat.core::pcf(
          spar = .5,
          method = "b") # effectively constrain g(0) = 0
      data.frame(
        radius = pcf$r,
        pcf = pcf$pcf)
    }, error = function(e){
      print("Error")
      data.frame()
    })
  }) %>%
  dplyr::ungroup()


plot <- ggplot2::ggplot(
  data = ripleys_k %>%
    dplyr::filter(time_point != "Uninfected")) +
  ggplot2::theme_bw() +
  ggplot2::geom_line(
    mapping = ggplot2::aes(
      x = radius,
      y = pcf,
      group = as.factor(paste(row, column, Image_Metadata_FieldID))),
    alpha = .05) +
  ggplot2::scale_x_sqrt(
    "Cell Neighbor Distance (Pixels)",
    breaks= c(1, 20, 50, 100, 200, 400, 600),
    expand = c(0, 0)) +
  ggplot2::scale_y_continuous(
    "Pair Correlation Function") +
  ggplot2::coord_cartesian(
    xlim = c(0, 800),
    ylim = c(-1, 10)) +
  ggplot2::scale_color_discrete("Field ID") +
  ggplot2::facet_wrap(facets = dplyr::vars(time_point))

ggplot2::ggsave(
  filename = "product/figures/TS2/spatial/pair_correlation_function_TS2PL1_by_field_facet_wrap_no_Uninfected_2022020418.pdf",
  plot = plot,
  width = 8,
  height = 4)


# very consistent location of the peak
ripleys_k %>%
  dplyr::group_by(plate_id, row, column, time_point, Image_Metadata_FieldID) %>%
  dplyr::arrange(desc(pcf)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(time_point) %>%
  dplyr::summarize(
    radius_mean = mean(radius),
    radius_sd = sd(radius))








