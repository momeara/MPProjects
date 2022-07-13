
library(tidyverse)


field_width_pixels <- 2560
field_height_pixels <- 2160
well_width_pixels <- 3 * field_width_pixels
well_height_pixels <- 3 * field_width_pixels

compute_well_coordinates <- function(
    cell_features,
    center_x_feature = Cells_AreaShape_Center_X,
    center_y_feature = Cells_AreaShape_Center_Y){
  cell_features %>%
    dplyr::mutate(
      field_x = {{center_x_feature}},
      field_y = {{center_y_feature}},
      field_row = dplyr::case_when(
        Image_Metadata_FieldID %in% c("0007", "0008", "0009") ~ 1,
        Image_Metadata_FieldID %in% c("0004", "0005", "0006") ~ 2,
        Image_Metadata_FieldID %in% c("0001", "0002", "0003") ~ 3),
      field_column = dplyr::case_when(
        Image_Metadata_FieldID %in% c("0007", "0004", "0001") ~ 1,
        Image_Metadata_FieldID %in% c("0008", "0005", "0002") ~ 2,
        Image_Metadata_FieldID %in% c("0009", "0006", "0003") ~ 3),
      well_x = (field_column-1)*field_width_pixels + {{center_x_feature}},
      well_y = (3-field_row)*field_height_pixels + {{center_y_feature}},
      time_point = factor(
        time_point,
          levels = c("Uninfected", "8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours"),
          labels = c("Uninfected", "8 hours", "18 hours", "24 hours", "30 hours", "36 hours", "48 hours"))) %>%
    dplyr::bind_cols(
      cell_features %>%
        dplyr::group_by(
          Image_Metadata_PlateID,
          Image_Metadata_WellID) %>%
        dplyr::transmute(well_index = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Image_Metadata_PlateID, -Image_Metadata_WellID))
}

compute_neighbor_graph_by_well <- function(
  cell_features,
  max_neighbor_distance_pixels = 100,
  verbose = FALSE){

  cell_features %>%
    dplyr::group_by(
      plate_id,
      row,
      column,
      is_control,
      time_point) %>%
    dplyr::do({
      data <- .
      if(verbose){
        cat(
          "plate_id: ", data$plate_id[1], " ",
          "row: ", data$row[1], " ",
          "column: ", data$column[1], " ",
          "n cells: ", nrow(data), "\n",
          sep = "")
      }
      data %>%
        dplyr::transmute(
          cell1_well_index = well_index,
          cell1_well_x = well_x,
          cell1_well_y = well_y,
          cell2_well_index = well_index,
          cell2_well_x = well_x,
          cell2_well_y = well_y) %>%
        tidyr::expand(
          tidyr::nesting(
            cell1_well_index,
            cell1_well_x,
            cell1_well_y),
          tidyr::nesting(
            cell2_well_index,
            cell2_well_x,
            cell2_well_y)) %>%
        dplyr::filter(
          cell1_well_x > max_neighbor_distance_pixels,
          cell1_well_x < well_width_pixels - max_neighbor_distance_pixels,
          cell1_well_y > max_neighbor_distance_pixels,
          cell1_well_y < well_height_pixels - max_neighbor_distance_pixels) %>%
        dplyr::mutate(
          distance = sqrt(
            (cell2_well_x - cell1_well_x)^2 +
            (cell2_well_y - cell1_well_y)^2)) %>%
        dplyr::filter(
          distance > 0,
          distance < max_neighbor_distance_pixels)
    })
}
