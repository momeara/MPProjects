

library(tidyverse)



make_scatter_boards <- function(
    features,
    feature_x,
    feature_y,
    feature_color,
    scales,
    output_dir,
    geom_point_args = NULL) {

    features %>%
        dplyr::group_by(plate_id) %>%        
        dplyr::do({
            data <- .
            data <- data %>%
                dplyr::mutate(
                    drug_1_label = reorder(drug_1_label, drug_1_concentration),
                    drug_2_label = reorder(drug_2_label, drug_2_concentration))
            plot <- ggplot2::ggplot() +
                ggplot2::theme_bw() +
                ggplot2::theme(legend.position = "bottom") +
                ggplot2::geom_point(
                    data = data,
                    mapping = ggplot2::aes(
                        x =  !!sym(feature_x),
                        y = !!sym(feature_y),
                        color = !!sym(feature_color)),
                    alpha = .5,
                    size = .4) +
                ggplot2::ggtitle(
                    label = paste0(
                        feature_y, " by ", feature_x),
                    subtitle = paste0("Plate: ", data$plate_id[1])) +
                facet_grid(
                    rows = dplyr::vars(drug_1_label),
                    cols = dplyr::vars(drug_2_label)) +
                scales

            if (!dir.exists(paste0(output_dir, "/", data$plate_id[1]))) {
                cat("Creating output directory '", paste0(output_dir, "/", data$plate_id[1]), "\n", sep = "")
                dir.create(paste0(output_dir, "/", data$plate_id[1]), recursive = TRUE)
            }
            output_basename <- paste0(output_dir, "/", data$plate_id[1], "/scatter_board_", feature_x, "_", feature_y, "_", feature_color)
            cat("Saving plots '", output_basename, ".[pdf/png] ...\n", sep = "")
            ggplot2::ggsave(
                filename = paste0(output_basename, ".pdf"),
                plot = plot,
                width = 14,
                height = 14)
            ggplot2::ggsave(
                filename = paste0(output_basename, ".png"),
                plot = plot,
                width = 14,
                height = 14)
            data.frame()
        })
}
